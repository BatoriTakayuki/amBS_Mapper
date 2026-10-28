"""Persistent serial jobs run as subprocesses, so a running analysis can be cancelled."""
from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import uuid

from .pipeline import validate_explore, validate_request
from .progress import PREFIX

ACTIVE = ("queued", "running")


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=1), encoding="utf-8")
    temporary.replace(path)


class Jobs:
    def __init__(self, data):
        self.data = Path(data)
        self.root = self.data / "results"
        self.root.mkdir(parents=True, exist_ok=True)
        (self.data / "input").mkdir(parents=True, exist_ok=True)
        self.lock = threading.Lock()
        self.wake = threading.Event()
        self.process = None
        self.current = None
        self.jobs = {}
        for file in self.root.glob("*/job.json"):
            try:
                job = json.loads(file.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if job.get("status") == "running":
                job.update(status="failed", error="サーバーが停止したため中断されました。もう一度実行してください。",
                           finished=now())
                atomic_json(file, job)
            self.jobs[job["id"]] = job
        self.thread = threading.Thread(target=self._loop, daemon=True)
        self.thread.start()

    def _save(self, job):
        atomic_json(self.root / job["id"] / "job.json", job)

    def list(self):
        with self.lock:
            return sorted((dict(job) for job in self.jobs.values()), key=lambda job: job["created"], reverse=True)

    def get(self, job_id):
        with self.lock:
            job = self.jobs.get(job_id)
            return dict(job) if job else None

    def submit(self, payload):
        mode = payload.get("mode", "analyze") if isinstance(payload, dict) else None
        if mode == "setup":
            request = dict(mode="setup", title="参照ゲノムの準備")
        elif mode in ("analyze", "check"):
            request = validate_request(payload)
        elif mode == "explore":
            request = validate_explore(payload)
        else:
            raise ValueError("modeが不正です。")
        job_id = uuid.uuid4().hex
        directory = self.root / job_id
        directory.mkdir(parents=True)
        atomic_json(directory / "request.json", request)
        job = dict(id=job_id, mode=request["mode"], title=request.get("title") or "", status="queued",
                   created=now(), started=None, finished=None, error=None, progress={})
        with self.lock:
            self.jobs[job_id] = job
            self._save(job)
        self.wake.set()
        return dict(job)

    def cancel(self, job_id):
        with self.lock:
            job = self.jobs.get(job_id)
            if not job:
                raise KeyError(job_id)
            if job["status"] == "queued":
                job.update(status="cancelled", finished=now())
                self._save(job)
            elif job["status"] == "running" and self.current == job_id and self.process:
                job["cancel_requested"] = True
                try:
                    os.killpg(self.process.pid, signal.SIGTERM)
                except (ProcessLookupError, PermissionError):
                    pass
            return dict(job)

    def delete(self, job_id):
        with self.lock:
            job = self.jobs.get(job_id)
            if not job:
                raise KeyError(job_id)
            if job["status"] in ACTIVE:
                raise ValueError("実行中・待機中のジョブは削除できません。先に中止してください。")
            directory = (self.root / job_id).resolve()
            if directory.parent != self.root.resolve():
                raise ValueError("削除対象が不正です。")
            shutil.rmtree(directory)
            del self.jobs[job_id]

    def _next(self):
        with self.lock:
            queued = sorted((job for job in self.jobs.values() if job["status"] == "queued"),
                            key=lambda job: job["created"])
            if not queued:
                return None
            job = queued[0]
            job.update(status="running", started=now())
            self.current = job["id"]
            self._save(job)
            return dict(job)

    def _loop(self):
        while True:
            job = self._next()
            if not job:
                self.wake.wait(2)
                self.wake.clear()
                continue
            try:
                self._run(job)
            except Exception as error:  # keep the queue alive whatever happens in one job
                with self.lock:
                    self.jobs[job["id"]].update(status="failed", error=f"内部エラー: {error}", finished=now())
                    self._save(self.jobs[job["id"]])
            finally:
                with self.lock:
                    self.current, self.process = None, None

    def _run(self, job):
        directory = self.root / job["id"]
        if job["mode"] == "setup":
            command = [sys.executable, "-m", "ambs_mapper", "setup", "--data", str(self.data)]
        else:
            command = [sys.executable, "-m", "ambs_mapper", "run", str(directory), "--data", str(self.data)]
        last_error = None
        with open(directory / "pipeline.log", "a", encoding="utf-8") as log:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                       encoding="utf-8", errors="replace", start_new_session=True)
            with self.lock:
                self.process = process
            for line in process.stdout:
                if line.startswith(PREFIX):
                    try:
                        progress = json.loads(line[len(PREFIX):])
                    except ValueError:
                        continue
                    with self.lock:
                        current = self.jobs[job["id"]]
                        current["progress"] = {**current.get("progress", {}), **progress}
                        if progress.get("message"):
                            current["message"] = progress["message"]
                    if progress.get("message"):
                        log.write(f"[{now()[11:19]}] {progress['message']}\n")
                        log.flush()
                    continue
                log.write(line)
                log.flush()
                if line.startswith("エラー: "):
                    last_error = line[len("エラー: "):].strip()
            code = process.wait()
        with self.lock:
            current = self.jobs[job["id"]]
            if current.pop("cancel_requested", False):
                current.update(status="cancelled", finished=now())
            elif code == 0:
                current.update(status="completed", finished=now())
            else:
                current.update(status="failed", finished=now(),
                               error=last_error or f"解析が異常終了しました（終了コード{code}）。ログを確認してください。")
            self._save(current)

"""Local browser application, implemented with the Python standard library."""
from __future__ import annotations

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import io
import json
import mimetypes
import os
from pathlib import Path
import re
import shutil
import signal
import threading
from urllib.parse import quote, urlsplit
import zipfile

from . import __version__, igv_view
from .caller import DEFAULT_PARAMS
from .jobs import Jobs
from .pipeline import list_inputs
from .bsmap import index_ready
from .reference import reference_ready
from .variants import strain_list

STATIC = Path(__file__).with_name("static")
MAX_BODY = 1024 * 1024
JOB = r"([a-f0-9]{32})"
INLINE = {".png": "image/png", ".json": "application/json; charset=utf-8"}


def handler(jobs):
    data = jobs.data
    resources = data / "resources"

    class Handler(BaseHTTPRequestHandler):
        server_version = "AmBSMapper"

        def setup(self):
            super().setup()
            self.connection.settimeout(60)

        def log_message(self, format, *args):
            pass

        def _headers(self, status, content_type, length, filename=None):
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(length))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Referrer-Policy", "same-origin")
            # igv.js sets inline styles and draws through blob: images
            self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; "
                             "style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; "
                             "object-src 'none'; base-uri 'none'; frame-ancestors 'none'")
            if filename:
                self.send_header("Content-Disposition", f"attachment; filename*=UTF-8''{quote(filename)}")
            self.end_headers()

        def send(self, status, value, content_type="application/json; charset=utf-8", filename=None):
            body = value if isinstance(value, bytes) else json.dumps(value, ensure_ascii=False).encode("utf-8")
            self._headers(status, content_type, len(body), filename)
            try:
                self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                pass

        def _file(self, file, base, download):
            file = file.resolve()
            if not file.is_relative_to(base.resolve()) or not file.is_file():
                return self.send(404, {"error": "ファイルが見つかりません。"})
            inline = INLINE.get(file.suffix.lower())
            content_type = inline if inline and not download else "application/octet-stream"
            try:
                with file.open("rb") as handle:
                    self._headers(200, content_type, os.fstat(handle.fileno()).st_size,
                                  None if inline and not download else file.name)
                    shutil.copyfileobj(handle, self.wfile, length=1024 * 1024)
            except (BrokenPipeError, ConnectionResetError):
                pass

        def do_GET(self):
            try:
                self._get()
            except (OSError, ValueError, RuntimeError) as error:
                self.send(500, {"error": str(error)})

        def _get(self):
            url = urlsplit(self.path)
            path = url.path
            if path in ("/", "/app.js", "/style.css", "/igv-panel.js", "/vendor/igv-3.0.0.min.js"):
                file = STATIC / ("index.html" if path == "/" else path[1:])
                mime = mimetypes.guess_type(file.name)[0] or "text/plain"
                return self.send(200, file.read_bytes(), mime + "; charset=utf-8")
            if path == "/api/health":
                tools = ("bcftools", "samtools")
                return self.send(200, dict(app="ambs-mapper", version=__version__, logical_cpus=os.cpu_count(),
                                           tools={name: bool(shutil.which(name)) for name in tools},
                                           reference_ready=reference_ready(resources) and index_ready(resources),
                                           disk_free_gb=round(shutil.disk_usage(data).free / 1024 ** 3, 1),
                                           defaults=DEFAULT_PARAMS))
            if path == "/api/inputs":
                files = list_inputs(data / "input")
                folders = sorted({str(Path(item["path"]).parent.as_posix()) for item in files} - {"."})
                return self.send(200, dict(files=files, folders=folders))
            if path == "/api/strains":
                try:
                    return self.send(200, dict(strains=strain_list(resources)))
                except RuntimeError as error:
                    return self.send(503, {"error": str(error)})
            if path == "/api/jobs":
                return self.send(200, jobs.list())
            igv_match = re.fullmatch(rf"/api/jobs/{JOB}/igv/([^/]+)", path)
            if igv_match:
                return self._igv(*igv_match.groups())
            match = re.fullmatch(rf"/api/jobs/{JOB}(?:/(request|log|results|zip|files/(.+)))?", path)
            if not match:
                return self.send(404, {"error": "Not found"})
            job_id, action, name = match.groups()
            job = jobs.get(job_id)
            if not job:
                return self.send(404, {"error": "ジョブが見つかりません。"})
            directory = jobs.root / job_id
            if not action:
                return self.send(200, job)
            if action == "request":
                return self.send(200, json.loads((directory / "request.json").read_text(encoding="utf-8")))
            if action == "log":
                file = directory / "pipeline.log"
                if not file.is_file():
                    return self.send(200, {"text": ""})
                with file.open("rb") as handle:
                    handle.seek(max(0, file.stat().st_size - 60000))
                    return self.send(200, {"text": handle.read().decode("utf-8", errors="replace")})
            if action == "results":
                file = directory / "output" / "results.json"
                if job["status"] != "completed" or not file.is_file():
                    return self.send(409, {"error": "結果はまだありません。"})
                return self.send(200, file.read_bytes())
            if action == "zip":
                if job["status"] != "completed":
                    return self.send(409, {"error": "結果はまだありません。"})
                buffer = io.BytesIO()
                with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
                    for file in sorted((directory / "output").rglob("*")):
                        # igv/ holds viewer copies of the per-sample BAMs
                        if file.is_file() and file.relative_to(directory / "output").parts[0] != "igv":
                            archive.write(file, file.relative_to(directory / "output").as_posix())
                    archive.write(directory / "request.json", "request.json")
                    if (directory / "pipeline.log").is_file():
                        archive.write(directory / "pipeline.log", "pipeline.log")
                title = re.sub(r"[^0-9A-Za-z_\-ぁ-んァ-ヶ一-龠]+", "_", job.get("title") or "") or job_id[:8]
                return self.send(200, buffer.getvalue(), "application/zip", f"ambs_{title}.zip")
            download = "download=1" in (url.query or "")
            return self._file(directory / "output" / name, directory / "output", download)

        def _igv(self, job_id, asset):
            """Embedded IGV: its configuration, or a file served with byte ranges."""
            job = jobs.get(job_id)
            output = jobs.root / job_id / "output"
            if not job or job["status"] != "completed" or not (output / "igv" / "tracks.json").is_file():
                return self.send(404, {"error": "この解析にはIGV表示用のデータがありません（v0.5.0より前の解析は再解析してください）。"})
            if asset == "config":
                config, alignments, amplicons = igv_view.configuration(output, job_id)
                return self.send(200, dict(config=config, alignments=alignments, amplicons=amplicons))
            file = igv_view.assets(output, resources).get(asset)
            if file is None or not file.is_file():
                return self.send(404, {"error": "IGVのファイルが見つかりません。"})
            with file.open("rb") as handle:
                size = os.fstat(handle.fileno()).st_size
                start, end = 0, size - 1
                partial = self.headers.get("Range")
                if partial:
                    try:
                        start, end = igv_view.byte_range(partial, size)
                    except ValueError:
                        self.send_response(416)
                        self.send_header("Content-Range", f"bytes */{size}")
                        self.send_header("Content-Length", "0")
                        self.end_headers()
                        return
                self.send_response(206 if partial else 200)
                self.send_header("Content-Type", "application/octet-stream")
                self.send_header("Accept-Ranges", "bytes")
                self.send_header("Content-Length", str(end - start + 1))
                self.send_header("Cache-Control", "no-store")
                if partial:
                    self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
                self.end_headers()
                handle.seek(start)
                remaining = end - start + 1
                try:
                    while remaining:
                        chunk = handle.read(min(1024 * 1024, remaining))
                        if not chunk:
                            break
                        self.wfile.write(chunk)
                        remaining -= len(chunk)
                except (BrokenPipeError, ConnectionResetError):
                    pass

        def do_POST(self):
            origin = self.headers.get("Origin")
            if origin and urlsplit(origin).netloc != self.headers.get("Host"):
                return self.send(403, {"error": "同じ画面から操作してください。"})
            if self.headers.get("Sec-Fetch-Site") == "cross-site":
                return self.send(403, {"error": "同じ画面から操作してください。"})
            path = urlsplit(self.path).path
            try:
                length = int(self.headers.get("Content-Length", "0"))
                if not 0 <= length <= MAX_BODY:
                    return self.send(413, {"error": "入力が大きすぎます。"})
                payload = json.loads(self.rfile.read(length) or b"{}")
                if path == "/api/jobs":
                    return self.send(202, jobs.submit(payload))
                match = re.fullmatch(rf"/api/jobs/{JOB}/(cancel|delete)", path)
                if not match:
                    return self.send(404, {"error": "Not found"})
                job_id, action = match.groups()
                if action == "cancel":
                    return self.send(200, jobs.cancel(job_id))
                jobs.delete(job_id)
                return self.send(200, {"deleted": job_id})
            except KeyError:
                return self.send(404, {"error": "ジョブが見つかりません。"})
            except (ValueError, TypeError) as error:
                return self.send(400, {"error": str(error)})
            except OSError as error:
                return self.send(500, {"error": f"保存できませんでした: {error}"})

    return Handler


def serve(host, port, data):
    jobs = Jobs(data)
    server = ThreadingHTTPServer((host, port), handler(jobs))
    server.daemon_threads = True
    if threading.current_thread() is threading.main_thread():
        def terminate(signum, frame):
            raise KeyboardInterrupt
        signal.signal(signal.SIGTERM, terminate)
    print(f"amBS Mapper: http://{host}:{server.server_port}", flush=True)
    try:
        server.serve_forever(poll_interval=0.2)
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()

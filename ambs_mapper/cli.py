"""Command line entry point: web server, analysis runner and reference setup."""
import argparse
import os
from pathlib import Path
import sys

DATA = Path(os.environ.get("AMBS_DATA", "./data"))


def main(argv=None):
    parser = argparse.ArgumentParser(prog="ambs_mapper", description="Amplicon bisulfite-seq mapper")
    commands = parser.add_subparsers(dest="command", required=True)
    web = commands.add_parser("web", help="ブラウザ画面を起動")
    web.add_argument("--host", default="127.0.0.1")
    web.add_argument("--port", type=int, default=8501)
    run = commands.add_parser("run", help="ジョブディレクトリのrequest.jsonを解析")
    run.add_argument("job_dir")
    setup = commands.add_parser("setup", help="参照ゲノムを準備")
    for command in (web, run, setup):
        command.add_argument("--data", default=str(DATA))
    args = parser.parse_args(argv)
    data = Path(args.data)
    resources = data / "resources"
    if args.command == "web":
        from .webapp import serve
        serve(args.host, args.port, data)
    elif args.command == "setup":
        import os
        from .bsmap import build_index
        from .reference import setup_reference
        from .variants import strain_list
        setup_reference(resources)
        build_index(resources, threads=min(16, os.cpu_count() or 1))
        strain_list(resources)
    else:
        from .pipeline import run as run_pipeline
        try:
            run_pipeline(args.job_dir, data / "input", resources)
        except (ValueError, RuntimeError) as error:
            print("エラー: " + str(error), flush=True)
            sys.exit(2)

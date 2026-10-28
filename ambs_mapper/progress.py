"""Machine-readable progress lines on stdout, parsed by the web job runner."""
import json
import sys

PREFIX = "@@progress "


def emit(stage, **values):
    print(PREFIX + json.dumps(dict(stage=stage, **values), ensure_ascii=False), file=sys.stdout, flush=True)


def log(message):
    print(message, file=sys.stdout, flush=True)

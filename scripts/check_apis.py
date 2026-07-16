#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
import urllib.request
import urllib.error

def get_json(url: str) -> dict:
    try:
        with urllib.request.urlopen(url.rstrip("/") + "/health", timeout=15) as r:
            return json.loads(r.read().decode("utf-8"))
    except Exception as exc:
        raise SystemExit(f"Health check failed for {url}: {exc}") from exc

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rnaseq-url", required=True)
    p.add_argument("--network-url", required=True)
    p.add_argument("--llm-url", required=True)
    a = p.parse_args()
    payload = {
        "rnaseq": get_json(a.rnaseq_url),
        "network": get_json(a.network_url),
        "llm_triage": get_json(a.llm_url),
    }
    print(json.dumps(payload, indent=2))

if __name__ == "__main__":
    main()

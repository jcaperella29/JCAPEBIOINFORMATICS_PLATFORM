#!/usr/bin/env python3
import json
import sys
import urllib.request

result = {}

for base in sys.argv[1:]:
    try:
        with urllib.request.urlopen(base.rstrip("/") + "/health", timeout=10) as r:
            result[base] = {"ok": True, "status": r.status}
    except Exception as e:
        print(f"API not reachable at {base}: {e}", file=sys.stderr)
        sys.exit(2)

print(json.dumps(result, indent=2))

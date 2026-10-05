#!/usr/bin/env python3
import argparse, json, sys, urllib.request


def check(name, base):
    url = base.rstrip('/') + '/health'
    try:
        with urllib.request.urlopen(url, timeout=10) as r:
            return {'name': name, 'url': base, 'ok': True, 'status': r.status}
    except Exception as e:
        print(f'{name} API not reachable at {base}: {e}', file=sys.stderr)
        raise SystemExit(2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rnaseq-url', required=True)
    ap.add_argument('--network-url', required=True)
    ap.add_argument('--llm-url', required=True)
    a = ap.parse_args()
    result = {
        'rnaseq': check('RNA-seq', a.rnaseq_url),
        'network': check('Network', a.network_url),
        'llm': check('LLM', a.llm_url),
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

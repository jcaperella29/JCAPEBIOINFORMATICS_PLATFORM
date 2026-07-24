#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path
from typing import List, Optional

from dotenv import load_dotenv
from openai import OpenAI


def eprint(*args):
    print(*args, file=sys.stderr)


def discover_md_files(playbook_dir: Path) -> List[Path]:
    if not playbook_dir.exists() or not playbook_dir.is_dir():
        raise FileNotFoundError(f"Playbook directory not found: {playbook_dir}")
    files = sorted(playbook_dir.glob("*.md"))
    if not files:
        raise FileNotFoundError(f"No .md files found in: {playbook_dir}")
    return files


def upload_files(client: OpenAI, files: List[Path]) -> List[str]:
    file_ids: List[str] = []
    for p in files:
        with p.open("rb") as f:
            # purpose must be "assistants" for file_search/vector store usage
            up = client.files.create(file=f, purpose="assistants")
        file_ids.append(up.id)
        print(f"Uploaded: {p.name} -> file_id={up.id}")
    return file_ids


def ensure_vector_store(client: OpenAI, name: str, vector_store_id: Optional[str]) -> str:
    if vector_store_id:
        # verify it exists
        vs = client.vector_stores.retrieve(vector_store_id)
        print(f"Using existing vector store: {vs.id} (name={getattr(vs, 'name', 'n/a')})")
        return vs.id

    vs = client.vector_stores.create(name=name)
    print(f"Created vector store: {vs.id} (name={name})")
    return vs.id


def attach_files(client: OpenAI, vector_store_id: str, file_ids: List[str]) -> None:
    for fid in file_ids:
        # attaches file into vector store for file_search
        res = client.vector_stores.files.create(vector_store_id=vector_store_id, file_id=fid)
        print(f"Attached: file_id={fid} -> vs_file_id={res.id}")


def wait_until_completed(client: OpenAI, vector_store_id: str, *, timeout_s: int = 600, poll_s: int = 2) -> None:
    """
    Poll until all vector store files are 'completed' or timeout.
    """
    t0 = time.time()
    while True:
        lst = client.vector_stores.files.list(vector_store_id=vector_store_id)

        statuses = [getattr(x, "status", None) for x in getattr(lst, "data", [])]
        # Statuses can include: in_progress, completed, failed (per API)
        in_prog = sum(1 for s in statuses if s == "in_progress")
        failed = [s for s in statuses if s == "failed"]
        completed = sum(1 for s in statuses if s == "completed")

        print(f"Status: completed={completed} in_progress={in_prog} failed={len(failed)}")

        if failed:
            raise RuntimeError("One or more files failed to index in the vector store.")
        if in_prog == 0:
            return
        if time.time() - t0 > timeout_s:
            raise TimeoutError(f"Timed out waiting for vector store indexing after {timeout_s}s.")
        time.sleep(poll_s)


def main():
    ap = argparse.ArgumentParser(description="Index playbook markdown files into an OpenAI vector store.")
    ap.add_argument(
        "--playbook-dir",
        required=False,
        default=None,
        help="Path to playbook directory containing *.md files.",
    )
    ap.add_argument(
        "--vector-store-id",
        required=False,
        default=None,
        help="Reuse an existing vector store id instead of creating a new one.",
    )
    ap.add_argument(
        "--vector-store-name",
        required=False,
        default="enrichment_playbook",
        help="Name for the vector store (if creating).",
    )
    ap.add_argument(
        "--timeout-s",
        required=False,
        type=int,
        default=600,
        help="Timeout for indexing/polling.",
    )
    args = ap.parse_args()

    # Load environment variables (explicit path if you want)
    # If you already solved dotenv path issues, you can keep it simple:
    load_dotenv(override=True)

    # Make sure the key exists
    if not os.environ.get("OPENAI_API_KEY"):
        raise RuntimeError(
            "OPENAI_API_KEY not found in environment. "
            "Ensure your .env is in the working directory or set it explicitly."
        )

    # Default playbook path based on your screenshot:
    # Windows: C:\Users\jcape\Downloads\Enrichment_results_analyis_app\playbook
    # WSL:     /mnt/c/Users/jcape/Downloads/Enrichment_results_analyis_app/playbook
    if args.playbook_dir:
        playbook_dir = Path(args.playbook_dir).expanduser().resolve()
    else:
        playbook_dir = Path("/mnt/c/Users/jcape/Downloads/Enrichment_results_analyis_app/playbook").resolve()

    files = discover_md_files(playbook_dir)
    print(f"Found {len(files)} markdown files in {playbook_dir}")

    client = OpenAI()

    file_ids = upload_files(client, files)
    vs_id = ensure_vector_store(client, name=args.vector_store_name, vector_store_id=args.vector_store_id)

    attach_files(client, vs_id, file_ids)

    print("Waiting for vector store indexing to complete...")
    wait_until_completed(client, vs_id, timeout_s=args.timeout_s)

    print("\n✅ Done.")
    print(f"VECTOR_STORE_ID={vs_id}")
    print("Add this to your .env:")
    print(f"VECTOR_STORE_ID={vs_id}")


if __name__ == "__main__":
    main()

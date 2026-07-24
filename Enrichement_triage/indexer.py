#!/usr/bin/env python3
"""
Create/update an OpenAI Vector Store for file_search RAG.

Usage examples:
  python rag/indexer.py --name "enrich-knowledge" --docs ./rag_docs
  python rag/indexer.py --vector-store-id vs_123 --docs ./rag_docs

Env:
  OPENAI_API_KEY=...
"""

from __future__ import annotations
import os
import time
import argparse
from pathlib import Path
from typing import List, Optional, Tuple

from openai import OpenAI

SUPPORTED_EXT = {".pdf", ".txt", ".md", ".html", ".htm", ".docx"}


def iter_files(docs_dir: Path) -> List[Path]:
    files = []
    for p in docs_dir.rglob("*"):
        if p.is_file() and p.suffix.lower() in SUPPORTED_EXT:
            files.append(p)
    return sorted(files)


def create_or_get_vector_store(client: OpenAI, name: str, vector_store_id: Optional[str]) -> str:
    if vector_store_id:
        # Trust user input; we’ll error naturally if invalid
        return vector_store_id

    vs = client.vector_stores.create(name=name)  # POST /v1/vector_stores :contentReference[oaicite:1]{index=1}
    return vs.id


def upload_file(client: OpenAI, path: Path) -> str:
    # Upload a "File" object first, then attach it to vector store.
    # This is the canonical flow described in the File Search guide. :contentReference[oaicite:2]{index=2}
    with path.open("rb") as f:
        up = client.files.create(file=f, purpose="assistants")
    return up.id


def attach_file_to_vector_store(client: OpenAI, vector_store_id: str, file_id: str):
    # POST /v1/vector_stores/{vector_store_id}/files :contentReference[oaicite:3]{index=3}
    return client.vector_stores.files.create(vector_store_id=vector_store_id, file_id=file_id)


def poll_vector_store_ready(client: OpenAI, vector_store_id: str, *, sleep_s: float = 2.0, timeout_s: int = 1800):
    """
    Wait until vector store has finished processing files (no in_progress).
    The docs recommend polling until files are out of in_progress. :contentReference[oaicite:4]{index=4}
    """
    t0 = time.time()
    while True:
        if time.time() - t0 > timeout_s:
            raise TimeoutError(f"Vector store not ready after {timeout_s}s: {vector_store_id}")

        # list files + check status
        page = client.vector_stores.files.list(vector_store_id=vector_store_id, limit=100)
        items = list(page.data)

        # if more than 100 files, paginate
        while getattr(page, "has_more", False):
            page = client.vector_stores.files.list(
                vector_store_id=vector_store_id,
                limit=100,
                after=page.data[-1].id,
            )
            items.extend(page.data)

        statuses = [it.status for it in items]  # status is typically "in_progress"/"completed"/"failed"
        n_total = len(statuses)
        n_inprog = sum(1 for s in statuses if s == "in_progress")
        n_failed = sum(1 for s in statuses if s == "failed")
        n_done = sum(1 for s in statuses if s == "completed")

        print(f"[poll] total={n_total} completed={n_done} in_progress={n_inprog} failed={n_failed}")

        if n_failed > 0:
            # surface failed file ids for debugging
            failed_ids = [it.id for it in items if it.status == "failed"]
            raise RuntimeError(f"Some vector store files failed: {failed_ids}")

        if n_total > 0 and n_inprog == 0:
            return  # ready

        time.sleep(sleep_s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docs", required=True, help="Directory containing docs to index (pdf/txt/md/html/docx).")
    ap.add_argument("--name", default="enrichment-rag", help="Vector store name (used only if creating new).")
    ap.add_argument("--vector-store-id", default=None, help="Existing vector store id to reuse.")
    ap.add_argument("--no-poll", action="store_true", help="Do not poll for completion.")
    args = ap.parse_args()

    docs_dir = Path(args.docs).resolve()
    if not docs_dir.exists() or not docs_dir.is_dir():
        raise FileNotFoundError(f"--docs must be a directory: {docs_dir}")

    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

    files = iter_files(docs_dir)
    if not files:
        raise RuntimeError(f"No supported docs found in {docs_dir}. Supported: {sorted(SUPPORTED_EXT)}")

    vs_id = create_or_get_vector_store(client, args.name, args.vector_store_id)
    print(f"[vector_store] id={vs_id}")

    # Upload + attach
    for p in files:
        print(f"[upload] {p}")
        file_id = upload_file(client, p)
        attach = attach_file_to_vector_store(client, vs_id, file_id)
        print(f"[attach] file_id={file_id} -> vs_file_id={attach.id}")

    if not args.no_poll:
        poll_vector_store_ready(client, vs_id)

    print("\n✅ Done.")
    print(f"Set this env var for your app:\n  OPENAI_VECTOR_STORE_ID={vs_id}")


if __name__ == "__main__":
    main()

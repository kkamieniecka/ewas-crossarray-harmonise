#!/usr/bin/env python3
"""Fetch the GSE237561 IDATs named in a sample sheet, verifying each one.

The series tarball (GSE237561_RAW.tar, ~1.3 GB) stalls repeatedly on this
connection, so files are pulled per sample from the GEO supplementary paths
the sheet already records in grn_url / red_url.

Completeness is checked with gzip.decompress rather than by file size: a
truncated IDAT is not rejected by minfi's reader, it surfaces much later as
`argument is of length zero` inside readIDAT, which is an expensive way to
discover a bad download.

Usage:
  python3 tools/fetch_idats.py --sheet data/sample_sheets/sample_sheet_GSE237561.csv \
      --out-dir data/idat --jobs 6
"""
from __future__ import annotations

import argparse, csv, gzip, os, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor


def ok(path: str) -> bool:
    """A file counts as present only if it is a complete gzip stream."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    try:
        with open(path, "rb") as fh:
            d = gzip.GzipFile(fileobj=fh)
            while d.read(1 << 20):
                pass
        return True
    except (OSError, EOFError, gzip.BadGzipFile):
        return False


def fetch(url: str, out_dir: str, tries: int = 4) -> tuple[str, bool]:
    url = url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
    dst = os.path.join(out_dir, os.path.basename(url))
    if ok(dst):
        return dst, True
    tmp = dst + ".part"
    for _ in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as fh:
                while chunk := r.read(1 << 20):
                    fh.write(chunk)
            os.replace(tmp, dst)
            if ok(dst):
                return dst, True
        except Exception:
            pass
    if os.path.exists(tmp):
        os.remove(tmp)
    return dst, False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    urls = []
    with open(a.sheet) as fh:
        for row in csv.DictReader(fh):
            urls += [row["grn_url"], row["red_url"]]
    todo = [u for u in urls if not ok(os.path.join(a.out_dir, os.path.basename(u)))]
    print(f"urls={len(urls)} already_complete={len(urls) - len(todo)} todo={len(todo)}", flush=True)

    bad = []
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for i, (dst, good) in enumerate(ex.map(lambda u: fetch(u, a.out_dir), todo), 1):
            if not good:
                bad.append(dst)
            if i % 20 == 0:
                print(f"  {i}/{len(todo)}", flush=True)

    have = sum(1 for u in urls if ok(os.path.join(a.out_dir, os.path.basename(u))))
    print(f"have={have}/{len(urls)} failed={len(bad)}", flush=True)
    for b in bad:
        print("FAIL", b, flush=True)
    return 0 if have == len(urls) else 1


if __name__ == "__main__":
    sys.exit(main())

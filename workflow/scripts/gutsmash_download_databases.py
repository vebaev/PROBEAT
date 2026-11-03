#!/usr/bin/env python3
"""Run gutSMASH's download_databases with HTTPS URLs.

The pinned gutSMASH source downloads Pfam over ftp:// and Resfams over http://,
which fails on HTTPS-only networks. The same files are served over HTTPS and the
upstream checksums are still verified, so only the URL scheme is rewritten.
"""
import sys

import antismash.download_databases as downloads


def https(url: str) -> str:
    for scheme in ("ftp://", "http://"):
        if url.startswith(scheme):
            return "https://" + url[len(scheme):]
    return url


for name in ("PFAM27_URL", "PFAM_LATEST_URL", "RESFAM_URL", "CLUSTERBLAST_URL"):
    if hasattr(downloads, name):
        setattr(downloads, name, https(getattr(downloads, name)))

sys.argv[0] = "download_databases.py"
downloads._main()

#!/usr/bin/env python3
import csv
import json
import sys
from pathlib import Path


session, output = map(Path, sys.argv[1:3])
source = session / "00.OverviewGeneTables.json"
rows = []
if source.exists():
    data = json.loads(source.read_text(errors="replace"))
    rows = data.get("ResultsTable") or []

fields = ["AOI", "start", "end", "class", "fullname"]
with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fields})

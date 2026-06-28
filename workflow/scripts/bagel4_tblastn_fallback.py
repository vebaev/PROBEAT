#!/usr/bin/env python3
import csv
import subprocess
import sys
from pathlib import Path


assembly, database_dir, tblastn, output = map(Path, sys.argv[1:5])
context = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
databases = [
    ("Class I", database_dir / "db_proteins/bagel4_class1bacteriocin_db.fasta"),
    ("Class II", database_dir / "db_proteins/bagel4_class2bacteriocin_db.fasta"),
    ("Class III", database_dir / "db_proteins/bagel4_class3bacteriocin_db.fasta"),
]

hits = []
for bacteriocin_class, database in databases:
    command = [
        str(tblastn), "-query", str(database), "-subject", str(assembly),
        "-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore",
        "-evalue", "1e-5", "-max_target_seqs", "1", "-num_threads", "4",
    ]
    result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 11:
            continue
        query, contig = fields[0], fields[1]
        identity, length, query_length = float(fields[2]), int(fields[3]), int(fields[4])
        start, end = sorted((int(fields[7]), int(fields[8])))
        coverage = 100 * length / query_length if query_length else 0
        if identity < 30 or coverage < 50:
            continue
        hits.append({
            "contig": contig, "start": start, "end": end, "class": bacteriocin_class,
            "name": query, "identity": identity, "coverage": coverage,
            "evalue": fields[9], "bitscore": float(fields[10]),
        })

hits.sort(key=lambda row: (row["contig"], row["start"], -row["bitscore"]))
clusters = []
for hit in hits:
    if clusters and clusters[-1]["contig"] == hit["contig"] and hit["start"] <= clusters[-1]["end"] + context:
        cluster = clusters[-1]
        cluster["end"] = max(cluster["end"], hit["end"] + context)
        cluster["hits"].append(hit)
    else:
        clusters.append({
            "contig": hit["contig"],
            "start": max(1, hit["start"] - context),
            "end": hit["end"] + context,
            "hits": [hit],
        })

fields = [
    "AOI", "start", "end", "class", "fullname", "identity", "coverage", "evalue",
    "hit_start", "hit_end", "source",
]
with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for index, cluster in enumerate(clusters, 1):
        best = max(cluster["hits"], key=lambda row: row["bitscore"])
        writer.writerow({
            "AOI": f"{cluster['contig']}.AOI_{index:02d}",
            "start": cluster["start"],
            "end": cluster["end"],
            "class": best["class"],
            "fullname": best["name"],
            "identity": f"{best['identity']:.2f}",
            "coverage": f"{best['coverage']:.2f}",
            "evalue": best["evalue"],
            "hit_start": best["start"],
            "hit_end": best["end"],
            "source": "BAGEL4 tblastn fallback",
        })

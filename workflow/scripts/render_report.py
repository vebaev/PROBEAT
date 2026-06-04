#!/usr/bin/env python3
import csv
import datetime as dt
import html
import os
import sys
from pathlib import Path

import yaml


def read_samples(path):
    if not Path(path).exists():
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def exists(path):
    return Path(path).exists()


def status(path):
    return "done" if exists(path) else "missing"


def write_minimal_pdf(path, title, lines):
    text = [title, ""] + lines
    stream_lines = ["BT", "/F1 11 Tf", "50 790 Td"]
    for i, line in enumerate(text[:60]):
        safe = line.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
        if i:
            stream_lines.append("0 -15 Td")
        stream_lines.append(f"({safe}) Tj")
    stream_lines.append("ET")
    stream = "\n".join(stream_lines).encode("latin-1", "replace")
    objects = []
    objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    objects.append(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
    objects.append(b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>")
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    objects.append(b"<< /Length %d >>\nstream\n" % len(stream) + stream + b"\nendstream")
    out = [b"%PDF-1.4\n"]
    offsets = [0]
    for i, obj in enumerate(objects, start=1):
        offsets.append(sum(len(x) for x in out))
        out.append(b"%d 0 obj\n" % i + obj + b"\nendobj\n")
    xref = sum(len(x) for x in out)
    out.append(b"xref\n0 %d\n0000000000 65535 f \n" % (len(objects) + 1))
    for off in offsets[1:]:
        out.append(b"%010d 00000 n \n" % off)
    out.append(b"trailer << /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n" % (len(objects) + 1, xref))
    Path(path).write_bytes(b"".join(out))


def main():
    config_path, samples_path, html_out, pdf_out = sys.argv[1:5]
    config = yaml.safe_load(Path(config_path).read_text())
    samples = read_samples(samples_path)
    title = config.get("params", {}).get("report_title", "Bacterial Genome Analysis Report")
    now = dt.datetime.now().isoformat(timespec="seconds")

    rows = []
    for sample in samples:
        name = sample["sample"]
        checks = {
            "trimmed reads": f"data/intermediate/trimmed/{name}.done",
            "assembly": f"data/results/assemblies/{name}/assembly.fasta",
            "QUAST": f"data/results/quast/{name}/report.tsv",
            "CheckM2": f"data/results/checkm2/{name}/quality_report.tsv",
            "GTDB-Tk": f"data/results/gtdbtk/{name}/gtdbtk.bac120.summary.tsv",
            "Bakta": f"data/results/bakta/{name}/{name}.gff3",
            "AMRFinderPlus": f"data/results/amrfinder/{name}.tsv",
            "BAGEL4": f"data/results/bagel4/{name}/bagel4.done",
            "antiSMASH": f"data/results/antismash/{name}/index.html",
            "MOB-suite": f"data/results/mobsuite/{name}/mobtyper_results.txt",
        }
        rows.append((name, sample.get("layout", ""), checks))

    css = "body{font-family:Arial,sans-serif;margin:32px;line-height:1.4} table{border-collapse:collapse;width:100%} th,td{border:1px solid #ccc;padding:6px;text-align:left} .done{color:#067a2f}.missing{color:#9a3412}"
    parts = [f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title><style>{css}</style></head><body>"]
    parts.append(f"<h1>{html.escape(title)}</h1><p>Generated: {html.escape(now)}</p>")
    parts.append("<h2>Samples</h2><table><tr><th>Sample</th><th>Layout</th><th>Step</th><th>Status</th><th>Expected output</th></tr>")
    for name, layout, checks in rows:
        for step, path in checks.items():
            st = status(path)
            parts.append(f"<tr><td>{html.escape(name)}</td><td>{html.escape(layout)}</td><td>{html.escape(step)}</td><td class='{st}'>{st}</td><td>{html.escape(path)}</td></tr>")
    parts.append("</table>")
    parts.append("<h2>MultiQC</h2><p>data/results/multiqc/multiqc_report.html</p>")
    parts.append("</body></html>")
    Path(html_out).write_text("\n".join(parts))

    lines = [f"Generated: {now}", f"Samples: {len(samples)}"]
    for name, layout, checks in rows:
        done = sum(1 for p in checks.values() if exists(p))
        lines.append(f"{name} ({layout}): {done}/{len(checks)} steps complete")
    write_minimal_pdf(pdf_out, title, lines)


if __name__ == "__main__":
    main()

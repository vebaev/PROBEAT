#!/usr/bin/env python3
"""Merge CGC summary with substrate predictions, convert synteny PDFs to PNG."""
import csv, sys, os, subprocess, shutil
from pathlib import Path
from subprocess import run


def main():
    cgc_file = Path(sys.argv[1])       # cgc_standard_out_summary.tsv
    substrate_file = Path(sys.argv[2]) # substrate_prediction.tsv
    synteny_dir = Path(sys.argv[3])     # synteny_pdf/
    out_table = Path(sys.argv[4])       # merged TSV

    # Read substrate predictions
    subs = {}  # CGC# -> {PULID, substrate, bitscore}
    with open(substrate_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            cols = line.split("\t")
            if len(cols) >= 4:
                cgcid = cols[0]  # "4|CGC9"
                cgc_num = cgcid.split("|")[-1]  # "CGC9"
                subs[cgc_num] = {
                    "PULID": cols[1],
                    "substrate": cols[2],
                    "bitscore": cols[3],
                }

    # Read CGC summary and merge
    merged = []
    with open(cgc_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cgc = row.get("CGC#", "")
            if cgc in subs:
                sub = subs[cgc]
                merged.append({
                    "CGC#": cgc,
                    "Contig": row.get("Contig ID", ""),
                    "Start": row.get("Cluster Start", ""),
                    "End": row.get("Cluster End", ""),
                    "Genes": row.get("Genes", ""),
                    "CAZymes": row.get("CAZymes", ""),
                    "Length": row.get("Length (bp)", ""),
                    "PULID": sub["PULID"],
                    "Substrate": sub["substrate"],
                    "Bitscore": sub["bitscore"],
                })

    # Write merged TSV
    if merged:
        with open(out_table, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=merged[0].keys(), delimiter="\t")
            w.writeheader()
            w.writerows(merged)

    # Convert synteny PDFs to PNG (only those for merged CGCs)
    cgc_nums = {r["CGC#"] for r in merged}
    png_dir = synteny_dir.parent / "synteny_png"
    png_dir.mkdir(exist_ok=True)
    pdftoppm = shutil.which("pdftoppm")
    sibling_pdftoppm = Path(sys.executable).with_name("pdftoppm")
    if not pdftoppm and sibling_pdftoppm.exists():
        pdftoppm = str(sibling_pdftoppm)
    python = shutil.which("python") or shutil.which("python3")
    if not pdftoppm:
        print("pdftoppm unavailable; trying Python PDF conversion fallback.")
    for pdf in sorted(synteny_dir.glob("*.pdf")):
        # File naming: "10_CGC18_PUL0744-syntenic.pdf"
        parts = pdf.stem.split("_")
        for p in parts:
            if p.startswith("CGC") and p in cgc_nums:
                png = png_dir / (pdf.stem + ".png")
                if not png.exists():
                    if pdftoppm:
                        run([pdftoppm, "-png", "-r", "150", "-singlefile",
                             str(pdf), str(png.with_suffix(""))],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        # pdftoppm appends -1.png or .png depending on version
                        actual = png_dir / (pdf.stem + "-1.png")
                        if actual.exists():
                            actual.rename(png)
                    if not png.exists() and python:
                        script = (
                            "import sys\n"
                            "from pathlib import Path\n"
                            "pdf,out=Path(sys.argv[1]),Path(sys.argv[2])\n"
                            "try:\n"
                            " import fitz\n"
                            " doc=fitz.open(pdf); page=doc.load_page(0); pix=page.get_pixmap(matrix=fitz.Matrix(2,2), alpha=False); pix.save(out)\n"
                            "except Exception:\n"
                            " from pdf2image import convert_from_path\n"
                            " imgs=convert_from_path(str(pdf), dpi=150, first_page=1, last_page=1); imgs[0].save(out)\n"
                        )
                        run([python, "-c", script, str(pdf), str(png)],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                break

    print("Merged: %d CGCs -> %s" % (len(merged), out_table))
    print("PNGs: %d in %s" % (len(list(png_dir.glob("*.png"))), png_dir))


if __name__ == "__main__":
    main()

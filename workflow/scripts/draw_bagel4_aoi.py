#!/usr/bin/env python3
import csv
import html
import sys
from pathlib import Path


bakta_tsv, clusters_tsv, output_dir = map(Path, sys.argv[1:4])
assembly_fasta = Path(sys.argv[4]) if len(sys.argv) > 4 else None
bakta_fasta = Path(sys.argv[5]) if len(sys.argv) > 5 else None
output_dir.mkdir(parents=True, exist_ok=True)


def color(product):
    text = product.lower()
    if any(word in text for word in ("bacteriocin", "plantaricin", "lantibiotic", "riipp")):
        return "#57b957"
    if any(word in text for word in ("transport", "permease", "abc")):
        return "#ef5350"
    if any(word in text for word in ("regulator", "transcription", "response regulator")):
        return "#f4d35e"
    if any(word in text for word in ("protease", "peptidase")):
        return "#c85adf"
    if any(word in text for word in ("immunity", "resistance")):
        return "#f28cb1"
    return "#93b8b8" if product and "hypothetical" not in text else "#d8dee5"


def parse_int(value):
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def overlaps(left_start, left_end, right_start, right_end):
    return left_start is not None and left_end is not None and right_start is not None and right_end is not None and left_end >= right_start and left_start <= right_end


def fasta_lengths(path):
    if not path or not path.exists():
        return {}
    lengths = {}
    name = None
    size = 0
    with path.open(errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = size
                name = line[1:].split()[0]
                size = 0
            else:
                size += len(line)
    if name is not None:
        lengths[name] = size
    return lengths


def contig_aliases(assembly_lengths, bakta_lengths):
    aliases = {}
    bakta_by_length = {}
    for name, size in bakta_lengths.items():
        bakta_by_length.setdefault(size, []).append(name)
    for name, size in assembly_lengths.items():
        matches = bakta_by_length.get(size, [])
        if len(matches) == 1 and matches[0] != name:
            aliases[name] = matches[0]
    return aliases


assembly_to_bakta = contig_aliases(fasta_lengths(assembly_fasta), fasta_lengths(bakta_fasta))


features = []
with bakta_tsv.open(errors="replace") as handle:
    for row in csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"):
        if len(row) < 8:
            continue
        contig_name = row[0]
        if contig_name.startswith("contig_"):
            contig_name = contig_name[len("contig_"):]
        features.append({
            "contig": contig_name, "type": row[1],
            "start": int(row[2]), "end": int(row[3]), "strand": row[4],
            "locus": row[5], "gene": row[6], "product": row[7],
        })
by_contig = {}
for feature in features:
    by_contig.setdefault(feature["contig"], []).append(feature)
    by_contig.setdefault(f"contig_{feature['contig']}", []).append(feature)

with clusters_tsv.open(errors="replace") as handle:
    clusters = list(csv.DictReader(handle, delimiter="\t"))

for cluster in clusters:
    aoi = cluster["AOI"]
    contig = aoi.split(".AOI_")[0] if ".AOI_" in aoi else aoi.split(".")[0]
    start, end = int(cluster["start"]), int(cluster["end"])
    hit_start = parse_int(cluster.get("hit_start"))
    hit_end = parse_int(cluster.get("hit_end"))
    if hit_start is not None and hit_end is not None and hit_start > hit_end:
        hit_start, hit_end = hit_end, hit_start
    stripped_contig = contig[len("contig_"):] if contig.startswith("contig_") else contig
    mapped_contig = assembly_to_bakta.get(contig, contig)
    stripped_mapped = mapped_contig[len("contig_"):] if mapped_contig.startswith("contig_") else mapped_contig
    candidate_contigs = [contig, stripped_contig, mapped_contig, stripped_mapped]
    genes = []
    seen = set()
    for candidate in candidate_contigs:
        if candidate in seen:
            continue
        seen.add(candidate)
        contig_features = by_contig.get(candidate) or []
        genes = [
            row for row in contig_features
            if row["type"].lower() == "cds" and row["end"] >= start and row["start"] <= end
        ]
        if genes:
            break
    if not genes:
        genes = [
            {"contig": contig, "type": "cds", "start": start, "end": end, "strand": "+",
             "locus": cluster.get("AOI", ""), "gene": cluster.get("fullname", "") or cluster.get("class", "") or "AOI",
             "product": cluster.get("fullname", "") or "BAGEL4 AOI"}
        ]
    width, height, margin, baseline = 1400, 430, 45, 260
    span = max(1, end - start)
    scale = (width - 2 * margin) / span
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{margin}' y='35' font-family='Arial' font-size='24' font-weight='bold'>{html.escape(cluster['AOI'])}</text>",
        f"<text x='{margin}' y='62' font-family='Arial' font-size='14'>Contig {html.escape(contig)}: {start:,}-{end:,}; {html.escape(cluster.get('fullname', ''))}</text>",
        f"<line x1='{margin}' y1='{baseline}' x2='{width-margin}' y2='{baseline}' stroke='#263238' stroke-width='2'/>",
    ]
    for index, gene in enumerate(genes):
        x1 = margin + (max(start, gene["start"]) - start) * scale
        x2 = margin + (min(end, gene["end"]) - start) * scale
        arrow = min(12, max(4, (x2 - x1) / 3))
        if gene["strand"] == "-":
            points = f"{x2},{baseline-18} {x1+arrow},{baseline-18} {x1},{baseline} {x1+arrow},{baseline+18} {x2},{baseline+18}"
        else:
            points = f"{x1},{baseline-18} {x2-arrow},{baseline-18} {x2},{baseline} {x2-arrow},{baseline+18} {x1},{baseline+18}"
        label = gene["gene"] or gene["locus"] or f"orf{index+1:04d}"
        product = gene["product"] or "No function determined"
        fill = "#57b957" if overlaps(gene["start"], gene["end"], hit_start, hit_end) else color(product)
        parts.append(f"<polygon points='{points}' fill='{fill}' stroke='#37474f' stroke-width='1'/>")
        tx = (x1 + x2) / 2
        parts.append(f"<text x='{tx:.1f}' y='{baseline-27}' transform='rotate(-55 {tx:.1f},{baseline-27})' font-family='Arial' font-size='10'>{html.escape(label)}</text>")
    if len(genes) == 1 and genes[0]["start"] == start and genes[0]["end"] == end:
        parts.append(f"<text x='{margin}' y='{baseline+58}' font-family='Arial' font-size='13' fill='#8a4b00'>Bakta CDS features were not mapped inside this AOI; showing BAGEL4 AOI interval.</text>")
    legend = [("#d8dee5", "No function determined"), ("#93b8b8", "Annotated gene"), ("#57b957", "Bacteriocin / RiPP"), ("#ef5350", "Transport"), ("#f4d35e", "Regulation"), ("#f28cb1", "Immunity"), ("#c85adf", "Protease")]
    for index, (fill, label) in enumerate(legend):
        x, y = width - 260, 25 + index * 22
        parts.append(f"<rect x='{x}' y='{y}' width='14' height='14' fill='{fill}' stroke='#37474f'/>")
        parts.append(f"<text x='{x+22}' y='{y+12}' font-family='Arial' font-size='12'>{label}</text>")
    parts.append("</svg>")
    safe_name = cluster["AOI"].replace("/", "_")
    (output_dir / f"{safe_name}.svg").write_text("\n".join(parts))

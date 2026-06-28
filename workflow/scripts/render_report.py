#!/usr/bin/env python3
import csv
import datetime as dt
import html
import json
import os
import re
import struct
import subprocess
import sys
import zlib
from collections import Counter
from pathlib import Path

import yaml


WORKFLOW_VERSION = os.environ.get("PROBEAT_WORKFLOW_VERSION", "").strip()


def workflow_version():
    global WORKFLOW_VERSION
    if WORKFLOW_VERSION:
        return WORKFLOW_VERSION
    try:
        result = subprocess.run(
            ["git", "describe", "--tags", "--abbrev=0"],
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode == 0 and result.stdout.strip():
            WORKFLOW_VERSION = result.stdout.strip()
            return WORKFLOW_VERSION
    except OSError:
        pass
    WORKFLOW_VERSION = "1.0"
    return WORKFLOW_VERSION


def read_samples(path):
    if not Path(path).exists():
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def exists(path):
    return Path(path).exists()


def read_json_file(path):
    path = Path(path)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except (OSError, json.JSONDecodeError):
        return {}


def gutsmash_status_record(sample):
    return read_json_file(f"data/results/{sample}/gutsmash/.probeat_gutsmash_status.json")


def status(path):
    path = Path(path)
    if path.name == ".probeat_gutsmash_status.json":
        record = read_json_file(path)
        if record.get("status") == "completed":
            return "done"
        if record:
            return "soft-failed"
        return "missing"
    return "done" if path.exists() else "missing"


def fmt_int(value):
    if value in (None, ""):
        return "NA"
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value)


def fmt_float(value, digits=2):
    if value in (None, ""):
        return "NA"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def fmt_percent(value):
    if value in (None, ""):
        return "NA"
    try:
        val = float(value)
        if val <= 1:
            val *= 100
        return f"{val:.2f}%"
    except (TypeError, ValueError):
        return str(value)


def read_fastp_summary(sample):
    path = Path(f"data/results/{sample}/qc/fastp/{sample}.json")
    if not path.exists():
        path = Path(f"data/results/{sample}/qc/porechop/{sample}.json")
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    summary = data.get("summary", {})
    filtering = data.get("filtering_result", {})
    duplication = data.get("duplication", {})
    insert_size = data.get("insert_size", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    return {
        "raw_reads": before.get("total_reads"),
        "raw_bases": before.get("total_bases"),
        "raw_q30": before.get("q30_rate"),
        "raw_gc": before.get("gc_content"),
        "filtered_reads": after.get("total_reads"),
        "filtered_bases": after.get("total_bases"),
        "filtered_q30": after.get("q30_rate"),
        "filtered_gc": after.get("gc_content"),
        "passed_filter_reads": filtering.get("passed_filter_reads"),
        "low_quality_reads": filtering.get("low_quality_reads"),
        "too_many_n_reads": filtering.get("too_many_N_reads"),
        "adapter_trimmed_reads": data.get("adapter_cutting", {}).get("adapter_trimmed_reads"),
        "adapter_trimmed_bases": data.get("adapter_cutting", {}).get("adapter_trimmed_bases"),
        "duplication_rate": duplication.get("rate"),
        "insert_size_peak": insert_size.get("peak"),
    }


def read_quast(sample):
    path = Path(f"data/results/{sample}/quast/report.tsv")
    if not path.exists():
        return {}
    values = {}
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                values[row[0]] = row[1]
    return {
        "contigs": values.get("# contigs"),
        "largest_contig": values.get("Largest contig"),
        "total_length": values.get("Total length"),
        "gc": values.get("GC (%)"),
        "n50": values.get("N50"),
        "n90": values.get("N90"),
        "l50": values.get("L50"),
        "l90": values.get("L90"),
        "ns_per_100kbp": values.get("# N's per 100 kbp"),
    }


def count_bakta_features(sample):
    path = Path(f"data/results/{sample}/bakta/{sample}.gff3")
    counts = {"CDS": 0, "tRNA": 0, "rRNA": 0, "tmRNA": 0, "ncRNA": 0, "CRISPR": 0}
    if not path.exists():
        return counts
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            feature = fields[2]
            if feature in counts:
                counts[feature] += 1
            elif feature == "crispr-repeat":
                counts["CRISPR"] += 0
    return counts


def count_amr_hits(sample):
    path = Path(f"data/results/{sample}/amrfinder/{sample}.tsv")
    if not path.exists():
        return None
    with path.open() as handle:
        rows = [line for line in handle if line.strip()]
    return max(len(rows) - 1, 0)


def count_tsv_hits(path):
    path = Path(path)
    if not path.exists():
        return None
    with path.open(errors="replace") as handle:
        rows = [
            line for line in handle
            if line.strip() and not line.startswith("#") and not line.lower().startswith("file\t")
        ]
    return len(rows)


def read_bagel4(sample):
    done = Path(f"data/results/{sample}/bagel4/bagel4.done")
    if not done.exists():
        return {"status": "not run"}
    result_files = list(Path(f"data/results/{sample}/bagel4").glob("*.tsv"))
    if not result_files:
        return {"status": "not configured", "summary": "not configured"}
    hits = sum(count_tsv_hits(path) or 0 for path in result_files)
    return {"status": "done", "hits": hits, "summary": f"{fmt_int(hits)} candidate hits"}


def bagel4_rows(sample):
    return read_tsv_rows(f"data/results/{sample}/bagel4/bacteriocin_clusters.tsv")


def bagel4_images(sample, suffixes=(".png", ".svg")):
    bagel_dir = Path(f"data/results/{sample}/bagel4")
    if not bagel_dir.exists():
        return []
    return sorted(path for path in bagel_dir.rglob("*") if path.suffix.lower() in suffixes)


def bagel4_preferred_images(sample):
    """Return one BAGEL4 graphic per AOI, preferring PNG over SVG for reports."""
    images = {}
    for path in bagel4_images(sample):
        key = path.with_suffix("")
        if key not in images or path.suffix.lower() == ".png":
            images[key] = path
    return [images[key] for key in sorted(images)]


def count_dbcan_hits(sample):
    path = Path(f"data/results/{sample}/dbcan/overview.tsv")
    if not path.exists():
        path = Path(f"data/results/{sample}/dbcan/overview.txt")
    if not path.exists():
        return None
    with path.open(newline="", errors="replace") as handle:
        return sum(
            1 for row in csv.DictReader(handle, delimiter="\t")
            if re.search(r"\b(?:GH|GT|PL|CE|AA|CBM)\d+", row.get("Recommend Results", "") or "")
        )

def read_tsv_rows(path):
    path = Path(path)
    if not path.exists():
        return []
    with path.open(newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_gutsmash_index_html(sample):
    """Parse gutSMASH index.html to extract region table rows and clusterblast SVG paths."""
    html_path = Path(f"data/results/{sample}/gutsmash/index.html")
    if not html_path.exists():
        return [], []
    text = html_path.read_text(errors="replace")
    rows = []
    # Extract rows from the single-record-tables tbody
    tbody_match = re.search(r'<div id="single-record-tables">.*?</div>\s*</div>\s*<div id="compact-record-table"', text, re.DOTALL)
    if not tbody_match:
        tbody_match = re.search(r'<div id="single-record-tables">.*?</div>\s*</div>', text, re.DOTALL)
    if tbody_match:
        tbody = tbody_match.group(0)
        for tr in re.findall(r'<tr class="linked-row[^"]*"[^>]*>(.*?)</tr>', tbody, re.DOTALL):
            cells = re.findall(r'<td[^>]*>(.*?)</td>', tr, re.DOTALL)
            if len(cells) >= 8:
                region = re.sub(r'<[^>]+>', '', cells[0]).replace('&nbsp;', ' ').replace('&nbsp', ' ').strip()
                type_ = re.sub(r'<[^>]+>', '', cells[1]).strip()
                class_ = re.sub(r'<[^>]+>', '', cells[2]).strip()
                from_ = re.sub(r'<[^>]+>', '', cells[3]).strip()
                to_ = re.sub(r'<[^>]+>', '', cells[4]).strip()
                cluster = re.sub(r'<[^>]+>', '', cells[5]).strip()
                cluster_type = re.sub(r'<[^>]+>', '', cells[6]).strip()
                similarity = re.sub(r'<[^>]+>', '', cells[7]).strip()
                rows.append({
                    "region": region,
                    "type": type_,
                    "class": class_,
                    "from": from_,
                    "to": to_,
                    "cluster": cluster,
                    "cluster_type": cluster_type,
                    "similarity": similarity,
                })
    # Find all knownclusterblast *_all.svg files
    svg_dir = Path(f"data/results/{sample}/gutsmash/svg")
    svg_files = sorted(svg_dir.glob("knownclusterblast_*_all.svg")) if svg_dir.exists() else []
    return rows, svg_files


def grouped_gut_brain_rows(sample):
    """Combine detected genes by Metabolite/Trait and deduplicate shared text."""
    grouped = {}
    for row in read_tsv_rows(f"data/results/{sample}/gut_brain_markers.tsv"):
        if row.get("bakta_locus", "-") == "-":
            continue
        traits = [
            " ".join(value.split())
            for value in re.split(r"\s*\|\s*", row.get("metabolite_or_trait", ""))
            if value.strip()
        ]
        for trait in traits:
            group = grouped.setdefault(trait, {
                "metabolite_or_trait": trait,
                "marker": [],
                "bakta_locus": [],
                "ko_or_ec_hints": [],
                "marker_logic": [],
                "interpretation": [],
            })
            fields = {
                "marker": (row.get("marker", ""), r"\s*;\s*"),
                "bakta_locus": (row.get("bakta_locus", ""), r"\s*,\s*"),
                "ko_or_ec_hints": (row.get("ko_or_ec_hints", ""), r"\s*;\s*"),
                "marker_logic": (row.get("marker_logic", ""), r"\s*\|\s*"),
                "interpretation": (row.get("interpretation", ""), r"\s*\|\s*"),
            }
            for field, (value, separator) in fields.items():
                for item in re.split(separator, value or ""):
                    item = " ".join(item.split())
                    if item and item != "-" and item not in group[field]:
                        group[field].append(item)

    rows = []
    for group in grouped.values():
        rows.append({
            "metabolite_or_trait": group["metabolite_or_trait"],
            "marker": "\n".join(group["marker"]),
            "bakta_locus": "\n".join(group["bakta_locus"]),
            "ko_or_ec_hints": "\n".join(group["ko_or_ec_hints"]),
            "marker_logic": "\n".join(group["marker_logic"]),
            "interpretation": "\n".join(group["interpretation"]),
        })
    return rows


def dbcan_class_counts(sample):
    counts = {key: 0 for key in ("GH", "GT", "PL", "CE", "AA", "CBM")}
    path = Path(f"data/results/{sample}/dbcan/overview.tsv")
    if not path.exists():
        path = Path(f"data/results/{sample}/dbcan/overview.txt")
    if not path.exists():
        return counts
    for row in read_tsv_rows(path):
        families = set(re.findall(
            r"\b(GH|GT|PL|CE|AA|CBM)\d+",
            row.get("Recommend Results", "") or "",
        ))
        for family in families:
            counts[family] += 1
    return counts


def dbcan_recommended_gene_count(sample):
    path = Path(f"data/results/{sample}/dbcan/overview.tsv")
    if not path.exists():
        path = Path(f"data/results/{sample}/dbcan/overview.txt")
    if not path.exists():
        return 0
    return sum(
        1 for row in read_tsv_rows(path)
        if re.search(r"\b(?:GH|GT|PL|CE|AA|CBM)\d+", row.get("Recommend Results", "") or "")
    )


def dbcan_detected_families(sample):
    families = {}
    overview_path = Path(f"data/results/{sample}/dbcan/overview.tsv")
    if not overview_path.exists():
        overview_path = Path(f"data/results/{sample}/dbcan/overview.txt")
    if not overview_path.exists():
        return families
    for row in read_tsv_rows(overview_path):
        gene_id = row.get("Gene ID", "")
        for family in re.findall(
            r"\b(?:GH|GT|PL|CE|AA|CBM)\d+",
            row.get("Recommend Results", "") or "",
        ):
            families.setdefault(family, set()).add(gene_id)
    return families


def dbcan_to_bakta_loci(sample):
    bakta = {}
    gff_path = Path(f"data/results/{sample}/bakta/{sample}.gff3")
    if gff_path.exists():
        with gff_path.open(errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "CDS":
                    continue
                attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
                bakta.setdefault((fields[0], fields[6]), []).append(
                    (int(fields[3]), int(fields[4]), attrs.get("locus_tag", attrs.get("ID", "")))
                )
    mapping = {}
    prodigal_path = Path(f"data/results/{sample}/dbcan/prodigal.gff")
    if prodigal_path.exists():
        with prodigal_path.open(errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "CDS":
                    continue
                attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
                gene_id = attrs.get("ID", "")
                start, stop = int(fields[3]), int(fields[4])
                overlaps = [
                    (max(0, min(stop, b_stop) - max(start, b_start) + 1), locus)
                    for b_start, b_stop, locus in bakta.get((f"contig_{fields[0]}", fields[6]), [])
                ]
                if overlaps:
                    overlap, locus = max(overlaps)
                    if overlap > 0:
                        mapping[gene_id] = locus
    return mapping


def probiotic_cazy_matches(sample):
    path = Path("data/databases/probiotic_cazy.txt")
    if not path.exists():
        return []
    detected = dbcan_detected_families(sample)
    locus_map = dbcan_to_bakta_loci(sample)
    matches = []
    for row in read_tsv_rows(path):
        requested = set(re.findall(r"\b(?:GH|GT|PL|CE|AA|CBM)\d+", row.get("Main dbCAN / CAZy families", "")))
        found = sorted(requested & set(detected), key=lambda value: (value[:3], int(re.search(r"\d+", value).group())))
        if found:
            genes = sorted({gene for family in found for gene in detected[family] if gene})
            loci = sorted({locus_map.get(g, g) for g in genes} if locus_map else genes)
            matches.append((
                row.get("Prebiotic substrate / group", ""),
                ", ".join(found),
                ", ".join(loci),
                row.get("Typical enzymes / functions", ""),
            ))
    return matches


PDF_PALETTES = {
    "blue": ((0.08, 0.30, 0.52), (0.88, 0.94, 0.98), (0.96, 0.98, 1.00)),
    "green": ((0.05, 0.42, 0.25), (0.88, 0.96, 0.91), (0.96, 0.99, 0.97)),
    "red": ((0.62, 0.13, 0.16), (0.98, 0.90, 0.90), (1.00, 0.97, 0.97)),
    "orange": ((0.68, 0.34, 0.06), (0.99, 0.93, 0.84), (1.00, 0.98, 0.94)),
    "purple": ((0.36, 0.22, 0.56), (0.93, 0.90, 0.97), (0.98, 0.97, 1.00)),
}


def pdf_palette(title):
    name = str(title).lower()
    if any(key in name for key in ("amr", "vfdb", "virulence", "resistance", "abricate", "card", "megares")):
        return PDF_PALETTES["red"]
    if any(key in name for key in ("mobrecon", "mob-recon", "genomad", "plasmid", "virus", "provirus")):
        return PDF_PALETTES["orange"]
    if any(key in name for key in ("pgmdb", "probiotic", "gut-brain", "bagel", "dbcan", "cazy")):
        return PDF_PALETTES["green"]
    if any(key in name for key in ("kegg", "kofam", "eggnog", "cog", "gutsmash", "antismash", "genomespot", "go ")):
        return PDF_PALETTES["purple"]
    return PDF_PALETTES["blue"]


def draw_pdf_table(pdf, title, subtitle, headers, rows, widths, row_height=18, new_page=True, start_y=800):
    def ascii_pdf_text(value):
        replacements = {
            "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta",
            "±": "+/-", "–": "-", "—": "-", "→": "->", "↔": "<->",
        }
        text = str(value)
        for source, target in replacements.items():
            text = text.replace(source, target)
        return text

    def wrap_cell(value, width):
        text = ascii_pdf_text(value)
        max_chars = max(8, int((width - 8) / 3.7))
        lines = []
        for paragraph in text.splitlines() or [""]:
            current = ""
            for token in re.split(r"(\s+|,\s*)", paragraph):
                if not token:
                    continue
                candidate = current + token
                if len(candidate) > max_chars and current.strip():
                    lines.append(current.strip())
                    current = token.lstrip()
                else:
                    current = candidate
            if current.strip():
                lines.append(current.strip())
        return lines or [""]

    if new_page:
        pdf.add_page()
    accent, header_fill, alternate_fill = pdf_palette(title)
    pdf.rect(40, start_y - 3, 5, 18 if new_page else 14, fill=accent)
    pdf.text(52, start_y, title, size=15 if new_page else 11, color=accent, bold=True)
    pdf.text(40, start_y - 20, subtitle, size=9, color=(0.38, 0.45, 0.55))
    y = start_y - 50
    x0 = 40
    for header, width in zip(headers, widths):
        pdf.rect(x0, y, width, 18, fill=header_fill, stroke=accent, line_width=0.4)
        pdf.text(x0 + 4, y + 5, header, size=7, color=accent, bold=True)
        x0 += width
    y -= row_height
    for row_index, row in enumerate(rows):
        wrapped = [wrap_cell(value, width) for value, width in zip(row, widths)]
        current_height = max(row_height, 10 + max(len(lines) for lines in wrapped) * 9)
        if y - current_height < 45:
            pdf.add_page()
            y = 795
            x0 = 40
            for header, width in zip(headers, widths):
                pdf.rect(x0, y, width, 18, fill=header_fill, stroke=accent, line_width=0.4)
                pdf.text(x0 + 4, y + 5, header, size=7, color=accent, bold=True)
                x0 += width
            y -= 22
        x0 = 40
        for lines, width in zip(wrapped, widths):
            cell_y = y - current_height + row_height
            fill = alternate_fill if row_index % 2 else (1, 1, 1)
            pdf.rect(x0, cell_y, width, current_height, fill=fill, stroke=(0.82, 0.86, 0.9), line_width=0.3)
            text_y = y + row_height - 10
            for line in lines:
                numeric_cell = bool(re.fullmatch(r"[-+]?[\d,.]+(?:%|x| bp| C| pH)?", line.strip()))
                text_x = x0 + width - 5 - len(line) * 3.25 if numeric_cell else x0 + 4
                pdf.text(max(x0 + 4, text_x), text_y, line, size=6.5, color=(0.16, 0.21, 0.28))
                text_y -= 9
            x0 += width
        y -= current_height
    return y

def read_rasusa(sample, layout):
    path = Path(f"data/results/{sample}/qc/rasusa/{sample}_{layout}.tsv")
    if not path.exists():
        return {"status": "not run"}
    with path.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"), {})
    return {
        "status": row.get("action", "done"),
        "input_reads": row.get("input_reads"),
        "output_reads": row.get("output_reads"),
        "target_reads": row.get("target_reads"),
    }


def count_keggdecoder_functions(sample):
    path = Path(f"data/results/{sample}/keggdecoder/pathway_completeness.tsv")
    if not path.exists():
        return None
    with path.open(errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")[1:]
        values = handle.readline().rstrip("\n").split("\t")[1:]
    numeric_values = []
    for value in values:
        if not value.strip():
            continue
        try:
            numeric_values.append(float(value))
        except ValueError:
            continue
    return sum(1 for value in numeric_values if value > 0), len(header)

def read_keggdecoder_functions(sample):
    path = Path(f"data/results/{sample}/keggdecoder/pathway_completeness.tsv")
    if not path.exists():
        return []
    with path.open(errors="replace") as handle:
        names = handle.readline().rstrip("\n").split("\t")[1:]
        values = handle.readline().rstrip("\n").split("\t")[1:]
    rows = []
    for name, value in zip(names, values):
        if not value.strip():
            continue
        try:
            numeric_value = float(value)
        except ValueError:
            continue
        if numeric_value > 0:
            rows.append((name, numeric_value))
    return rows


def keggdecoder_category(name):
    text = name.lower()
    groups = [
        ("Vitamins & cofactors", ("thiamin", "riboflavin", "cobalamin", "vitamin")),
        ("Amino acid biosynthesis", ("histidine", "arginine", "lysine", "serine", "threonine", "asparagine", "glutamine", "cysteine", "glycine", "proline", "alanine", "valine", "methionine", "phenylalanine", "isoleucine", "leucine", "tryptophan", "tyrosine", "aspartate", "glutamate")),
        ("Transport, stress & interaction", ("transporter", "biofilm", "competence", "adhesion", "secretion", "arsenic", "flagellum", "chemotaxis")),
        ("Carbohydrates & fermentation", ("glycolysis", "gluconeogenesis", "entner", "mixed acid", "starch", "glycogen", "pectin", "chitin", "cellul", "galacturon", "amylase", "glucosidase", "pullulanase")),
        ("Energy & central metabolism", ("atpase", "cytochrome", "tca", "anaplerotic", "hydrogenase", "oxidoreductase", "photosystem", "rubisco", "cbb cycle", "rtca", "wood-ljungdahl", "hydroxypropionate", "hydroxybutyrate")),
        ("Storage & other metabolism", ("polyhydroxybutyrate", "polyphosphate", "mevalonate", "mep-doxp", "carotenoid", "retinal", "serine pathway")),
    ]
    for category, terms in groups:
        if any(term in text for term in terms):
            return category
    return "Other detected functions"


def keggdecoder_heatmap_color(value):
    value = max(0.0, min(1.0, value))
    return (1 - 0.72 * value, 0.94 - 0.28 * value, 0.98 - 0.18 * value)

def probiotic_relevant_keggdecoder(functions):
    terms = (
        "lactate", "acetate", "glycolysis", "starch", "glycogen",
        "thiamin", "riboflavin", "cobalamin", "vitamin",
        "biofilm", "adhesion", "competence", "transporter",
        "arsenic", "copper", "fe-mn", "atpase", "cytochrome bd",
        "arginine", "cysteine", "glutamate", "glutamine", "lysine",
        "methionine", "tryptophan", "polyphosphate",
    )
    return [(name, value) for name, value in functions if any(term in name.lower() for term in terms)]


def read_pgmdb(sample):
    summary_path = Path(f"data/results/{sample}/pgmdb/summary.tsv")
    hits_path = Path(f"data/results/{sample}/pgmdb/hits.tsv")
    if not summary_path.exists():
        return {"status": "not run"}
    categories = []
    total_hits = unique_markers = 0
    with summary_path.open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["metric"] == "total_hits":
                total_hits = int(row["value"])
                unique_markers = int(row["unique_markers"])
            else:
                categories.append({
                    "category": row["metric"],
                    "hits": int(row["value"]),
                    "unique_markers": int(row["unique_markers"]),
                })
    return {
        "status": "done",
        "total_hits": total_hits,
        "unique_markers": unique_markers,
        "categories": categories,
        "summary": "; ".join(
            f"{item['category']} {item['unique_markers']}" for item in categories[:4]
        ) or "no curated marker hits",
        "hits_path": str(hits_path),
    }


def split_pgmdb_functions(value):
    pieces = [item.strip() for item in re.split(r"\s*/\s*|\s*;\s*", value or "") if item.strip()]
    seen = []
    for piece in pieces:
        normalized = normalize_pgmdb_function(piece)
        if normalized not in seen:
            seen.append(normalized)
    return seen


def normalize_pgmdb_function(value):
    cleaned = " ".join(str(value or "").split())
    if cleaned and cleaned == cleaned.lower():
        cleaned = cleaned.capitalize()
    return cleaned


def pgmdb_sankey_edges(sample):
    edges = []
    seen = set()
    for row in read_tsv_rows(f"data/results/{sample}/pgmdb/hits.tsv"):
        gene = (row.get("gene") or row.get("marker") or "").strip()
        locus = (row.get("locus_tag") or "").strip()
        if not gene or not locus:
            continue
        left = f"{gene} ({locus})"
        for function in split_pgmdb_functions(row.get("function", "Other")):
            key = (left, function)
            if key not in seen:
                edges.append({"gene_label": left, "function": function})
                seen.add(key)
    return edges


def pgmdb_sankey_data(sample):
    edges = pgmdb_sankey_edges(sample)
    left_labels = sorted({edge["gene_label"] for edge in edges}, key=str.lower)
    right_labels = sorted({edge["function"] for edge in edges}, key=str.lower)
    return {"edges": edges, "left_labels": left_labels, "right_labels": right_labels}


def pgmdb_sankey_plotly(sample, div_id):
    data = pgmdb_sankey_data(sample)
    edges = data["edges"]
    if not edges:
        return ""
    left_labels = data["left_labels"]
    right_labels = data["right_labels"]
    labels = left_labels + right_labels
    left_index = {label: idx for idx, label in enumerate(left_labels)}
    right_index = {label: len(left_labels) + idx for idx, label in enumerate(right_labels)}
    function_colors = [
        "#7c9cf5", "#62b4ff", "#6cd39c", "#f7c948",
        "#f48fb1", "#f97066", "#9d8df1", "#4fd1c5",
    ]
    function_color_map = {
        label: function_colors[index % len(function_colors)]
        for index, label in enumerate(right_labels)
    }
    node_colors = ["#edf2f7"] * len(left_labels) + [function_color_map[label] for label in right_labels]
    sources = [left_index[edge["gene_label"]] for edge in edges]
    targets = [right_index[edge["function"]] for edge in edges]
    values = [1] * len(edges)
    link_colors = [function_color_map[edge["function"]] for edge in edges]
    payload = {
        "labels": labels,
        "sources": sources,
        "targets": targets,
        "values": values,
        "node_colors": node_colors,
        "link_colors": link_colors,
        "height": max(520, 120 + max(len(left_labels), len(right_labels)) * 24),
    }
    payload_json = json.dumps(payload, ensure_ascii=False)
    return (
        f"<div id='{html.escape(div_id)}' class='plotly-sankey'></div>"
        "<script>"
        f"(function(){{"
        f"var data={payload_json};"
        f"var trace={{type:'sankey',orientation:'h',arrangement:'snap',"
        f"node:{{pad:16,thickness:18,line:{{color:'#94a3b8',width:1}},label:data.labels,color:data.node_colors,"
        f"hovertemplate:'%{{label}}<extra></extra>'}},"
        f"link:{{source:data.sources,target:data.targets,value:data.values,color:data.link_colors,"
        f"hovertemplate:'%{{source.label}} -> %{{target.label}}<extra></extra>'}}}};"
        f"var layout={{height:data.height,margin:{{l:18,r:18,t:8,b:8}},paper_bgcolor:'#ffffff',plot_bgcolor:'#ffffff',"
        f"font:{{family:'Arial, sans-serif',size:13,color:'#0f172a'}}}};"
        f"Plotly.newPlot('{div_id}',[trace],layout,{{displayModeBar:false,responsive:true}});"
        f"}})();"
        "</script>"
    )


def pgmdb_sankey_svg(sample):
    data = pgmdb_sankey_data(sample)
    edges = data["edges"]
    left_labels = data["left_labels"]
    right_labels = data["right_labels"]
    if not edges:
        return None

    row_h = 28
    top = 54
    bottom = 38
    left_x = 24
    right_x = 505
    box_w = 170
    box_h = 18
    width = 700
    height = max(240, top + bottom + max(len(left_labels), len(right_labels)) * row_h)
    colors = ["#c7d2fe", "#bfdbfe", "#bbf7d0", "#fde68a", "#fbcfe8", "#fecaca", "#ddd6fe", "#a7f3d0"]
    function_colors = {label: colors[index % len(colors)] for index, label in enumerate(right_labels)}
    left_y = {label: top + index * row_h for index, label in enumerate(left_labels)}
    right_y = {label: top + index * row_h for index, label in enumerate(right_labels)}

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}' role='img' aria-label='PGMDB Sankey diagram for {html.escape(sample)}'>",
        "<style>",
        ".title{font:700 16px Arial,sans-serif;fill:#0f2742}",
        ".subtitle{font:12px Arial,sans-serif;fill:#64748b}",
        ".node-label{font:12px Arial,sans-serif;fill:#0f172a}",
        ".node-box{stroke:#94a3b8;stroke-width:1}",
        ".link{fill:none;stroke-linecap:round;opacity:.56}",
        "</style>",
        f"<rect x='0' y='0' width='{width}' height='{height}' fill='#ffffff'/>",
        f"<text class='title' x='24' y='24'>PGMDB probiotic genes to functional categories: {html.escape(sample)}</text>",
        "<text class='subtitle' x='24' y='42'>gene (Bakta locus tag) -&gt; probiotic function</text>",
    ]

    for edge in edges:
        y1 = left_y[edge["gene_label"]] + box_h / 2
        y2 = right_y[edge["function"]] + box_h / 2
        x1 = left_x + box_w
        x2 = right_x
        c1 = x1 + 90
        c2 = x2 - 90
        color = function_colors[edge["function"]]
        parts.append(
            f"<path class='link' d='M{x1},{y1} C{c1},{y1} {c2},{y2} {x2},{y2}' stroke='{color}' stroke-width='6'/>"
        )

    for label in left_labels:
        y = left_y[label]
        parts.append(f"<rect class='node-box' x='{left_x}' y='{y}' width='{box_w}' height='{box_h}' rx='4' fill='#f8fafc'/>")
        parts.append(f"<text class='node-label' x='{left_x + 8}' y='{y + 13}'>{html.escape(label)}</text>")
    for label in right_labels:
        y = right_y[label]
        parts.append(f"<rect class='node-box' x='{right_x}' y='{y}' width='{box_w}' height='{box_h}' rx='4' fill='{function_colors[label]}'/>")
        parts.append(f"<text class='node-label' x='{right_x + 8}' y='{y + 13}'>{html.escape(label)}</text>")

    parts.append("</svg>")
    return "\n".join(parts)


def write_pgmdb_sankey_assets(sample, output_dir):
    svg = pgmdb_sankey_svg(sample)
    if not svg:
        return None
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    svg_path = output_dir / f"{sample}_pgmdb_sankey.svg"
    png_path = output_dir / f"{sample}_pgmdb_sankey.png"
    svg_path.write_text(svg)
    converter = Path("mamba_envs/snakemake-conda-helper-2026-06-04/bin/rsvg-convert")
    if converter.exists():
        subprocess.run([str(converter), "-w", "1400", "-o", str(png_path), str(svg_path)], check=False)
    return {"svg": svg_path, "png": png_path if png_path.exists() else None}


def read_busco(sample):
    path = Path(f"data/results/{sample}/busco/busco_summary.txt")
    if not path.exists():
        return {"status": "not run"}
    text = path.read_text(errors="replace")
    match = re.search(
        r"C:(?P<complete>[\d.]+)%\[S:(?P<single>[\d.]+)%,D:(?P<duplicated>[\d.]+)%\],F:(?P<fragmented>[\d.]+)%,M:(?P<missing>[\d.]+)%,n:(?P<n>\d+)",
        text,
    )
    if not match:
        return {"status": "done"}
    values = match.groupdict()
    values["status"] = "done"
    values["summary"] = (
        f"Complete {values['complete']}%; duplicated {values['duplicated']}%; "
        f"fragmented {values['fragmented']}%; missing {values['missing']}%"
    )
    return values


def parse_mash_species(reference):
    label = Path(reference).name.rsplit(".", 1)[0]
    organism = label.split("-")[-1].replace("_", " ").strip()
    words = organism.split()
    if len(words) >= 2:
        return " ".join(words[:2])
    return organism or reference


def mash_hash_percent(hashes):
    if not hashes or "/" not in hashes:
        return None
    numerator, denominator = hashes.split("/", 1)
    try:
        return (int(numerator) / int(denominator)) * 100
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def read_mash(sample):
    path = Path(f"data/results/{sample}/mash/top_hits.tsv")
    if not path.exists():
        return {"status": "not run"}
    top_species = []
    seen = set()
    first = None
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5:
                species = parse_mash_species(fields[0])
                percent = mash_hash_percent(fields[4])
                hit = {
                    "species": species,
                    "percent": percent,
                    "percent_text": f"{percent:.1f}%" if percent is not None else "NA",
                    "reference": fields[0],
                    "distance": fields[2],
                    "pvalue": fields[3],
                    "hashes": fields[4],
                }
                if first is None:
                    first = hit
                if species not in seen:
                    top_species.append(hit)
                    seen.add(species)
                if len(top_species) >= 3:
                    break
    if first:
        summary = "; ".join(f"{hit['species']} {hit['percent_text']}" for hit in top_species)
        return {
                    "status": "done",
                    "summary": summary,
                    "reference": first["reference"],
                    "species": first["species"],
                    "distance": first["distance"],
                    "pvalue": first["pvalue"],
                    "hashes": first["hashes"],
                    "top_species": top_species,
                }
    return {"status": "done", "summary": "no hits", "reference": "no hits", "top_species": []}


def read_sourmash(sample):
    path = Path(f"data/results/{sample}/sourmash/gather.with-lineages.csv")
    if not path.exists():
        return {"status": "not run"}
    hits = []
    seen = set()
    with path.open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle):
            species = row.get("species", "").replace("s__", "").strip()
            if not species:
                lineage = row.get("lineage", "")
                species = next((part[3:] for part in lineage.split(";") if part.startswith("s__")), "")
            if not species or species in seen:
                continue
            value = row.get("f_match") or row.get("f_unique_to_query") or row.get("match") or ""
            try:
                percent = float(value) * 100
                percent_text = f"{percent:.1f}%"
            except (TypeError, ValueError):
                percent_text = "NA"
            hits.append({"species": species, "percent_text": percent_text})
            seen.add(species)
            if len(hits) >= 3:
                break
    summary = "; ".join(f"{hit['species']} {hit['percent_text']}" for hit in hits) or "no hits"
    return {"status": "done", "summary": summary, "top_species": hits}


def load_species_references(path):
    references = {}
    epithets = {}
    source = Path(path)
    if not source.exists():
        return references
    with source.open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = row.get("organism_name", "").strip()
            words = name.split()
            species = " ".join(words[:2]).lower() if len(words) >= 2 else name.lower()
            if species:
                references[species] = row
                epithet = species.split()[-1]
                epithets.setdefault(epithet, []).append(row)
    for epithet, rows in epithets.items():
        if len(rows) == 1:
            references[f"__epithet__:{epithet}"] = rows[0]
    return references


def reference_for_summary(summary, references):
    candidates = summary["sourmash"].get("top_species", []) + summary["mash"].get("top_species", [])
    for hit in candidates:
        species = hit.get("species", "").lower().replace("s__", "").strip()
        if species in references:
            return references[species]
        epithet = species.split()[-1] if species else ""
        if f"__epithet__:{epithet}" in references:
            return references[f"__epithet__:{epithet}"]
    return {}


def numeric(value):
    try:
        return float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return None


def qc_assessment(summary, reference):
    fastp = summary["fastp"]
    quast = summary["quast"]
    busco = summary["busco"]
    genome_size = numeric(quast.get("total_length"))
    ref_size = numeric(reference.get("genome_size"))
    gc = numeric(quast.get("gc"))
    ref_gc = numeric(reference.get("gc_percent"))
    coverage = numeric(fastp.get("filtered_bases"))
    coverage = coverage / genome_size if coverage is not None and genome_size else None
    size_diff = ((genome_size - ref_size) / ref_size * 100) if genome_size and ref_size else None
    gc_diff = abs(gc - ref_gc) if gc is not None and ref_gc is not None else None
    duplicated = numeric(busco.get("duplicated"))
    return [
        ("Q30 after QC", ">90%", fmt_percent(fastp.get("filtered_q30")), "pass" if numeric(fastp.get("filtered_q30")) and numeric(fastp.get("filtered_q30")) >= 0.9 else "caution"),
        ("Estimated coverage", "50-200x optimal", f"{coverage:.1f}x" if coverage is not None else "NA", "pass" if coverage is not None and 50 <= coverage <= 200 else "caution"),
        ("Contigs", "<50 good; <100 acceptable", fmt_int(quast.get("contigs")), "pass" if numeric(quast.get("contigs")) is not None and numeric(quast.get("contigs")) < 50 else "caution"),
        ("N50", ">100 kb minimum; >500 kb good", fmt_int(quast.get("n50")), "pass" if numeric(quast.get("n50")) is not None and numeric(quast.get("n50")) >= 100000 else "caution"),
        ("Genome size vs reference", "close to selected species", f"{size_diff:+.1f}%" if size_diff is not None else "NA", "pass" if size_diff is not None and abs(size_diff) <= 10 else "caution"),
        ("GC vs reference", "within +/-1-2%", f"{gc_diff:.2f} percentage points" if gc_diff is not None else "NA", "pass" if gc_diff is not None and gc_diff <= 2 else "caution"),
        ("BUSCO duplicated", "low expected", f"{duplicated:.1f}%" if duplicated is not None else "NA", "pass" if duplicated is not None and duplicated <= 10 else "red flag"),
    ]


def assessment_label(value):
    return f"⚠️ {value}" if value == "red flag" else value


def pdf_problem_suffix(summary, metric):
    assessments = {name: status for name, _, _, status in summary.get("qc_assessment", [])}
    return " (!)" if assessments.get(metric) in {"caution", "red flag"} else ""


def sample_summary(sample):
    name = sample["sample"]
    return {
        "sample": name,
        "layout": sample.get("layout", ""),
        "fastp": read_fastp_summary(name),
        "rasusa": read_rasusa(name, sample.get("layout", "")),
        "quast": read_quast(name),
        "bakta": count_bakta_features(name),
        "amr_hits": count_amr_hits(name),
        "plasmid_hits": count_tsv_hits(f"data/results/{name}/plasmidfinder/{name}.tsv"),
        "vfdb_hits": count_tsv_hits(f"data/results/{name}/vfdb/{name}.tsv"),
        "kofam_hits": count_tsv_hits(f"data/results/{name}/kofamscan/{name}.tsv"),
        "keggdecoder": count_keggdecoder_functions(name),
        "dbcan_hits": count_dbcan_hits(name),
        "pgmdb": read_pgmdb(name),
        "busco": read_busco(name),
        "mash": read_mash(name),
        "sourmash": read_sourmash(name),
    }


def compact_analysis_summary(sample):
    result = {}

    fastani = []
    for line in Path(f"data/results/{sample}/fastani/{sample}_fastani.tsv").read_text(errors="replace").splitlines() if Path(f"data/results/{sample}/fastani/{sample}_fastani.tsv").exists() else []:
        fields = line.split("\t")
        if len(fields) >= 3:
            try:
                fastani.append((float(fields[2]), Path(fields[1]).name.split("_genomic")[0]))
            except ValueError:
                pass
    references = {}
    for row in read_tsv_rows(f"data/results/{sample}/fastani/refs/references.tsv"):
        accession = row.get("assembly_accession") or row.get("Assembly Accession") or ""
        if not accession and len(row) == 1:
            continue
        references[accession] = row
    if fastani:
        ani, accession = max(fastani)
        ref = references.get(accession, {})
        label = ref.get("species") or ref.get("label") or accession
        strain = ref.get("strain") or ref.get("Strain") or ""
        result["fastani"] = f"{ani:.2f}% to {label}{f' strain {strain}' if strain and strain not in label else ''}; {'PASS' if ani > 95 else 'caution'}"

    mob_contigs = [
        row for row in read_tsv_rows(f"data/results/{sample}/mobsuite/contig_report.txt")
        if row.get("molecule_type", "").lower() == "plasmid"
    ]
    mob_plasmids = read_tsv_rows(f"data/results/{sample}/mobsuite/mobtyper_results.txt")
    result["mob"] = f"{len(mob_contigs)} plasmid contigs; {len(mob_plasmids)} reconstructed plasmids; {sum(r.get('predicted_mobility') == 'mobilizable' for r in mob_plasmids)} mobilizable"

    ge_plasmids = read_tsv_rows(f"data/results/{sample}/genomad/assembly_summary/assembly_plasmid_summary.tsv")
    ge_viruses = read_tsv_rows(f"data/results/{sample}/genomad/assembly_summary/assembly_virus_summary.tsv")
    result["genomad"] = f"{len(ge_plasmids)} plasmid contigs; {len(ge_viruses)} virus/provirus contigs"

    gutsmash, _ = parse_gutsmash_index_html(sample)
    gutsmash_status = gutsmash_status_record(sample)
    if gutsmash_status and gutsmash_status.get("status") != "completed":
        error = gutsmash_status.get("error", "").strip() or "runtime unavailable"
        result["gutsmash"] = f"unavailable; {error}"
    else:
        top_cluster = gutsmash[0].get("cluster", "") if gutsmash else ""
        result["gutsmash"] = f"{len(gutsmash)} metabolic cluster{'s' if len(gutsmash) != 1 else ''}" + (f"; top match {top_cluster}" if top_cluster else "")

    genome_rows = {row.get("target", ""): row for row in read_tsv_rows(f"data/results/{sample}/genomespot/{sample}.predictions.tsv")}
    oxygen = genome_rows.get("oxygen", {}).get("value", "NA")
    temp = numeric(genome_rows.get("temperature_optimum", {}).get("value"))
    ph = numeric(genome_rows.get("ph_optimum", {}).get("value"))
    result["genomespot"] = f"oxygen {oxygen}; optimum {temp:.1f} C; pH {ph:.1f}" if temp is not None and ph is not None else f"oxygen {oxygen}"

    antismash = read_tsv_rows(f"data/results/{sample}/antismash/{sample}_regions.tsv")
    products = sorted({row.get("Products", "") for row in antismash if row.get("Products")})
    result["antismash"] = f"{len(antismash)} biosynthetic clusters" + (f"; {', '.join(products[:3])}" if products else "")

    gut_rows = read_tsv_rows(f"data/results/{sample}/gut_brain_markers.tsv")
    gut_markers = {
        marker.strip()
        for row in gut_rows for marker in row.get("marker", "").split(";")
        if marker.strip() and marker.strip() != "-"
    }
    trait_counts = Counter(
        trait.strip()
        for row in gut_rows for trait in row.get("metabolite_or_trait", "").split("|")
        if trait.strip()
    )
    result["gut_brain"] = f"{len(gut_markers)} primary markers; top traits " + ", ".join(name for name, _ in trait_counts.most_common(2))

    go_rows = read_tsv_rows(f"data/results/{sample}/goprobioticprofiler/probiotic_go_categories.tsv")
    significant = [row for row in go_rows if row.get("Evidence_Tier") == "Significant"]
    supportive = [row for row in go_rows if row.get("Evidence_Tier") == "Supportive"]
    categories = [row.get("Probiotic_Category", "") for row in significant + supportive if row.get("Probiotic_Category")]
    result["go"] = (
        f"{len(significant)} significant, {len(supportive)} supportive probiotic GO categories"
        + (f"; {', '.join(categories[:3])}" if categories else "")
    )
    return result


def workflow_steps():
    return [
        ("Raw FASTQ", "Input reads"),
        ("fastp / Porechop", "Adapter/QC"),
        ("Rasusa", "Read subsampling"),
        ("FastQC", "Read QC"),
        ("Unicycler / Flye", "Assembly"),
        ("Medaka", "ONT polishing"),
        ("QUAST", "Assembly QC"),
        ("Mash", "Assembly ID"),
        ("sourmash", "GTDB ID"),
        ("BUSCO", "Completeness"),
        ("Bakta light", "Annotation"),
        ("AMRFinderPlus", "AMR scan"),
        ("PlasmidFinder", "Plasmids"),
        ("VFDB", "Virulence"),
        ("KofamScan", "KO functions"),
        ("KEGGDecoder", "Pathway completeness"),
        ("dbCAN", "CAZymes"),
        ("PGMDB", "Probiotic genes"),
        ("antiSMASH", "BGC scan"),
        ("MultiQC", "Combined QC"),
        ("Final report", "HTML/PDF"),
    ]


def workflow_layout():
    return {"rows": 7, "box_w": 180, "box_h": 48, "gap_x": 62, "gap_y": 28, "x0": 34, "y0": 28}


def workflow_svg():
    steps = workflow_steps()
    layout = workflow_layout()
    rows = layout["rows"]
    box_w = layout["box_w"]
    box_h = layout["box_h"]
    gap_x = layout["gap_x"]
    gap_y = layout["gap_y"]
    x0 = layout["x0"]
    y0 = layout["y0"]
    cols = (len(steps) + rows - 1) // rows
    width = x0 * 2 + cols * box_w + (cols - 1) * gap_x
    height = y0 * 2 + rows * box_h + (rows - 1) * gap_y
    parts = [
        f"<svg class='workflow' viewBox='0 0 {width} {height}' role='img' aria-label='Workflow diagram'>",
        "<defs><marker id='arrow' markerWidth='10' markerHeight='10' refX='7' refY='3' orient='auto' markerUnits='strokeWidth'><path d='M0,0 L0,6 L8,3 z' fill='#334155'/></marker></defs>",
        f"<rect x='0' y='0' width='{width}' height='{height}' fill='#f8fafc'/>",
    ]
    for i, (name, desc) in enumerate(steps):
        col = i // rows
        row = i % rows
        x = x0 + col * (box_w + gap_x)
        y = y0 + row * (box_h + gap_y)
        fill = "#e0f2fe" if i < 2 else "#ecfdf5" if i < 7 else "#fef3c7" if i < 13 else "#ede9fe"
        parts.append(f"<rect x='{x}' y='{y}' width='{box_w}' height='{box_h}' rx='6' fill='{fill}' stroke='#334155' stroke-width='1.2'/>")
        parts.append(f"<text x='{x + 12}' y='{y + 22}' font-family='Arial' font-size='13' font-weight='700' fill='#0f172a'>{html.escape(name)}</text>")
        parts.append(f"<text x='{x + 12}' y='{y + 38}' font-family='Arial' font-size='11' fill='#475569'>{html.escape(desc)}</text>")
        if i < len(steps) - 1:
            if row < rows - 1 and i + 1 < len(steps):
                parts.append(f"<line x1='{x + box_w / 2}' y1='{y + box_h}' x2='{x + box_w / 2}' y2='{y + box_h + gap_y - 8}' stroke='#334155' stroke-width='1.4' marker-end='url(#arrow)'/>")
            else:
                nx = x0 + (col + 1) * (box_w + gap_x)
                ny = y0
                parts.append(f"<path d='M{x + box_w},{y + box_h / 2} C{x + box_w + 28},{y + box_h / 2} {nx - 28},{ny + box_h / 2} {nx},{ny + box_h / 2}' fill='none' stroke='#334155' stroke-width='1.4' marker-end='url(#arrow)'/>")
    parts.append("</svg>")
    return "\n".join(parts)


def pdf_escape(text):
    return str(text).replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


class PdfCanvas:
    width = 595
    height = 842

    def __init__(self):
        self.pages = []
        self.page_titles = []
        self.current = []
        self.current_title = ""
        self.images = []

    def add_page(self):
        if self.current:
            self.pages.append("\n".join(self.current))
            self.page_titles.append(self.current_title)
        self.current = []
        self.current_title = ""

    def finish(self):
        if self.current:
            self.pages.append("\n".join(self.current))
            self.page_titles.append(self.current_title)
            self.current = []

    def rect(self, x, y, w, h, fill=None, stroke=None, line_width=1):
        cmds = ["q"]
        if fill:
            cmds.append(f"{fill[0]} {fill[1]} {fill[2]} rg")
        if stroke:
            cmds.append(f"{stroke[0]} {stroke[1]} {stroke[2]} RG")
            cmds.append(f"{line_width} w")
        op = "B" if fill and stroke else "f" if fill else "S"
        cmds.append(f"{x} {y} {w} {h} re {op}")
        cmds.append("Q")
        self.current.append("\n".join(cmds))

    def line(self, x1, y1, x2, y2, color=(0.2, 0.25, 0.33), width=1):
        self.current.append(f"q {color[0]} {color[1]} {color[2]} RG {width} w {x1} {y1} m {x2} {y2} l S Q")

    def warning(self, x, y, size=10):
        self.current.append(
            f"q 1 0.72 0 rg 0.55 0.38 0 RG 0.7 w "
            f"{x} {y} m {x + size} {y} l {x + size / 2} {y + size} l h B Q"
        )
        self.text(x + size / 2 - 1.2, y + 2.2, "!", size=max(5, size - 3), color=(0.2, 0.15, 0.05), bold=True)

    def text(self, x, y, text, size=10, color=(0.08, 0.1, 0.16), bold=False):
        if bold and size >= 14 and y >= 770 and not self.current_title:
            self.current_title = str(text)
        if bold and size >= 12 and y >= 770 and x <= 42 and color == (0.05, 0.16, 0.28):
            color = pdf_palette(text)[0]
            self.rect(x, y - 3, 5, 18 if size >= 14 else 14, fill=color)
            x += 12
        font = "F2" if bold else "F1"
        safe = pdf_escape(text)
        self.current.append(f"BT /{font} {size} Tf {color[0]} {color[1]} {color[2]} rg {x} {y} Td ({safe}) Tj ET")

    def warning_text(self, x, y, text, size=9, color=(0.16, 0.21, 0.28)):
        parts = str(text).split(" (!)")
        cursor = x
        for index, part in enumerate(parts):
            if part:
                self.text(cursor, y, part, size=size, color=color)
                cursor += len(part) * size * 0.52
            if index < len(parts) - 1:
                self.text(cursor + 3, y, "(!)", size=size, color=(0.78, 0.05, 0.05), bold=True)
                cursor += size * 2.5

    def image(self, x, y, w, h, image):
        if image is None:
            return
        name = f"Im{len(self.images) + 1}"
        self.images.append((name, image))
        self.current.append(f"q {w} 0 0 {h} {x} {y} cm /{name} Do Q")

    def chart_image(self, x, y, w, h, image):
        if image is None:
            return
        self.rect(x - 6, y - 6, w + 12, h + 12, fill=(1, 1, 1), stroke=(0.78, 0.83, 0.88), line_width=0.5)
        self.image(x, y, w, h, image)

    def wrapped_text(self, x, y, text, max_chars=70, size=9, color=(0.2, 0.25, 0.33), bold=False, leading=12):
        lines = []
        words = str(text).split()
        current = ""
        for word in words:
            candidate = f"{current} {word}".strip()
            if len(candidate) > max_chars and current:
                lines.append(current)
                current = word
            else:
                current = candidate
        if current:
            lines.append(current)
        for line in lines:
            self.text(x, y, line, size=size, color=color, bold=bold)
            y -= leading
        return y


def add_results_tree_page(pdf, summaries):
    if not summaries:
        return
    pdf.add_page()
    pdf.text(40, 800, "Results Directory Structure", size=15, color=(0.05, 0.16, 0.28), bold=True)
    for sm in summaries:
        sample = sm["sample"]
        pdf.text(40, 780, f"Sample {sample} - main output files from each analysis tool", size=9, color=(0.38, 0.45, 0.55))
        tree_file = Path(f"data/results/{sample}/tree.txt")
        if not tree_file.exists():
            continue
        lines = tree_file.read_text(encoding="utf-8").strip().split("\n")
        yy = 755
        step = 10
        detail_size = 5.2
        program_size = 6.2
        for line in lines:
            if yy < 48:
                pdf.add_page()
                pdf.text(40, 800, "Results Directory Structure - continued", size=12, color=(0.05, 0.16, 0.28), bold=True)
                yy = 775
            if line.startswith("data/results/"):
                pdf.rect(40, yy - 3, 515, step, fill=(0.91, 0.95, 0.99), stroke=(0.75, 0.82, 0.9))
                pdf.text(48, yy, line, size=max(5.2, program_size), color=(0.05, 0.16, 0.28), bold=True)
            elif line.startswith("+-- "):
                successful = "✅" in line
                program = line.replace("+-- ✅ ", "").replace("+-- ❌ ", "")
                pdf.rect(40, yy - 3, 515, step, fill=(0.91, 0.98, 0.94) if successful else (0.96, 0.96, 0.96))
                pdf.text(48, yy, "OK" if successful else "--", size=program_size, color=(0.03, 0.48, 0.18) if successful else (0.48, 0.5, 0.54), bold=True)
                pdf.text(70, yy, program, size=program_size, color=(0.05, 0.16, 0.28), bold=True)
            else:
                pdf.text(52, yy, line.replace("✅", "OK").replace("❌", ""), size=detail_size, color=(0.28, 0.34, 0.42))
            yy -= step
        break


def add_toc_and_page_numbers(pdf):
    """Insert a contents page after the overview and number every page."""
    titles = list(pdf.page_titles)
    if len(titles) > 1:
        titles[1] = "Overview and Workflow"
    toc = PdfCanvas()
    toc.add_page()
    toc.text(40, 800, "Contents", size=18, color=(0.05, 0.16, 0.28), bold=True)
    y = 770
    for index, title in enumerate(titles):
        if index == 0:
            continue
        if not title:
            continue
        page_number = index + 1 if index == 1 else index + 2
        label = title[:76]
        accent, _, _ = pdf_palette(title)
        toc.rect(40, y + 1, 4, 4, fill=accent)
        toc.text(48, y, label, size=8, color=(0.16, 0.21, 0.28))
        line_y = y + 2
        line_start = 48 + len(label) * 4.0 + 4
        line_end = 520
        if line_start < line_end:
            toc.current.append(
                f"q 0.72 0.76 0.82 RG 0.4 w [2 3] 0 d {line_start} {line_y} m {line_end} {line_y} l S Q"
            )
        toc.text(535, y, str(page_number), size=8, color=(0.05, 0.16, 0.28), bold=True)
        y -= 13
        if y < 48:
            break
    toc.finish()
    pdf.pages.insert(2, toc.pages[0])
    pdf.page_titles.insert(2, "Contents")
    total = len(pdf.pages)
    for index in range(total):
        if index == 0:
            footer = (
                f"q 0.72 0.76 0.82 RG 0.4 w 40 28 m 555 28 l S Q\n"
                f"BT /F2 8 Tf 0.12 0.15 0.32 rg 257 15 Td (workflow version {pdf_escape(workflow_version())}) Tj ET"
            )
        else:
            footer = (
                f"q 0.72 0.76 0.82 RG 0.4 w 40 28 m 555 28 l S Q\n"
                f"BT /F1 7 Tf 0.38 0.45 0.55 rg 40 15 Td (University of Plovdiv, Dept. Molecular Biology \\251) Tj ET\n"
                f"BT /F1 7 Tf 0.38 0.45 0.55 rg 520 15 Td (Page {index + 1} of {total}) Tj ET"
            )
        pdf.pages[index] += "\n" + footer


def svg_to_png(svg_path):
    svg_path = Path(svg_path)
    if not svg_path.exists():
        return None
    png_path = svg_path.with_suffix(".report.png")
    converter = Path("mamba_envs/snakemake-conda-helper-2026-06-04/bin/rsvg-convert")
    if not png_path.exists() or png_path.stat().st_mtime < svg_path.stat().st_mtime:
        subprocess.run([str(converter), "-o", str(png_path), str(svg_path)], check=True)
    return png_path


def draw_genomad_pdf(pdf, sample):
    plasmid_path = Path(f"data/results/{sample}/genomad/assembly_summary/assembly_plasmid_summary.tsv")
    virus_path = Path(f"data/results/{sample}/genomad/assembly_summary/assembly_virus_summary.tsv")
    plasmids = read_tsv_rows(plasmid_path) if plasmid_path.exists() else []
    viruses = read_tsv_rows(virus_path) if virus_path.exists() else []
    if not plasmid_path.exists() and not virus_path.exists():
        return
    y = draw_pdf_table(
        pdf, f"geNomad Plasmid and Virus Summary: {sample}",
        f"Plasmid sequences identified by geNomad; {len(plasmids)} contig(s) classified as plasmids",
        ("Contig", "Length", "Topology", "Plasmid score", "FDR", "Hallmarks", "Marker enrichment"),
        [(r.get("seq_name", ""), r.get("length", ""), r.get("topology", ""), r.get("plasmid_score", ""),
          r.get("fdr", ""), r.get("n_hallmarks", ""), r.get("marker_enrichment", "")) for r in plasmids],
        (105, 50, 75, 65, 45, 65, 110),
        row_height=18,
    )
    draw_pdf_table(
        pdf, "geNomad Virus Summary",
        f"Virus sequences identified by geNomad; {len(viruses)} contig(s) classified as viruses",
        ("Contig", "Length", "Topology", "Virus score", "FDR", "Hallmarks", "Marker enrichment"),
        [(r.get("seq_name", ""), r.get("length", ""), r.get("topology", ""), r.get("virus_score", ""),
          r.get("fdr", ""), r.get("n_hallmarks", ""), r.get("marker_enrichment", "")) for r in viruses],
        (105, 50, 75, 65, 45, 65, 110),
        row_height=18,
        new_page=False,
        start_y=y - 25,
    )


def draw_mobrecon_pdf(pdf, sample):
    contig_path = Path(f"data/results/{sample}/mobsuite/contig_report.txt")
    plasmid_path = Path(f"data/results/{sample}/mobsuite/mobtyper_results.txt")
    contigs = [
        row for row in read_tsv_rows(contig_path)
        if row.get("molecule_type", "").strip().lower() == "plasmid"
    ] if contig_path.exists() else []
    plasmids = read_tsv_rows(plasmid_path) if plasmid_path.exists() else []
    if not contig_path.exists() and not plasmid_path.exists():
        return
    y = draw_pdf_table(
        pdf, f"MOB-recon Contig Classification: {sample}",
        f"Only contigs classified as plasmids are shown; {len(contigs)} plasmid contig(s)",
        ("Contig", "Length", "Rep types", "Relaxases", "Mobility", "Nearest neighbor"),
        [(r.get("contig_id", ""), r.get("size", ""), r.get("rep_type(s)", ""), r.get("relaxase_type(s)", ""),
          r.get("predicted_mobility", ""), r.get("mash_neighbor_identification", "")) for r in contigs],
        (95, 55, 80, 80, 80, 125),
        row_height=20,
    )
    draw_pdf_table(
        pdf, "MOB-recon Reconstructed Plasmids",
        f"Reconstructed plasmid groups and predicted host range; {len(plasmids)} group(s)",
        ("Plasmid", "Contigs", "Length", "Rep types", "Mobility", "Nearest neighbor", "Host range"),
        [(r.get("sample_id", ""), r.get("num_contigs", ""), r.get("size", ""), r.get("rep_type(s)", ""),
          r.get("predicted_mobility", ""), r.get("mash_neighbor_identification", ""),
          r.get("predicted_host_range_overall_name", "")) for r in plasmids],
        (65, 50, 50, 70, 70, 110, 100),
        row_height=20,
        new_page=False,
        start_y=y - 25,
    )


def draw_bagel_antismash_pdf(pdf, sample):
    bagel_rows = bagel4_rows(sample)
    if bagel_rows:
        bagel_y = draw_pdf_table(
            pdf, f"BAGEL4 Bacteriocin and RiPP Clusters: {sample}",
            "Areas of Interest detected by BAGEL4",
            ("AOI", "Start", "End", "Class", "Sequence"),
            [(row.get("AOI", ""), row.get("start", ""), row.get("end", ""), row.get("class", ""), row.get("fullname", "")) for row in bagel_rows],
            (120, 65, 65, 115, 150),
            row_height=22,
        )
        image_y = bagel_y - 30
        for image_path in bagel4_images(sample, (".png",)):
            image = load_png_for_pdf(image_path)
            if image:
                width = 500
                height = width * image["height"] / image["width"]
                if height > 270:
                    scale = 270 / height
                    width *= scale
                    height *= scale
                if image_y - height < 55:
                    pdf.add_page()
                    pdf.text(40, 800, f"BAGEL4 Cluster Graphics: {sample}", size=15, color=(0.05, 0.16, 0.28), bold=True)
                    image_y = 760
                pdf.text(48, image_y + 8, image_path.stem, size=8, color=(0.25, 0.31, 0.39), bold=True)
                pdf.chart_image(40 + (515 - width) / 2, image_y - height, width, height, image)
                image_y -= height + 34

    as_regions_path = Path(f"data/results/{sample}/antismash/{sample}_regions.tsv")
    if as_regions_path.exists():
        as_regions = read_tsv_rows(as_regions_path)
        if as_regions:
            draw_pdf_table(
                pdf, f"antiSMASH Biosynthetic Gene Clusters: {sample}",
                "Predicted biosynthetic gene clusters with product type and gene count",
                ("Region", "Location", "Length", "Products", "Genes"),
                [(r.get("Region", ""), r.get("Location", ""), r.get("Length", ""),
                  r.get("Products", ""), r.get("Genes", "")) for r in as_regions],
                (50, 100, 65, 160, 50),
                row_height=22,
            )


def write_pdf(path, title, summaries, now):
    pdf = PdfCanvas()
    report_dir = Path(path).parent
    header_blue = (31 / 255, 38 / 255, 82 / 255)

    # Dedicated cover page.
    pdf.add_page()
    pdf.rect(0, 760, 595, 82, fill=header_blue)
    pdf.text(40, 812, title, size=13, color=(1, 1, 1), bold=True)
    pdf.text(40, 790, f"Generated: {now}", size=10, color=(0.8, 0.9, 1))
    pdf.image(428, 771, 142, 62, load_png_for_pdf("data/logo.png", header_blue))
    cover_logo = load_png_for_pdf("data/probeat_logo2.png")
    if cover_logo:
        logo_size = 420
        pdf.image((595 - logo_size) / 2, 190, logo_size, logo_size, cover_logo)

    # Existing overview page.
    pdf.add_page()
    pdf.rect(0, 760, 595, 82, fill=header_blue)
    pdf.text(40, 812, title, size=13, color=(1, 1, 1), bold=True)
    pdf.text(40, 790, f"Generated: {now}", size=10, color=(0.8, 0.9, 1))
    pdf.image(428, 771, 142, 62, load_png_for_pdf("data/logo.png", header_blue))

    if summaries:
        first = summaries[0]
        fastp = first["fastp"]
        quast = first["quast"]
        sourmash_species = first["sourmash"].get("top_species", [])
        mash_species = first["mash"].get("top_species", [])
        species = (
            sourmash_species[0].get("species")
            if sourmash_species
            else mash_species[0].get("species") if mash_species else "Not identified"
        )
        cards = [
            ("Raw reads", fmt_int(fastp.get("raw_reads"))),
            ("Contigs (>= 500 bp)", fmt_int(quast.get("contigs"))),
            ("Genome size", f"{fmt_int(quast.get('total_length'))} bp"),
            ("GC", fmt_percent(quast.get("gc"))),
            ("Species", species),
            ("N50", fmt_int(quast.get("n50"))),
        ]
        x = 40
        y = 706
        for i, (label, value) in enumerate(cards):
            col = i % 3
            row = i // 3
            cx = x + col * 178
            cy = y - row * 58
            pdf.rect(cx, cy, 160, 43, fill=(0.96, 0.98, 1), stroke=(0.75, 0.82, 0.9))
            value_size = 9 if label == "Species" else 13
            if label == "Species" and len(str(value)) > 28:
                value_size = 7
            pdf.text(cx + 10, cy + 26, value, size=value_size, color=(0.05, 0.16, 0.28), bold=True)
            pdf.text(cx + 10, cy + 10, label, size=7, color=(0.38, 0.45, 0.55))

    pdf.text(40, 582, "Workflow", size=14, color=(0.05, 0.16, 0.28), bold=True)
    pdf.image(62, 20, 470, 524, load_png_for_pdf("data/workflow.png"))

    pdf.add_page()
    pdf.text(40, 790, "Per-sample Summary", size=15, color=(0.05, 0.16, 0.28), bold=True)
    y = 760
    for summary in summaries:
        fastp = summary["fastp"]
        quast = summary["quast"]
        bakta = summary["bakta"]
        busco = summary["busco"]
        mash = summary["mash"]
        sourmash = summary["sourmash"]
        pgmdb = summary["pgmdb"]
        rasusa = summary["rasusa"]
        extra = compact_analysis_summary(summary["sample"])
        duplicated = busco.get("duplicated")
        if busco.get("status") == "done":
            busco_lines = [
                f"BUSCO: Complete {busco.get('complete', 'NA')}%; single-copy {busco.get('single', 'NA')}%",
                f"BUSCO duplicated: {duplicated or 'NA'}%; "
                f"fragmented {busco.get('fragmented', 'NA')}%; missing {busco.get('missing', 'NA')}%",
            ]
        else:
            busco_lines = [f"BUSCO: {busco.get('status', 'not run')}"]
        sourmash_summary = sourmash.get("summary", sourmash.get("status", "not run"))
        sourmash_species = [item.strip() for item in sourmash_summary.split(";") if item.strip()]
        pg_categories = ", ".join(
            f"{cat['category']} {cat['unique_markers']}"
            for cat in (pgmdb.get("categories") or [])
            if cat["category"].lower() in ("acid stress", "bile resistance", "exopolysaccharide production", "vitamin biosynthesis")
        )
        lines = [
            (f"Sample {summary['sample']} ({summary['layout']})", 0),
            (f"Raw reads/bases: {fmt_int(fastp.get('raw_reads'))} / {fmt_int(fastp.get('raw_bases'))}", 0),
            (f"After adapter/QC reads/bases: {fmt_int(fastp.get('filtered_reads'))} / {fmt_int(fastp.get('filtered_bases'))}; Q30 {fmt_percent(fastp.get('filtered_q30'))}", 0),
            (f"Adapter-trimmed reads/bases: {fmt_int(fastp.get('adapter_trimmed_reads'))} / {fmt_int(fastp.get('adapter_trimmed_bases'))}", 0),
            (f"Rasusa: {rasusa.get('status', 'not run')}; {fmt_int(rasusa.get('input_reads'))} -> {fmt_int(rasusa.get('output_reads'))} reads", 0),
            (f"Coverage: {summary['qc_assessment'][1][2]}", 0),
            (f"Assembly: {fmt_int(quast.get('contigs'))} contigs; {fmt_int(quast.get('total_length'))} bp; GC {fmt_percent(quast.get('gc'))}; N50 {fmt_int(quast.get('n50'))}", 0),
            (f"Taxonomy: sourmash {sourmash_species[0] if sourmash_species else 'not run'}; Mash {mash.get('species', mash.get('status', 'not run'))}", 0),
            (f"FastANI: {extra.get('fastani', 'not run')}", 0),
            (f"Bakta: CDS {fmt_int(bakta.get('CDS'))}; tRNA {fmt_int(bakta.get('tRNA'))}; rRNA {fmt_int(bakta.get('rRNA'))}; tmRNA {fmt_int(bakta.get('tmRNA'))}; CRISPR {fmt_int(bakta.get('CRISPR'))}", 0),
            (f"AMRFinderPlus AMR hits: {fmt_int(summary['amr_hits'])}", 0),
            (f"PlasmidFinder replicon hits: {fmt_int(summary['plasmid_hits'])}", 0),
            (f"VFDB virulence-factor hits: {fmt_int(summary['vfdb_hits'])}", 0),
            (f"KofamScan KO annotations: {fmt_int(summary['kofam_hits'])}", 0),
            (f"KEGGDecoder functions present: {fmt_int(summary['keggdecoder'][0])} / {fmt_int(summary['keggdecoder'][1])}" if summary["keggdecoder"] else "KEGGDecoder: not run", 0),
            (f"dbCAN CAZymes: {fmt_int(summary['dbcan_hits'])}", 0),
            (f"PGMDB: {fmt_int(pgmdb.get('unique_markers'))} unique; {pg_categories}" if pgmdb.get("status") == "done" else "PGMDB: not run", 0),
            (f"MOB-recon: {extra.get('mob', 'not run')}", 0),
            (f"geNomad: {extra.get('genomad', 'not run')}", 0),
            (f"gutSMASH: {extra.get('gutsmash', 'not run')[:92]}", 0),
            (f"GenomeSPOT: {extra.get('genomespot', 'not run')}", 0),
            (f"antiSMASH: {extra.get('antismash', 'not run')[:92]}", 0),
            (f"Gut-Brain: {extra.get('gut_brain', 'not run')}", 0),
            (f"GO Probiotic: {extra.get('go', 'not run')}", 0),
            (f"BAGEL4 bacteriocin/RiPP clusters: {fmt_int(len(bagel4_rows(summary['sample'])))}", 0),
        ] + [(line, 0) for line in busco_lines]
        section_specs = [
            ("QC & Assembly", PDF_PALETTES["blue"][0], (
                "Raw reads", "After adapter", "Adapter-trimmed", "Rasusa", "Coverage", "Assembly", "BUSCO",
            )),
            ("Identification & Annotation", PDF_PALETTES["orange"][0], (
                "Taxonomy", "FastANI", "Bakta", "MOB-recon", "geNomad",
            )),
            ("Safety", PDF_PALETTES["red"][0], (
                "AMRFinderPlus", "PlasmidFinder", "VFDB",
            )),
            ("Functional & Probiotic Potential", PDF_PALETTES["green"][0], (
                "KofamScan", "KEGGDecoder", "dbCAN", "PGMDB", "gutSMASH", "GenomeSPOT",
                "antiSMASH", "Gut-Brain", "GO Probiotic", "BAGEL4",
            )),
        ]
        grouped_lines = []
        for section_name, section_color, prefixes in section_specs:
            section_lines = [item for item in lines[1:] if item[0].startswith(prefixes)]
            if section_lines:
                grouped_lines.append((f"__SECTION__{section_name}", section_color))
                grouped_lines.extend(section_lines)
        lines = [lines[0]] + grouped_lines
        rendered_lines = []
        for line, indent in lines[1:]:
            if line.startswith("__SECTION__"):
                rendered_lines.append((line, indent, None))
                continue
            label, separator, value = line.partition(":")
            if separator:
                max_chars = 86 if label == "GO Probiotic" else 104
                wrapped = []
                current = ""
                for word in value.strip().split():
                    candidate = f"{current} {word}".strip()
                    if len(candidate) > max_chars and current:
                        wrapped.append(current)
                        current = word
                    else:
                        current = candidate
                wrapped.append(current or "")
                rendered_lines.append((line, indent, wrapped))
            else:
                rendered_lines.append((line, indent, [line]))
        line_height = 14
        box_height = 34 + sum(line_height if line.startswith("__SECTION__") else line_height * max(1, len(wrapped or [])) for line, _, wrapped in rendered_lines)
        if y - box_height < 45:
            pdf.add_page()
            y = 790
        pdf.rect(40, y - box_height + 10, 515, box_height, fill=(0.99, 0.99, 0.98), stroke=(0.82, 0.86, 0.9))
        pdf.text(55, y - 4, lines[0][0], size=12, color=(0.05, 0.16, 0.28), bold=True)
        yy = y - 24
        for line, indent, wrapped_value in rendered_lines:
            if line.startswith("__SECTION__"):
                label = line.replace("__SECTION__", "")
                pdf.rect(50, yy - 4, 495, 13, fill=(0.94, 0.96, 0.98))
                pdf.rect(50, yy - 4, 4, 13, fill=indent)
                pdf.text(61, yy, label, size=8, color=indent, bold=True)
                yy -= 14
                continue
            label, separator, value = line.partition(":")
            if separator:
                pdf.text(55 + indent, yy, f"{label}:", size=9, color=(0.10, 0.16, 0.24), bold=True)
                value_x = 55 + indent + len(label + ":") * 5.8 + 8
                for offset, wrapped_line in enumerate(wrapped_value or [""]):
                    pdf.warning_text(value_x if offset == 0 else 55 + indent + 14, yy, wrapped_line, size=9, color=(0.16, 0.21, 0.28))
                    if offset < len(wrapped_value or []) - 1:
                        yy -= 14
            else:
                pdf.warning_text(55 + indent, yy, line, size=9, color=(0.16, 0.21, 0.28))
            yy -= 14
        y -= box_height + 18

    # Flye Assembly Report (ONT samples only)
    for summary in summaries:
        if summary.get("layout") != "ONT":
            continue
        sample = summary["sample"]
        flye_info_path = Path(f"data/results/{sample}/flye/assembly_info.txt")
        if not flye_info_path.exists():
            continue
        flye_rows = []
        for line in flye_info_path.read_text(errors="replace").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 7:
                flye_rows.append(fields)
        if not flye_rows:
            continue
        flye_headers = ("Contig", "Length", "Coverage", "Circular", "Repeat", "Multiplicity", "Alt group")
        flye_widths = (90, 70, 60, 55, 55, 70, 55)
        flye_table_rows = []
        for fields in flye_rows:
            circ = "Yes" if fields[3] == "Y" else "No"
            repeat = "Yes" if fields[4] == "Y" else "No"
            flye_table_rows.append((fields[0], fmt_int(fields[1]), fields[2], circ, repeat, fields[5], fields[6]))
        draw_pdf_table(pdf, f"Flye Assembly Report: {sample}", "ONT long-read contig statistics", flye_headers, flye_table_rows, flye_widths)

    for summary in summaries:
        sample = summary["sample"]
        pdf.add_page()
        pdf.text(40, 800, f"Bakta Report: {sample}", size=15, color=(0.05, 0.16, 0.28), bold=True)
        bakta_text_path = Path(f"data/results/{sample}/bakta/{sample}.txt")
        bakta_lines = bakta_text_path.read_text(errors="replace").splitlines() if bakta_text_path.exists() else []
        yy = 775
        for line in bakta_lines:
            pdf.text(40, yy, line, size=8, color=(0.16, 0.21, 0.28), bold=line.endswith(":"))
            yy -= 12
        bakta_image = load_png_for_pdf(f"data/results/{sample}/bakta/{sample}.png")
        pdf.chart_image(285, 410, 270, 270, bakta_image)

        genome_map = load_png_for_pdf(f"data/results/{sample}/genome_map/{sample}_genome_map.png")
        if genome_map:
            pdf.add_page()
            pdf.text(40, 800, f"Integrated Bacterial Genome Map: {sample}", size=15, color=(0.05, 0.16, 0.28), bold=True)
            pdf.text(
                40, 780,
                "Bakta CDS, ABRicate, antiSMASH, BAGEL4, dbCAN CAZymes, PGMDB probiotic genes, GC content and GC skew",
                size=8, color=(0.38, 0.45, 0.55),
            )
            map_width = 490
            map_height = map_width * genome_map["height"] / genome_map["width"]
            pdf.chart_image(52, 742 - map_height, map_width, map_height, genome_map)

        # FastANI page
        heatmap_png = Path(f"data/results/{sample}/fastani/{sample}_ani_heatmap.png")
        if heatmap_png.exists():
            sourmash_csv = Path(f"data/results/{sample}/sourmash/gather.with-lineages.csv")
            species = ""
            if sourmash_csv.exists():
                import csv as _csv
                with sourmash_csv.open() as f:
                    for row in _csv.DictReader(f):
                        lineage = row.get("lineage", "")
                        parts = lineage.split(";")
                        for p in parts:
                            if p.startswith("s__"):
                                species = p[3:]
                        break
            pdf.add_page()
            pdf.text(40, 800, f"FastANI Species Confirmation: {sample}", size=15, color=(0.05, 0.16, 0.28), bold=True)
            pdf.text(40, 780, f"Assumed species: {species or 'unknown'} (from sourmash GTDB)", size=9, color=(0.38, 0.45, 0.55))
            image = load_png_for_pdf(heatmap_png)
            if image:
                width = 515
                height = width * image["height"] / image["width"]
                pdf.chart_image(40, 750 - height, width, height, image)
                ani_path = Path(f"data/results/{sample}/fastani/{sample}_fastani.tsv")
                ani_values = []
                if ani_path.exists():
                    for line in ani_path.read_text(errors="replace").splitlines():
                        fields = line.split("\t")
                        if len(fields) >= 3:
                            try:
                                ani_values.append(float(fields[2]))
                            except ValueError:
                                pass
                if ani_values and max(ani_values) > 95:
                    pdf.text(
                        40, 725 - height,
                        f"PASS: Best FastANI {max(ani_values):.2f}% > 95%.",
                        size=10, color=(0.02, 0.48, 0.20), bold=True,
                    )
                    pdf.text(
                        40, 708 - height,
                        f"The assembled strain is correctly assigned to {species}.",
                        size=10, color=(0.02, 0.48, 0.20), bold=True,
                    )

        draw_genomad_pdf(pdf, sample)
        draw_mobrecon_pdf(pdf, sample)

        amr_rows = read_tsv_rows(f"data/results/{sample}/amrfinder/{sample}.tsv")
        vfdb_rows = read_tsv_rows(f"data/results/{sample}/vfdb/{sample}.tsv")
        amr_exists = Path(f"data/results/{sample}/amrfinder/{sample}.tsv").exists()
        vfdb_exists = Path(f"data/results/{sample}/vfdb/{sample}.tsv").exists()

        if (amr_exists and not amr_rows) and (vfdb_exists and not vfdb_rows):
            # Both zero hits: one page
            pdf.add_page()
            pdf.text(40, 800, f"AMRFinderPlus & VFDB: {sample}", size=15, color=(0.05, 0.16, 0.28), bold=True)
            pdf.text(40, 770, "No antimicrobial resistance genes detected.", size=10, color=(0.38, 0.45, 0.55))
            pdf.text(40, 748, "No virulence-factor database matches detected.", size=10, color=(0.38, 0.45, 0.55))
        else:
            if amr_rows:
                table_rows = [
                    (row.get("Element symbol", ""), row.get("Element name", ""), row.get("Class", ""), row.get("% Identity to reference", ""))
                    for row in amr_rows
                ]
                draw_pdf_table(pdf, f"AMRFinderPlus Hits: {sample}", "Detected antimicrobial resistance-associated elements", ("Gene", "Element", "Class", "Identity"), table_rows, (85, 230, 145, 65))
            elif amr_exists:
                pdf.add_page()
                pdf.text(40, 800, f"AMRFinderPlus: {sample}", size=15, color=(0.05, 0.16, 0.28), bold=True)
                pdf.text(40, 770, "No antimicrobial resistance genes detected.", size=10, color=(0.38, 0.45, 0.55))
            if vfdb_rows:
                table_rows = [
                    (row.get("SEQUENCE", ""), row.get("GENE", ""), row.get("PRODUCT", ""), row.get("%COVERAGE", ""), row.get("%IDENTITY", ""))
                    for row in vfdb_rows
                ]
                draw_pdf_table(pdf, f"VFDB Hits: {sample}", "Detected virulence-factor database matches", ("Contig", "Gene", "Product", "Coverage", "Identity"), table_rows, (90, 65, 230, 65, 65))
            elif vfdb_exists:
                pdf.add_page()
                pdf.text(40, 800, f"VFDB: {sample}", size=15, color=(0.05, 0.16, 0.28), bold=True)
                pdf.text(40, 770, "No virulence-factor database matches detected.", size=10, color=(0.38, 0.45, 0.55))

        card_path = Path(f"data/results/{sample}/abricate/{sample}_card.tsv")
        megares_path = Path(f"data/results/{sample}/abricate/{sample}_megares.tsv")
        card_rows = [r for r in (read_tsv_rows(card_path) if card_path.exists() else []) if r.get("SEQUENCE")]
        megares_rows = [r for r in (read_tsv_rows(megares_path) if megares_path.exists() else []) if r.get("SEQUENCE")]

        if (card_path.exists() and not card_rows) and (megares_path.exists() and not megares_rows):
            pdf.add_page()
            pdf.text(40, 800, f"ABRicate CARD & MEGARES: {sample}", size=15, color=(0.05, 0.16, 0.28), bold=True)
            pdf.text(40, 770, "No CARD antimicrobial resistance genes detected.", size=10, color=(0.38, 0.45, 0.55))
            pdf.text(40, 748, "No MEGARES antimicrobial resistance genes detected.", size=10, color=(0.38, 0.45, 0.55))
        else:
            for db_name, rows in [("CARD", card_rows), ("MEGARES", megares_rows)]:
                if rows:
                    headers = ("Contig", "Start", "End", "Gene", "%Coverage", "%Identity")
                    widths = (90, 40, 40, 120, 65, 65)
                    table_rows = [
                        (r.get("SEQUENCE", ""), r.get("START", ""), r.get("END", ""),
                         r.get("GENE", ""), r.get("%COVERAGE", ""), r.get("%IDENTITY", ""))
                        for r in rows
                    ]
                    if db_name == "CARD":
                        headers = ("Contig", "Start", "End", "Gene", "%Coverage", "%Identity", "Resistance")
                        widths = (80, 38, 38, 105, 55, 55, 144)
                        table_rows = [
                            (r.get("SEQUENCE", ""), r.get("START", ""), r.get("END", ""),
                             r.get("GENE", ""), r.get("%COVERAGE", ""), r.get("%IDENTITY", ""),
                             r.get("RESISTANCE", ""))
                            for r in rows
                        ]
                    draw_pdf_table(
                        pdf, f"ABRicate {db_name} Hits: {sample}",
                        f"Antimicrobial resistance genes detected by {db_name} database",
                        headers,
                        table_rows,
                        widths,
                        row_height=18,
                    )

        draw_bagel_antismash_pdf(pdf, sample)

        classes = dbcan_class_counts(sample)
        class_names = {
            "GH": "Glycoside Hydrolases", "GT": "Glycosyltransferases",
            "PL": "Polysaccharide Lyases", "CE": "Carbohydrate Esterases",
            "AA": "Auxiliary Activities", "CBM": "Carbohydrate-Binding Modules",
        }
        draw_pdf_table(
            pdf, f"dbCAN CAZyme Summary: {sample}",
            f"CAZymes supported by dbCAN consensus Recommend Results; {fmt_int(dbcan_recommended_gene_count(sample))} recommended genes",
            ("Class", "Description", "Consensus genes"),
            [(key, class_names[key], classes[key]) for key in class_names],
            (70, 340, 105),
            row_height=24,
        )
        prebiotic_matches = probiotic_cazy_matches(sample)
        if prebiotic_matches:
            draw_pdf_table(
                pdf, f"dbCAN Prebiotic Carbohydrate Utilization: {sample}",
                "Detected CAZy families matched against data/databases/probiotic_cazy.txt",
                ("Prebiotic substrate / group", "Families", "Bakta locus tags", "Typical enzymes / functions"),
                prebiotic_matches,
                (155, 70, 130, 160),
                row_height=24,
            )

        cgc_table = Path(f"data/results/{sample}/dbcan/cgc/cgc_with_substrates.tsv")
        if not cgc_table.exists():
            cgc_table = Path(f"data/results/{sample}/dbcan/cgc_with_substrates.tsv")
        if cgc_table.exists():
            cgc_rows = read_tsv_rows(cgc_table)
            if cgc_rows:
                draw_pdf_table(
                    pdf, f"dbCAN CGC Substrate Predictions: {sample}",
                    "CAZyme Gene Clusters with predicted substrates from dbCAN-PUL",
                    ("CGC#", "Location", "Length", "CAZymes", "PULID", "Substrate", "Bitscore"),
                    [(r.get("CGC#", ""), f"{r.get('Contig','')}:{r.get('Start','')}-{r.get('End','')}",
                      f"{r.get('Length','')} bp", r.get("CAZymes", ""),
                      r.get("PULID", ""), r.get("Substrate", ""), r.get("Bitscore", ""))
                     for r in cgc_rows],
                    (40, 95, 55, 55, 65, 120, 55),
                    row_height=18,
                )
            # 2-column CGC images
            pngs = sorted(Path(f"data/results/{sample}/dbcan/cgc/synteny_png").glob("*.png"))
            if not pngs:
                pngs = sorted(Path(f"data/results/{sample}/dbcan/synteny_png").glob("*.png"))
            if pngs:
                pdf.add_page()
                pdf.text(40, 800, f"dbCAN CGC Synteny Maps: {sample}", size=15, color=(0.05, 0.16, 0.28), bold=True)
                img_w = 470
                yy = 770
                for png in pngs:
                    image = load_png_for_pdf(png)
                    if image:
                        img_h = img_w * image["height"] / image["width"]
                        if yy - img_h < 40:
                            pdf.add_page()
                            yy = 770
                        pdf.chart_image(62, yy - img_h, img_w, img_h, image)
                        yy -= img_h + 5

        pg_summary = read_tsv_rows(f"data/results/{sample}/pgmdb/summary.tsv")
        pg_hits = read_tsv_rows(f"data/results/{sample}/pgmdb/hits.tsv")
        draw_pdf_table(
            pdf, f"PGMDB Probiotic Gene Summary: {sample}",
            "Curated probiotic-potential marker matches; presence does not demonstrate activity",
            ("Probiotic function", "Hits", "Unique markers"),
            [(row.get("metric", ""), row.get("value", ""), row.get("unique_markers", "")) for row in pg_summary],
            (330, 85, 100),
            row_height=20,
        )
        pgmdb_potential = load_png_for_pdf(
            Path(f"data/results/{sample}/pgmdb/pgmdb_functional_potential.png")
        )
        if pgmdb_potential:
            pdf.add_page()
            title = f"PGMDB Probiotic Gene Coverage: {sample}"
            accent, _, _ = pdf_palette(title)
            pdf.text(40, 800, title, size=15, color=accent, bold=True)
            pdf.text(
                40, 780,
                "Only categories represented by at least 3 unique markers in PGMDB are shown",
                size=8, color=(0.38, 0.45, 0.55),
            )
            width = 515
            height = width * pgmdb_potential["height"] / pgmdb_potential["width"]
            if height > 720:
                scale = 720 / height
                width *= scale
                height *= scale
            pdf.chart_image(40 + (515 - width) / 2, 750 - height, width, height, pgmdb_potential)
        grouped = {}
        for row in pg_hits:
            key = (row.get("function", "Other"), row.get("marker") or row.get("gene", ""))
            grouped.setdefault(key, []).append(row)
        gene_rows = []
        prev_func = None
        for (function, gene), gene_rows_hits in sorted(grouped.items()):
            loci = "\n".join(sorted(set(r.get("locus_tag", "") for r in gene_rows_hits if r.get("locus_tag"))))
            evidence = "\n".join(sorted({
                r.get("product", "")
                for r in gene_rows_hits
                if r.get("product")
            }))
            gene_rows.append((
                function if function != prev_func else "",
                gene,
                loci,
                evidence,
            ))
            prev_func = function
        draw_pdf_table(
            pdf, f"PGMDB Probiotic Genes: {sample}",
            "PGMDB genes searched first in Bakta, then in eggNOG when absent from Bakta",
            ("Probiotic function", "PGMDB gene / aliases", "Bakta locus tags", "Annotation"),
            gene_rows,
            (135, 95, 105, 180),
            row_height=18,
            )
        # sankey diagram disabled

        gba_path = Path(f"data/results/{sample}/gut_brain_markers.tsv")
        if gba_path.exists():
            gba_rows = grouped_gut_brain_rows(sample)
            if gba_rows:
                draw_pdf_table(
                    pdf, f"Gut-Brain Axis Markers: {sample}",
                    "Detected genes grouped by Metabolite/Trait; repeated logic and interpretation shown once",
                    ("Metabolite/Trait", "Markers", "Bakta loci", "KO/EC hints", "Marker logic", "Interpretation"),
                    [(r["metabolite_or_trait"], r["marker"], r["bakta_locus"],
                      r["ko_or_ec_hints"], r["marker_logic"], r["interpretation"]) for r in gba_rows],
                    (75, 70, 105, 65, 90, 110),
                    row_height=16,
                )

    for summary in summaries:
        detected_functions = read_keggdecoder_functions(summary["sample"])
        categorized = {}
        for name, value in detected_functions:
            categorized.setdefault(keggdecoder_category(name), []).append((name, value))
        pdf.add_page()
        pdf.text(40, 800, "KEGGDecoder Functional Completeness", size=15, color=(0.05, 0.16, 0.28), bold=True)
        pdf.text(40, 780, f"Sample {summary['sample']} - all {len(detected_functions)} detected functions; completeness from 0 to 100%", size=9, color=(0.38, 0.45, 0.55))
        yy = 748
        relevant = sorted(probiotic_relevant_keggdecoder(detected_functions), key=lambda item: (-item[1], item[0].lower()))
        pdf.text(40, yy, "Probiotic-relevant detected functions", size=11, color=(0.05, 0.16, 0.28), bold=True)
        yy -= 20
        for name, value in relevant:
            if yy < 55:
                pdf.add_page()
                yy = 795
            pdf.text(48, yy + 3, name[:72], size=8, color=(0.16, 0.21, 0.28))
            pdf.rect(420, yy, 100, 12, fill=(0.94, 0.96, 0.98), stroke=(0.82, 0.86, 0.9), line_width=0.4)
            pdf.rect(420, yy, max(1, 100 * value), 12, fill=keggdecoder_heatmap_color(value))
            pdf.text(526, yy + 3, f"{value * 100:.0f}%", size=8, color=(0.16, 0.21, 0.28), bold=True)
            yy -= 18
        yy -= 12
        for category, functions in categorized.items():
            functions.sort(key=lambda item: (-item[1], item[0].lower()))
            required_height = 26 + len(functions) * 18
            if yy - required_height < 45:
                pdf.add_page()
                yy = 795
            pdf.text(40, yy, category, size=11, color=(0.05, 0.16, 0.28), bold=True)
            yy -= 20
            for name, value in functions:
                pdf.text(48, yy + 3, name[:72], size=8, color=(0.16, 0.21, 0.28))
                pdf.rect(420, yy, 100, 12, fill=(0.94, 0.96, 0.98), stroke=(0.82, 0.86, 0.9), line_width=0.4)
                pdf.rect(420, yy, max(1, 100 * value), 12, fill=keggdecoder_heatmap_color(value))
                pdf.text(526, yy + 3, f"{value * 100:.0f}%", size=8, color=(0.16, 0.21, 0.28), bold=True)
                yy -= 18
            yy -= 12

    # KEGGaNOG charts at the end
    for summary in summaries:
        kegganog_png = f"data/results/{summary['sample']}/kegganog/{summary['sample']}_pathways_filtered.png"
        if Path(kegganog_png).exists():
            pdf.add_page()
            pdf.text(40, 800, "KEGGaNOG Pathway Completeness (from eggNOG)", size=15, color=(0.05, 0.16, 0.28), bold=True)
            pdf.text(40, 780, f"Sample {summary['sample']} — non-zero KEGG pathways from eggNOG-mapper annotations", size=9, color=(0.38, 0.45, 0.55))
            image = load_png_for_pdf(kegganog_png)
            if image:
                width = 515
                height = width * image["height"] / image["width"]
                available = 750
                if height > available:
                    scale = available / height
                    width *= scale
                    height *= scale
                pdf.chart_image(40, 760 - height, width, height, image)

    # Transport table at the end
    for summary in summaries:
        tsv_path = Path(f"data/results/{summary['sample']}/eggnog/{summary['sample']}_transport.tsv")
        if tsv_path.exists():
            rows = read_tsv_rows(tsv_path)
            draw_pdf_table(
                pdf, f"Nutrient Uptake Potential (KEGG TC) - {summary['sample']}",
                "Transporter classification from eggNOG-mapper; probiotic relevance annotation",
                ("Category", "TC class", "Transporter family", "Probiotic relevance", "Unique TCs"),
                [(row.get("Category", ""), row.get("TC class", ""), row.get("Transporter family", ""),
                  row.get("Probiotic relevance", ""), row.get("Unique TCs", "")) for row in rows],
                (105, 55, 140, 160, 55),
                row_height=20,
            )
        ec_path = Path(f"data/results/{summary['sample']}/eggnog/{summary['sample']}_ec.tsv")
        if ec_path.exists():
            rows = read_tsv_rows(ec_path)
            draw_pdf_table(
                pdf, f"EC Enzyme Classification - {summary['sample']}",
                "Enzyme Commission classes from eggNOG-mapper annotations",
                ("EC class", "Subcategory", "EC prefix", "Description", "Unique ECs"),
                [(row.get("EC class", ""), row.get("Subcategory", ""), row.get("EC prefix", ""),
                  row.get("Description", ""), row.get("Unique ECs", "")) for row in rows],
                (55, 105, 80, 170, 55),
                row_height=18,
            )
        cog_png = Path(f"data/results/{summary['sample']}/eggnog/{summary['sample']}_cog_radar.png")
        if cog_png.exists():
            pdf.add_page()
            pdf.text(40, 800, "COG Functional Categories", size=15, color=(0.05, 0.16, 0.28), bold=True)
            pdf.text(40, 780, f"Sample {summary['sample']} - COG category distribution from eggNOG-mapper", size=9, color=(0.38, 0.45, 0.55))
            image = load_png_for_pdf(cog_png)
            if image:
                width = 440
                height = width * image["height"] / image["width"]
                if height > 740:
                    scale = 740 / height
                    width *= scale
                    height *= scale
                pdf.chart_image(77, 740 - height, width, height, image)

    # GenomeSPOT growth predictions at the very end
    for summary in summaries:
        pred_path = Path(f"data/results/{summary['sample']}/genomespot/{summary['sample']}.predictions.tsv")
        if pred_path.exists():
            preds = read_tsv_rows(pred_path)
            if preds:
                draw_pdf_table(
                    pdf, f"GenomeSPOT Growth Condition Predictions: {summary['sample']}",
                    "Predicted oxygen tolerance, temperature, salinity and pH from genome features",
                    ("Condition", "Predicted value", "Error", "Units"),
                    [(r.get("target", "").replace("_", " ").title(),
                      r.get("value", "") if r.get("target") == "oxygen" else f"{float(r.get('value', 0)):.2f}" if r.get("value") else "",
                      r.get("error", "") if r.get("target") == "oxygen" else f"{float(r.get('error', 0)):.2f}" if r.get("error") else "",
                      r.get("units", "")) for r in preds],
                    (120, 110, 110, 100),
                    row_height=20,
                )

        gutsmash_rows, gutsmash_svgs = parse_gutsmash_index_html(summary["sample"])
        if gutsmash_rows:
            table_y = draw_pdf_table(
                pdf, f"gutSMASH Metabolic Clusters: {summary['sample']}",
                "Identified primary metabolite regions from gutSMASH index.html",
                ("Region", "Type", "Class", "From", "To", "Most similar known cluster", "Similarity"),
                [
                    (r["region"], r["type"], r["class"], r["from"], r["to"],
                     f"{r['cluster']} ({r['cluster_type']})" if r["cluster_type"] else r["cluster"],
                     r["similarity"])
                    for r in gutsmash_rows
                ],
                (65, 90, 50, 55, 55, 130, 55),
                row_height=20,
            )
            known_cluster = svg_to_png(f"data/results/{summary['sample']}/gutsmash/svg/knownclusterblast_r10c1_all.svg")
            image = load_png_for_pdf(known_cluster) if known_cluster else None
            if image:
                width = 515
                height = width * image["height"] / image["width"]
                available = max(70, table_y - 48)
                if height > available:
                    scale = available / height
                    width *= scale
                    height *= scale
                pdf.chart_image(40, max(36, table_y - height - 18), width, height, image)
        else:
            gutsmash_regions = sorted(Path(f"data/results/{summary['sample']}/gutsmash").glob("*.region*.gbk"))
            if gutsmash_regions:
                draw_pdf_table(
                    pdf, f"gutSMASH Metabolic Clusters: {summary['sample']}",
                    "Gut-associated metabolic regions detected by gutSMASH (fallback to GenBank files)",
                    ("Detected region", "Result directory"),
                    [(path.name, str(path.parent)) for path in gutsmash_regions],
                    (180, 335),
                    row_height=20,
                )
            else:
                gutsmash_status = gutsmash_status_record(summary["sample"])
                if gutsmash_status and gutsmash_status.get("status") != "completed":
                    error = gutsmash_status.get("error", "").strip() or "runtime unavailable"
                    pdf.add_page()
                    pdf.text(40, 800, f"gutSMASH Metabolic Clusters: {summary['sample']}", size=15, color=(0.05, 0.16, 0.28), bold=True)
                    pdf.text(40, 776, "Unavailable in this run", size=11, color=(0.57, 0.25, 0.05), bold=True)
                    pdf.wrapped_text(40, 752, error, max_chars=90, size=9, color=(0.35, 0.39, 0.45), leading=14)

        goprobioticprofiler_dir = Path(f"data/results/{summary['sample']}/goprobioticprofiler")
        if goprobioticprofiler_dir.exists():
            go_enrichment_path = goprobioticprofiler_dir / "go_enrichment.tsv"
            if go_enrichment_path.exists():
                go_rows = []
                with go_enrichment_path.open() as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    for row in reader:
                        fdr = float(row.get("FDR", 1))
                        if fdr < 0.05:
                            go_rows.append((row["GO_ID"], row["GO_Name"], row.get("Namespace", ""), row.get("Fold_Enrichment", ""), f"{fdr:.4f}"))
                if go_rows:
                    draw_pdf_table(
                        pdf, f"GO Enrichment: {summary['sample']}",
                        "Significantly enriched GO terms (FDR < 0.05)",
                        ("GO ID", "Name", "Namespace", "Fold Enrichment", "FDR"),
                        go_rows,
                        (75, 200, 80, 70, 70),
                        row_height=18,
                    )
            for image_name, title, subtitle in [
                (
                    "general_go_bubble.png",
                    f"Genome GO Functional Landscape: {summary['sample']}",
                    "Bubble size represents gene count; position represents genome coverage; color represents GO branch",
                ),
                (
                    "probiotic_go_bubble.png",
                    f"Probiotic GO Categories: {summary['sample']}",
                    "Category-level analysis: position is fold enrichment, bubble size is unique PGMDB-supported genes, and color is evidence tier",
                ),
            ]:
                image = load_png_for_pdf(goprobioticprofiler_dir / image_name)
                if image:
                    pdf.add_page()
                    pdf.text(40, 800, title, size=15, color=(0.05, 0.16, 0.28), bold=True)
                    pdf.text(40, 780, subtitle, size=8, color=(0.38, 0.45, 0.55))
                    width = 515
                    height = width * image["height"] / image["width"]
                    if height > 720:
                        scale = 720 / height
                        width *= scale
                        height *= scale
                    pdf.chart_image(40, 750 - height, width, height, image)

    add_results_tree_page(pdf, summaries)
    pdf.finish()
    add_toc_and_page_numbers(pdf)
    write_pdf_objects(path, pdf.pages, pdf.images)


def draw_workflow_pdf(pdf, x0, y0):
    colors = {
        "data": (0.89, 0.92, 0.95), "qc": (0.86, 0.92, 1.0),
        "assembly": (0.86, 0.98, 0.91), "taxonomy": (0.81, 0.97, 0.98),
        "annotation": (1.0, 0.95, 0.78), "safety": (1.0, 0.89, 0.89),
        "functional": (0.93, 0.91, 1.0), "probiotic": (0.99, 0.90, 0.95),
        "report": (0.12, 0.15, 0.32),
    }
    nodes = {
        "FASTQ": (40, 500, "RAW FASTQ", "input reads", "data"),
        "FASTP": (145, 500, "FASTP", "adapter + QC", "qc"),
        "RASUSA": (250, 500, "RASUSA", "subsampling", "qc"),
        "UNICYCLER": (355, 500, "UNICYCLER", "assembly", "assembly"),
        "FASTQC": (145, 450, "FASTQC", "raw-read QC", "qc"),
        "QUAST": (40, 400, "QUAST", "assembly QC", "qc"),
        "BUSCO": (145, 400, "BUSCO", "completeness", "qc"),
        "MASH": (250, 400, "MASH", "RefSeq ID", "taxonomy"),
        "SOURMASH": (355, 400, "SOURMASH", "GTDB ID", "taxonomy"),
        "BAKTA": (460, 400, "BAKTA LIGHT", "annotation", "annotation"),
        "AMR": (40, 350, "AMRFINDERPLUS", "AMR genes", "safety"),
        "VFDB": (145, 350, "VFDB", "virulence", "safety"),
        "PLASMID": (250, 350, "PLASMIDFINDER", "replicons", "safety"),
        "ANTISMASH": (355, 350, "ANTISMASH", "BGCs", "functional"),
        "KOFAM": (40, 300, "KOFAMSCAN", "KO functions", "functional"),
        "KEGG": (145, 300, "KEGGDECODER", "pathways", "functional"),
        "DBCAN": (250, 300, "DBCAN", "CAZymes", "functional"),
        "CAZY": (355, 300, "PROBIOTIC_CAZY", "prebiotic use", "probiotic"),
        "PGMDB": (460, 300, "PGMDB", "probiotic genes", "probiotic"),
        "REPORT": (460, 500, "FINAL REPORT", "PDF + HTML", "report"),
    }
    edges = [
        ("FASTQ", "FASTP"), ("FASTQ", "FASTQC"), ("FASTP", "RASUSA"), ("RASUSA", "UNICYCLER"),
        ("UNICYCLER", "QUAST"), ("UNICYCLER", "BUSCO"), ("UNICYCLER", "MASH"),
        ("UNICYCLER", "SOURMASH"), ("UNICYCLER", "BAKTA"), ("UNICYCLER", "AMR"),
        ("UNICYCLER", "VFDB"), ("UNICYCLER", "PLASMID"), ("UNICYCLER", "ANTISMASH"),
        ("BAKTA", "KOFAM"), ("KOFAM", "KEGG"), ("BAKTA", "DBCAN"), ("DBCAN", "CAZY"),
        ("BAKTA", "PGMDB"), ("KEGG", "REPORT"), ("CAZY", "REPORT"), ("PGMDB", "REPORT"),
    ]
    box_w, box_h = 86, 30
    for source, target in edges:
        sx, sy = nodes[source][0], nodes[source][1]
        tx, ty = nodes[target][0], nodes[target][1]
        pdf.line(sx + box_w / 2, sy + box_h / 2, tx + box_w / 2, ty + box_h / 2, color=(0.45, 0.5, 0.58), width=0.45)
    for x, y, name, desc, kind in nodes.values():
        fill = colors[kind]
        text_color = (1, 1, 1) if kind == "report" else (0.05, 0.16, 0.28)
        pdf.rect(x, y, box_w, box_h, fill=fill, stroke=(0.2, 0.25, 0.33), line_width=0.65)
        pdf.text(x + 4, y + 18, name[:18], size=5.5, color=text_color, bold=True)
        pdf.text(x + 4, y + 7, desc[:20], size=4.5, color=text_color)


def load_png_for_pdf(path, background=(1, 1, 1)):
    png = Path(path)
    if not png.exists():
        return None
    data = png.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    pos = 8
    width = height = color_type = bit_depth = None
    compressed = []
    while pos < len(data):
        length = struct.unpack(">I", data[pos:pos + 4])[0]
        chunk_type = data[pos + 4:pos + 8]
        chunk_data = data[pos + 8:pos + 8 + length]
        pos += 12 + length
        if chunk_type == b"IHDR":
            width, height, bit_depth, color_type, compression, filter_method, interlace = struct.unpack(">IIBBBBB", chunk_data)
            if bit_depth != 8 or color_type not in {2, 6} or compression != 0 or filter_method != 0 or interlace != 0:
                return None
        elif chunk_type == b"IDAT":
            compressed.append(chunk_data)
        elif chunk_type == b"IEND":
            break
    if width is None or height is None:
        return None

    channels = 4 if color_type == 6 else 3
    stride = width * channels
    raw = zlib.decompress(b"".join(compressed))
    rows = []
    idx = 0
    previous = bytearray(stride)
    for _ in range(height):
        filter_type = raw[idx]
        idx += 1
        row = bytearray(raw[idx:idx + stride])
        idx += stride
        for i in range(stride):
            left = row[i - channels] if i >= channels else 0
            up = previous[i]
            up_left = previous[i - channels] if i >= channels else 0
            if filter_type == 1:
                row[i] = (row[i] + left) & 0xFF
            elif filter_type == 2:
                row[i] = (row[i] + up) & 0xFF
            elif filter_type == 3:
                row[i] = (row[i] + ((left + up) // 2)) & 0xFF
            elif filter_type == 4:
                row[i] = (row[i] + paeth(left, up, up_left)) & 0xFF
            elif filter_type != 0:
                return None
        rows.append(row)
        previous = row

    bg = tuple(int(max(0, min(1, c)) * 255) for c in background)
    rgb = bytearray()
    for row in rows:
        for i in range(0, len(row), channels):
            r, g, b = row[i], row[i + 1], row[i + 2]
            if channels == 4:
                a = row[i + 3] / 255
                r = int(r * a + bg[0] * (1 - a))
                g = int(g * a + bg[1] * (1 - a))
                b = int(b * a + bg[2] * (1 - a))
            rgb.extend((r, g, b))
    return {"width": width, "height": height, "data": bytes(rgb)}


def paeth(a, b, c):
    p = a + b - c
    pa = abs(p - a)
    pb = abs(p - b)
    pc = abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    if pb <= pc:
        return b
    return c


def write_pdf_objects(path, page_streams, images=None):
    images = images or []
    objects = []
    objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    kids = []
    first_page_obj = 3
    for i in range(len(page_streams)):
        kids.append(f"{first_page_obj + i * 2} 0 R")
    objects.append(f"<< /Type /Pages /Kids [{' '.join(kids)}] /Count {len(page_streams)} >>".encode())
    font_regular_obj = first_page_obj + len(page_streams) * 2
    font_bold_obj = font_regular_obj + 1
    image_first_obj = font_bold_obj + 1
    image_resource = ""
    if images:
        image_resource = " /XObject << " + " ".join(f"/{name} {image_first_obj + idx} 0 R" for idx, (name, _) in enumerate(images)) + " >>"
    for i, stream_text in enumerate(page_streams):
        page_obj = first_page_obj + i * 2
        content_obj = page_obj + 1
        resources = f"<< /Font << /F1 {font_regular_obj} 0 R /F2 {font_bold_obj} 0 R >>{image_resource} >>"
        objects.append(f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Resources {resources} /Contents {content_obj} 0 R >>".encode())
        stream = stream_text.encode("latin-1", "replace")
        objects.append(b"<< /Length %d >>\nstream\n" % len(stream) + stream + b"\nendstream")
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>")
    for _, image in images:
        stream = zlib.compress(image["data"])
        objects.append(
            f"<< /Type /XObject /Subtype /Image /Width {image['width']} /Height {image['height']} "
            f"/ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length {len(stream)} >>\nstream\n".encode()
            + stream
            + b"\nendstream"
        )

    out = [b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n"]
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
    selected_samples = set(sys.argv[5:])
    if selected_samples:
        samples = [sample for sample in samples if sample["sample"] in selected_samples]
        missing = selected_samples - {sample["sample"] for sample in samples}
        if missing:
            raise SystemExit(f"Unknown sample(s): {', '.join(sorted(missing))}")
    title = config.get("params", {}).get("report_title", "PRObiotic Bacterial gEnome Analysis Toolkit (PROBEAT)")
    params = config.get("params", {})
    now = dt.datetime.now().isoformat(timespec="seconds")

    rows = []
    summaries = []
    references = load_species_references(config.get("databases", {}).get("ncbi_refseq_species", ""))
    for sample in samples:
        name = sample["sample"]
        checks = {
            "trimmed reads": f"data/intermediate/trimmed/{name}.done",
            "assembly": (
                f"data/results/{name}/medaka/consensus.fasta"
                if sample.get("layout") == "ONT"
                else f"data/results/{name}/assembly/assembly.fasta"
            ),
            "QUAST": f"data/results/{name}/quast/report.tsv",
            "Bakta": f"data/results/{name}/bakta/{name}.gff3",
            "AMRFinderPlus": f"data/results/{name}/amrfinder/{name}.tsv",
            "antiSMASH": f"data/results/{name}/antismash/index.html",
        }
        if sample.get("layout") == "ONT":
            checks["Porechop"] = f"data/intermediate/trimmed/{name}_ONT.fastq.gz"
            checks["Flye"] = f"data/results/{name}/flye/assembly.fasta"
            checks["Medaka"] = f"data/results/{name}/medaka/consensus.fasta"
        if params.get("run_mash", False):
            checks["Mash"] = f"data/results/{name}/mash/top_hits.tsv"
        if params.get("run_sourmash", False):
            checks["sourmash"] = f"data/results/{name}/sourmash/gather.with-lineages.csv"
        if params.get("run_busco", False):
            checks["BUSCO"] = f"data/results/{name}/busco/busco_summary.txt"
        if params.get("run_plasmidfinder", False):
            checks["PlasmidFinder"] = f"data/results/{name}/plasmidfinder/{name}.tsv"
        if params.get("run_vfdb", False):
            checks["VFDB"] = f"data/results/{name}/vfdb/{name}.tsv"
        if params.get("run_kofamscan", False):
            checks["KofamScan"] = f"data/results/{name}/kofamscan/{name}.tsv"
        checks["Rasusa"] = f"data/results/{name}/qc/rasusa/{name}_{sample.get('layout', '')}.tsv"
        if params.get("run_keggdecoder", False):
            checks["KEGGDecoder"] = f"data/results/{name}/keggdecoder/pathway_completeness.tsv"
        if params.get("run_dbcan", False):
            checks["dbCAN"] = f"data/results/{name}/dbcan/overview.tsv"
        if params.get("run_pgmdb", False):
            checks["PGMDB probiotic markers"] = f"data/results/{name}/pgmdb/summary.tsv"
        checks["ABRicate CARD"] = f"data/results/{name}/abricate/{name}_card.tsv"
        checks["ABRicate MEGARes"] = f"data/results/{name}/abricate/{name}_megares.tsv"
        if params.get("run_genomad", False):
            checks["geNomad plasmids"] = f"data/results/{name}/genomad/assembly_summary/assembly_plasmid_summary.tsv"
            checks["geNomad viruses"] = f"data/results/{name}/genomad/assembly_summary/assembly_virus_summary.tsv"
        if params.get("run_eggnog", False):
            checks["eggNOG-mapper"] = f"data/results/{name}/eggnog/{name}.emapper.annotations"
        if params.get("run_gtdbtk", False):
            checks["GTDB-Tk"] = f"data/results/{name}/gtdbtk/gtdbtk.bac120.summary.tsv"
        if params.get("run_mobsuite", False):
            checks["MOB-recon"] = f"data/results/{name}/mobsuite/contig_report.txt"
        if params.get("run_kegganog", False):
            checks["KEGGaNOG"] = f"data/results/{name}/kegganog/{name}_pathways.tsv"
        if params.get("run_genomespot", False):
            checks["GenomeSPOT"] = f"data/results/{name}/genomespot/{name}.predictions.tsv"
        if params.get("run_gutsmash", False):
            checks["gutSMASH"] = f"data/results/{name}/gutsmash/.probeat_gutsmash_status.json"
        if params.get("run_fastmlst", False):
            checks["FastMLST"] = f"data/results/{name}/fastmlst/{name}_mlst.tsv"
        if params.get("run_fastani", False):
            checks["FastANI"] = f"data/results/{name}/fastani/{name}_ani_heatmap.png"
        if params.get("run_goprobioticprofiler", False):
            checks["GO Probiotic"] = f"data/results/{name}/goprobioticprofiler/probiotic_go_summary.json"
        if params.get("run_pycirclize", False):
            checks["pyCirclize genome map"] = f"data/results/{name}/genome_map/{name}_genome_map.png"
        rows.append((name, sample.get("layout", ""), checks))
        summary = sample_summary(sample)
        summary["reference"] = reference_for_summary(summary, references)
        summary["qc_assessment"] = qc_assessment(summary, summary["reference"])
        summaries.append(summary)

    css = (
        "body{font-family:Arial,sans-serif;margin:0;background:#f8fafc;color:#0f172a;line-height:1.45}"
        ".hero{background:#1f2652;color:white;padding:28px 42px;display:flex;justify-content:space-between;gap:24px;align-items:center}.hero p{color:#cfe5ff;margin:6px 0 0}"
        ".brand-logo{width:180px;max-width:30vw;height:auto;display:block}"
        "main{padding:28px 42px}.panel{background:white;border:1px solid #d9e2ec;padding:18px;margin:0 0 22px}"
        "h1{margin:0;font-size:30px}h2{margin:0 0 12px;font-size:20px;color:#0f2742}"
        "table{border-collapse:collapse;width:100%;margin:12px 0 24px;background:white;table-layout:fixed}th,td{border:1px solid #d9e2ec;padding:7px;text-align:left;font-size:13px;overflow-wrap:anywhere;word-break:break-word}"
        "th{background:#eef4fb;color:#0f2742}.done{color:#067a2f;font-weight:700}.missing{color:#9a3412;font-weight:700}.soft-failed{color:#92400e;font-weight:700}.note{color:#64748b}"
        ".workflow{width:100%;max-width:1040px;border:1px solid #d9e2ec;background:#f8fafc}.cards{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin:14px 0 22px}"
        ".card{background:white;border:1px solid #d9e2ec;padding:13px}.card b{display:block;font-size:20px;color:#0f2742}.card span{color:#64748b;font-size:12px}"
        ".plotly-sankey{width:100%;min-height:520px;border:1px solid #d9e2ec;background:white}"
    )
    html_dir = Path(html_out).parent
    logo_src = os.path.relpath("data/logo.png", html_dir)
    workflow_src = os.path.relpath("data/workflow.png", html_dir)
    parts = [
        f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title>"
        "<script src='https://cdn.plot.ly/plotly-2.35.2.min.js'></script>"
        f"<style>{css}</style></head><body>"
    ]
    parts.append(
        f"<section class='hero'><div><h1>{html.escape(title)}</h1><p>Generated: {html.escape(now)}</p></div>"
        f"<img class='brand-logo' src='{html.escape(logo_src)}' alt='Paisii Hilendarski University of Plovdiv logo'>"
        "</section><main>"
    )
    if summaries:
        first = summaries[0]
        fastp = first["fastp"]
        quast = first["quast"]
        cards = [
            ("Raw reads", fmt_int(fastp.get("raw_reads"))),
            ("After QC reads", fmt_int(fastp.get("filtered_reads"))),
            ("Contigs (>= 500 bp)", fmt_int(quast.get("contigs"))),
            ("N50", fmt_int(quast.get("n50"))),
        ]
        parts.append("<section class='cards'>")
        for label, value in cards:
            parts.append(f"<div class='card'><b>{html.escape(value)}</b><span>{html.escape(label)}</span></div>")
        parts.append("</section>")
    parts.append("<section class='panel'><h2>Workflow</h2>")
    parts.append(f"<img class='workflow' src='{html.escape(workflow_src)}' alt='Workflow diagram'>")
    parts.append("</section>")
    parts.append("<section class='panel'><h2>Run Summary</h2>")
    parts.append("<table><tr><th>Sample</th><th>Raw reads</th><th>Raw bases</th><th>After QC reads</th><th>After QC bases</th><th>Q30 after QC</th><th>Adapter-trimmed reads</th><th>Contigs</th><th>Genome size</th><th>GC</th><th>N50</th><th>Mash top hit</th><th>sourmash</th><th>Bakta CDS</th><th>AMR hits</th><th>Plasmids</th><th>VFDB</th><th>KOfam</th><th>KEGGDecoder</th><th>dbCAN</th><th>PGMDB</th><th>BUSCO</th></tr>")
    for summary in summaries:
        fastp = summary["fastp"]
        quast = summary["quast"]
        busco = summary["busco"]
        mash = summary["mash"]
        sourmash = summary["sourmash"]
        pgmdb = summary["pgmdb"]
        rasusa = summary["rasusa"]
        busco_text = busco.get("summary", busco.get("status", "not run"))
        parts.append(
            "<tr>"
            f"<td>{html.escape(summary['sample'])}</td>"
            f"<td>{fmt_int(fastp.get('raw_reads'))}</td>"
            f"<td>{fmt_int(fastp.get('raw_bases'))}</td>"
            f"<td>{fmt_int(fastp.get('filtered_reads'))}</td>"
            f"<td>{fmt_int(fastp.get('filtered_bases'))}</td>"
            f"<td>{fmt_percent(fastp.get('filtered_q30'))}</td>"
            f"<td>{fmt_int(fastp.get('adapter_trimmed_reads'))}</td>"
            f"<td>{fmt_int(quast.get('contigs'))}</td>"
            f"<td>{fmt_int(quast.get('total_length'))}</td>"
            f"<td>{fmt_percent(quast.get('gc'))}</td>"
            f"<td>{fmt_int(quast.get('n50'))}</td>"
            f"<td>{html.escape(mash.get('summary', mash.get('status', 'not run')))}</td>"
            f"<td>{html.escape(sourmash.get('summary', sourmash.get('status', 'not run')))}</td>"
            f"<td>{fmt_int(summary['bakta'].get('CDS'))}</td>"
            f"<td>{fmt_int(summary['amr_hits'])}</td>"
            f"<td>{fmt_int(summary['plasmid_hits'])}</td>"
            f"<td>{fmt_int(summary['vfdb_hits'])}</td>"
            f"<td>{fmt_int(summary['kofam_hits'])}</td>"
            f"<td>{fmt_int(summary['keggdecoder'][0]) if summary['keggdecoder'] else 'not run'}</td>"
            f"<td>{fmt_int(summary['dbcan_hits'])}</td>"
            f"<td>{fmt_int(pgmdb.get('unique_markers'))}</td>"
            f"<td>{html.escape(busco_text)}</td>"
            "</tr>"
        )
    parts.append("</table>")
    parts.append("</section>")

    for summary in summaries:
        fastp = summary["fastp"]
        quast = summary["quast"]
        bakta = summary["bakta"]
        busco = summary["busco"]
        mash = summary["mash"]
        sourmash = summary["sourmash"]
        pgmdb = summary["pgmdb"]
        # Flye Assembly Report (ONT only)
        if summary.get("layout") == "ONT":
            flye_info_path = Path(f"data/results/{summary['sample']}/flye/assembly_info.txt")
            if flye_info_path.exists():
                flye_rows = []
                for line in flye_info_path.read_text(errors="replace").splitlines():
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    if len(fields) >= 7:
                        flye_rows.append(fields)
                if flye_rows:
                    parts.append(f"<section class='panel'><h2>Flye Assembly Report: {html.escape(summary['sample'])}</h2>")
                    parts.append("<p class='note'>ONT long-read contig statistics from Flye assembler.</p>")
                    parts.append("<table><tr><th>Contig</th><th>Length</th><th>Coverage</th><th>Circular</th><th>Repeat</th><th>Multiplicity</th><th>Alt group</th></tr>")
                    for fields in flye_rows:
                        circ = "Yes" if fields[3] == "Y" else "No"
                        repeat = "Yes" if fields[4] == "Y" else "No"
                        parts.append(f"<tr><td>{html.escape(fields[0])}</td><td>{fmt_int(fields[1])}</td><td>{html.escape(fields[2])}</td><td>{html.escape(circ)}</td><td>{html.escape(repeat)}</td><td>{html.escape(fields[5])}</td><td>{html.escape(fields[6])}</td></tr>")
                    parts.append("</table></section>")
        parts.append(f"<section class='panel'><h2>Sample {html.escape(summary['sample'])}</h2>")
        genome_map_path = Path(f"data/results/{summary['sample']}/genome_map/{summary['sample']}_genome_map.png")
        if genome_map_path.exists():
            parts.append("<h2>Integrated Bacterial Genome Map</h2>")
            parts.append(
                "<p class='note'>Bakta CDS, ABRicate, antiSMASH, BAGEL4, dbCAN CAZymes, "
                "PGMDB probiotic genes, GC content and GC skew.</p>"
            )
            parts.append(f"<img class='workflow' src='{html.escape(str(genome_map_path))}' alt='Integrated bacterial genome map'>")
        parts.append("<table><tr><th>Category</th><th>Metric</th><th>Value</th></tr>")
        details = [
            ("Raw data", "reads", fmt_int(fastp.get("raw_reads"))),
            ("Raw data", "bases", fmt_int(fastp.get("raw_bases"))),
            ("Raw data", "Q30", fmt_percent(fastp.get("raw_q30"))),
            ("After adapter/QC", "reads", fmt_int(fastp.get("filtered_reads"))),
            ("After adapter/QC", "bases", fmt_int(fastp.get("filtered_bases"))),
            ("After adapter/QC", "Q30", fmt_percent(fastp.get("filtered_q30"))),
            ("After adapter/QC", "low-quality reads removed", fmt_int(fastp.get("low_quality_reads"))),
            ("After adapter/QC", "too-many-N reads removed", fmt_int(fastp.get("too_many_n_reads"))),
            ("After adapter/QC", "adapter-trimmed reads", fmt_int(fastp.get("adapter_trimmed_reads"))),
            ("Rasusa", "action", rasusa.get("status", "not run")),
            ("Rasusa", "reads used for assembly", fmt_int(rasusa.get("output_reads"))),
            ("Assembly", "contigs", fmt_int(quast.get("contigs"))),
            ("Assembly", "largest contig", fmt_int(quast.get("largest_contig"))),
            ("Assembly", "total length", fmt_int(quast.get("total_length"))),
            ("Assembly", "GC", fmt_percent(quast.get("gc"))),
            ("Assembly", "N50", fmt_int(quast.get("n50"))),
            ("Mash", "top species", mash.get("summary", mash.get("status", "not run"))),
            ("Mash", "best distance", mash.get("distance", "NA")),
            ("sourmash", "top species", sourmash.get("summary", sourmash.get("status", "not run"))),
            ("Bakta", "CDS", fmt_int(bakta.get("CDS"))),
            ("Bakta", "tRNA", fmt_int(bakta.get("tRNA"))),
            ("Bakta", "rRNA", fmt_int(bakta.get("rRNA"))),
            ("Bakta", "tmRNA", fmt_int(bakta.get("tmRNA"))),
            ("Bakta", "ncRNA", fmt_int(bakta.get("ncRNA"))),
            ("Bakta", "CRISPR arrays", fmt_int(bakta.get("CRISPR"))),
            ("AMRFinderPlus", "AMR hits", fmt_int(summary["amr_hits"])),
            ("PlasmidFinder", "replicon hits", fmt_int(summary["plasmid_hits"])),
            ("VFDB", "virulence-factor hits", fmt_int(summary["vfdb_hits"])),
            ("KofamScan", "KO annotations", fmt_int(summary["kofam_hits"])),
            ("KEGGDecoder", "functions present", fmt_int(summary["keggdecoder"][0]) if summary["keggdecoder"] else "not run"),
            ("KEGGDecoder", "functions assessed", fmt_int(summary["keggdecoder"][1]) if summary["keggdecoder"] else "not run"),
            ("dbCAN", "CAZyme hits", fmt_int(summary["dbcan_hits"])),
            ("PGMDB", "total marker hits", fmt_int(pgmdb.get("total_hits"))),
            ("PGMDB", "unique curated markers", fmt_int(pgmdb.get("unique_markers"))),
            ("PGMDB", "top functional categories", pgmdb.get("summary", pgmdb.get("status", "not run"))),
            ("BAGEL4", "bacteriocin/RiPP clusters", fmt_int(len(bagel4_rows(summary["sample"])))),
            ("BUSCO", "status", busco.get("status", "not run")),
            ("BUSCO", "complete", fmt_percent(busco.get("complete"))),
            ("BUSCO", "single-copy", fmt_percent(busco.get("single"))),
            ("BUSCO", "duplicated", fmt_percent(busco.get("duplicated"))),
            ("BUSCO", "fragmented", fmt_percent(busco.get("fragmented"))),
            ("BUSCO", "missing", fmt_percent(busco.get("missing"))),
        ]
        for category, metric, value in details:
            parts.append(f"<tr><td>{html.escape(category)}</td><td>{html.escape(metric)}</td><td>{html.escape(str(value))}</td></tr>")
        parts.append("</table>")
        keggdecoder_svg = Path(f"data/results/{summary['sample']}/keggdecoder/pathway_completeness.svg")
        if keggdecoder_svg.exists():
            relative = os.path.relpath(keggdecoder_svg, Path(html_out).parent)
            parts.append("<h2>KEGGDecoder Functional Completeness</h2>")
            parts.append(f"<img class='workflow' src='{html.escape(relative)}' alt='KEGGDecoder pathway completeness heatmap'>")
        reference = summary.get("reference", {})
        parts.append("<h2>Isolate Genome QC Assessment</h2>")
        parts.append(
            f"<p class='note'>Selected NCBI reference: {html.escape(reference.get('organism_name', 'not available'))}; "
            f"assembly {html.escape(reference.get('assembly_accession', 'NA'))}; "
            f"genome size {html.escape(reference.get('genome_size', 'NA'))}; GC {html.escape(reference.get('gc_percent', 'NA'))}%.</p>"
        )
        parts.append("<table><tr><th>Metric</th><th>Good practice</th><th>Observed</th><th>Assessment</th></tr>")
        for metric, good_practice, observed, assessment in summary.get("qc_assessment", []):
            parts.append(
                f"<tr><td>{html.escape(metric)}</td><td>{html.escape(good_practice)}</td>"
                f"<td>{html.escape(observed)}</td><td>{html.escape(assessment_label(assessment))}</td></tr>"
            )
        parts.append("</table>")
        parts.append("<h2>PGMDB Probiotic Marker Categories</h2>")
        parts.append("<p class='note'>Curated genomic potential markers searched first by Bakta gene symbol and then, when absent from Bakta, by eggNOG Preferred_name. Presence alone does not demonstrate probiotic activity.</p>")
        parts.append("<table><tr><th>Functional category</th><th>Matched loci</th><th>Unique curated markers</th></tr>")
        for item in pgmdb.get("categories", []):
            parts.append(
                f"<tr><td>{html.escape(item['category'])}</td><td>{fmt_int(item['hits'])}</td>"
                f"<td>{fmt_int(item['unique_markers'])}</td></tr>"
            )
        if not pgmdb.get("categories"):
            parts.append("<tr><td colspan='3'>No PGMDB result available.</td></tr>")
        parts.append("</table>")
        pgmdb_potential = Path(
            f"data/results/{summary['sample']}/pgmdb/pgmdb_functional_potential.png"
        )
        if pgmdb_potential.exists():
            parts.append("<h2>PGMDB Probiotic Gene Coverage</h2>")
            parts.append(
                "<p class='note'>Only categories represented by at least 3 unique markers "
                "in PGMDB are shown.</p>"
            )
            parts.append(
                f"<img src='{html.escape(str(pgmdb_potential))}' "
                "style='max-width:100%;height:auto' />"
            )
        parts.append("<h2>PGMDB Probiotic Genes</h2>")
        parts.append("<table><tr><th>Probiotic function</th><th>PGMDB gene / aliases</th><th>Bakta locus tags</th><th>Annotation</th></tr>")
        pg_hits = read_tsv_rows(f"data/results/{summary['sample']}/pgmdb/hits.tsv")
        grouped_hits = {}
        for row in pg_hits:
            key = (normalize_pgmdb_function(row.get("function", "Other")), row.get("marker") or row.get("gene", ""))
            grouped_hits.setdefault(key, []).append(row)
        if grouped_hits:
            prev_func = None
            for (function, gene), gene_rows_hits in sorted(grouped_hits.items()):
                loci = sorted(set(r.get("locus_tag", "") for r in gene_rows_hits if r.get("locus_tag")))
                evidence = sorted({
                    r.get("product", "")
                    for r in gene_rows_hits
                    if r.get("product")
                })
                parts.append(
                    f"<tr><td>{html.escape(function if function != prev_func else '')}</td>"
                    f"<td>{html.escape(gene)}</td>"
                    f"<td style='white-space:normal;overflow-wrap:anywhere'>{'<br>'.join(html.escape(value) for value in loci)}</td>"
                    f"<td style='white-space:normal;overflow-wrap:anywhere'>{'<br>'.join(html.escape(value) for value in evidence)}</td></tr>"
                )
                prev_func = function
        else:
            parts.append("<tr><td colspan='4'>No PGMDB gene-level hits available.</td></tr>")
        parts.append("</table>")
        # sankey diagram disabled

        card_path = Path(f"data/results/{summary['sample']}/abricate/{summary['sample']}_card.tsv")
        megares_path = Path(f"data/results/{summary['sample']}/abricate/{summary['sample']}_megares.tsv")
        amrfinder_path = Path(f"data/results/{summary['sample']}/amrfinder/{summary['sample']}.tsv")
        if amrfinder_path.exists():
            amrfinder_rows = read_tsv_rows(amrfinder_path)
            if amrfinder_rows:
                parts.append("<h2>AMRFinderPlus Hits</h2>")
                parts.append("<table><tr><th>Gene</th><th>Element</th><th>Class</th><th>Identity</th></tr>")
                for row in amrfinder_rows:
                    parts.append(
                        "<tr>"
                        f"<td>{html.escape(row.get('Element symbol', ''))}</td>"
                        f"<td>{html.escape(row.get('Element name', ''))}</td>"
                        f"<td>{html.escape(row.get('Class', ''))}</td>"
                        f"<td>{html.escape(row.get('% Identity to reference', ''))}</td>"
                        "</tr>"
                    )
                parts.append("</table>")
        if card_path.exists() and megares_path.exists():
            card_rows = [row for row in read_tsv_rows(card_path) if row.get("SEQUENCE")]
            megares_rows = [row for row in read_tsv_rows(megares_path) if row.get("SEQUENCE")]
            parts.append("<h2>ABRicate CARD &amp; MEGARes</h2>")
            parts.append("<table><tr><th>Database</th><th>Detected hits</th></tr>")
            parts.append(f"<tr><td>CARD</td><td>{len(card_rows)}</td></tr>")
            parts.append(f"<tr><td>MEGARes</td><td>{len(megares_rows)}</td></tr>")
            parts.append("</table>")
            if card_rows:
                parts.append("<h2>ABRicate CARD Hits</h2>")
                parts.append("<table><tr><th>Contig</th><th>Start</th><th>End</th><th>Gene</th><th>%Coverage</th><th>%Identity</th><th>Resistance</th></tr>")
                for row in card_rows:
                    parts.append(
                        "<tr>"
                        f"<td>{html.escape(row.get('SEQUENCE', ''))}</td>"
                        f"<td>{html.escape(row.get('START', ''))}</td>"
                        f"<td>{html.escape(row.get('END', ''))}</td>"
                        f"<td>{html.escape(row.get('GENE', ''))}</td>"
                        f"<td>{html.escape(row.get('%COVERAGE', ''))}</td>"
                        f"<td>{html.escape(row.get('%IDENTITY', ''))}</td>"
                        f"<td>{html.escape(row.get('RESISTANCE', ''))}</td>"
                        "</tr>"
                    )
                parts.append("</table>")
            if megares_rows:
                parts.append("<h2>ABRicate MEGARes Hits</h2>")
                parts.append("<table><tr><th>Contig</th><th>Start</th><th>End</th><th>Gene</th><th>%Coverage</th><th>%Identity</th></tr>")
                for row in megares_rows:
                    parts.append(
                        "<tr>"
                        f"<td>{html.escape(row.get('SEQUENCE', ''))}</td>"
                        f"<td>{html.escape(row.get('START', ''))}</td>"
                        f"<td>{html.escape(row.get('END', ''))}</td>"
                        f"<td>{html.escape(row.get('GENE', ''))}</td>"
                        f"<td>{html.escape(row.get('%COVERAGE', ''))}</td>"
                        f"<td>{html.escape(row.get('%IDENTITY', ''))}</td>"
                        "</tr>"
                    )
                parts.append("</table>")

        gba_path = Path(f"data/results/{summary['sample']}/gut_brain_markers.tsv")
        if gba_path.exists():
            parts.append("<h2>Gut-Brain Axis Markers</h2>")
            parts.append("<p class='note'>Detected genes grouped by Metabolite/Trait; repeated marker logic and interpretation are shown once.</p>")
            parts.append("<table><tr><th>Metabolite/Trait</th><th>Markers</th><th>Bakta loci</th><th>KO/EC hints</th><th>Marker logic</th><th>Interpretation</th></tr>")
            for row in grouped_gut_brain_rows(summary["sample"]):
                cells = {
                    key: html.escape(row.get(key, "")).replace("\n", "<br>")
                    for key in row
                }
                parts.append(
                    f"<tr><td>{cells['metabolite_or_trait']}</td>"
                    f"<td>{cells['marker']}</td>"
                    f"<td>{cells['bakta_locus']}</td>"
                    f"<td>{cells['ko_or_ec_hints']}</td>"
                    f"<td>{cells['marker_logic']}</td>"
                    f"<td>{cells['interpretation']}</td></tr>"
                )
            parts.append("</table>")

        eggnog_path = Path(f"data/results/{summary['sample']}/eggnog/{summary['sample']}.emapper.annotations")
        if eggnog_path.exists():
            annotation_count = sum(1 for line in eggnog_path.open(errors="replace") if line.strip() and not line.startswith("#"))
            parts.append("<h2>eggNOG Functional Annotation</h2>")
            parts.append(f"<p>{annotation_count} proteins received eggNOG-mapper annotations.</p>")
            transport_path = Path(f"data/results/{summary['sample']}/eggnog/{summary['sample']}_transport.tsv")
            if transport_path.exists():
                parts.append("<table><tr><th>Transport category</th><th>TC class</th><th>Family</th><th>Probiotic relevance</th><th>Unique TCs</th></tr>")
                for row in read_tsv_rows(transport_path):
                    parts.append(
                        f"<tr><td>{html.escape(row.get('Category', ''))}</td><td>{html.escape(row.get('TC class', ''))}</td>"
                        f"<td>{html.escape(row.get('Transporter family', ''))}</td><td>{html.escape(row.get('Probiotic relevance', ''))}</td>"
                        f"<td>{html.escape(row.get('Unique TCs', ''))}</td></tr>"
                    )
                parts.append("</table>")
            cog_path = Path(f"data/results/{summary['sample']}/eggnog/{summary['sample']}_cog_radar.png")
            if cog_path.exists():
                relative = os.path.relpath(cog_path, Path(html_out).parent)
                parts.append(f"<img class='workflow' src='{html.escape(relative)}' alt='eggNOG COG functional category radar'>")

        bagel_rows = bagel4_rows(summary["sample"])
        if bagel_rows:
            parts.append("<h2>BAGEL4 Bacteriocin and RiPP Clusters</h2>")
            parts.append("<table><tr><th>AOI</th><th>Start</th><th>End</th><th>Class</th><th>Sequence</th></tr>")
            for row in bagel_rows:
                parts.append(
                    f"<tr><td>{html.escape(row.get('AOI', ''))}</td><td>{html.escape(row.get('start', ''))}</td>"
                    f"<td>{html.escape(row.get('end', ''))}</td><td>{html.escape(row.get('class', ''))}</td>"
                    f"<td>{html.escape(row.get('fullname', ''))}</td></tr>"
                )
            parts.append("</table>")
            for image_path in bagel4_preferred_images(summary["sample"]):
                relative = os.path.relpath(image_path, Path(html_out).parent)
                parts.append(f"<img class='workflow' src='{html.escape(relative)}' alt='BAGEL4 cluster graphic'>")
        antismash_path = Path(f"data/results/{summary['sample']}/antismash/{summary['sample']}_regions.tsv")
        if antismash_path.exists():
            parts.append("<h2>antiSMASH Biosynthetic Gene Clusters</h2>")
            parts.append("<table><tr><th>Region</th><th>Location</th><th>Length</th><th>Products</th><th>Genes</th></tr>")
            for row in read_tsv_rows(antismash_path):
                parts.append(
                    f"<tr><td>{html.escape(row.get('Region', ''))}</td><td>{html.escape(row.get('Location', ''))}</td>"
                    f"<td>{html.escape(row.get('Length', ''))}</td><td>{html.escape(row.get('Products', ''))}</td>"
                    f"<td>{html.escape(row.get('Genes', ''))}</td></tr>"
                )
            parts.append("</table>")
        genomad_plasmid_path = Path(f"data/results/{summary['sample']}/genomad/assembly_summary/assembly_plasmid_summary.tsv")
        genomad_virus_path = Path(f"data/results/{summary['sample']}/genomad/assembly_summary/assembly_virus_summary.tsv")
        if genomad_plasmid_path.exists():
            genomad_plasmid_rows = read_tsv_rows(genomad_plasmid_path)
            parts.append(f"<h2>geNomad Plasmid Summary</h2>")
            parts.append(f"<p>{len(genomad_plasmid_rows)} contig(s) classified as plasmids</p>")
            parts.append("<table><tr><th>Contig</th><th>Length</th><th>Topology</th><th>Plasmid score</th><th>FDR</th><th>Hallmarks</th><th>Marker enrichment</th></tr>")
            for row in genomad_plasmid_rows:
                parts.append(
                    f"<tr><td>{html.escape(row.get('seq_name', ''))}</td><td>{html.escape(row.get('length', ''))}</td>"
                    f"<td>{html.escape(row.get('topology', ''))}</td><td>{html.escape(row.get('plasmid_score', ''))}</td>"
                    f"<td>{html.escape(row.get('fdr', ''))}</td><td>{html.escape(row.get('n_hallmarks', ''))}</td>"
                    f"<td>{html.escape(row.get('marker_enrichment', ''))}</td></tr>"
                )
            parts.append("</table>")
        if genomad_virus_path.exists():
            genomad_virus_rows = read_tsv_rows(genomad_virus_path)
            parts.append("<h2>geNomad Virus Summary</h2>")
            parts.append(f"<p>{len(genomad_virus_rows)} contig(s) classified as viruses</p>")
            parts.append("<table><tr><th>Contig</th><th>Length</th><th>Topology</th><th>Virus score</th><th>FDR</th><th>Hallmarks</th><th>Marker enrichment</th></tr>")
            for row in genomad_virus_rows:
                parts.append(
                    f"<tr><td>{html.escape(row.get('seq_name', ''))}</td><td>{html.escape(row.get('length', ''))}</td>"
                    f"<td>{html.escape(row.get('topology', ''))}</td><td>{html.escape(row.get('virus_score', ''))}</td>"
                    f"<td>{html.escape(row.get('fdr', ''))}</td><td>{html.escape(row.get('n_hallmarks', ''))}</td>"
                    f"<td>{html.escape(row.get('marker_enrichment', ''))}</td></tr>"
                )
            parts.append("</table>")

        fastani_path = Path(f"data/results/{summary['sample']}/fastani/{summary['sample']}_fastani.tsv")
        if fastani_path.exists():
            ani_rows = []
            for line in fastani_path.read_text(errors="replace").splitlines():
                fields = line.split("\t")
                if len(fields) >= 5:
                    ani_rows.append(fields)
            parts.append("<h2>FastANI Species Confirmation</h2>")
            parts.append("<table><tr><th>Query</th><th>Reference</th><th>ANI</th><th>Aligned fragments</th><th>Total fragments</th></tr>")
            for row in ani_rows:
                parts.append(f"<tr><td>{html.escape(Path(row[0]).name)}</td><td>{html.escape(Path(row[1]).name)}</td><td>{html.escape(row[2])}%</td><td>{html.escape(row[3])}</td><td>{html.escape(row[4])}</td></tr>")
            parts.append("</table>")
            if ani_rows:
                best_ani = max(float(row[2]) for row in ani_rows)
                if best_ani > 95:
                    parts.append(
                        f"<p class='done'>PASS: The assembled strain is correctly assigned to the "
                        f"reference species (best FastANI {best_ani:.2f}% &gt; 95%).</p>"
                    )

        mob_path = Path(f"data/results/{summary['sample']}/mobsuite/contig_report.txt")
        if mob_path.exists():
            mob_rows = [
                row for row in read_tsv_rows(mob_path)
                if row.get("molecule_type", "").strip().lower() == "plasmid"
            ]
            parts.append("<h2>MOB-recon Contig Classification</h2>")
            parts.append("<p>Only contigs classified as plasmids are shown.</p>")
            parts.append("<table><tr><th>Contig</th><th>Type</th><th>Length</th><th>Rep types</th><th>Relaxases</th><th>Mobility</th><th>Nearest neighbor</th></tr>")
            for row in mob_rows:
                parts.append(f"<tr><td>{html.escape(row.get('contig_id', ''))}</td><td>{html.escape(row.get('molecule_type', ''))}</td><td>{html.escape(row.get('size', ''))}</td><td>{html.escape(row.get('rep_type(s)', ''))}</td><td>{html.escape(row.get('relaxase_type(s)', ''))}</td><td>{html.escape(row.get('predicted_mobility', ''))}</td><td>{html.escape(row.get('mash_neighbor_identification', ''))}</td></tr>")
            parts.append("</table>")
        plasmid_path = Path(f"data/results/{summary['sample']}/mobsuite/mobtyper_results.txt")
        if plasmid_path.exists():
            plasmid_rows = read_tsv_rows(plasmid_path)
            parts.append("<h2>MOB-recon Reconstructed Plasmids</h2>")
            parts.append("<table><tr><th>Plasmid</th><th>Contigs</th><th>Length</th><th>Rep types</th><th>Mobility</th><th>Nearest neighbor</th><th>Host range</th></tr>")
            for row in plasmid_rows:
                parts.append(f"<tr><td>{html.escape(row.get('sample_id', ''))}</td><td>{html.escape(row.get('num_contigs', ''))}</td><td>{html.escape(row.get('size', ''))}</td><td>{html.escape(row.get('rep_type(s)', ''))}</td><td>{html.escape(row.get('predicted_mobility', ''))}</td><td>{html.escape(row.get('mash_neighbor_identification', ''))}</td><td>{html.escape(row.get('predicted_host_range_overall_name', ''))}</td></tr>")
            parts.append("</table>")

        gutsmash_rows, gutsmash_svgs = parse_gutsmash_index_html(summary["sample"])
        if gutsmash_rows:
            parts.append("<h2>gutSMASH Metabolic Clusters</h2>")
            parts.append("<table><tr><th>Region</th><th>Type</th><th>Class</th><th>From</th><th>To</th><th>Most similar known cluster</th><th>Similarity</th></tr>")
            for r in gutsmash_rows:
                cluster = html.escape(r["cluster"])
                if r["cluster_type"]:
                    cluster += f" ({html.escape(r['cluster_type'])})"
                parts.append(
                    f"<tr><td>{html.escape(r['region'])}</td><td>{html.escape(r['type'])}</td><td>{html.escape(r['class'])}</td>"
                    f"<td>{html.escape(r['from'])}</td><td>{html.escape(r['to'])}</td>"
                    f"<td>{cluster}</td><td>{html.escape(r['similarity'])}</td></tr>"
                )
            parts.append("</table>")
            if gutsmash_svgs:
                parts.append("<h3>gutSMASH ClusterBlast Images</h3>")
                for svg in gutsmash_svgs:
                    svg_content = svg.read_text(errors="replace")
                    parts.append(f"<div class='gutsmash-svg'><h4>{html.escape(svg.name)}</h4>")
                    parts.append(svg_content)
                    parts.append("</div>")
        elif sorted(Path(f"data/results/{summary['sample']}/gutsmash").glob("*.region*.gbk")):
            parts.append("<h2>gutSMASH Metabolic Clusters</h2>")
            _gs_count = len(sorted(Path(f"data/results/{summary['sample']}/gutsmash").glob("*.region*.gbk")))
            parts.append(f"<p>{_gs_count} gut-associated metabolic region(s) detected.</p>")
        else:
            gutsmash_status = gutsmash_status_record(summary["sample"])
            if gutsmash_status and gutsmash_status.get("status") != "completed":
                error = gutsmash_status.get("error", "").strip() or "runtime unavailable"
                parts.append("<h2>gutSMASH Metabolic Clusters</h2>")
                parts.append(f"<p class='note'>Unavailable in this run: {html.escape(error)}</p>")

        goprobioticprofiler_dir = Path(f"data/results/{summary['sample']}/goprobioticprofiler")
        if goprobioticprofiler_dir.exists():
            parts.append("<h2>GO Probiotic</h2>")
            # GO enrichment table
            go_enrichment_path = goprobioticprofiler_dir / "go_enrichment.tsv"
            if go_enrichment_path.exists():
                parts.append("<h3>GO Enrichment</h3>")
                parts.append("<table><tr><th>GO ID</th><th>Name</th><th>Namespace</th><th>Fold Enrichment</th><th>FDR</th></tr>")
                with go_enrichment_path.open() as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    for row in reader:
                        fdr = float(row.get("FDR", 1))
                        fdr_str = f"{fdr:.4f}" if fdr < 0.05 else f"{fdr:.4f}"
                        sig_class = "significant" if fdr < 0.05 else ""
                        parts.append(
                            f"<tr class='{sig_class}'><td>{html.escape(row['GO_ID'])}</td>"
                            f"<td>{html.escape(row['GO_Name'])}</td>"
                            f"<td>{html.escape(row.get('Namespace', ''))}</td>"
                            f"<td>{html.escape(row.get('Fold_Enrichment', ''))}</td>"
                            f"<td>{html.escape(fdr_str)}</td></tr>"
                        )
                parts.append("</table>")
            # Images
            for img_name in ["general_go_bubble.png", "probiotic_go_bubble.png"]:
                img_path = goprobioticprofiler_dir / img_name
                if img_path.exists():
                    parts.append(f"<h3>{html.escape(img_name.replace('_', ' ').replace('.png', '').title())}</h3>")
                    parts.append(f"<img src='{img_path}' style='max-width:100%;' />")

        if busco.get("status") != "done":
            parts.append("<p class='note'>BUSCO result is not available in this run.</p>")
        parts.append("</section>")

    parts.append("<section class='panel'><h2>Step Status</h2><table><tr><th>Sample</th><th>Layout</th><th>Step</th><th>Status</th><th>Expected output</th></tr>")
    for name, layout, checks in rows:
        for step, path in checks.items():
            st = status(path)
            parts.append(f"<tr><td>{html.escape(name)}</td><td>{html.escape(layout)}</td><td>{html.escape(step)}</td><td class='{st}'>{st}</td><td>{html.escape(path)}</td></tr>")
    parts.append("</table>")
    multiqc_path = Path("data/results/multiqc/multiqc_report.html")
    if multiqc_path.exists():
        parts.append("<h2>MultiQC</h2><p>data/results/multiqc/multiqc_report.html</p>")
    parts.append("</section>")
    parts.append("</main></body></html>")
    Path(html_out).write_text("\n".join(parts))

    write_pdf(pdf_out, title, summaries, now)


if __name__ == "__main__":
    main()

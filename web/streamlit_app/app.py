"""Streamlit dashboard for submitting and monitoring PROBEAT analyses."""

from __future__ import annotations

import csv
import base64
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import zipfile

import redis
import streamlit as st


UPLOADS = Path(os.environ.get("PROBEAT_UPLOADS", "/app/data/uploads"))
RESULTS = Path(os.environ.get("PROBEAT_RESULTS", "/app/data/results"))
STATE = Path(os.environ.get("PROBEAT_STATE", "/app/data/state"))
DATABASES = Path(os.environ.get("PROBEAT_DATABASES", "/app/data/databases"))
ROOT = Path(os.environ.get("PROBEAT_ROOT", Path(__file__).resolve().parents[2]))
REDIS_URL = os.environ.get("REDIS_URL", "redis://redis:6379/0")
QUEUE = "probeat:jobs"
JOB_STATE = "probeat:job-state"
HEADER_LOGO = ROOT / "data" / "logo.png"
SIDEBAR_LOGO = ROOT / "data" / "probeat_logo2.png"
ZIP_DIR = STATE / "downloads"
MAX_UPLOAD_BYTES = 1024 * 1024 * 1024
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

RULE_LABELS = {
    "porechop_ont": "Porechop",
    "rasusa_ont": "Rasusa",
    "flye_ont": "Flye",
    "medaka_ont": "Medaka",
    "fastp_pe": "fastp",
    "fastp_se": "fastp",
    "fastqc": "FastQC",
    "unicycler": "Unicycler",
    "quast": "QUAST",
    "mash": "Mash",
    "sourmash_gtdb": "sourmash",
    "busco": "BUSCO",
    "bakta": "Bakta",
    "amrfinderplus": "AMRFinderPlus",
    "abricate_card": "ABRicate CARD",
    "abricate_megares": "ABRicate MEGARes",
    "plasmidfinder": "PlasmidFinder",
    "vfdb": "VFDB",
    "genomad": "geNomad",
    "mobsuite": "MOB-recon",
    "eggnog": "eggNOG-mapper",
    "kofamscan": "KofamScan",
    "keggdecoder": "KEGGDecoder",
    "dbcan": "dbCAN",
    "pgmdb": "PGMDB",
    "goprobioticprofiler": "GO Probiotic",
    "gut_brain_markers": "Gut-Brain markers",
    "kegganog": "KEGGaNOG",
    "genomespot": "GenomeSPOT",
    "gutsmash": "gutSMASH",
    "bagel4": "BAGEL4",
    "antismash": "antiSMASH",
    "fastani_references": "FastANI",
    "fastani": "FastANI",
    "pycirclize": "pyCirclize",
    "multiqc": "MultiQC",
    "final_report": "Report",
}


st.set_page_config(page_title="PROBEAT", page_icon="P", layout="wide")


st.markdown(
    """
    <style>
    :root {
        --probeat-blue: #1f2652;
        --probeat-gold: #d8b65c;
        --probeat-ink: #10233d;
        --probeat-muted: #65748b;
        --probeat-line: #d9e1ec;
        --probeat-panel: #f7f9fc;
    }
    .block-container { padding-top: 2.6rem; padding-bottom: 2.5rem; }
    .probeat-header {
        background: var(--probeat-blue);
        color: white;
        border-radius: 10px;
        padding: 18px 22px 16px 22px;
        margin: 18px 0 18px 0;
        border: 1px solid rgba(216, 182, 92, .45);
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 24px;
        overflow: visible;
    }
    .probeat-header-text { flex: 1 1 auto; min-width: 0; }
    .probeat-header h1 {
        margin: 0;
        font-size: 26px;
        line-height: 1.35;
        letter-spacing: 0;
        color: white;
        padding-top: 2px;
    }
    .probeat-header p {
        margin: 8px 0 0 0;
        color: #dce5f2;
        font-size: 14px;
    }
    .probeat-header-logo {
        flex: 0 0 auto;
        max-width: 112px;
        max-height: 64px;
        border-radius: 8px;
        object-fit: contain;
        display: block;
        background: white;
    }
    .sidebar-logo {
        width: 100%;
        border-radius: 14px;
        display: block;
        margin: 0 0 16px 0;
        border: 1px solid rgba(216, 182, 92, .45);
        background: white;
    }
    @media (max-width: 800px) {
        .probeat-header { align-items: flex-start; }
        .probeat-header-logo { max-width: 88px; max-height: 58px; }
    }
    .metric-panel {
        border: 1px solid var(--probeat-line);
        background: white;
        border-radius: 10px;
        padding: 14px 16px;
        min-height: 92px;
    }
    .metric-panel .label {
        color: var(--probeat-muted);
        font-size: 12px;
        text-transform: uppercase;
        font-weight: 700;
    }
    .metric-panel .value {
        color: var(--probeat-ink);
        font-size: 28px;
        font-weight: 800;
        margin-top: 4px;
    }
    .small-note { color: var(--probeat-muted); font-size: 13px; }
    .status-pill {
        display: inline-block;
        border-radius: 999px;
        padding: 3px 10px;
        font-size: 12px;
        font-weight: 800;
        border: 1px solid transparent;
    }
    .status-completed { background: #e6f7ef; color: #107044; border-color: #9ed8bd; }
    .status-running { background: #eaf2ff; color: #1d58a5; border-color: #a9c7f4; }
    .status-queued { background: #fff5d8; color: #866000; border-color: #e2c56d; }
    .status-failed { background: #fde8e8; color: #a32626; border-color: #efaaaa; }
    .status-cancelled { background: #eef1f5; color: #52606d; border-color: #c9d1db; }
    .sidebar-run-card {
        border: 1px solid #c7d3e3;
        background: #f7f9fc;
        border-radius: 10px;
        padding: 12px;
        margin: 12px 0;
    }
    .spinner-dot {
        width: 18px;
        height: 18px;
        border: 3px solid #c7d3e3;
        border-top-color: var(--probeat-blue);
        border-radius: 50%;
        display: inline-block;
        vertical-align: middle;
        margin-right: 8px;
        animation: probeat-spin 1s linear infinite;
    }
    @keyframes probeat-spin { to { transform: rotate(360deg); } }
    div[data-testid="stDataFrame"] { border: 1px solid var(--probeat-line); border-radius: 8px; }
    </style>
    """,
    unsafe_allow_html=True,
)


@st.cache_data(show_spinner=False)
def image_data_uri(path: str) -> str:
    image_path = Path(path)
    if not image_path.exists():
        return ""
    encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def render_header(title: str, subtitle: str) -> None:
    logo = image_data_uri(str(HEADER_LOGO))
    image = f"<img class='probeat-header-logo' src='{logo}' alt='PROBEAT logo'>" if logo else ""
    st.markdown(
        f"""
        <div class="probeat-header">
            <div class="probeat-header-text">
                <h1>{title}</h1>
                <p>{subtitle}</p>
            </div>
            {image}
        </div>
        """,
        unsafe_allow_html=True,
    )


@st.cache_resource(show_spinner=False)
def redis_client() -> redis.Redis | None:
    try:
        client = redis.Redis.from_url(REDIS_URL, decode_responses=True)
        client.ping()
        return client
    except redis.RedisError:
        return None


def now() -> float:
    return time.time()


def fmt_time(value: float | int | str | None) -> str:
    if not value:
        return "-"
    try:
        return time.strftime("%Y-%m-%d %H:%M", time.localtime(float(value)))
    except (TypeError, ValueError):
        return "-"


def fmt_duration(start: float | int | str | None, end: float | int | str | None = None) -> str:
    if not start:
        return "-"
    try:
        seconds = int((float(end) if end else now()) - float(start))
    except (TypeError, ValueError):
        return "-"
    if seconds < 60:
        return f"{seconds}s"
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes}m {seconds}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes}m"


def status_badge(status: str) -> str:
    clean = status or "unknown"
    css = clean if clean in {"completed", "running", "queued", "failed", "cancelled"} else "cancelled"
    return f"<span class='status-pill status-{css}'>{clean.upper()}</span>"


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    except OSError:
        return []


def read_quast(sample: str) -> dict[str, str]:
    rows = read_tsv(RESULTS / sample / "quast" / "report.tsv")
    if not rows:
        return {}
    row = rows[0]
    return {
        "contigs": row.get("# contigs", row.get("# contigs (>= 0 bp)", "")),
        "genome_size": row.get("Total length", row.get("Total length (>= 0 bp)", "")),
        "gc": row.get("GC (%)", ""),
        "n50": row.get("N50", ""),
    }


def read_fastp(sample: str) -> dict[str, str]:
    data = read_json(RESULTS / sample / "qc" / "fastp" / f"{sample}.json")
    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    return {
        "raw_reads": before.get("total_reads", ""),
        "clean_reads": after.get("total_reads", ""),
        "q30": after.get("q30_rate", ""),
    }


def read_busco(sample: str) -> str:
    path = RESULTS / sample / "busco" / "busco_summary.txt"
    if not path.exists():
        return ""
    try:
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.strip().startswith("C:"):
                return line.strip()
    except OSError:
        return ""
    return ""


def read_fastani_species(sample: str) -> str:
    rows = read_tsv(RESULTS / sample / "fastani" / f"{sample}_fastani.tsv")
    if not rows:
        return ""
    first = rows[0]
    for key in ("species", "reference_species", "reference", "query_species"):
        if first.get(key):
            return first[key]
    return next(iter(first.values()), "")


def read_pgmdb(sample: str) -> dict[str, str]:
    rows = read_tsv(RESULTS / sample / "pgmdb" / "summary.tsv")
    values = {row.get("metric", ""): row.get("value", "") for row in rows}
    return {
        "unique_markers": values.get("unique_markers", values.get("Unique markers", "")),
        "categories": values.get("top_categories", values.get("Top categories", "")),
    }


def result_report_paths(sample: str) -> tuple[Path, Path]:
    report_dir = RESULTS / sample / "report"
    return report_dir / f"{sample}_report.pdf", report_dir / f"{sample}_report.html"


def sample_metrics(sample: str) -> dict[str, str]:
    quast = read_quast(sample)
    fastp = read_fastp(sample)
    pgmdb = read_pgmdb(sample)
    pdf, html = result_report_paths(sample)
    return {
        "sample": sample,
        "raw_reads": str(fastp.get("raw_reads", "") or "-"),
        "clean_reads": str(fastp.get("clean_reads", "") or "-"),
        "contigs": str(quast.get("contigs", "") or "-"),
        "genome_size": str(quast.get("genome_size", "") or "-"),
        "gc": str(quast.get("gc", "") or "-"),
        "n50": str(quast.get("n50", "") or "-"),
        "species": read_fastani_species(sample) or "-",
        "busco": read_busco(sample) or "-",
        "pgmdb": str(pgmdb.get("unique_markers", "") or "-"),
        "pdf": "ready" if pdf.exists() else "missing",
        "html": "ready" if html.exists() else "missing",
    }


def redis_jobs() -> list[dict]:
    client = redis_client()
    if not client:
        return []
    try:
        return [json.loads(value) for value in client.hgetall(JOB_STATE).values()]
    except (redis.RedisError, json.JSONDecodeError):
        return []


def state_jobs() -> list[dict]:
    records = []
    for path in sorted((STATE / "jobs").glob("*.json")):
        data = read_json(path)
        if data:
            records.append(data)
    return records


def discovered_result_jobs() -> list[dict]:
    records = []
    if not RESULTS.exists():
        return records
    for path in sorted(RESULTS.iterdir()):
        if not path.is_dir() or path.name == "report":
            continue
        pdf, html = result_report_paths(path.name)
        if not pdf.exists() and not html.exists():
            continue
        records.append(
            {
                "sample": path.name,
                "layout": "unknown",
                "profile": "unknown",
                "status": "completed",
                "created_at": path.stat().st_mtime,
                "finished_at": max(pdf.stat().st_mtime if pdf.exists() else 0, html.stat().st_mtime if html.exists() else 0),
                "source": "results",
            }
        )
    return records


def jobs() -> list[dict]:
    by_sample: dict[str, dict] = {}
    for item in discovered_result_jobs() + state_jobs() + redis_jobs():
        sample = item.get("sample")
        if not sample:
            continue
        current = by_sample.get(sample, {})
        current.update(item)
        by_sample[sample] = current
    return sorted(by_sample.values(), key=lambda item: item.get("updated_at") or item.get("created_at") or 0, reverse=True)


def job_for(sample: str) -> dict | None:
    for item in jobs():
        if item.get("sample") == sample:
            return item
    return None


def save_upload(sample_dir: Path, uploaded) -> str:
    target = sample_dir / Path(uploaded.name).name
    uploaded.seek(0)
    with target.open("wb") as handle:
        while chunk := uploaded.read(1024 * 1024):
            handle.write(chunk)
    uploaded.seek(0)
    if target.name.lower().endswith((".fastq.gz", ".fq.gz")):
        return str(decompress_fastq_gz(target))
    return str(target)


def decompress_fastq_gz(path: Path) -> Path:
    lower_name = path.name.lower()
    if lower_name.endswith(".fastq.gz"):
        output = path.with_name(path.name[:-3])
    elif lower_name.endswith(".fq.gz"):
        output = path.with_name(path.name[:-3])
    else:
        return path
    tmp_output = output.with_suffix(output.suffix + ".tmp")
    try:
        with gzip.open(path, "rb") as source, tmp_output.open("wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
        tmp_output.replace(output)
    except (OSError, EOFError) as exc:
        tmp_output.unlink(missing_ok=True)
        raise RuntimeError(f"Could not decompress {path.name}. Please upload a valid gzip FASTQ file.") from exc
    return output


def valid_fastq_name(name: str) -> bool:
    return name.lower().endswith(FASTQ_SUFFIXES)


def sample_name_from_upload(uploaded) -> str:
    if uploaded is None:
        return ""
    name = Path(uploaded.name).name
    lower_name = name.lower()
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if lower_name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")
    return safe


def upload_error(uploaded, label: str) -> str | None:
    if uploaded is None:
        return f"{label} is required."
    if not valid_fastq_name(uploaded.name):
        return f"{label} must be .fastq, .fq, .fastq.gz, or .fq.gz."
    if getattr(uploaded, "size", 0) > MAX_UPLOAD_BYTES:
        return f"{label} is larger than 1 GB."
    return None


def sample_step_paths(sample: str, job: dict | None = None) -> list[tuple[str, list[Path]]]:
    layout = (job or {}).get("layout") or "PE"
    steps: list[tuple[str, list[Path]]] = []
    if layout == "ONT":
        steps.extend(
            [
                ("Porechop", list((RESULTS / sample / "qc" / "porechop").glob("*"))),
                ("Rasusa", list((RESULTS / sample / "qc" / "rasusa").glob("*"))),
                ("Flye", [RESULTS / sample / "flye" / "assembly.fasta"]),
                ("Medaka", [RESULTS / sample / "medaka" / "consensus.fasta"]),
            ]
        )
    elif layout == "HYBRID":
        steps.extend(
            [
                ("fastp", [RESULTS / sample / "qc" / "fastp" / f"{sample}.json"]),
                ("FastQC", list((RESULTS / sample / "qc" / "fastqc").glob("*_fastqc.zip"))),
                ("Porechop", list((RESULTS / sample / "qc" / "porechop").glob("*"))),
                ("Rasusa", list((RESULTS / sample / "qc" / "rasusa").glob("*"))),
                ("Unicycler", [RESULTS / sample / "assembly" / "assembly.fasta"]),
            ]
        )
    else:
        steps.extend(
            [
                ("fastp", [RESULTS / sample / "qc" / "fastp" / f"{sample}.json"]),
                ("FastQC", list((RESULTS / sample / "qc" / "fastqc").glob("*_fastqc.zip"))),
                ("Rasusa", list((RESULTS / sample / "qc" / "rasusa").glob("*"))),
                ("Unicycler", [RESULTS / sample / "assembly" / "assembly.fasta"]),
            ]
        )
    steps.extend(
        [
            ("QUAST", [RESULTS / sample / "quast" / "report.tsv"]),
            ("Mash", [RESULTS / sample / "mash" / "top_hits.tsv"]),
            ("sourmash", [RESULTS / sample / "sourmash" / "gather.with-lineages.csv"]),
            ("BUSCO", [RESULTS / sample / "busco" / "busco_summary.txt"]),
            ("Bakta", [RESULTS / sample / "bakta" / f"{sample}.gff3"]),
            ("AMRFinderPlus", [RESULTS / sample / "amrfinder" / f"{sample}.tsv"]),
            ("ABRicate CARD", [RESULTS / sample / "abricate" / f"{sample}_card.tsv"]),
            ("ABRicate MEGARes", [RESULTS / sample / "abricate" / f"{sample}_megares.tsv"]),
            ("PlasmidFinder", [RESULTS / sample / "plasmidfinder" / f"{sample}.tsv"]),
            ("VFDB", [RESULTS / sample / "vfdb" / f"{sample}.tsv"]),
            ("geNomad", [RESULTS / sample / "genomad" / "assembly_summary" / "assembly_plasmid_summary.tsv"]),
            ("MOB-recon", [RESULTS / sample / "mobsuite" / "contig_report.txt"]),
            ("eggNOG-mapper", [RESULTS / sample / "eggnog" / f"{sample}.emapper.annotations"]),
            ("KofamScan", [RESULTS / sample / "kofamscan" / f"{sample}.tsv"]),
            ("KEGGDecoder", [RESULTS / sample / "keggdecoder" / "pathway_completeness.tsv"]),
            ("dbCAN", [RESULTS / sample / "dbcan" / "overview.tsv"]),
            ("PGMDB", [RESULTS / sample / "pgmdb" / "summary.tsv"]),
            ("GO Probiotic", [RESULTS / sample / "goprobioticprofiler" / "probiotic_go_summary.json"]),
            ("Gut-Brain markers", [RESULTS / sample / "gut_brain_markers.tsv"]),
            ("KEGGaNOG", [RESULTS / sample / "kegganog" / f"{sample}_pathways.tsv"]),
            ("GenomeSPOT", [RESULTS / sample / "genomespot" / f"{sample}.predictions.tsv"]),
            ("gutSMASH", [RESULTS / sample / "gutsmash" / ".probeat_gutsmash_status.json"]),
            ("BAGEL4", [RESULTS / sample / "bagel4" / "bacteriocin_clusters.tsv"]),
            ("antiSMASH", [RESULTS / sample / "antismash" / "index.html"]),
            ("FastANI", [RESULTS / sample / "fastani" / f"{sample}_ani_heatmap.png"]),
            ("pyCirclize", [RESULTS / sample / "genome_map" / f"{sample}_genome_map.png"]),
            ("Report", [RESULTS / sample / "report" / f"{sample}_report.pdf"]),
        ]
    )
    return steps


def step_done(paths: list[Path]) -> bool:
    return any(path.exists() and path.stat().st_size > 0 for path in paths)


def rule_label(rule: str) -> str:
    return RULE_LABELS.get(rule, rule.replace("_", " ").strip().title())


def failed_rule(sample: str) -> str:
    text = log_tail(probeat_log_path(sample), limit=320)
    matches = re.findall(r"Error in rule ([A-Za-z0-9_]+):", text)
    return matches[-1] if matches else ""


def active_rule(sample: str) -> str:
    text = log_tail(probeat_log_path(sample), limit=420)
    started = re.findall(r"\n(?:localrule|rule) ([A-Za-z0-9_]+):[\s\S]*?\n\s+jobid: (\d+)", text)
    finished = set(re.findall(r"Finished jobid: (\d+)", text))
    active = [rule for rule, jobid in started if jobid not in finished]
    if active:
        return active[-1]
    matches = re.findall(r"(?:localrule|rule) ([A-Za-z0-9_]+):", text)
    return matches[-1] if matches else ""


def active_environment_setup(sample: str) -> str:
    text = log_tail(probeat_log_path(sample), limit=420)
    created = re.findall(r"Environment for /app/workflow/envs/([A-Za-z0-9_.-]+)\.yaml created", text)
    creating = re.findall(r"Creating conda environment workflow/envs/([A-Za-z0-9_.-]+)\.yaml", text)
    if not creating:
        return ""
    if len(creating) > len(created):
        return creating[-1]
    if creating[-1] not in created:
        return creating[-1]
    return ""


def derived_progress(sample: str, job: dict) -> dict:
    steps = sample_step_paths(sample, job)
    completed = sum(step_done(paths) for _, paths in steps)
    total = len(steps)
    first_missing = next((name for name, paths in steps if not step_done(paths)), "")
    percent = int(round((completed / total) * 100)) if total else 0
    if job.get("status") == "completed":
        percent = 100
        first_missing = "Complete"
    status = job.get("status")
    worker_step = job.get("current_step") or ""
    if status == "completed":
        current_step = "Complete"
    elif status == "failed":
        rule = failed_rule(sample)
        current_step = f"Failed: {rule_label(rule)}" if rule else f"Failed: {first_missing or worker_step or '-'}"
    elif status == "running":
        environment = active_environment_setup(sample)
        if environment:
            current_step = f"Preparing tool environments: {environment}"
        else:
            rule = active_rule(sample) or (worker_step if worker_step != "queued" else "")
            current_step = rule_label(rule) if rule else (first_missing or "Running")
    elif status == "queued":
        current_step = "Queued"
    else:
        current_step = first_missing or worker_step or "-"
    return {
        "completed_steps": completed,
        "total_steps": total,
        "progress_percent": percent,
        "current_step": current_step,
        "steps": [{"Step": name, "Status": "ready" if step_done(paths) else "pending"} for name, paths in steps],
    }


def probeat_log_path(sample: str) -> Path:
    return RESULTS / sample / "run" / "probeat.log"


def failure_summary(sample: str) -> str:
    text = log_tail(probeat_log_path(sample), limit=260)
    if not text:
        return ""
    patterns = [
        r"Error in rule [\s\S]*?(?=\nComplete log|\nWorkflowError|\Z)",
        r"Missing [^\n]+",
        r"Exception:[^\n]+",
        r"WorkflowError:[\s\S]*?(?=\n\S|\Z)",
        r"RuleException:[\s\S]*?(?=\n\[|\Z)",
    ]
    hits: list[str] = []
    for pattern in patterns:
        for match in re.finditer(pattern, text):
            value = match.group(0).strip()
            if value and value not in hits:
                hits.append(value)
    if hits:
        return "\n\n".join(hits[-3:])
    failed_lines = [line for line in text.splitlines() if re.search(r"error|exception|missing|failed", line, re.I)]
    return "\n".join(failed_lines[-12:])


def delete_sample_data(sample: str) -> None:
    for path in [RESULTS / sample, UPLOADS / sample, STATE / "jobs" / f"{sample}.json"]:
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    client = redis_client()
    if client:
        try:
            client.hdel(JOB_STATE, sample)
            client.delete(f"probeat:cancel:{sample}")
        except redis.RedisError:
            pass


def queue_job(job: dict) -> None:
    client = redis_client()
    if not client:
        raise RuntimeError("Redis is not available.")
    STATE.joinpath("jobs").mkdir(parents=True, exist_ok=True)
    try:
        STATE.joinpath("jobs", f"{job['sample']}.json").write_text(json.dumps(job, indent=2) + "\n", encoding="utf-8")
    except OSError as exc:
        raise RuntimeError(f"State path is not writable for sample {job['sample']}: {exc}") from exc
    client.hset(JOB_STATE, job["sample"], json.dumps(job))
    client.lpush(QUEUE, json.dumps(job))


def requeue_job(job: dict) -> None:
    job["status"] = "queued"
    job["created_at"] = now()
    job["updated_at"] = now()
    queue_job(job)


def clean_gutsmash_outputs(sample: str) -> None:
    shutil.rmtree(RESULTS / sample / "gutsmash", ignore_errors=True)
    log_dir = ROOT / "logs" / "gutsmash"
    for path in log_dir.glob(f"*{sample}*"):
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()


def run_command(command: list[str], timeout: int = 900) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=ROOT, text=True, capture_output=True, timeout=timeout, check=False)


def render_metric(label: str, value: str | int, note: str = "") -> None:
    st.markdown(
        f"""
        <div class="metric-panel">
            <div class="label">{label}</div>
            <div class="value">{value}</div>
            <div class="small-note">{note}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def full_profile_for_layout(layout: str) -> str:
    if layout == "ONT":
        return "ont"
    if layout == "HYBRID":
        return "full"
    return "probiotic"


def available_read_layouts() -> list[str]:
    layouts = ["PE", "SE"]
    ont_markers = {
        "porechop": [ROOT / "mamba_envs" / "porechop-0.2.4" / "bin" / "porechop"],
        "flye": [ROOT / "mamba_envs" / "flye-2.9.6" / "bin" / "flye"],
        "medaka": [
            ROOT / "mamba_envs" / "medaka" / "bin" / "medaka_consensus",
            ROOT / "mamba_envs" / "medaka" / "bin" / "medaka",
        ],
    }
    if all(any(path.exists() for path in candidates) for candidates in ont_markers.values()):
        layouts.extend(["ONT", "HYBRID"])
    return layouts


def workflow_name(layout: str) -> str:
    if layout == "ONT":
        return "Complete ONT workflow"
    if layout == "HYBRID":
        return "Complete hybrid workflow"
    return "Complete Illumina workflow"


def render_dashboard() -> None:
    render_header("Dashboard", "Current analyses and downloadable results.")
    records = jobs()
    if not records:
        st.info("No analyses yet. Start a new analysis from the left menu.")
        return

    active = [item for item in records if item.get("status") in {"queued", "running"}]
    selected_default = active[0]["sample"] if active else records[0]["sample"]
    selected = st.selectbox("Sample", [item["sample"] for item in records], index=[item["sample"] for item in records].index(selected_default))
    selected_job = job_for(selected) or {}
    progress = derived_progress(selected, selected_job)

    cols = st.columns([1, 1, 1, 1])
    with cols[0]:
        render_metric("Status", selected_job.get("status", "unknown").upper(), workflow_name(selected_job.get("layout", "PE")))
    with cols[1]:
        render_metric("Current Step", progress["current_step"], f"{progress['completed_steps']} of {progress['total_steps']} tools")
    with cols[2]:
        render_metric("Progress", f"{progress['progress_percent']}%", "workflow")
    with cols[3]:
        render_metric("Runtime", fmt_duration(selected_job.get("started_at"), selected_job.get("finished_at")), "wall time")
    st.progress(progress["progress_percent"] / 100, text=f"{progress['current_step']} - {progress['progress_percent']}%")

    actions = st.columns([1, 1, 1, 1.35])
    if actions[0].button("Refresh", width="stretch"):
        st.rerun()
    client = redis_client()
    if selected_job.get("status") in {"queued", "running"} and actions[1].button("Cancel", width="stretch"):
        if client:
            if selected_job.get("status") == "queued":
                selected_job["status"] = "cancelled"
                client.hset(JOB_STATE, selected, json.dumps(selected_job))
            else:
                client.set(f"probeat:cancel:{selected}", "1", ex=3600)
            st.rerun()
    rerunnable = selected_job.get("status") in {"failed", "failed_preflight", "completed", "cancelled"}
    if rerunnable and actions[2].button("Run again", width="stretch"):
        if client:
            requeue_job(selected_job)
            st.rerun()
    if rerunnable and actions[3].button("Clean gutSMASH + Run again", width="stretch"):
        if client:
            clean_gutsmash_outputs(selected)
            requeue_job(selected_job)
            st.rerun()

    st.subheader("Downloads")
    for item in records:
        sample = item.get("sample", "")
        pdf, _ = result_report_paths(sample)
        row = st.columns([1.3, 1, 1, 1])
        row[0].markdown(
            f"**{sample}** {status_badge(item.get('status', 'unknown'))}<br>"
            f"<span class='small-note'>{fmt_time(item.get('updated_at') or item.get('finished_at') or item.get('created_at'))}</span>",
            unsafe_allow_html=True,
        )
        if pdf.exists():
            row[1].download_button("PDF", pdf.read_bytes(), pdf.name, "application/pdf", key=f"pdf-{sample}", width="stretch")
        else:
            row[1].button("PDF", disabled=True, key=f"pdf-missing-{sample}", width="stretch")
        if zip_is_current(sample):
            zip_path = zip_result_path(sample)
            with zip_path.open("rb") as zip_handle:
                row[2].download_button("ZIP", zip_handle, zip_path.name, "application/zip", key=f"zip-{sample}", width="stretch")
        elif (RESULTS / sample).exists():
            if row[2].button("Create ZIP", key=f"zip-create-{sample}", width="stretch"):
                with st.spinner(f"Creating ZIP archive for {sample}..."):
                    try:
                        build_zip_results(sample)
                    except OSError as exc:
                        st.error(f"Could not create ZIP archive: {exc}")
                        continue
                st.rerun()
        else:
            row[2].button("ZIP", disabled=True, key=f"zip-missing-{sample}", width="stretch")
        if row[3].button("Delete", key=f"delete-{sample}", width="stretch", disabled=item.get("status") in {"queued", "running"}):
            delete_sample_data(sample)
            st.rerun()

    with st.expander("Recent log", expanded=False):
        st.caption(f"Showing probeat.log for {selected}")
        summary = failure_summary(selected) if selected_job.get("status") == "failed" else ""
        if summary:
            st.error(summary)
        tail = log_tail(probeat_log_path(selected))
        if tail:
            st.code(tail, language="text")
        else:
            st.info("No probeat.log text is available yet.")

    if selected_job.get("status") in {"queued", "running"}:
        time.sleep(5)
        st.rerun()


def render_new_analysis() -> None:
    render_header("New Analysis", "Upload reads and start the complete PROBEAT workflow.")
    if not redis_client():
        st.error("Redis is not reachable. Start Docker Compose services before submitting new analyses.")
        return
    sample = st.text_input("Sample name", help="Letters, numbers, dots, underscores, and hyphens are allowed.")
    layouts = available_read_layouts()
    layout = st.segmented_control("Read type", layouts, default=layouts[0])
    profile = full_profile_for_layout(layout)
    st.caption(f"{workflow_name(layout)} will run automatically. Profiles are intentionally hidden.")
    col_a, col_b = st.columns(2)
    cores = col_a.number_input("CPU cores", min_value=1, max_value=128, value=16)
    memory_gb = col_b.number_input("Memory limit (GB)", min_value=4, max_value=1024, value=64)
    first_label = "R1 FASTQ" if layout in {"PE", "HYBRID"} else "FASTQ"
    first = st.file_uploader(first_label, type=["fastq", "fq", "gz"], key=f"reads-{layout}-first")
    second = st.file_uploader("R2 FASTQ", type=["fastq", "fq", "gz"], key=f"reads-{layout}-second") if layout in {"PE", "HYBRID"} else None
    long_reads = st.file_uploader("Long reads FASTQ", type=["fastq", "fq", "gz"], key="reads-hybrid-long") if layout == "HYBRID" else None
    st.caption("Accepted files: .fastq, .fq, .fastq.gz, .fq.gz. Maximum size: 1 GB per file.")
    check_databases = st.checkbox("Check databases before start", value=True)
    submitted = st.button("Start analysis", type="primary", width="stretch")

    with st.expander("Database check", expanded=False):
        st.write("Checks whether the databases required for the complete workflow are available.")
        db_layout = st.segmented_control("Workflow", ["PE / SE", "ONT", "HYBRID"], default="PE / SE", key="db-layout-check")
        db_profile = full_profile_for_layout("PE" if db_layout == "PE / SE" else db_layout)
        if st.button("Check databases", width="stretch"):
            with st.spinner("Checking databases..."):
                result = run_command([str(ROOT / "probeat"), "databases", "verify", "--profile", db_profile])
            st.success("All required databases were found.") if result.returncode == 0 else st.error("Some databases are missing.")
            st.code(result.stdout + result.stderr, language="text")

    if not submitted:
        return
    sample = sample.strip() or sample_name_from_upload(first)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sample or ""):
        st.error("Use a sample name containing only letters, numbers, dots, underscores, or hyphens, or upload a FASTQ file with a valid filename.")
        return
    errors = [upload_error(first, first_label)]
    if layout in {"PE", "HYBRID"}:
        errors.append(upload_error(second, "R2 FASTQ"))
    if layout == "HYBRID":
        errors.append(upload_error(long_reads, "Long reads FASTQ"))
    errors = [error for error in errors if error]
    if errors:
        for error in errors:
            st.error(error)
        return
    existing = job_for(sample)
    if existing and existing.get("status") in {"queued", "running"}:
        st.error("A queued or running job already exists for this sample.")
        return
    if check_databases:
        with st.spinner("Checking databases..."):
            preflight = run_command([str(ROOT / "probeat"), "databases", "verify", "--profile", profile])
        if preflight.returncode:
            st.error("Required databases are missing.")
            st.code(preflight.stdout + preflight.stderr, language="text")
            return
    sample_dir = UPLOADS / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    try:
        reads = [save_upload(sample_dir, first)]
        if second:
            reads.append(save_upload(sample_dir, second))
        if long_reads:
            reads.append(save_upload(sample_dir, long_reads))
    except RuntimeError as exc:
        st.error(str(exc))
        return
    job = {
        "sample": sample,
        "layout": layout,
        "profile": profile,
        "cores": int(cores),
        "memory_gb": int(memory_gb),
        "reads": reads,
        "status": "queued",
        "created_at": now(),
        "source": "streamlit",
    }
    try:
        queue_job(job)
    except RuntimeError as exc:
        st.error(str(exc))
        return
    st.success(f"{sample} was added to the analysis queue.")
    st.info("Open the Dashboard page to follow the current step and download results.")


def log_tail(path: Path, limit: int = 160) -> str:
    if not path.exists():
        return ""
    try:
        return "\n".join(path.read_text(encoding="utf-8", errors="replace").splitlines()[-limit:])
    except OSError:
        return ""


def zip_result_path(sample: str) -> Path:
    return ZIP_DIR / f"{sample}_results.zip"


def newest_result_mtime(sample: str) -> float:
    base = RESULTS / sample
    newest = base.stat().st_mtime if base.exists() else 0.0
    if not base.exists():
        return newest
    for path in base.rglob("*"):
        if path.is_file():
            try:
                newest = max(newest, path.stat().st_mtime)
            except OSError:
                continue
    return newest


def zip_is_current(sample: str) -> bool:
    zip_path = zip_result_path(sample)
    if not zip_path.exists():
        return False
    try:
        if not zipfile.is_zipfile(zip_path):
            return False
        return zip_path.stat().st_mtime >= newest_result_mtime(sample)
    except OSError:
        return False


def build_zip_results(sample: str) -> Path:
    base = RESULTS / sample
    if not base.exists():
        raise FileNotFoundError(f"Missing results directory: {base}")
    ZIP_DIR.mkdir(parents=True, exist_ok=True)
    final_zip = zip_result_path(sample)
    tmp_zip = ZIP_DIR / f".{sample}.{os.getpid()}.{int(time.time() * 1000)}.zip.tmp"
    try:
        with zipfile.ZipFile(tmp_zip, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in base.rglob("*"):
                if path.is_file():
                    archive.write(path, path.relative_to(base.parent))
        tmp_zip.replace(final_zip)
    finally:
        if tmp_zip.exists():
            tmp_zip.unlink()
    return final_zip


PAGES = {
    "New Analysis": render_new_analysis,
    "Dashboard": render_dashboard,
}

active = [item for item in jobs() if item.get("status") in {"queued", "running"}]
page_names = list(PAGES)
sidebar_logo = image_data_uri(str(SIDEBAR_LOGO))
if sidebar_logo:
    st.sidebar.markdown(
        f"<img class='sidebar-logo' src='{sidebar_logo}' alt='PROBEAT workflow logo'>",
        unsafe_allow_html=True,
    )
page = st.sidebar.radio("PROBEAT", page_names, index=1 if active else 0)
if active:
    first = active[0]
    active_progress = derived_progress(first.get("sample", ""), first)
    st.sidebar.markdown(
        f"""
        <div class="sidebar-run-card">
            <span class="spinner-dot"></span><strong>Workflow running</strong><br>
            <span class="small-note">{first.get('sample', 'sample')} - {active_progress['current_step']}</span>
        </div>
        """,
        unsafe_allow_html=True,
    )
st.sidebar.caption(f"Workflow version {os.environ.get('PROBEAT_WORKFLOW_VERSION', '1.0')}")
st.sidebar.caption(f"Port: {os.environ.get('PROBEAT_WEB_PORT', '9797')}")
PAGES[page]()

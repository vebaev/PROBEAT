#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace="${PROBEAT_ROOT:-$(cd "$script_dir/.." && pwd)}"
cd "$workspace"

usage() {
    cat <<'EOF'
Usage:
  workflow/run_sample.sh SAMPLE PE  R1_FASTQ R2_FASTQ [SNAKEMAKE_OPTIONS...]
  workflow/run_sample.sh SAMPLE SE  FASTQ             [SNAKEMAKE_OPTIONS...]
  workflow/run_sample.sh SAMPLE ONT FASTQ             [SNAKEMAKE_OPTIONS...]
  workflow/run_sample.sh SAMPLE HYBRID R1_FASTQ R2_FASTQ LONG_FASTQ [SNAKEMAKE_OPTIONS...]

Examples:
  workflow/run_sample.sh PE01 PE /data/PE01_R1.fastq.gz /data/PE01_R2.fastq.gz --cores 16
  workflow/run_sample.sh SE01 SE /data/SE01.fastq.gz --cores 16
  workflow/run_sample.sh ONT01 ONT /data/ONT01.fastq.gz --cores 16
  workflow/run_sample.sh HY01 HYBRID /data/HY01_R1.fastq.gz /data/HY01_R2.fastq.gz /data/HY01_ONT.fastq.gz --cores 16
  workflow/run_sample.sh ONT01 ONT /data/ONT01.fastq.gz --dry-run --cores 16
EOF
}

if [[ $# -lt 3 ]]; then
    usage >&2
    exit 2
fi

sample="$1"
layout="${2^^}"
r1="$3"
shift 3
long_reads=""

if [[ ! "$sample" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
    echo "Invalid sample name: $sample" >&2
    echo "Use letters, numbers, dots, underscores, or hyphens." >&2
    exit 2
fi

case "$layout" in
    PE)
        if [[ $# -lt 1 ]]; then
            echo "PE layout requires both R1 and R2 FASTQ files." >&2
            usage >&2
            exit 2
        fi
        r2="$1"
        shift
        ;;
    HYBRID)
        if [[ $# -lt 2 ]]; then
            echo "HYBRID layout requires R1, R2, and long-read FASTQ files." >&2
            usage >&2
            exit 2
        fi
        r2="$1"
        long_reads="$2"
        shift 2
        ;;
    SE|ONT)
        r2=""
        ;;
    *)
        echo "Unsupported layout: $layout. Use PE, SE, ONT, or HYBRID." >&2
        exit 2
        ;;
esac

if [[ ! -f "$r1" ]]; then
    echo "FASTQ file does not exist: $r1" >&2
    exit 2
fi
if [[ "$layout" == "PE" && ! -f "$r2" ]]; then
    echo "FASTQ file does not exist: $r2" >&2
    exit 2
fi
if [[ "$layout" == "HYBRID" ]]; then
    if [[ ! -f "$r2" ]]; then
        echo "FASTQ file does not exist: $r2" >&2
        exit 2
    fi
    if [[ ! -f "$long_reads" ]]; then
        echo "FASTQ file does not exist: $long_reads" >&2
        exit 2
    fi
fi

r1="$(realpath "$r1")"
if [[ -n "$r2" ]]; then
    r2="$(realpath "$r2")"
fi
if [[ -n "$long_reads" ]]; then
    long_reads="$(realpath "$long_reads")"
fi

run_dir="data/results/$sample/run"
sample_file="$run_dir/samples.tsv"
mkdir -p "$run_dir"
printf 'sample\tlayout\tr1\tr2\tlong_reads\n%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$layout" "$r1" "$r2" "$long_reads" > "$sample_file"

target="data/results/$sample/report/${sample}_report.pdf"
extra_args=("$@")

has_rerun_triggers=false
for arg in "${extra_args[@]}"; do
    if [[ "$arg" == "--rerun-triggers" ]]; then
        has_rerun_triggers=true
        break
    fi
done
if [[ "$has_rerun_triggers" == false ]]; then
    extra_args+=(--rerun-triggers mtime)
fi

echo "Sample: $sample"
echo "Layout: $layout"
echo "Input: $r1${r2:+, $r2}${long_reads:+, $long_reads}"
echo "Target: $target"

config_args=("samples_file=$sample_file")
if [[ -n "${PROBEAT_ANALYSIS_PROFILE:-}" ]]; then
    config_args+=("analysis_profile=$PROBEAT_ANALYSIS_PROFILE")
fi

workflow/run_snakemake.sh "$target" \
    --config "${config_args[@]}" \
    "${extra_args[@]}"

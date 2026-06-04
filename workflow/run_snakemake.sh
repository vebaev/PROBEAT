#!/usr/bin/env bash
set -euo pipefail

workspace="/home/vebaev/codex"
cd "$workspace"

mkdir -p logs

timestamp="$(date --iso-8601=seconds)"
command=(snakemake --snakefile workflow/Snakefile --profile workflow/profiles/local "$@")

{
    printf '[%s] cwd=%s command=' "$timestamp" "$workspace"
    printf '%q ' "${command[@]}"
    printf '\n'
} >> logs/commands.log

"${command[@]}" 2>&1 | tee -a "logs/snakemake-${timestamp}.log"

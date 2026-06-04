#!/usr/bin/env bash
set -euo pipefail

workspace="/home/vebaev/codex"
helper_env="$workspace/mamba_envs/snakemake-helper-2026-06-04"
cd "$workspace"

export XDG_CACHE_HOME="$workspace/.cache"
export APPTAINER_CACHEDIR="$workspace/.cache/apptainer"
export APPTAINER_TMPDIR="$workspace/.cache/apptainer-tmp"
export TMPDIR="$workspace/.tmp"

mkdir -p logs "$XDG_CACHE_HOME" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$TMPDIR"

timestamp="$(date --iso-8601=seconds)"
command=("$helper_env/bin/snakemake" --snakefile workflow/Snakefile --profile workflow/profiles/local "$@")

{
    printf '[%s] cwd=%s command=' "$timestamp" "$workspace"
    printf '%q ' "${command[@]}"
    printf '\n'
} >> logs/commands.log

"${command[@]}" 2>&1 | tee -a "logs/snakemake-${timestamp}.log"

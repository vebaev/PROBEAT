#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace="${PROBEAT_ROOT:-$(cd "$script_dir/.." && pwd)}"
helper_env="${PROBEAT_HELPER_ENV:-$workspace/mamba_envs/snakemake-conda-helper-2026-06-04}"
cd "$workspace"

if [[ ! -x "$helper_env/bin/snakemake" ]]; then
    echo "Snakemake helper environment is missing: $helper_env" >&2
    echo "Run: ./probeat setup --profile core" >&2
    exit 2
fi

export PATH="$helper_env/bin:$PATH"
export PROBEAT_ROOT="$workspace"
export PROBEAT_PYTHON="$helper_env/bin/python"
export LANG=C
export LC_ALL=C
export XDG_CACHE_HOME="$workspace/.cache"
export MPLCONFIGDIR="$workspace/.cache/matplotlib"
export CONDARC="$workspace/workflow/config/condarc"
export CONDA_PKGS_DIRS="$workspace/micromamba-root/pkgs"
export CONDA_ENVS_PATH="$workspace/snakemake_conda_envs"
export APPTAINER_CACHEDIR="$workspace/.cache/apptainer"
export APPTAINER_TMPDIR="$workspace/.cache/apptainer-tmp"
export TMPDIR="$workspace/.tmp"

mkdir -p logs "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$CONDA_PKGS_DIRS" "$CONDA_ENVS_PATH" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$TMPDIR"

timestamp="$(date --iso-8601=seconds)"
command=("$helper_env/bin/snakemake" --snakefile workflow/Snakefile --profile workflow/profiles/local "$@")

{
    printf '[%s] cwd=%s command=' "$timestamp" "$workspace"
    printf '%q ' "${command[@]}"
    printf '\n'
} >> logs/commands.log

"${command[@]}" 2>&1 | tee -a "logs/snakemake-${timestamp}.log"

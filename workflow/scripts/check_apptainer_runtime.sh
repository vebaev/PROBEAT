#!/usr/bin/env bash
set -euo pipefail

workspace="/home/vebaev/codex"
apptainer="$workspace/mamba_envs/snakemake-helper-2026-06-04/bin/apptainer"
image="$workspace/app_containers/fastp_1.3.3.sif"

if [[ ! -s "$image" ]]; then
    echo "Missing test image: $image" >&2
    exit 2
fi

"$apptainer" exec --cleanenv --bind "$workspace:$workspace" "$image" fastp --version

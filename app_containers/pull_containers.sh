#!/usr/bin/env bash
set -euo pipefail

workspace="/home/vebaev/codex"
manifest="$workspace/app_containers/containers.tsv"
apptainer="$workspace/mamba_envs/snakemake-helper-2026-06-04/bin/apptainer"

export APPTAINER_CACHEDIR="$workspace/.cache/apptainer"
export APPTAINER_TMPDIR="$workspace/.cache/apptainer-tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$workspace/logs/containers"

tail -n +2 "$manifest" | while IFS=$'\t' read -r name image uri status; do
    if [[ "$uri" == "manual" ]]; then
        printf '[SKIP] %s: %s\n' "$name" "$status" | tee -a "$workspace/logs/containers/pull.log"
        continue
    fi

    if [[ -s "$workspace/$image" ]]; then
        printf '[OK] %s already exists: %s\n' "$name" "$image" | tee -a "$workspace/logs/containers/pull.log"
        continue
    fi

    printf '[PULL] %s <- %s\n' "$name" "$uri" | tee -a "$workspace/logs/containers/pull.log"
    "$apptainer" pull "$workspace/$image" "$uri" 2>&1 | tee -a "$workspace/logs/containers/${name}.log"
done

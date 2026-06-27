#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
output="$root/manifests/conda-lock"
mkdir -p "$output"

if ! command -v conda-lock >/dev/null 2>&1; then
    echo "conda-lock is required. Install it in a dedicated helper environment." >&2
    exit 2
fi

for environment in "$root"/workflow/envs/*.yaml; do
    name="$(basename "$environment" .yaml)"
    echo "==> $name"
    conda-lock lock \
        --file "$environment" \
        --platform linux-64 \
        --lockfile "$output/$name.conda-lock.yaml"
done

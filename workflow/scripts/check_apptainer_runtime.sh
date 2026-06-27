#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace="${PROBEAT_ROOT:-$(cd "$script_dir/../.." && pwd)}"
apptainer="${APPTAINER_BIN:-$workspace/mamba_envs/snakemake-helper-2026-06-04/bin/apptainer}"
image="$workspace/app_containers/fastp_1.3.3.sif"

if [[ ! -s "$image" ]]; then
    echo "Missing test image: $image" >&2
    exit 2
fi

if [[ ! -e /dev/fuse ]]; then
    echo "Missing /dev/fuse. Apptainer cannot execute these containers on this server without FUSE or system Apptainer support." >&2
    exit 3
fi

"$apptainer" exec --cleanenv --bind "$workspace:$workspace" "$image" fastp --version

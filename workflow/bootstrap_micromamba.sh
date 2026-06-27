#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target="$root/.tools/bin"
mkdir -p "$target"

if [[ -x "$target/micromamba" ]]; then
    echo "Micromamba is already installed: $target/micromamba"
    exit 0
fi

archive="$root/.tools/micromamba.tar.bz2"
curl -L --fail --output "$archive" https://micro.mamba.pm/api/micromamba/linux-64/latest
tar -xjf "$archive" -C "$root/.tools" bin/micromamba
chmod +x "$target/micromamba"
echo "Installed Micromamba: $target/micromamba"

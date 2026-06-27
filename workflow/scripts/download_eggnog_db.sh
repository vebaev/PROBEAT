#!/usr/bin/env bash
set -euo pipefail

DB_DIR="${1:-data/databases/eggnog}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EGGNOG_ENV="$(dirname "$SCRIPT_DIR")/../mamba_envs/eggnog"

BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"

mkdir -p "$DB_DIR"

echo "=== eggNOG v5 database downloader ==="
echo "Target: $DB_DIR"
echo "Base URL: $BASE_URL"
echo ""

ALL_OK=true
for f in eggnog.db eggnog.taxa.db eggnog_proteins.dmnd; do
    if [[ -s "$DB_DIR/$f" ]]; then
        echo "[OK] $f ($(du -h "$DB_DIR/$f" | cut -f1))"
    else
        echo "[MISSING] $f"
        ALL_OK=false
    fi
done

if $ALL_OK; then
    echo ""
    echo "All eggNOG database files are present."
    exit 0
fi

echo ""
echo "=== Checking eggNOG download server ==="

for url in \
    "$BASE_URL/eggnog.db.gz" \
    "$BASE_URL/eggnog.taxa.tar.gz" \
    "$BASE_URL/eggnog_proteins.dmnd.gz"
do
    echo "Testing: $url"
    if ! wget -4 -q --spider "$url"; then
        echo ""
        echo "ERROR: Cannot reach:"
        echo "  $url"
        echo ""
        echo "Try manually:"
        echo "  wget -4 $url"
        exit 1
    fi
done

echo ""
echo "=== Checking eggNOG environment ==="

if [[ ! -d "$EGGNOG_ENV" ]]; then
    echo "Error: eggNOG conda environment not found at:"
    echo "  $EGGNOG_ENV" >&2
    exit 1
fi

DOWNLOAD_PY="$EGGNOG_ENV/bin/download_eggnog_data.py"

if [[ ! -f "$DOWNLOAD_PY" ]]; then
    echo "Error: download_eggnog_data.py not found at:"
    echo "  $DOWNLOAD_PY" >&2
    exit 1
fi

echo "Downloader: $DOWNLOAD_PY"

if grep -q "eggnogdb.embl.de" "$DOWNLOAD_PY"; then
    echo "Patching old eggNOG domain in downloader..."
    cp "$DOWNLOAD_PY" "$DOWNLOAD_PY.bak"
    sed -i 's|eggnogdb.embl.de|eggnog5.embl.de|g' "$DOWNLOAD_PY"
fi

echo ""
echo "=== Starting automatic download ==="

yes | "$DOWNLOAD_PY" --data_dir "$DB_DIR"

echo ""
echo "=== Verifying downloaded files ==="

for f in eggnog.db eggnog.taxa.db eggnog_proteins.dmnd; do
    if [[ -s "$DB_DIR/$f" ]]; then
        echo "[OK] $f ($(du -h "$DB_DIR/$f" | cut -f1))"
    else
        echo "[ERROR] Missing after download: $f"
        exit 1
    fi
done

echo ""
echo "=== eggNOG database installed successfully ==="
ls -lh "$DB_DIR"
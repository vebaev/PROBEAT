#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$root"

export DOCKER_CONFIG="$root/.docker"
mkdir -p "$DOCKER_CONFIG"

docker build --build-arg INSTALL_PROFILE=core -t probeat:docker-smoke .
PROBEAT_IMAGE=probeat:docker-smoke PROBEAT_INSTALL_PROFILE=core \
    docker compose up -d --no-build redis streamlit worker
trap 'docker compose down' EXIT
timeout 180 bash -c 'until curl -fsS http://localhost:9797/_stcore/health; do sleep 5; done'

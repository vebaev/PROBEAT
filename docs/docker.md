# Docker Deployment

The documented deployment uses the full PROBEAT image, Streamlit, Redis, one
worker, and external storage volumes.

## Build

```bash
cp .env.example .env
docker build --build-arg INSTALL_PROFILE=full -t probeat:full .
```

## Databases

```bash
docker-compose run --rm database-init
docker-compose run --rm database-manager
```

## Start Services

```bash
docker-compose up -d redis worker streamlit
```

The interface is available at:

```text
http://localhost:9797
```

## Storage

The default `.env.example` maps runtime data to:

```text
./storage/uploads
./storage/databases
./storage/results
./storage/logs
./storage/state
```

These directories are ignored by Git.

## Docker Config

On servers where `$HOME` is read-only, set Docker's config directory to a
writable project-local path before building:

```bash
export DOCKER_CONFIG="$PWD/.docker"
mkdir -p "$DOCKER_CONFIG"
```

# Installation

PROBEAT is documented for Docker deployment with the full profile.

## Full Docker Installation

```bash
git clone https://github.com/vebaev/PROBEAT.git
cd PROBEAT
cp .env.example .env
docker build --build-arg INSTALL_PROFILE=full -t probeat:full .
docker-compose run --rm database-init
docker-compose run --rm database-manager
docker-compose up -d redis worker streamlit
```

Open:

```text
http://localhost:9797
```

The database download step stores large biological databases under
`./storage/databases`. Runtime outputs are written under `./storage/results`.

## Validation

```bash
docker-compose run --rm worker ./probeat doctor --profile full --strict
docker-compose run --rm worker workflow/run_sample.sh SMOKE_SE SE tests/data/smoke.fastq --dry-run --cores 2
```

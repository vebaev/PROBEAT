# PROBEAT

PROBEAT (PRObiotic Bacterial gEnome Analysis Toolkit) is a Dockerized
Snakemake workflow for bacterial isolate genome analysis with a probiotic
interpretation layer, Streamlit web UI, Redis worker, and HTML/PDF reporting.

This repository documents the **full** Docker profile only. The full profile
builds the complete runtime image for Illumina PE/SE, Oxford Nanopore ONT, and
HYBRID analyses.

## Requirements

- Docker Engine with Docker Compose (`docker-compose`)
- Internet access for the first image build and database download
- Sufficient disk space for the full image, workflow outputs, and downloaded
  biological databases
- Linux server or workstation recommended for production runs

The databases are not stored in Git and are not baked into the image. They are
downloaded into `./storage/databases` by the `database-init` Compose service.

## Quickstart

Clone the repository:

```bash
git clone https://github.com/vebaev/PROBEAT.git
cd PROBEAT
```

Create the local environment file:

```bash
cp .env.example .env
```

Build the full Docker image:

```bash
docker build --build-arg INSTALL_PROFILE=full -t probeat:full .
```

Download all required databases:

```bash
docker-compose run --rm database-init
```

Verify the downloaded databases:

```bash
docker-compose run --rm database-manager
```

Start the application stack:

```bash
docker-compose up -d redis worker streamlit
```

Check that the web UI is healthy:

```bash
curl http://localhost:9797/_stcore/health
```

Open the interface:

```text
http://localhost:9797
```

## Validate Installation

Run the full runtime doctor:

```bash
docker-compose run --rm worker ./probeat doctor --profile full --strict
```

Run the gutSMASH runtime self-test:

```bash
docker-compose run --rm worker \
  /app/mamba_envs/snakemake-conda-helper-2026-06-04/bin/python \
  /app/workflow/scripts/run_gutsmash_wrapper.py --self-test
```

Run a Snakemake dry-run with the bundled smoke FASTQ:

```bash
docker-compose run --rm worker \
  workflow/run_sample.sh SMOKE_SE SE tests/data/smoke.fastq \
  --dry-run --cores 2
```

## Data Locations

The default `.env.example` keeps runtime data under `./storage`:

- `./storage/uploads`: uploaded FASTQ files
- `./storage/results`: workflow outputs and reports
- `./storage/databases`: downloaded biological databases
- `./storage/logs`: workflow and service logs
- `./storage/state`: job state files

These directories are intentionally ignored by Git.

## Running Samples

Use the Streamlit UI for routine operation. For command-line runs inside the
worker container:

```bash
docker-compose run --rm worker \
  workflow/run_sample.sh SAMPLE PE /app/data/uploads/SAMPLE_R1.fastq.gz /app/data/uploads/SAMPLE_R2.fastq.gz --cores 16

docker-compose run --rm worker \
  workflow/run_sample.sh SAMPLE SE /app/data/uploads/SAMPLE.fastq.gz --cores 16

docker-compose run --rm worker \
  workflow/run_sample.sh SAMPLE ONT /app/data/uploads/SAMPLE.fastq.gz --cores 16
```

Results are written to:

```text
storage/results/<sample>/
```

The final report is written to:

```text
storage/results/<sample>/report/<sample>_report.html
storage/results/<sample>/report/<sample>_report.pdf
```

## Stop Services

```bash
docker-compose down
```

This stops the services but keeps `./storage` data.

## Troubleshooting

Check container status:

```bash
docker-compose ps
```

Check Streamlit logs:

```bash
docker logs --tail=120 probeat-streamlit-1
```

Check worker logs:

```bash
docker logs --tail=120 probeat-worker-1
```

Re-run database verification:

```bash
docker-compose run --rm database-manager
```

If a workflow step fails, inspect:

```text
storage/results/<sample>/run/probeat.log
storage/logs/
```

## Repository Hygiene

Do not commit runtime data, downloaded databases, FASTQ/SRA inputs, Docker
caches, Mamba environments, generated reports, or local agent handoff files.
The repository is intended to contain source, configuration, tests, and
documentation only.

## License

See [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff).

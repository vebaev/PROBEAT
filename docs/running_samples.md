# Running Samples

Use the Streamlit UI at `http://localhost:9797` for routine runs.

For command-line execution, run samples inside the worker container:

```bash
docker-compose run --rm worker \
  workflow/run_sample.sh SAMPLE PE /app/data/uploads/SAMPLE_R1.fastq.gz /app/data/uploads/SAMPLE_R2.fastq.gz --cores 16

docker-compose run --rm worker \
  workflow/run_sample.sh SAMPLE SE /app/data/uploads/SAMPLE.fastq.gz --cores 16

docker-compose run --rm worker \
  workflow/run_sample.sh SAMPLE ONT /app/data/uploads/SAMPLE.fastq.gz --cores 16
```

Run a dry-run smoke test:

```bash
docker-compose run --rm worker \
  workflow/run_sample.sh SMOKE_SE SE tests/data/smoke.fastq --dry-run --cores 2
```

Outputs are written to:

```text
storage/results/<sample>/
```

The final reports are:

```text
storage/results/<sample>/report/<sample>_report.html
storage/results/<sample>/report/<sample>_report.pdf
```

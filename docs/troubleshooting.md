# Troubleshooting

Start with:

```bash
./probeat doctor --profile full
./probeat status SAMPLE
docker-compose logs --tail 200 worker
```

Common failures are missing databases, unavailable NCBI/DNS access, insufficient
disk space, incompatible ONT Medaka models, or stale Snakemake incomplete
markers. Check active processes before unlocking Snakemake or clearing markers.

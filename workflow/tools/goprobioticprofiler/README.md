# GOProbioticProfiler

Genome-wide and probiotic-focused GO bubble profiler for bacterial genome analysis.

## Goal

Build a custom GO enrichment and GO-Slim profiler for the PROBEAT workflow.

## Input

1. `sample.emapper.annotations` – eggNOG-mapper annotations
2. `pgmdb_hits.tsv` – PGMDB probiotic marker results
3. `go-basic.obo` – GO ontology

## Output

1. `general_go_profile.tsv` and `general_go_bubble.png` – genome-wide GO landscape
2. `probiotic_go_terms.tsv` – deduplicated curated probiotic GO term evidence
3. `probiotic_go_categories.tsv` and `probiotic_go_bubble.png` – category coverage,
   unique gene support, enrichment, and Significant/Supportive/Detected evidence tiers
4. Existing enrichment, GO-Slim, profile, and JSON outputs for compatibility

## Installation

```bash
conda env create -f environment.yml
conda activate goprobioticprofiler
pip install -e .
```

## Usage

```bash
goprobioticprofiler \
  --eggnog sample.emapper.annotations \
  --pgmdb pgmdb_hits.tsv \
  --obo go-basic.obo \
  --outdir results \
  --verbose
```

## Integration

See `workflow/tools/goprobioticprofiler/Snakefile` for Snakemake integration.

## Performance

- Optimized for bacterial genomes with 2,000–10,000 CDS
- Uses GOATOOLS for fast GO enrichment
- Supports batch processing via Snakemake

## License

MIT

# Data Analysis Environment

This workspace is organized for reproducible data analysis with Snakemake as the main workflow manager.

## Operating Rules

1. Do not use `sudo`.
2. Do not modify files outside `/home/vebaev/codex/`.
3. Install software only through Apptainer containers or Mamba environments.
4. Do not update existing environments.
5. Create a new environment for each new tool.
6. Run all analyses through Snakemake.
7. Record commands and logs for every analysis.
8. Never delete raw FASTQ files.
9. Track workflow files with Git.
10. Use Snakemake as the primary workflow manager.
11. Use Mamba only for Snakemake and helper tools.
12. Use Apptainer for all bioinformatics applications.

## Layout

- `workflow/`: Snakemake workflow, configuration, environments, and profiles.
- `data/raw/`: immutable raw input data, including FASTQ files.
- `data/intermediate/`: generated intermediate files.
- `data/results/`: final machine-readable results.
- `app_containers/`: Apptainer image files and definitions.
- `logs/`: Snakemake and rule logs.
- `reports/`: human-readable reports.

## Standard Run Pattern

Create or activate a dedicated Mamba environment for Snakemake, then run:

```bash
workflow/run_snakemake.sh
```

Bioinformatics tools should be called from Apptainer containers inside Snakemake rules.

The local helper environment is expected at:

```bash
/home/vebaev/codex/mamba_envs/snakemake-helper-2026-06-04
```

It contains Snakemake, Mamba, Apptainer, Git, Python, and helper libraries only.

## Bacterial Genome Workflow

Add Illumina samples to:

```bash
workflow/config/samples.tsv
```

Format:

```text
sample	layout	r1	r2
sample_a	PE	data/raw/sample_a_R1.fastq.gz	data/raw/sample_a_R2.fastq.gz
sample_b	SE	data/raw/sample_b.fastq.gz	
```

Then run:

```bash
workflow/run_snakemake.sh
```

The workflow includes fastp, FastQC, Unicycler, QUAST, CheckM2, GTDB-Tk, Bakta, AMRFinderPlus, BAGEL4 placeholder, antiSMASH, MOB-suite, MultiQC, and final HTML/PDF reports.

Database-dependent tools expect databases under `data/databases/`; see `data/databases/README.md`.

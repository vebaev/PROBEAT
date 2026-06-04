# Apptainer Containers

Place Apptainer definition files and documented image build instructions here.

Do not install bioinformatics applications directly into the system or into the Snakemake helper environment. Each bioinformatics tool should have its own Apptainer image and be called from a Snakemake rule with `container:`.

Local `.sif` images are ignored by Git because they can be large. Keep reproducible definition files or build notes under version control.

## Pull configured containers

```bash
app_containers/pull_containers.sh
```

The configured public containers are listed in `app_containers/containers.tsv`.

BAGEL4 is marked as manual because no public Biocontainers image is configured here. To enable it, provide `app_containers/bagel4.manual.sif` and set `params.run_bagel4: true` in `workflow/config/config.yaml`.

## Server runtime note

On this server, Apptainer image execution currently fails because `/dev/fuse` is unavailable to the current user:

```text
fuse: device not found, try 'modprobe fuse' first
```

The containers have been pulled successfully, but running them requires a working Apptainer runtime with FUSE or a setuid/system Apptainer configuration. Do not use `sudo` from this workspace to change that. See `reports/server_runtime_and_storage_plan.md`.

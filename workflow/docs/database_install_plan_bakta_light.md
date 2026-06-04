# Database Install Plan: Bakta Light

Date: 2026-06-04

## Current blocker

Apptainer containers are downloaded, but cannot run because `/dev/fuse` is not available to user `vebaev`.

An attempted one-time `sudo` fix failed before running any FUSE command:

```text
sudo: /etc/sudo.conf is owned by uid 65534, should be 0
sudo: /usr/bin/sudo must be owned by uid 0 and have the setuid bit set
```

This means `sudo` itself is broken on the server. A system administrator must repair `sudo` and/or FUSE outside this workspace.

Suggested administrator-side commands:

```bash
chown root:root /usr/bin/sudo /etc/sudo.conf
chmod 4755 /usr/bin/sudo
modprobe fuse
test -e /dev/fuse || mknod -m 666 /dev/fuse c 10 229
chmod 666 /dev/fuse
```

After this, verify as `vebaev`:

```bash
workflow/scripts/check_apptainer_runtime.sh
```

## Database layout

Use only project-local database paths:

```text
data/databases/
├── checkm2/
├── gtdbtk/
├── bakta/
├── amrfinderplus/
├── antismash/
├── mob_suite/
└── bagel4/
```

Do not use `/tmp` for database downloads or extraction.

## Recommended installation order

1. Fix Apptainer/FUSE runtime.
2. Install small databases first:
   - AMRFinderPlus
   - MOB-suite
   - CheckM2
3. Install Bakta light.
4. Install antiSMASH databases or switch to a full antiSMASH container strategy.
5. Install GTDB-Tk last because it is the largest and most storage-sensitive.
6. Run a test sample through the workflow.

## Space plan

Use Bakta light by default.

| Component | Planning size |
| --- | ---: |
| Existing Apptainer SIF images | 5.3 GB |
| Existing Apptainer/cache data | 11 GB |
| GTDB-Tk R232 | 110-130 GB |
| Bakta light | 5-10 GB |
| CheckM2 | 5-10 GB |
| AMRFinderPlus | 1 GB |
| antiSMASH | 10-20 GB |
| MOB-suite | 2-5 GB |
| Intermediate/results reserve | 50-100 GB |

Recommended reserved space: 220-300 GB.

Current `/home/vebaev/codex` free space was about 3.0 TB, so storage is sufficient.

## Commands after FUSE is fixed

Pull containers again only if any image is missing:

```bash
app_containers/pull_containers.sh
```

Then install databases through containerized tools or documented project-local downloads.

Expected workflow run:

```bash
workflow/run_snakemake.sh
```

For cautious first execution:

```bash
workflow/run_snakemake.sh --cores 24 --resources mem_mb=96000
```

For normal server execution:

```bash
workflow/run_snakemake.sh
```

The local profile allows up to 48 CPU threads and 178 GB RAM total.

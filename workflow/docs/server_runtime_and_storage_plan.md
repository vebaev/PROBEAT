# Server Runtime and Storage Plan

Date: 2026-06-04

## Apptainer / FUSE status

Apptainer containers were pulled successfully into `app_containers/`, but execution currently fails on this server:

```text
fuse: device not found, try 'modprobe fuse' first
```

Observed locally:

- `/dev/fuse` is not available to user `vebaev`.
- `fuse`, `fusectl`, `squashfs`, and `overlay` appear in `/proc/filesystems`.
- `fuse-overlayfs`, `squashfuse`, and `squashfuse_ll` binaries exist.
- SIF execution, `--unsquash`, and sandbox execution all fail because Apptainer still needs FUSE/overlay support.

This cannot be fixed safely from `/home/vebaev/codex/` without `sudo` or system configuration changes.

Administrator-side options:

1. Provide a usable `/dev/fuse` device for user jobs.
2. Ensure user namespaces and FUSE mounts are allowed for Apptainer.
3. Install/configure a system Apptainer/Singularity runtime that can execute SIF images on this server.

Do not use `sudo` from this workspace.

## Current local space

- `/home/vebaev/codex`: 23 TB filesystem, about 3.0 TB free at the time of check.
- `/tmp`: about 8.5 GB free, too small for database work.
- Current pulled containers: about 5.3 GB.
- Current Apptainer/cache data: about 11 GB.

Use `/home/vebaev/codex/.tmp` and `data/databases/`, not `/tmp`, for large work.

## Database size estimates

Recommended planning target: reserve at least 250-300 GB for databases and 50-100 GB for intermediate/results for a modest batch of bacterial genomes.

| Tool | Database / data | Checked size | Planning size |
| --- | --- | ---: | ---: |
| GTDB-Tk 2.7.x | R232 reference data | 60.8 GB compressed; ~100 GB unarchived | 110-130 GB |
| Bakta 1.12 | full DB | 30 GB zipped; 84 GB unzipped | 90-100 GB |
| Bakta 1.12 | light DB | 1.3 GB zipped; 3.9 GB unzipped | 5-10 GB |
| CheckM2 1.1 | Zenodo reference DB | 1.735 GB compressed | 5-10 GB |
| AMRFinderPlus | latest NCBI DB listing | ~123 MB visible files | 1 GB |
| antiSMASH 8 | official standalone-full image includes DB | ~9 GB Docker download | 10-20 GB |
| MOB-suite 3.1.9 | Zenodo data.tar.gz | 473 MB compressed | 2-5 GB after setup |
| BAGEL4 | no configured public CLI container/DB | unknown | manual |

If Bakta light is acceptable, a practical database budget is about 140-170 GB. If Bakta full is used, plan for about 230-270 GB before analysis outputs.

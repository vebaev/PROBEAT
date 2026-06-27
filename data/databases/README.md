# PROBEAT Databases

Large biological databases are intentionally not tracked by Git.

For the documented full Docker installation, download all required databases
with:

```bash
docker-compose run --rm database-init
```

Verify the installed database set with:

```bash
docker-compose run --rm database-manager
```

By default, Docker Compose stores downloaded databases under:

```text
./storage/databases
```

The small curated PROBEAT interpretation files are tracked in this directory:

- `CURATED_DATABASES.md`
- `PGMDB.tab`
- `gut_brain_axis_marker_database.tsv`
- `gut_brain_axis_marker_database.modules.tsv`
- `probiotic_cazy.txt`

All other database directories and large downloaded files should remain
untracked.

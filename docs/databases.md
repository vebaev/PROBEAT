# Databases

Database requirements are declared in `manifests/databases.yaml`.

```bash
./probeat databases list
./probeat databases install
./probeat databases verify
./probeat databases versions
```

Large databases are not committed or included in Docker images. Database
installation always downloads all declared PROBEAT databases. Installation must
download into staging storage and activate data only after verification.
Curated PROBEAT tables require provenance and license review before release;
see `data/databases/CURATED_DATABASES.md`.

Downloads are cached under `/app/downloads` (host: `PROBEAT_DOWNLOADS_HOST`,
default `./storage/downloads`), and targets that are already installed are
skipped, so `./probeat databases install` can be rerun safely. When one target
fails, the remaining targets still run and the failed ones are listed at the end.

## eggNOG manual install

`download_eggnog_data.py` (eggNOG-mapper 2.1.15) downloads from
`http://eggnog5.embl.de`, which serves HTTP only. On networks that allow only
HTTPS the eggNOG target fails. Copy these files from `emapperdb-5.0.2` into
`storage/databases/eggnog/` (container path `/app/data/databases/eggnog/`):

- `eggnog.db`
- `eggnog.taxa.db`
- `eggnog_proteins.dmnd`

Then run `./probeat databases verify`; the installer skips eggNOG once these
files are present.

## geNomad database version

`workflow/envs/genomad.yaml` pins geNomad 1.12, which reads database v1.8 or
newer only ([compatibility table](https://portal.nersc.gov/genomad/__data__/releases.txt)).
If an older database (for example v1.7) is installed, rerun
`workflow/scripts/install_databases.sh genomad`; it replaces the old database.

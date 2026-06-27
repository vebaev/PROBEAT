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

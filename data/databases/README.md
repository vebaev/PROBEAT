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


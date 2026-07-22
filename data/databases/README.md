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

## Database Sources

Most workflow tools manage their required databases or models through their own
downloaders, including Bakta, antiSMASH, dbCAN, geNomad, eggNOG-mapper, BUSCO,
AMRFinderPlus, MOB-suite, and Medaka.

PROBEAT also downloads and uses the following external database resources:

### KOfamScan

- https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
- https://www.genome.jp/ftp/db/kofam/ko_list.gz

### Mash

- https://gembox.cbcb.umd.edu/mash/RefSeqSketchesDefaults.msh.gz

### sourmash

- https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs220/gtdb-reps-rs220-k31.dna.zip
- https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs220/gtdb-rs220.lineages.csv

### GO

- https://current.geneontology.org/ontology/go-basic.obo

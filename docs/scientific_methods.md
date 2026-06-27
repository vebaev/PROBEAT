# Scientific methods and reproducibility

PROBEAT records the workflow version, Git commit, input checksums, selected
profile, parameters, tool versions, and database versions for each run.
Tool/database declarations are maintained in `manifests/`.

The workflow supports Illumina PE/SE and ONT isolate reads. It performs read
QC and normalization, assembly, quality/taxonomic assessment, Bakta annotation,
safety/mobile-element analysis, functional/metabolic analysis, probiotic
marker analysis, and integrated reporting.

PGMDB, GO categories, gut-brain markers, and probiotic CAZyme outputs are
descriptive genomic-potential screens. They are not experimental validation.
Changes that alter biological results must pass and document the BV7/DPL5
scientific regression comparison.

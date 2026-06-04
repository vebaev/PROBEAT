# Apptainer Containers

Place Apptainer definition files and documented image build instructions here.

Do not install bioinformatics applications directly into the system or into the Snakemake helper environment. Each bioinformatics tool should have its own Apptainer image and be called from a Snakemake rule with `container:`.

Local `.sif` images are ignored by Git because they can be large. Keep reproducible definition files or build notes under version control.

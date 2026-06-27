# Output structure

All sample outputs are grouped below `data/results/<sample>/`.

- `run/`: immutable run metadata and CLI log.
- `qc/`, `assembly/`, `flye/`, `medaka/`: read and assembly processing.
- analysis-specific directories: machine-readable tables and figures.
- `report/`: canonical HTML/PDF report and timestamped report copies.

Raw reads, databases, generated outputs, and runtime state are excluded from
Git.

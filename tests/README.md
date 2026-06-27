# PROBEAT tests

- `unit/`: fast tests for CLI and scientific helper logic.
- `smoke/`: PE, SE, and ONT DAG checks using tiny FASTQ inputs.
- `regression/`: comparisons against the recorded BV7 and DPL5 baseline.
- `data/`: deliberately tiny, non-biological test inputs.

Run:

```bash
./probeat test smoke
tests/smoke/run.sh
python tests/regression/check_baseline.py
tests/docker/smoke.sh
```

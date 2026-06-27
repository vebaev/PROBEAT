# Development

```bash
./probeat test smoke
python tests/regression/check_baseline.py
docker-compose config
```

Keep helper scripts focused, argument-driven, non-destructive, and free of
absolute paths. Add unit tests for scientific parsing and interpretation logic.
Do not update the regression baseline merely to make a failing test pass;
review and explain the biological change first.

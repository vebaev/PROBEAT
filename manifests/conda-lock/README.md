# Conda lock files

Release lock files are generated for Linux x86-64 with:

```bash
workflow/scripts/lock_environments.sh
```

They must be regenerated and reviewed before a tagged release. Lock files are
not fabricated when dependency resolution or channel access is unavailable.

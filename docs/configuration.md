# Configuration

The documented deployment uses `.env.example` with the full Docker profile.

Copy it before running Compose:

```bash
cp .env.example .env
```

Default values:

```env
PROBEAT_IMAGE=probeat:full
PROBEAT_INSTALL_PROFILE=full
PROBEAT_WEB_PORT=9797
```

Runtime storage defaults to `./storage`:

```env
PROBEAT_UPLOADS_HOST=./storage/uploads
PROBEAT_DATABASES_HOST=./storage/databases
PROBEAT_RESULTS_HOST=./storage/results
PROBEAT_LOGS_HOST=./storage/logs
PROBEAT_STATE_HOST=./storage/state
```

The active workflow configuration is `workflow/config/config.yaml`. Public
workflow files should not contain workstation-specific absolute paths.

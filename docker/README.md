# PROBEAT Docker

PROBEAT is documented for the full Docker profile.

Build the full image:

```bash
docker build --build-arg INSTALL_PROFILE=full -t probeat:full .
```

Download databases into the mounted `./storage/databases` directory:

```bash
docker-compose run --rm database-init
```

Verify databases:

```bash
docker-compose run --rm database-manager
```

Start Redis, worker, and Streamlit:

```bash
docker-compose up -d redis worker streamlit
```

The Streamlit UI listens on:

```text
http://localhost:9797
```

Runtime data and downloaded databases stay outside the image under `./storage`.

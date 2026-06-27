FROM mambaorg/micromamba:2.0.5

ARG INSTALL_PROFILE=probiotic
ARG GUTSMASH_COMMIT=1651919f7e93d09a995ea080a3c029be84ee48c6
ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root
WORKDIR /app
RUN chown $MAMBA_USER:$MAMBA_USER /app \
    && mkdir -p /opt/probeat \
    && chown $MAMBA_USER:$MAMBA_USER /opt/probeat

# Copy patch file early for Layer 1
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/diamond_tmpdir.patch /app/docker/diamond_tmpdir.patch

# ── Layer 1: gutSMASH source (rarely changes, ~320MB) ──────────────
# Clone at pinned commit so the build is self-sufficient even when the
# host build context does not contain workflow/gutsmash_src (it is
# git-ignored and absent on fresh checkouts).  The entire gutsmash_src
# directory is excluded from the Docker build context via .dockerignore
# so neither the explicit COPY nor the final `COPY . /app` can clobber
# this clone with a partial / stale local copy.
USER $MAMBA_USER
RUN micromamba install -y -n base git curl patch \
    && git clone --quiet \
      https://github.com/victoriapascal/gutsmash.git \
      /tmp/gutsmash-clone \
    && git -C /tmp/gutsmash-clone checkout --quiet --detach "${GUTSMASH_COMMIT}" \
    && mkdir -p /app/workflow /opt/probeat \
    && cp -a /tmp/gutsmash-clone /opt/probeat/gutsmash_src \
    && ln -sfn /opt/probeat/gutsmash_src /app/workflow/gutsmash_src \
    && rm -rf /tmp/gutsmash-clone /opt/probeat/gutsmash_src/.git \
    && micromamba remove -y -n base git curl \
    && for f in $(grep -rln --include='*.py' 'list\[\|dict\[\|tuple\[\|set\[' /opt/probeat/gutsmash_src/antismash/ 2>/dev/null | grep -v __pycache__ | grep -v /test/); do \
       if [ -f "$f" ] && ! head -1 "$f" | grep -q 'from __future__ import annotations'; then \
         sed -i '1s/^/from __future__ import annotations\n/' "$f"; \
       fi; \
     done \
     && patch -p1 -i /app/docker/diamond_tmpdir.patch -d /opt/probeat/gutsmash_src

# ── Layer 2: env manifests + tool manifests (rarely change) ─────────
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs /app/workflow/envs
COPY --chown=$MAMBA_USER:$MAMBA_USER manifests /app/manifests

# ── Layer 3: entry points + micromamba wrapper (rarely change) ─────
COPY --chown=$MAMBA_USER:$MAMBA_USER probeat /app/probeat
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/probeat_cli.py /app/workflow/probeat_cli.py
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/run_snakemake.sh /app/workflow/run_snakemake.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/run_sample.sh /app/workflow/run_sample.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/install_databases.sh /app/workflow/scripts/install_databases.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/keggdecoder_kofam_fallback.py /app/workflow/scripts/keggdecoder_kofam_fallback.py
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/gutsmash_runtime.py /app/workflow/scripts/gutsmash_runtime.py
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/run_gutsmash_wrapper.py /app/workflow/scripts/run_gutsmash_wrapper.py
COPY --chown=$MAMBA_USER:$MAMBA_USER docker /app/docker
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/tools /app/workflow/tools
RUN chmod +x /app/probeat /app/workflow/run_snakemake.sh /app/workflow/run_sample.sh \
    /app/workflow/scripts/install_databases.sh \
    /app/workflow/scripts/keggdecoder_kofam_fallback.py \
    /app/workflow/scripts/run_gutsmash_wrapper.py \
    && mkdir -p /app/mamba_envs /app/snakemake_conda_envs /app/micromamba-root/pkgs \
       /app/data/databases /app/data/results /app/data/state /app/data/uploads /app/logs

# ── Layer 4: conda/micromamba compatibility shims ───────────────────
USER root
RUN cp /app/docker/micromamba_conda_wrapper.py /usr/local/bin/conda \
    && cp /app/docker/micromamba_conda_wrapper.py /usr/local/bin/mamba \
    && chmod 755 /usr/local/bin/conda /usr/local/bin/mamba \
    && mkdir -p /opt/conda/bin \
    && cp /app/docker/conda_activate_shim.sh /opt/conda/bin/activate \
    && chmod 755 /opt/conda/bin/activate

# ── Layer 5: environment setup (expensive — cached if Layers 1-4 unchanged) ─
USER $MAMBA_USER
RUN micromamba create -y \
      -p /app/mamba_envs/snakemake-conda-helper-2026-06-04 \
      -f /app/workflow/envs/snakemake-conda-helper.yaml \
    && /app/probeat setup --profile "${INSTALL_PROFILE}" \
    && if [ "${INSTALL_PROFILE}" = "full" ] || [ "${INSTALL_PROFILE}" = "probiotic" ] || [ "${INSTALL_PROFILE}" = "ont" ]; then \
         test -x /app/mamba_envs/gutsmash/bin/python; \
         test -s /opt/probeat/gutsmash_src/run_gutsmash.py; \
         test -s /opt/probeat/gutsmash_src/antismash/__main__.py; \
         PROBEAT_ROOT=/app PROBEAT_GUTSMASH_SRC=/opt/probeat/gutsmash_src \
           /app/mamba_envs/snakemake-conda-helper-2026-06-04/bin/python /app/workflow/scripts/run_gutsmash_wrapper.py --self-test; \
       fi \
    && micromamba clean --all --yes \
    && find /app/mamba_envs -type d \( -name "__pycache__" -o -name ".pytest_cache" \) -prune -exec rm -rf {} + \
    && find /app/mamba_envs -type f \( -name "*.pyc" -o -name "*.pyo" -o -name "*.a" -o -name "*.la" \) -delete \
    && find /app/mamba_envs -type d \( -path "*/share/doc" -o -path "*/share/man" -o -path "*/share/info" \) -prune -exec rm -rf {} + \
    && rm -rf /app/.cache /app/.conda /app/.mamba /app/micromamba-root/pkgs/*

# ── Layer 6: full source copy (changes often, but after expensive setup) ───
# gutsmash_src is excluded from the build context via .dockerignore so
# the clone from Layer 1 survives this COPY intact.
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
RUN chmod +x /app/probeat /app/workflow/run_snakemake.sh /app/workflow/run_sample.sh \
    /app/workflow/scripts/install_databases.sh \
    /app/workflow/scripts/keggdecoder_kofam_fallback.py \
    /app/workflow/scripts/run_gutsmash_wrapper.py \
    && if [ ! -e /app/workflow/gutsmash_src ]; then \
         ln -sfn /opt/probeat/gutsmash_src /app/workflow/gutsmash_src; \
       fi

ENV PROBEAT_ROOT=/app \
    PROBEAT_GUTSMASH_SRC=/opt/probeat/gutsmash_src \
    PROBEAT_HELPER_ENV=/app/mamba_envs/snakemake-conda-helper-2026-06-04 \
    PROBEAT_DATABASES=/app/data/databases \
    PROBEAT_RESULTS=/app/data/results \
    PROBEAT_STATE=/app/data/state \
    PROBEAT_UPLOADS=/app/data/uploads \
    PROBEAT_SNAKEMAKE_CONDA_ENVS=/app/snakemake_conda_envs \
    PROBEAT_MICROMAMBA_PKGS=/app/micromamba-root/pkgs \
    PROBEAT_WORKFLOW_VERSION=1.0 \
    PATH=/app/mamba_envs/snakemake-conda-helper-2026-06-04/bin:$PATH

EXPOSE 9797
CMD ["streamlit", "run", "web/streamlit_app/app.py", "--server.address=0.0.0.0", "--server.port=9797", "--server.maxUploadSize=1024"]

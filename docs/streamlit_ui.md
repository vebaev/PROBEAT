# Streamlit UI

The UI submits jobs to Redis and never runs Snakemake inside the Streamlit
process. The worker calls the same `./probeat run` command used by local runs.

The default Docker Compose port is `9797`.

Pages:

- Dashboard: job and report overview, status counts, discovered completed
  reports, compact sample metrics, live activity for running jobs, and the
  per-sample progress/log panel. The dashboard header includes the PROBEAT
  workflow logo and the sidebar shows a rotating running indicator while a job
  is queued or running.
- New Analysis: upload PE, SE, or ONT reads, select profile, CPU/RAM limits,
  optionally run strict preflight, then submit the job to Redis.
- Results: completed-sample summary cards, HTML preview, PDF/HTML downloads,
  and ZIP download of the complete `results/<sample>/` directory.
- System Health: simplified Redis, storage, database, and `probeat doctor`
  checks with expandable technical details.

FASTQ uploads are streamed to disk in chunks under the configured uploads
volume. The UI also discovers already completed reports under the configured
results volume, so it remains useful after container restarts.

Log viewing is intentionally bounded to recent lines. Use the result directory
and worker log files for deeper debugging rather than adding unbounded live log
streams to the UI.

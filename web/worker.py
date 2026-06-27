#!/usr/bin/env python3
"""Redis-backed single-worker executor for PROBEAT jobs."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import time
import re
import sys

import redis


ROOT = Path(os.environ.get("PROBEAT_ROOT", "/app"))
SCRIPTS = ROOT / "workflow" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from gutsmash_runtime import profile_requires_gutsmash_name, runtime_errors as shared_gutsmash_runtime_errors  # noqa: E402

STATE = Path(os.environ.get("PROBEAT_STATE", "/data/state"))
REDIS_URL = os.environ.get("REDIS_URL", "redis://redis:6379/0")
QUEUE = "probeat:jobs"
WRITABLE_PATHS = {
    "results": Path(os.environ.get("PROBEAT_RESULTS", ROOT / "data" / "results")),
    "state": STATE,
    "uploads": Path(os.environ.get("PROBEAT_UPLOADS", ROOT / "data" / "uploads")),
    "logs": ROOT / "logs",
    "snakemake conda cache": Path(os.environ.get("PROBEAT_SNAKEMAKE_CONDA_ENVS", ROOT / "snakemake_conda_envs")),
    "micromamba package cache": Path(os.environ.get("PROBEAT_MICROMAMBA_PKGS", ROOT / "micromamba-root" / "pkgs")),
}


def update(client: redis.Redis, job: dict, **changes) -> dict:
    job.update(changes)
    job["updated_at"] = time.time()
    client.hset("probeat:job-state", job["sample"], json.dumps(job))
    path = STATE / "jobs" / f"{job['sample']}.json"
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(job, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    except OSError as exc:
        print(f"[WARN] Failed to mirror job state for {job['sample']} to {path}: {exc}", file=sys.stderr)
    return job


def recent_run_lines(sample: str) -> list[str]:
    path = ROOT / "data" / "results" / sample / "run" / "probeat.log"
    try:
        return path.read_text(encoding="utf-8", errors="replace").splitlines()[-160:]
    except OSError:
        return []


def environment_setup_step(lines: list[str]) -> str:
    created = [
        match.group(1)
        for line in lines
        if (match := re.search(r"Environment for /app/workflow/envs/([A-Za-z0-9_.-]+)\.yaml created", line))
    ]
    creating = [
        match.group(1)
        for line in lines
        if (match := re.search(r"Creating conda environment workflow/envs/([A-Za-z0-9_.-]+)\.yaml", line))
    ]
    if not creating:
        return ""
    if len(creating) > len(created) or creating[-1] not in created:
        return f"Preparing tool environments: {creating[-1]}"
    return ""


def writable_mount_errors() -> list[str]:
    errors: list[str] = []
    for label, path in WRITABLE_PATHS.items():
        try:
            path.mkdir(parents=True, exist_ok=True)
            probe = path / ".probeat_write_test"
            probe.write_text("ok\n", encoding="utf-8")
            probe.unlink()
        except OSError as exc:
            errors.append(f"{label}: {path} ({exc})")
    return errors


def gutsmash_runtime_errors(profile: str) -> list[str]:
    if not profile_requires_gutsmash_name(profile):
        return []
    return shared_gutsmash_runtime_errors(ROOT)


def main() -> None:
    client = redis.Redis.from_url(REDIS_URL, decode_responses=True)
    while True:
        try:
            queued = client.brpop(QUEUE, timeout=10)
        except redis.RedisError:
            time.sleep(5)
            client = redis.Redis.from_url(REDIS_URL, decode_responses=True)
            continue
        if not queued:
            continue
        _, payload = queued
        job = json.loads(payload)
        current = client.hget("probeat:job-state", job["sample"])
        if current and json.loads(current).get("status") == "cancelled":
            continue
        command = [
            str(ROOT / "probeat"),
            "run",
            job["sample"],
            job["layout"],
            *job["reads"],
            "--profile",
            job["profile"],
            "--cores",
            str(job["cores"]),
            "--memory-gb",
            str(job["memory_gb"]),
            "--no-doctor",
        ]
        update(client, job, status="running", command=command, started_at=time.time())
        log_path = STATE / "worker-logs" / f"{job['sample']}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        mount_errors = writable_mount_errors()
        gut_errors = gutsmash_runtime_errors(job.get("profile", ""))
        preflight_errors = []
        if mount_errors:
            preflight_errors.append("Writable Docker mount check failed:\n" + "\n".join(f"- {error}" for error in mount_errors))
        if gut_errors:
            preflight_errors.append("gutSMASH runtime check failed:\n" + "\n".join(f"- {error}" for error in gut_errors))
        if preflight_errors:
            message = "\n\n".join(preflight_errors)
            with log_path.open("a", encoding="utf-8") as log:
                log.write(message + "\n")
            update(
                client,
                job,
                status="failed_preflight",
                exit_code=2,
                finished_at=time.time(),
                error=message,
                log=str(log_path),
            )
            continue
        with log_path.open("a", encoding="utf-8") as log:
            process = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
            update(client, job, pid=process.pid, log=str(log_path))
            while process.poll() is None:
                if client.get(f"probeat:cancel:{job['sample']}") == "1":
                    process.terminate()
                    try:
                        process.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        process.kill()
                    client.delete(f"probeat:cancel:{job['sample']}")
                    update(client, job, status="cancelled")
                    break
                try:
                    recent = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-80:]
                    recent += recent_run_lines(job["sample"])
                    rules = [match.group(1) for line in recent if (match := re.search(r"rule ([A-Za-z0-9_]+):", line))]
                    progress = [match.groups() for line in recent if (match := re.search(r"(\d+) of (\d+) steps \((\d+)%\) done", line))]
                    changes = {}
                    environment = environment_setup_step(recent)
                    if environment:
                        changes["current_step"] = environment
                    elif rules:
                        changes["current_step"] = rules[-1]
                    if progress:
                        changes["completed_steps"], changes["total_steps"], changes["progress_percent"] = map(int, progress[-1])
                    if changes:
                        update(client, job, **changes)
                except OSError:
                    pass
                time.sleep(2)
            return_code = process.wait()
        if job.get("status") == "cancelled":
            continue
        update(
            client,
            job,
            status="completed" if return_code == 0 else "failed",
            exit_code=return_code,
            finished_at=time.time(),
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Stable command-line entry point for local, CI, and container PROBEAT runs."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

try:
    import yaml
except ImportError as exc:
    raise SystemExit("PyYAML is required. Run: ./probeat setup --profile core") from exc


ROOT = Path(__file__).resolve().parents[1]
MANIFESTS = ROOT / "manifests"
STATE = Path(os.environ.get("PROBEAT_STATE", ROOT / "data" / "state"))
RESULTS = Path(os.environ.get("PROBEAT_RESULTS", ROOT / "data" / "results"))
DATABASES = Path(os.environ.get("PROBEAT_DATABASES", ROOT / "data" / "databases"))
UPLOADS = Path(os.environ.get("PROBEAT_UPLOADS", ROOT / "data" / "uploads"))
SNAKEMAKE_CONDA_ENVS = Path(os.environ.get("PROBEAT_SNAKEMAKE_CONDA_ENVS", ROOT / "snakemake_conda_envs"))
MICROMAMBA_PKGS = Path(os.environ.get("PROBEAT_MICROMAMBA_PKGS", ROOT / "micromamba-root" / "pkgs"))
VERSION = "0.1.0-dev"
SCRIPTS = ROOT / "workflow" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from gutsmash_runtime import (  # noqa: E402
    profile_requires_gutsmash_config,
    runtime_errors as shared_gutsmash_runtime_errors,
    source_candidates as shared_gutsmash_source_candidates,
    valid_source as valid_gutsmash_source,
)


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def writable_directory(path: Path) -> tuple[bool, str]:
    try:
        path.mkdir(parents=True, exist_ok=True)
        probe = path / ".probeat_write_test"
        probe.write_text("ok\n", encoding="utf-8")
        probe.unlink()
        return True, ""
    except OSError as exc:
        return False, str(exc)


def git_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "--short=12", "HEAD"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return result.stdout.strip() or "unknown"


def profile_config(name: str) -> dict[str, Any]:
    release = load_yaml(MANIFESTS / "release.yaml")
    profiles = release.get("profiles", {})
    if name not in profiles:
        raise SystemExit(f"Unknown profile: {name}. Choose from: {', '.join(sorted(profiles))}")
    selected = dict(profiles[name])
    parent_name = selected.pop("inherits", None)
    if not parent_name:
        return selected
    parent = dict(profiles[parent_name])
    for key, value in selected.items():
        if isinstance(value, list) and isinstance(parent.get(key), list):
            parent[key] = list(dict.fromkeys(parent[key] + value))
        else:
            parent[key] = value
    return parent


GENERIC_EXECUTABLES = {"python", "python3", "perl", "bash", "sh", "java", "Rscript"}


def tool_path(tool: dict[str, Any]) -> str | None:
    candidates = tool.get("candidates", [])
    executable = tool.get("executable", tool["name"])
    for candidate in candidates:
        path = ROOT / candidate
        if path.is_file() and os.access(path, os.X_OK):
            return str(path)
    if candidates and executable in GENERIC_EXECUTABLES:
        # An interpreter from any other env would match; only the declared env is this tool.
        return None
    direct = ROOT / "mamba_envs" / tool["name"] / "bin" / executable
    if direct.is_file() and os.access(direct, os.X_OK):
        return str(direct)
    for path in (ROOT / "snakemake_conda_envs").glob(f"*/bin/{executable}"):
        if path.is_file() and os.access(path, os.X_OK):
            return str(path)
    return shutil.which(executable)


def database_path(item: dict[str, Any]) -> Path:
    return DATABASES / item["path"]


def database_valid(item: dict[str, Any]) -> tuple[bool, list[str]]:
    base = database_path(item)
    missing: list[str] = []
    if not base.exists():
        return False, [str(base)]
    for pattern in item.get("required_files", []):
        matches = list(base.glob(pattern))
        if not matches:
            missing.append(str(base / pattern))
            continue
        unreadable = [str(path) for path in matches if not os.access(path, os.R_OK)]
        if unreadable:
            missing.extend(f"unreadable: {path}" for path in unreadable)
    return not missing, missing


def profile_requires_gutsmash(profile: dict[str, Any]) -> bool:
    return "gutsmash" in profile.get("tools", []) or profile_requires_gutsmash_config(profile)


def gutsmash_source_candidates() -> list[Path]:
    return shared_gutsmash_source_candidates(ROOT)


def gutsmash_runtime_errors() -> list[str]:
    return shared_gutsmash_runtime_errors(ROOT)


def cmd_setup(args: argparse.Namespace) -> int:
    profile = profile_config(args.profile)
    bundled_micromamba = ROOT / ".tools" / "bin" / "micromamba"
    manager = shutil.which("micromamba") or shutil.which("mamba") or (str(bundled_micromamba) if bundled_micromamba.exists() else None)
    if not manager:
        raise SystemExit("micromamba or mamba is required. Run: workflow/bootstrap_micromamba.sh")
    env_root = ROOT / "mamba_envs"
    env_root.mkdir(exist_ok=True)
    target_aliases = {
        "snakemake-conda-helper": "snakemake-conda-helper-2026-06-04",
        "porechop": "porechop-0.2.4",
        "rasusa": "rasusa-2.1.0",
        "flye": "flye-2.9.6",
        "busco": "busco-6.0.0",
        "abricate": "abricate-1.0.1",
        "kofamscan": "kofamscan-1.3.0",
        "keggdecoder": "keggdecoder-1.3",
        "dbcan": "dbcan_v5",
    }
    for env_name in profile.get("environments", []):
        definition = ROOT / "workflow" / "envs" / f"{env_name}.yaml"
        target = env_root / target_aliases.get(env_name, env_name)
        if not definition.exists():
            raise SystemExit(f"Missing environment definition: {definition}")
        if target.exists() and not args.force:
            print(f"[OK] {env_name}: already exists")
            continue
        command = [manager, "create", "-y", "-p", str(target), "-f", str(definition)]
        print(f"[CREATE] {env_name}")
        for attempt in range(1, 4):
            try:
                subprocess.run(command, cwd=ROOT, check=True)
                break
            except subprocess.CalledProcessError:
                if attempt >= 3:
                    raise
                print(f"[WARN] {env_name}: creation failed; cleaning partial environment and retrying ({attempt}/3)")
                shutil.rmtree(target, ignore_errors=True)
                subprocess.run([manager, "clean", "--all", "--yes"], cwd=ROOT, check=False)
    goprobiotic_env = env_root / ("functional-python" if "functional-python" in profile.get("environments", []) else "goprobioticprofiler")
    if (
        ("functional-python" in profile.get("environments", []) or "goprobioticprofiler" in profile.get("environments", []))
        and (goprobiotic_env / "bin" / "pip").exists()
    ):
        subprocess.run(
            [str(goprobiotic_env / "bin" / "pip"), "install", "--no-deps", "-e", str(ROOT / "workflow" / "tools" / "goprobioticprofiler")],
            cwd=ROOT,
            check=True,
        )
    if "gutsmash" in profile.get("environments", []):
        source_candidates = [
            Path(os.environ["PROBEAT_GUTSMASH_SRC"]) if os.environ.get("PROBEAT_GUTSMASH_SRC") else None,
            ROOT / "workflow" / "gutsmash_src",
            Path("/opt/probeat/gutsmash_src"),
        ]
        commit = load_yaml(MANIFESTS / "tools.lock.yaml")["git_pins"]["gutsmash"]
        source = next(
            (
                candidate
                for candidate in source_candidates
                if candidate
                and (candidate / "run_gutsmash.py").exists()
                and (candidate / "run_gutsmash.py").stat().st_size > 0
                and (candidate / "antismash" / "__main__.py").exists()
                and (candidate / "antismash" / "__main__.py").stat().st_size > 0
            ),
            ROOT / "workflow" / "gutsmash_src",
        )
        required_source_files = [source / "run_gutsmash.py", source / "antismash" / "__main__.py"]
        bundled_ok = all(path.exists() and path.stat().st_size > 0 for path in required_source_files)
        if bundled_ok:
            print(f"[OK] gutsmash source: bundled copy present ({commit[:8]})")
        elif (source / ".git").exists():
            subprocess.run(["git", "-C", str(source), "checkout", "--detach", commit], check=True)
        elif not source.exists() or not any(source.iterdir()):
            source.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(["git", "clone", "-b", "gutsmash", "https://github.com/victoriapascal/gutsmash.git", str(source)], check=True)
            subprocess.run(["git", "-C", str(source), "checkout", "--detach", commit], check=True)
        else:
            missing = ", ".join(str(path.relative_to(ROOT)) for path in required_source_files if not path.exists())
            raise SystemExit(f"Invalid gutSMASH source at {source}; missing: {missing}")
    print(f"Profile '{args.profile}' is installed.")
    return 0


def cmd_databases(args: argparse.Namespace) -> int:
    manifest = load_yaml(MANIFESTS / "databases.yaml")
    entries = manifest.get("databases", [])
    selected = entries
    if args.database_command == "list":
        for item in entries:
            valid, _ = database_valid(item)
            print(f"{'[OK]' if valid else '[--]'} {item['id']:<24} {item.get('version', 'unknown'):<12} {item['path']}")
        return 0
    if args.database_command == "versions":
        for item in entries:
            valid, _ = database_valid(item)
            print(f"{item['id']}\t{item.get('version', 'unknown')}\t{'installed' if valid else 'missing'}")
        return 0
    if args.database_command == "verify":
        failed = 0
        for item in selected:
            valid, missing = database_valid(item)
            print(f"{'[OK]' if valid else '[FAIL]'} {item['id']}")
            for path in missing:
                print(f"       missing: {path}")
            failed += int(not valid)
        return 1 if failed else 0
    if args.database_command == "install":
        release = load_yaml(MANIFESTS / "release.yaml")
        targets = release.get("database_install_targets") or release.get("profiles", {}).get("full", {}).get("database_install_targets", [])
        if not targets:
            print("No automatic database targets are defined.")
            return 0
        print("Installing all declared PROBEAT databases.")
        command = [str(ROOT / "workflow" / "scripts" / "install_databases.sh"), *targets]
        return subprocess.run(command, cwd=ROOT).returncode
    raise SystemExit("Unsupported database command")


def cmd_doctor(args: argparse.Namespace) -> int:
    profile = profile_config(args.profile)
    tools = {item["name"]: item for item in load_yaml(MANIFESTS / "tools.lock.yaml").get("tools", [])}
    databases = {item["id"]: item for item in load_yaml(MANIFESTS / "databases.yaml").get("databases", [])}
    failures = 0
    print(f"PROBEAT doctor {VERSION} | profile={args.profile} | root={ROOT}")
    config_path = ROOT / "workflow" / "config" / "config.yaml"
    schema_path = ROOT / "workflow" / "schemas" / "config.schema.yaml"
    try:
        import jsonschema

        jsonschema.validate(load_yaml(config_path), load_yaml(schema_path))
        print(f"[OK] config schema: {config_path}")
    except ImportError:
        print("[WARN] jsonschema is unavailable; configuration schema was not checked")
    except Exception as exc:
        print(f"[FAIL] config schema: {exc}")
        failures += 1
    writable_paths = {
        "results": RESULTS,
        "state": STATE,
        "uploads": UPLOADS,
        "logs": ROOT / "logs",
        "snakemake conda cache": SNAKEMAKE_CONDA_ENVS,
        "micromamba package cache": MICROMAMBA_PKGS,
    }
    for label, directory in writable_paths.items():
        writable, error = writable_directory(directory)
        print(f"{'[OK]' if writable else '[FAIL]'} writable {label}: {directory}")
        if error:
            print(f"       {error}")
        failures += int(not writable)
    if not DATABASES.exists():
        try:
            DATABASES.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
    readable = DATABASES.is_dir() and os.access(DATABASES, os.R_OK | os.X_OK)
    print(f"{'[OK]' if readable else '[FAIL]'} readable: {DATABASES}")
    failures += int(not readable)
    free_gb = shutil.disk_usage(ROOT).free / 1024**3
    print(f"{'[OK]' if free_gb >= args.min_free_gb else '[WARN]'} free disk: {free_gb:.1f} GB")
    if profile_requires_gutsmash(profile):
        gut_errors = gutsmash_runtime_errors()
        print(f"{'[OK]' if not gut_errors else '[FAIL]'} gutSMASH runtime")
        for error in gut_errors:
            print(f"       {error}")
        failures += int(bool(gut_errors))
    for name in profile.get("tools", []):
        item = tools[name]
        path = tool_path(item)
        details: list[str] = []
        ok = path is not None
        if path:
            try:
                check_env = {
                    **os.environ,
                    "PATH": f"{Path(path).parent}:{os.environ.get('PATH', '')}",
                    "MPLCONFIGDIR": str(ROOT / ".cache" / "matplotlib"),
                }
                for key, value in item.get("check_env", {}).items():
                    check_env[key] = str(ROOT / value) if value.startswith("./") else value
                result = subprocess.run(
                    [path, *item.get("check_args", ["--version"])],
                    text=True,
                    capture_output=True,
                    timeout=15,
                    check=False,
                    cwd=ROOT,
                    env=check_env,
                )
                version_line = (result.stdout or result.stderr).strip().splitlines()
                if version_line:
                    details.append(f"version: {version_line[0][:160]}")
                if result.returncode and "check_args" in item:
                    details.append(f"check failed with exit code {result.returncode}")
                    ok = False
            except (OSError, subprocess.TimeoutExpired):
                details.append("version: could not query")
                if "check_args" in item:
                    ok = False
        print(f"{'[OK]' if ok else '[FAIL]'} tool: {name}{f' ({path})' if path else ''}")
        for detail in details:
            print(f"       {detail}")
        failures += int(not ok)
    for name in profile.get("databases", []):
        item = databases[name]
        valid, _ = database_valid(item)
        print(f"{'[OK]' if valid else '[FAIL]'} database: {name}")
        failures += int(not valid)
    if failures and args.strict:
        print(f"Doctor found {failures} blocking problem(s).")
        return 1
    print(f"Doctor completed with {failures} problem(s).")
    return 0


def validate_sample(sample: str, layout: str, reads: list[str]) -> list[Path]:
    if not sample or any(char not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-" for char in sample):
        raise SystemExit("Sample names may contain only letters, numbers, dots, underscores, and hyphens.")
    expected = 3 if layout == "HYBRID" else 2 if layout == "PE" else 1
    if len(reads) != expected:
        raise SystemExit(f"{layout} requires {expected} FASTQ input file(s).")
    paths = [Path(value).expanduser().resolve() for value in reads]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise SystemExit("Missing FASTQ input(s): " + ", ".join(missing))
    return paths


def cmd_run(args: argparse.Namespace) -> int:
    layout = args.layout.upper()
    reads = validate_sample(args.sample, layout, args.reads)
    job_path = STATE / "jobs" / f"{args.sample}.json"
    run_dir = RESULTS / args.sample / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    command = [str(ROOT / "workflow" / "run_sample.sh"), args.sample, layout, *map(str, reads), "--cores", str(args.cores)]
    if args.memory_gb:
        command.extend(["--resources", f"mem_mb={args.memory_gb * 1024}"])
    if args.dry_run:
        command.append("--dry-run")
    command.extend(args.snakemake_args)
    log_path = run_dir / "probeat.log"
    manifest = {
        "schema_version": 1,
        "run_id": run_id,
        "sample": args.sample,
        "layout": layout,
        "profile": args.profile,
        "workflow_version": VERSION,
        "git_commit": git_commit(),
        "configuration_sha256": checksum(ROOT / "workflow" / "config" / "config.yaml"),
        "tools_manifest_sha256": checksum(MANIFESTS / "tools.lock.yaml"),
        "databases_manifest_sha256": checksum(MANIFESTS / "databases.yaml"),
        "started_at": started,
        "inputs": [{"path": str(path), "sha256": checksum(path), "bytes": path.stat().st_size} for path in reads],
        "cores": args.cores,
        "memory_gb": args.memory_gb,
        "command": command,
        "log": str(log_path),
    }
    manifest_path = run_dir / "manifests" / f"{run_id}.json"
    write_json(manifest_path, manifest)
    if not (run_dir / "run_manifest.json").exists():
        write_json(run_dir / "run_manifest.json", manifest)
    state = {**manifest, "status": "running"}
    write_json(job_path, state)
    if not args.no_doctor:
        doctor_args = argparse.Namespace(profile=args.profile, strict=True, min_free_gb=1.0)
        if cmd_doctor(doctor_args):
            state["status"] = "failed_preflight"
            write_json(job_path, state)
            write_json(run_dir / "execution_status.json", state)
            write_json(run_dir / "manifests" / f"{run_id}.status.json", state)
            return 2
    run_env = {**os.environ, "PROBEAT_ANALYSIS_PROFILE": args.profile}
    with log_path.open("a", encoding="utf-8") as log:
        process = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, env=run_env)
    state["finished_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    state["exit_code"] = process.returncode
    state["status"] = "completed" if process.returncode == 0 else "failed"
    write_json(run_dir / "execution_status.json", state)
    write_json(run_dir / "manifests" / f"{run_id}.status.json", state)
    write_json(job_path, state)
    print(f"{args.sample}: {state['status']} (log: {log_path})")
    return process.returncode


def cmd_status(args: argparse.Namespace) -> int:
    path = STATE / "jobs" / f"{args.sample}.json"
    if not path.exists():
        report = RESULTS / args.sample / "report" / f"{args.sample}_report.pdf"
        print(json.dumps({"sample": args.sample, "status": "completed" if report.exists() else "unknown"}, indent=2))
        return 0
    print(path.read_text(encoding="utf-8"), end="")
    return 0


def cmd_test(args: argparse.Namespace) -> int:
    if args.suite != "smoke":
        raise SystemExit("Only the smoke suite is currently supported.")
    unit = [sys.executable, "-m", "unittest", "discover", "-s", "tests/unit", "-v"]
    if subprocess.run(unit, cwd=ROOT).returncode:
        return 1
    cases = [
        ["SMOKE_SE", "SE", "tests/data/smoke.fastq"],
        ["SMOKE_PE", "PE", "tests/data/smoke_R1.fastq", "tests/data/smoke_R2.fastq"],
        ["SMOKE_ONT", "ONT", "tests/data/smoke.fastq"],
        ["SMOKE_HYBRID", "HYBRID", "tests/data/smoke_R1.fastq", "tests/data/smoke_R2.fastq", "tests/data/smoke.fastq"],
    ]
    for case in cases:
        command = [str(ROOT / "probeat"), "run", *case, "--cores", "1", "--dry-run", "--no-doctor"]
        if subprocess.run(command, cwd=ROOT).returncode:
            return 1
    return 0


def parser() -> argparse.ArgumentParser:
    main = argparse.ArgumentParser(prog="probeat", description="PROBEAT reproducible bacterial genome analysis")
    main.add_argument("--version", action="version", version=VERSION)
    sub = main.add_subparsers(dest="command", required=True)

    setup = sub.add_parser("setup", help="Install isolated Micromamba environments")
    setup.add_argument("--profile", choices=("core", "full", "probiotic", "ont"), default="core")
    setup.add_argument("--force", action="store_true")
    setup.set_defaults(func=cmd_setup)

    db = sub.add_parser("databases", help="Install and verify workflow databases")
    db_sub = db.add_subparsers(dest="database_command", required=True)
    db_sub.add_parser("list").set_defaults(func=cmd_databases)
    db_sub.add_parser("versions").set_defaults(func=cmd_databases)
    for name in ("install", "verify"):
        command = db_sub.add_parser(name)
        command.add_argument("--profile", choices=("core", "full", "probiotic", "ont"), default=None, help=argparse.SUPPRESS)
        command.set_defaults(func=cmd_databases)

    doctor = sub.add_parser("doctor", help="Check tools, databases, storage, and permissions")
    doctor.add_argument("--profile", choices=("core", "full", "probiotic", "ont"), default="core")
    doctor.add_argument("--strict", action="store_true")
    doctor.add_argument("--min-free-gb", type=float, default=10.0)
    doctor.set_defaults(func=cmd_doctor)

    run = sub.add_parser("run", help="Run exactly one PE, SE, ONT, or HYBRID sample")
    run.add_argument("sample")
    run.add_argument("layout", choices=("PE", "SE", "ONT", "HYBRID"))
    run.add_argument("reads", nargs="+")
    run.add_argument("--profile", choices=("core", "full", "probiotic", "ont"), default="full")
    run.add_argument("--cores", type=int, default=max(1, os.cpu_count() or 1))
    run.add_argument("--memory-gb", type=int)
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--no-doctor", action="store_true")
    run.add_argument("--snakemake-args", nargs=argparse.REMAINDER, default=[])
    run.set_defaults(func=cmd_run)

    status = sub.add_parser("status", help="Show the stored state for one sample")
    status.add_argument("sample")
    status.set_defaults(func=cmd_status)

    test = sub.add_parser("test", help="Run repository tests")
    test.add_argument("suite", choices=("smoke",))
    test.set_defaults(func=cmd_test)
    return main


if __name__ == "__main__":
    os.chdir(ROOT)
    args = parser().parse_args()
    raise SystemExit(args.func(args))

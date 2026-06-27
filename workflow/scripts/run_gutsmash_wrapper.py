#!/usr/bin/env python3
"""Run gutSMASH through the PROBEAT runtime contract."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable, List

from gutsmash_runtime import resolve_python, resolve_source, runtime_environment, runtime_errors

HMMPRESS_SUFFIXES = ("h3f", "h3i", "h3m", "h3p")


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def hmm_sidecars(path: Path) -> List[Path]:
    return [path.with_name(f"{path.name}.{suffix}") for suffix in HMMPRESS_SUFFIXES]


def hmm_needs_press(path: Path) -> bool:
    if not nonempty(path):
        return False
    try:
        hmm_mtime = path.stat().st_mtime
        return any(not sidecar.is_file() or sidecar.stat().st_mtime < hmm_mtime for sidecar in hmm_sidecars(path))
    except OSError:
        return True


def resolve_env_executable(python_path: Path, name: str) -> Path | None:
    candidate = python_path.parent / name
    if candidate.is_file() and os.access(str(candidate), os.X_OK):
        return candidate
    resolved = shutil.which(name)
    return Path(resolved) if resolved else None


def press_hmm(path: Path, python_path: Path, env: dict[str, str], log_path: Path | None = None) -> None:
    if not nonempty(path) or not hmm_needs_press(path):
        return
    hmmpress = resolve_env_executable(python_path, "hmmpress")
    if not hmmpress:
        raise RuntimeError(f"Missing hmmpress executable needed to prepare {path}")
    result = subprocess.run(
        [str(hmmpress), "-f", str(path)],
        cwd=path.parent,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        message = (result.stderr or result.stdout or "no output").strip()
        raise RuntimeError(f"hmmpress failed for {path}: {message}")
    if log_path:
        write_log(log_path, [f"hmmpress prepared {path}"])


def prepare_bundled_hmms(source: Path, python_path: Path, env: dict[str, str], log_path: Path | None = None) -> None:
    """Prepare HMM databases bundled inside the staged gutSMASH source tree."""
    candidates = [
        source / "antismash" / "detection" / "genefunctions" / "gut_data" / "smcogs.hmm",
        source / "antismash" / "detection" / "genefunctions" / "data" / "smcogs.hmm",
    ]
    for candidate in candidates:
        press_hmm(candidate, python_path, env, log_path)


def verify_hmmscan_can_read(hmm_path: Path, python_path: Path, env: dict[str, str], work_dir: Path) -> None:
    if not nonempty(hmm_path):
        raise RuntimeError(f"Missing HMM database for hmmscan check: {hmm_path}")
    hmmscan = resolve_env_executable(python_path, "hmmscan")
    if not hmmscan:
        raise RuntimeError("Missing hmmscan executable needed to validate gutSMASH HMM databases")
    query = work_dir / "hmmscan_selftest.faa"
    query.write_text(">probe\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRK\n", encoding="utf-8")
    result = subprocess.run(
        [str(hmmscan), "--noali", str(hmm_path), str(query)],
        cwd=work_dir,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        message = (result.stderr or result.stdout or "no output").strip()
        raise RuntimeError(f"hmmscan could not read {hmm_path}: {message}")


def replace_symlink(source: Path, target: Path) -> None:
    if not source.exists():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        target.unlink()
    target.symlink_to(source.resolve())


def ensure_dir_or_replace(path: Path) -> None:
    if (path.exists() or path.is_symlink()) and not path.is_dir():
        path.unlink()
    path.mkdir(parents=True, exist_ok=True)


def prepare_databases(database: Path, source: Path) -> None:
    db_root = source / "antismash" / "databases"
    db_root.mkdir(parents=True, exist_ok=True)

    resfam = database / "resfam" / "resfam"
    if resfam.is_dir():
        target = db_root / "resfam"
        ensure_dir_or_replace(target)
        replace_symlink(resfam / "Resfams.hmm", target / "Resfams.hmm")
        for suffix in ("h3f", "h3i", "h3m", "h3p"):
            replace_symlink(resfam / f"Resfams.hmm.{suffix}", target / f"Resfams.hmm.{suffix}")

    clusterblast = database / "ClusterBlast"
    if clusterblast.is_dir():
        target = db_root / "clusterblast"
        ensure_dir_or_replace(target)
        replace_symlink(clusterblast / "clusters.txt", target / "clusters.txt")
        replace_symlink(clusterblast / "proteins.fasta", target / "proteins.fasta")
        replace_symlink(clusterblast / "proteins.dmnd", target / "proteins.dmnd")

    known = database / "KnownClusterBlast"
    if known.is_dir():
        target = db_root / "knownclusterblast"
        ensure_dir_or_replace(target)
        replace_symlink(known / "clusters.txt", target / "clusters.txt")
        replace_symlink(known / "proteins.fasta", target / "proteins.fasta")
        replace_symlink(known / "proteins.dmnd", target / "proteins.dmnd")

    smcogs = database / "smcogs"
    if smcogs.is_dir():
        target = source / "antismash" / "detection" / "genefunctions" / "gut_data"
        target.mkdir(parents=True, exist_ok=True)
        replace_symlink(smcogs / "smcogs.hmm", target / "smcogs.hmm")
        for suffix in HMMPRESS_SUFFIXES:
            replace_symlink(smcogs / f"smcogs.hmm.{suffix}", target / f"smcogs.hmm.{suffix}")


def tail(path: Path, limit: int = 80) -> str:
    try:
        return "\n".join(path.read_text(encoding="utf-8", errors="replace").splitlines()[-limit:])
    except OSError:
        return ""


def write_log(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line.rstrip() + "\n")


def normalize_index(output_dir: Path, log_path: Path) -> Path:
    expected = output_dir / "index.html"
    if nonempty(expected):
        return expected
    for candidate in sorted(output_dir.glob("**/index.html")):
        if candidate == expected:
            continue
        if nonempty(candidate):
            write_log(log_path, [f"gutSMASH wrote index.html at {candidate}; normalizing to {expected}"])
            shutil.copy2(str(candidate), str(expected))
            return expected
    return expected


def write_marker(output_dir: Path, sample: str, python_path: Path, source: Path, index: Path) -> None:
    marker = output_dir / ".probeat_gutsmash_done.json"
    payload = {
        "schema_version": 1,
        "status": "completed",
        "sample": sample,
        "python": str(python_path),
        "source": str(source),
        "index_html": str(index),
        "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    marker.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_status(
    output_dir: Path,
    sample: str,
    python_path: Path | None,
    source: Path | None,
    log_path: Path,
    status: str,
    error: str = "",
    index: Path | None = None,
) -> Path:
    status_path = output_dir / ".probeat_gutsmash_status.json"
    payload = {
        "schema_version": 1,
        "sample": sample,
        "status": status,
        "python": str(python_path) if python_path else "",
        "source": str(source) if source else "",
        "index_html": str(index) if index else "",
        "log": str(log_path),
        "error": error,
        "updated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    status_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return status_path


def fail(
    output_dir: Path,
    sample: str,
    log_path: Path,
    python_path: Path | None,
    source: Path | None,
    error: str,
    code: int,
    soft_fail: bool,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    done_marker = output_dir / ".probeat_gutsmash_done.json"
    if done_marker.exists():
        done_marker.unlink()
    write_status(output_dir, sample, python_path, source, log_path, "failed", error=error)
    if soft_fail:
        write_log(log_path, [f"gutSMASH soft-failed: {error}"])
        return 0
    return code


def self_test() -> int:
    errors = runtime_errors()
    if errors:
        print("[FAIL] gutSMASH runtime")
        for error in errors:
            print(f"       {error}")
        return 2
    python_path = resolve_python()
    source = resolve_source()
    if not python_path or not source:
        print("[FAIL] gutSMASH runtime")
        print("       Internal gutSMASH runtime resolution error after successful validation.")
        return 2
    try:
        with tempfile.TemporaryDirectory(prefix="probeat_gutsmash_selftest_") as tmpdir:
            staged_source = Path(tmpdir) / "source"
            shutil.copytree(source, staged_source, symlinks=True)
            env = runtime_environment(python_path, staged_source)
            prepare_bundled_hmms(staged_source, python_path, env)
            verify_hmmscan_can_read(
                staged_source / "antismash" / "detection" / "genefunctions" / "gut_data" / "smcogs.hmm",
                python_path,
                env,
                Path(tmpdir),
            )
    except Exception as exc:
        print("[FAIL] gutSMASH runtime")
        print(f"       gutSMASH HMM preparation check failed: {exc}")
        return 2
    print("[OK] gutSMASH runtime")
    print(f"       python: {python_path}")
    print(f"       source: {source}")
    print("       HMM preparation: ok")
    return 0


def run(args: argparse.Namespace) -> int:
    database = Path(args.database)
    if not database.is_dir():
        error = f"Missing gutSMASH database: {database}"
        print(error, file=sys.stderr)
        return fail(Path(args.output_dir), args.sample, Path(args.log), None, None, error, 2, args.soft_fail)

    errors = runtime_errors()
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        if args.log:
            write_log(Path(args.log), errors)
        return fail(Path(args.output_dir), args.sample, Path(args.log), resolve_python(), resolve_source(), "; ".join(errors), 2, args.soft_fail)

    python_path = resolve_python()
    source = resolve_source()
    if not python_path or not source:
        error = "Internal gutSMASH runtime resolution error after successful validation."
        print(error, file=sys.stderr)
        return fail(Path(args.output_dir), args.sample, Path(args.log), python_path, source, error, 2, args.soft_fail)

    output_dir = Path(args.output_dir)
    log_path = Path(args.log)
    shutil.rmtree(str(output_dir), ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=output_dir.parent, prefix=f".gutsmash_runtime_{args.sample}_") as tmpdir:
        staged_source = Path(tmpdir) / "source"
        shutil.copytree(source, staged_source, symlinks=True)
        prepare_databases(database, staged_source)

        env = runtime_environment(python_path, staged_source)
        try:
            prepare_bundled_hmms(staged_source, python_path, env, log_path)
        except Exception as exc:
            error = f"gutSMASH HMM preparation failed: {exc}"
            print(error, file=sys.stderr)
            return fail(output_dir, args.sample, log_path, python_path, source, error, 2, args.soft_fail)
        command: List[str] = [
            str(python_path),
            str(staged_source / "run_gutsmash.py"),
            "--genefinding-tool",
            "prodigal",
            "--minimal",
            "--cb-knownclusters",
            "--enable-genefunctions",
            "--output-dir",
            str(output_dir),
            "--cpus",
            str(args.cpus),
            str(args.input_fasta),
        ]
        with log_path.open("w", encoding="utf-8") as log_handle:
            result = subprocess.run(command, cwd=Path.cwd(), stdout=log_handle, stderr=subprocess.STDOUT, env=env, check=False)
        if result.returncode:
            error = "gutSMASH execution failed"
            print(f"{error}. Last log lines:", file=sys.stderr)
            print(tail(log_path), file=sys.stderr)
            return fail(output_dir, args.sample, log_path, python_path, source, error, result.returncode, args.soft_fail)

    index = normalize_index(output_dir, log_path)
    if not nonempty(index):
        error = "Missing gutSMASH index.html after real gutSMASH execution."
        print(f"{error} Last log lines:", file=sys.stderr)
        print(tail(log_path), file=sys.stderr)
        return fail(output_dir, args.sample, log_path, python_path, source, error, 2, args.soft_fail)

    write_marker(output_dir, args.sample, python_path, source, index)
    write_status(output_dir, args.sample, python_path, source, log_path, "completed", index=index)
    return 0


def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--sample", default="unknown")
    parser.add_argument("--input-fasta")
    parser.add_argument("--output-dir")
    parser.add_argument("--database")
    parser.add_argument("--log", default="logs/gutsmash/gutsmash.log")
    parser.add_argument("--cpus", default="1")
    parser.add_argument("--soft-fail", action="store_true")
    args = parser.parse_args(argv)
    if not args.self_test:
        missing = [name for name in ("input_fasta", "output_dir", "database") if not getattr(args, name)]
        if missing:
            parser.error("missing required arguments: " + ", ".join("--" + name.replace("_", "-") for name in missing))
    return args


def main(argv: List[str]) -> int:
    args = parse_args(argv)
    if args.self_test:
        return self_test()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

"""Shared gutSMASH runtime contract for PROBEAT."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
from typing import Dict, List, Optional


GUTSMASH_PROFILES = {"full", "probiotic", "ont"}


def default_root() -> Path:
    return Path(os.environ.get("PROBEAT_ROOT", Path(__file__).resolve().parents[2]))


def profile_requires_gutsmash_name(profile: str) -> bool:
    return profile in GUTSMASH_PROFILES


def profile_requires_gutsmash_config(profile: Dict[str, object]) -> bool:
    environments = profile.get("environments", [])
    return isinstance(environments, list) and "gutsmash" in environments


def runtime_python(root: Optional[Path] = None) -> Path:
    root = root or default_root()
    return root / "mamba_envs" / "gutsmash" / "bin" / "python"


def valid_source(path: Path) -> bool:
    try:
        return (
            (path / "run_gutsmash.py").is_file()
            and (path / "run_gutsmash.py").stat().st_size > 0
            and (path / "antismash" / "__main__.py").is_file()
            and (path / "antismash" / "__main__.py").stat().st_size > 0
        )
    except OSError:
        return False


def source_candidates(root: Optional[Path] = None) -> List[Path]:
    root = root or default_root()
    candidates: List[Path] = []
    env_source = os.environ.get("PROBEAT_GUTSMASH_SRC")
    if env_source:
        candidates.append(Path(env_source))
    candidates.extend(
        [
            Path("/opt/probeat/gutsmash_src"),
            Path("/app/workflow/gutsmash_src"),
            root / "workflow" / "gutsmash_src",
        ]
    )
    unique: List[Path] = []
    seen = set()
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            unique.append(candidate)
            seen.add(key)
    return unique


def resolve_source(root: Optional[Path] = None) -> Optional[Path]:
    for candidate in source_candidates(root):
        if valid_source(candidate):
            try:
                return candidate.resolve()
            except OSError:
                return candidate
    return None


def local_python_fallbacks(root: Optional[Path] = None) -> List[Path]:
    root = root or default_root()
    candidates: List[Path] = []
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.append(Path(conda_prefix) / "bin" / "python")
    candidates.extend(sorted((root / "snakemake_conda_envs").glob("*/bin/python")))
    candidates.extend(sorted((root / ".snakemake" / "conda").glob("*/bin/python")))
    for executable in ("python3", "python"):
        resolved = shutil.which(executable)
        if resolved:
            candidates.append(Path(resolved))
    return candidates


def allow_local_fallback() -> bool:
    return os.environ.get("PROBEAT_ALLOW_LOCAL_GUTSMASH_FALLBACK") == "1"


def resolve_python(root: Optional[Path] = None, allow_fallback: Optional[bool] = None) -> Optional[Path]:
    root = root or default_root()
    primary = runtime_python(root)
    if primary.is_file() and os.access(str(primary), os.X_OK):
        return primary
    if allow_fallback is None:
        allow_fallback = allow_local_fallback()
    if allow_fallback:
        for candidate in local_python_fallbacks(root):
            if candidate.is_file() and os.access(str(candidate), os.X_OK):
                return candidate
    return None


def runtime_environment(python_path: Path, source: Path) -> Dict[str, str]:
    env = dict(os.environ)
    env["PATH"] = f"{python_path.parent}:{env.get('PATH', '')}"
    old_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{source}{os.pathsep}{old_pythonpath}" if old_pythonpath else str(source)
    return env


def command_output(command: List[str], env: Dict[str, str]) -> Optional[str]:
    try:
        result = subprocess.run(
            command,
            cwd=default_root(),
            env=env,
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        return f"{command[0]}: {exc}"
    output = (result.stdout or "").strip()
    if not output:
        output = (result.stderr or "").strip()
    if result.returncode != 0:
        return f"{command[0]} exited {result.returncode}: {output or 'no output'}"
    return output


def compatibility_errors(python_path: Path, source: Path) -> List[str]:
    errors: List[str] = []
    env = runtime_environment(python_path, source)

    muscle_path = python_path.parent / "muscle"
    if not (muscle_path.is_file() and os.access(str(muscle_path), os.X_OK)):
        errors.append(f"Missing gutSMASH MUSCLE executable. Expected {muscle_path}.")
    else:
        muscle_version = command_output([str(muscle_path), "-version"], env)
        if muscle_version is None or not muscle_version.splitlines()[0].startswith("MUSCLE"):
            errors.append(
                "Incompatible gutSMASH MUSCLE executable. Expected classic MUSCLE output starting with "
                f"'MUSCLE', got: {muscle_version or 'no output'}"
            )

    help_check = command_output([str(python_path), str(source / "run_gutsmash.py"), "--help"], env)
    if help_check is None or help_check.startswith(str(python_path)):
        errors.append(
            "gutSMASH source/runtime compatibility check failed for run_gutsmash.py --help: "
            f"{help_check or 'no output'}"
        )
    return errors


def runtime_errors(root: Optional[Path] = None, allow_fallback: Optional[bool] = None) -> List[str]:
    root = root or default_root()
    errors: List[str] = []
    primary = runtime_python(root)
    if not (primary.is_file() and os.access(str(primary), os.X_OK)):
        if allow_fallback is None:
            allow_fallback = allow_local_fallback()
        fallback = resolve_python(root, allow_fallback=allow_fallback)
        if not fallback:
            if allow_fallback:
                errors.append(
                    "Missing gutSMASH Python executable. Searched "
                    f"{primary}, active conda env, snakemake_conda_envs, .snakemake/conda, python3, and python."
                )
            else:
                errors.append(
                    "Missing gutSMASH Python executable. Expected "
                    f"{primary}. Set PROBEAT_ALLOW_LOCAL_GUTSMASH_FALLBACK=1 only for local development fallback."
                )
    source = resolve_source(root)
    if not source:
        searched = ", ".join(str(candidate) for candidate in source_candidates(root))
        errors.append(f"Missing gutSMASH source. Searched {searched}.")
    python_path = resolve_python(root, allow_fallback=allow_fallback)
    if python_path and source:
        errors.extend(compatibility_errors(python_path, source))
    return errors

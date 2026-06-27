#!/app/mamba_envs/snakemake-conda-helper-2026-06-04/bin/python
"""Small conda/mamba compatibility wrapper around micromamba for Snakemake."""

from __future__ import annotations

import json
import os
import subprocess
import sys


def main() -> int:
    args = sys.argv[1:]
    if args in (["--version"], ["-V"]):
        print("conda 25.1.1")
        return 0
    if args == ["config", "--get", "channel_priority", "--json"]:
        print(json.dumps({"get": {}, "rc_path": "/opt/conda/.condarc", "success": True, "warnings": []}))
        return 0
    if args == ["info", "--json"]:
        result = subprocess.run(
            ["micromamba", "info", "--json"],
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode:
            sys.stderr.write(result.stderr)
            return result.returncode
        info = json.loads(result.stdout)
        prefix = info.get("base environment") or info.get("env location") or "/opt/conda"
        info.setdefault("conda_prefix", prefix)
        info.setdefault("conda_version", "25.1.1")
        info.setdefault("root_prefix", prefix)
        info.setdefault("envs_dirs", info.get("envs directories", [os.path.join(prefix, "envs")]))
        info.setdefault("pkgs_dirs", info.get("package cache", [os.path.join(prefix, "pkgs")]))
        print(json.dumps(info))
        return 0
    if args[:2] == ["env", "create"]:
        args = [arg for arg in args if arg != "--no-default-packages"]
    os.execvp("micromamba", ["micromamba", *args])
    return 127


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env bash
if [ "$#" -gt 0 ]; then
  export CONDA_PREFIX="$1"
  export PATH="$CONDA_PREFIX/bin:$PATH"
fi

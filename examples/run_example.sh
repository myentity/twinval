#!/usr/bin/env bash
# Smoke run of the whole pipeline on the bundled example data (about 4 minutes with R, seconds without).
set -euo pipefail
cd "$(dirname "$0")/.."
twinval validate --data examples/data
twinval run --data examples/data --out results/example --fast
echo "Open results/example/report.md"

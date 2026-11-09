#!/bin/bash
# Run ctx_finder on the synthetic example data and compare the calls with expected_calls.tsv.
# Runs without Slurm (the changepoint step runs in this process).
# Usage: bash run_example.sh [output_dir]     (default: ./example_output)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
PL=$HERE/../pipeline
OUT=${1:-example_output}

python "$PL/01_prepare.py" --manifest "$HERE/manifest.tsv" --out "$OUT" --run-changepoints-here
python "$PL/02_call.py" --manifest "$HERE/manifest.tsv" --out "$OUT"

python - "$OUT/ctx_finder_samples.tsv" "$HERE/expected_calls.tsv" <<'EOF'
import sys
import pandas as pd
got = pd.read_csv(sys.argv[1], sep="\t", dtype=str, keep_default_na=False)
exp = pd.read_csv(sys.argv[2], sep="\t", dtype=str, keep_default_na=False)
m = exp.merge(got, on="sample_id", suffixes=("_expected", ""))
m["match"] = (m["chromothripsis"] == m["chromothripsis_expected"]) & \
             (m["called_chromosomes"] == m["called_chromosomes_expected"])
print(m[["sample_id", "chromothripsis", "called_chromosomes", "match", "description"]].to_string(index=False))
print("All calls match the expected results." if m["match"].all() else "Some calls differ from the expected results.")
EOF

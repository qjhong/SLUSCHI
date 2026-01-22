#!/usr/bin/env bash
# script_v4_wrapper.sh
# Usage: ./script_v4_wrapper.sh [N]
# Runs the original csh pipeline, captures logs, postprocesses into pretty stdout and results.json

set -euo pipefail

# Where the original script lives (adjust if needed)
SLUSCHI_SRC="/home/qhong7/data/sluschi/sluschi_latest/src/mds_src"
ORIG_SCRIPT="${SLUSCHI_SRC}/script_v4.csh"

if [ ! -x "$ORIG_SCRIPT" ]; then
  echo "ERROR: original script not found or not executable: $ORIG_SCRIPT" >&2
  exit 2
fi

# job run local workspace (use current directory)
RUN_DIR="$(pwd)"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
STDOUT_RAW="${RUN_DIR}/stdout.log"
STDERR_RAW="${RUN_DIR}/stderr.log"
STDOUT_PRETTY="${RUN_DIR}/stdout_pretty.txt"
RESULTS_JSON="${RUN_DIR}/results.json"

# forward argument (n) if provided
ARG="${1:-0}"

echo "=== SLUSCHI RUN START: ${TIMESTAMP} ==="
echo "Running: ${ORIG_SCRIPT} ${ARG}"
echo ""

# Run the original C-shell script, capture raw logs
# Use /bin/csh to execute so environment is similar to before
/bin/csh -x "${ORIG_SCRIPT}" "${ARG}" > "${STDOUT_RAW}" 2> "${STDERR_RAW}" || true
# we use || true to allow postprocessing even if exit non-zero; you can change if you'd like

# Postprocess raw stdout (calls python script)
python3 "${SLUSCHI_SRC}/../api/postprocess_stdout.py" \
  --stdout "${STDOUT_RAW}" \
  --stderr "${STDERR_RAW}" \
  --pretty "${STDOUT_PRETTY}" \
  --json "${RESULTS_JSON}" \
  --timestamp "${TIMESTAMP}" \
  || {
    echo "WARNING: postprocessing script failed" >&2
}

# Print the pretty stdout to the real stdout so web wrapper can capture it
if [ -f "${STDOUT_PRETTY}" ]; then
  cat "${STDOUT_PRETTY}"
fi

# Print the delimited JSON block so UI can parse it from stdout
if [ -f "${RESULTS_JSON}" ]; then
  echo ""
  echo "===SLUSCHI_SUMMARY_JSON_START==="
  cat "${RESULTS_JSON}"
  echo ""
  echo "===SLUSCHI_SUMMARY_JSON_END==="
fi

# Final marker line
EXITCODE=0
echo ""
echo "=== SLUSCHI DONE: status=OK exit=${EXITCODE} ==="

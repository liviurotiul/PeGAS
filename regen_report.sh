#!/usr/bin/env bash
# Re-render the PeGAS R report from an existing run's output directory using the
# CURRENT dev template + CSS (src/pegas/report_template.Rmd, report_static.css) —
# no need to rerun the pipeline. It only reads the dataframe/*.tsv, results/ and
# fastqc/ that the run already produced.
#
# Usage:
#   ./regen_report.sh <pegas_output_dir> [output_html]
#
# If <output_html> is omitted it writes <output_dir>/report/report_r.html when
# that is writable, otherwise ~/.pegas/report_preview/<run>_report_r.html.
# Override the report env with PEGAS_REPORT_ENV.
set -euo pipefail

OUT_DIR="${1:?usage: regen_report.sh <pegas_output_dir> [output_html]}"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/src/pegas" && pwd)"
REPORT_ENV="${PEGAS_REPORT_ENV:-$HOME/.pegas/envs/report_r}"

if [[ ! -f "$OUT_DIR/dataframe/results.csv" ]]; then
  echo "[regen] $OUT_DIR/dataframe/results.csv not found — is this a completed PeGAS run?" >&2
  exit 1
fi

# Decide where to write the HTML (run dirs from Docker are often root-owned).
DEST="${2:-}"
if [[ -z "$DEST" ]]; then
  if [[ -w "$OUT_DIR/report" ]] || { [[ ! -e "$OUT_DIR/report" ]] && [[ -w "$OUT_DIR" ]]; }; then
    DEST="$OUT_DIR/report/report_r.html"
  else
    DEST="$HOME/.pegas/report_preview/$(basename "$OUT_DIR")_report_r.html"
    echo "[regen] $OUT_DIR/report is not writable; writing preview to a user location instead." >&2
  fi
fi
mkdir -p "$(dirname "$DEST")"

# Best-effort data dir (only used for the FastQC fallback / input table; the
# report prefers the manifest + output_dir, so a stale path is harmless).
DATA=""
if [[ -f "$OUT_DIR/pegas_config.json" ]]; then
  DATA="$(python3 -c "import json,sys;print(json.load(open('$OUT_DIR/pegas_config.json')).get('data',''))" 2>/dev/null || true)"
fi

echo "[regen] rendering with dev template: $SRC/report_template.Rmd"
conda run -p "$REPORT_ENV" --no-capture-output Rscript "$SRC/render_report.R" \
  --rmd "$SRC/report_template.Rmd" \
  --output "$DEST" \
  --dataframe_csv "$OUT_DIR/dataframe/results.csv" \
  --report_html "" \
  --data_dir "$DATA" \
  --output_dir "$OUT_DIR" \
  --pegas_version "dev-preview" \
  --pegas_install_dir "$SRC"

echo "[regen] done"
echo "[regen] open: file://$DEST"

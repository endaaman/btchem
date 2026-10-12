#!/usr/bin/env bash
# 予備解析一式を再生成する（repo root から: bash prelim/run_all.sh）
# 前提: data/seq/<GSE>/ に GEO の補足ファイルと SOFT が落としてあること（docs/prelim_de.md 参照）
set -euo pipefail
cd "$(dirname "$0")/.."
uv run python prelim/prepare.py
uv run python prelim/qc.py
uv run python prelim/de.py
uv run python prelim/volcano.py
uv run python prelim/enrich.py    # Enrichr のライブラリ取得にネットが要る
uv run python prelim/overlap.py

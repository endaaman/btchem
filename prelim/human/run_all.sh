#!/usr/bin/env bash
# ヒト骨 RNA-seq の予備解析一式（repo root から: bash prelim/human/run_all.sh）
# 前提: data/seq/{GSE141595,GSE141610,GSE72815,GSE276529}/ に GEO の補足ファイルと SOFT、
#       out/prelim/human/genesets/c2.cgp.v2023.2.Hs.symbols.gmt（SenMayo 用）。マウス側の out/prelim/{de,enrich} も使う。
set -euo pipefail
cd "$(dirname "$0")/../.."
uv run python -I prelim/human/prepare.py
uv run python -I prelim/human/qc_h.py
uv run python -I prelim/human/de_h.py
uv run python -I prelim/human/volcano_h.py
uv run python -I prelim/human/enrich_h.py   # Enrichr ライブラリ取得にネットが要る
uv run python -I prelim/human/compare_h.py

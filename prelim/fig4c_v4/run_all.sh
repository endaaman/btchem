#!/usr/bin/env bash
# 4C v4（案1 GSE304562 / 案2 GSE276391 の読み直し）。repo root から: bash prelim/fig4c_v4/run_all.sh
# data/seq/GSE304562/ に GEO の GSE304562_FPKM_allSamples.txt.gz と GSE304562_family.soft.gz が要る
# （https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304562/）。既存の DE / GSEA 表（out/prelim/{de,chemo,chemo2,human,invitro,aging_f,summary}）は読むだけ。
set -euo pipefail
cd "$(dirname "$0")/../.."
uv run python -I prelim/fig4c_v4/prepare.py
uv run python -I prelim/fig4c_v4/de.py
uv run python -I prelim/fig4c_v4/volcano.py
uv run python -I prelim/fig4c_v4/gsea.py
uv run python -I prelim/fig4c_v4/compare.py

#!/bin/bash
# in vitro ヒト骨系細胞 × DNA 傷害の予備解析。repo root から: bash prelim/invitro/run_all.sh
set -euo pipefail
uv run python -I prelim/invitro/prepare.py
uv run python -I prelim/invitro/de_iv.py
uv run python -I prelim/invitro/gsea_iv.py   # out/prelim/invitro/msigdb/ に human GMT が要る（docstring 参照）
uv run python -I prelim/invitro/figs_iv.py

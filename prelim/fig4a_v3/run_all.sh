#!/usr/bin/env bash
# Fig 4A v3（案1 GSE141595 デノスマブ / 案4 Kaya 2022）の再生成。repo root から:  bash prelim/fig4a_v3/run_all.sh
# 入力: data/seq/{GSE141595,GSE141610,GSE72815,GSE222752,GSE205237,GSE276391}/, data/NIHMS1776347-supplement-tS1.xlsx,
#       out/prelim/human/de/*.csv（既存のヒト DE）, out/prelim/invitro/msigdb/*.gmt, out/prelim/chemo2/msigdb/*.gmt
# GSEA は結果ファイルがあれば飛ばす（回し直すときは out/prelim/fig4a_v3/gsea/ と kaya/gsea_*.csv を消す）。所要 ~45 分。
set -euo pipefail
P=prelim/fig4a_v3
uv run python -I $P/pairing.py        # GSE141595 と GSE141610 が同じ人か
uv run python -I $P/composition.py    # 細胞種マーカーの組成
uv run python -I $P/dmab_de.py        # 共変量つき DE（感度分析）
uv run python -I $P/gsea4a.py         # Hallmark / GO BP / Reactome と要約モジュール形式
uv run python -I $P/dmab_sets.py      # デノスマブの中身
uv run python -I $P/reversal.py       # 閉経の変化をデノスマブが戻すか
uv run python -I $P/inflam_apop.py    # 炎症・アポトーシスの強さ
uv run python -I $P/summary_cols.py   # 要約モジュール図の列
uv run python -I $P/kaya.py           # 案4 Kaya の読み直し

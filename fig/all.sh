#!/usr/bin/env bash
# 本番図をすべて再生成する（入力は out/prelim/ の DE / GSEA 表。再計算しない）。
#   bash fig/all.sh   -> out/fig/{fig4,figS1,figS2,figS3}.{png,pdf,tiff}
#   + fig2 / fig3（臨床 xlsx が無ければ fig2_mock / fig3_mock）
set -euo pipefail
cd "$(dirname "$0")/.."
for f in fig4 figS1 figS2 figS3; do
    echo "== $f"
    uv run python "fig/$f.py"
done

# Figure 2 / 3（臨床 KM）。臨床 xlsx（fig/km.py の CLINICAL_XLSX）が無い間は --mock で fig2_mock / fig3_mock を出す。
xlsx=$(cd fig && uv run python -c 'import km; print(km.CLINICAL_XLSX)')
for f in fig2 fig3; do
    if [ -f "$xlsx" ]; then
        echo "== $f"
        uv run python "fig/$f.py"
    else
        echo "== $f (MOCK: $xlsx not found)"
        uv run python "fig/$f.py" --mock
    fi
done

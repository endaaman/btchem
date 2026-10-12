"""Figure 2: 単変量 KM — 骨・肺・肝の転移なし生存、周術期化療 (+) vs (−)。

パネル: A Bone / B Lung / C Liver（全 785 例。脳は載せない＝清水指示）
凡例: log-rank p, Cox 単変量 HR (化療 (+) vs (−)), 95% CI。書式は fig/km.py（1 本目 fig3 踏襲）。

実行:
  uv run python fig/fig2.py --mock   -> out/fig/fig2_mock.{png,pdf,tiff}, fig2_mock_stats.csv（擬似データ）
  uv run python fig/fig2.py          -> out/fig/fig2.{png,pdf,tiff}, fig2_stats.csv（臨床 xlsx 到着後）
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import km  # noqa: E402

PANELS = [('A', 'Bone', 'bone'), ('B', 'Lung', 'lung'), ('C', 'Liver', 'liver')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mock', action='store_true', help='1 本目の csv から擬似データで描く')
    ap.add_argument('--xlsx', default=km.CLINICAL_XLSX, help='臨床 xlsx（既定 km.CLINICAL_XLSX）')
    args = ap.parse_args()

    common.setup()
    df = km.load(args.mock, args.xlsx, endpoints=[p[2] for p in PANELS])
    name = 'fig2_mock' if args.mock else 'fig2'
    fig = plt.figure(figsize=(common.mm(km.FIG_W), common.mm(km.figure_height(args.mock))))
    rows = km.draw(fig, [(label, title, ep, df) for label, title, ep in PANELS], args.mock)
    km.write_stats(rows, name)
    common.save(fig, name)


if __name__ == '__main__':
    main()

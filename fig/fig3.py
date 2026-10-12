"""Figure 3: 骨転移なし生存の KM を閉経状態で層別、周術期化療 (+) vs (−)。

パネル: A Enrolled women / B Premenopausal women / C Postmenopausal women（骨 = BMFS）
凡例: log-rank p, Cox 単変量 HR (化療 (+) vs (−)), 95% CI。書式は fig/km.py（1 本目 fig3 踏襲）。

実行:
  uv run python fig/fig3.py --mock   -> out/fig/fig3_mock.{png,pdf,tiff}, fig3_mock_stats.csv（擬似データ）
  uv run python fig/fig3.py          -> out/fig/fig3.{png,pdf,tiff}, fig3_stats.csv（臨床 xlsx 到着後）
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import km  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mock', action='store_true', help='1 本目の csv から擬似データで描く')
    ap.add_argument('--xlsx', default=km.CLINICAL_XLSX, help='臨床 xlsx（既定 km.CLINICAL_XLSX）')
    args = ap.parse_args()

    common.setup()
    df = km.load(args.mock, args.xlsx, endpoints=['bone'])
    name = 'fig3_mock' if args.mock else 'fig3'
    panels = [
        ('A', 'Enrolled women', 'bone', df),
        ('B', 'Premenopausal women', 'bone', df[df.menopause == 0]),
        ('C', 'Postmenopausal women', 'bone', df[df.menopause == 1]),
    ]
    fig = plt.figure(figsize=(common.mm(km.FIG_W), common.mm(km.figure_height(args.mock))))
    rows = km.draw(fig, panels, args.mock)
    km.write_stats(rows, name)
    common.save(fig, name)


if __name__ == '__main__':
    main()

"""4C 案1/案2 の volcano（予備図）。描画は既存 prelim/volcano.py の plot() をそのまま使う
（基本閾値 |log2FC| ≥ 1 かつ padj < 0.05、注目遺伝子パネル prelim/genes.py を白抜き丸、padj < 0.1 にラベル）。

  uv run python -I prelim/fig4c_v4/volcano.py      # repo root から（de.py の後）

出力（out/prelim/fig4c_v4/volcano/）: Cort_MC3T3.png|pdf（単独 80 mm）, C4_options.png|pdf（案1 | 案2 を横並び 170 mm）
案2 の DE 表は既存 out/prelim/de/C_DEX.csv を読むだけ（再計算しない）。
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
import common  # noqa: E402
from volcano import plot  # noqa: E402

OUT = ROOT / 'out/prelim/fig4c_v4/volcano'
PANELS = [
    ('Cort_MC3T3', ROOT / 'out/prelim/fig4c_v4/de/Cort_MC3T3.csv',
     'Option 1: GSE304562\nMC3T3-E1 osteoblasts, corticosterone 24 h (3 vs 3)'),
    ('C_DEX', ROOT / 'out/prelim/de/C_DEX.csv',
     'Option 2: GSE276391\nmouse tibia ♂, dexamethasone 5 d (5 vs 5)'),
]


def load(p):
    return pd.read_csv(p).dropna(subset=['log2FoldChange', 'pvalue'])


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'{name}.{ext}', dpi=common.DPI, bbox_inches='tight', pad_inches=0.02)
        print('saved', OUT / f'{name}.{ext}')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    common.setup()
    lb, p, t = PANELS[0]
    fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(72)))
    plot(ax, load(p), t)
    ax.title.set_fontsize(7.5)
    save(fig, lb)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(common.mm(170), common.mm(72)))
    for ax, (lb, p, t), pl in zip(axes, PANELS, 'AB'):
        plot(ax, load(p), t)
        ax.title.set_fontsize(7.5)
        common.panel_label(ax, pl, dx=-0.18)
    fig.tight_layout()
    save(fig, 'C4_options')
    plt.close(fig)


if __name__ == '__main__':
    main()

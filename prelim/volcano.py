"""DE 結果から volcano と集計表を作る。

閾値（基本）: |log2FC| >= 1 かつ padj < 0.05。併記: padj < 0.1、|log2FC| >= 0.585 (1.5 倍)。
注目遺伝子パネル（prelim/genes.py）は有意かどうかに関係なく全点を色付けし、padj < 0.1 のものにラベル。

出力:
  out/prelim/volcano/<label>.png|pdf
  out/prelim/volcano/ABC.png|pdf           主 3 contrast を横並び
  out/prelim/deg_counts.csv                閾値ごとの DEG 数
  out/prelim/panel_table.csv               パネル遺伝子 x contrast の log2FC / padj（未検出は空欄）
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from contrasts import C, MAIN  # noqa: E402
from genes import INFLAM, APOP, ALL_GENES  # noqa: E402

DE = Path(common.OUT_DIR) / 'de'
OUT = Path(common.OUT_DIR)
LFC, PADJ = 1.0, 0.05
COL_UP, COL_DOWN, COL_NS = '#D98880', '#7FB3D5', '#CFCFCF'
COL_INFLAM, COL_APOP = '#B03A2E', '#1F4E79'
THRESHOLDS = [(0.05, 1.0), (0.05, 0.585), (0.05, 0.0), (0.1, 1.0), (0.1, 0.585), (0.1, 0.0)]


def load(label):
    d = pd.read_csv(DE / f'{label}.csv')
    return d.dropna(subset=['log2FoldChange', 'pvalue'])


def plot(ax, d, title, lfc=LFC, padj=PADJ):
    p = d.pvalue.clip(lower=1e-300)
    y = -np.log10(p)
    sig = (d.padj < padj) & (d.log2FoldChange.abs() >= lfc)
    up = sig & (d.log2FoldChange > 0)
    dn = sig & (d.log2FoldChange < 0)
    ns = ~sig
    ax.scatter(d.log2FoldChange[ns], y[ns], s=1.5, c=COL_NS, lw=0, rasterized=True)
    ax.scatter(d.log2FoldChange[up], y[up], s=2, c=COL_UP, lw=0, rasterized=True,
               label=f'Up ({up.sum()})')
    ax.scatter(d.log2FoldChange[dn], y[dn], s=2, c=COL_DOWN, lw=0, rasterized=True,
               label=f'Down ({dn.sum()})')
    for genes, col, name in [(INFLAM, COL_INFLAM, 'Inflam./OC'), (APOP, COL_APOP, 'Apopt./senesc.')]:
        m = d.gene.isin(genes)
        if name == 'Apopt./senesc.':
            m &= d.gene != 'Il6'   # Il6 は inflam 側の色で描く
        ax.scatter(d.log2FoldChange[m], y[m], s=9, facecolor='none', edgecolor=col, lw=0.6,
                   label=name)
        lab = m & (d.padj < 0.1)
        for gx, gy, gn in zip(d.log2FoldChange[lab], y[lab], d.gene[lab]):
            ax.annotate(gn, (gx, gy), xytext=(2, 2), textcoords='offset points',
                        fontsize=6, color=col, fontstyle='italic')
    # padj 閾値に相当する p 値（BH なので近似: 有意遺伝子の最大 p）
    if (d.padj < padj).any():
        ax.axhline(-np.log10(d.pvalue[d.padj < padj].max()), ls='--', lw=0.5, c='k')
    for v in (-lfc, lfc):
        ax.axvline(v, ls='--', lw=0.5, c='k')
    xm = np.nanpercentile(d.log2FoldChange.abs(), 99.9) * 1.1
    ax.set_xlim(-max(xm, 3), max(xm, 3))
    ax.set_xlabel('log$_2$ fold change')
    ax.set_ylabel('$-$log$_{10}$ P')
    ax.set_title(title)
    ax.legend(fontsize=5.5, frameon=False, loc='upper left', markerscale=1.5, handletextpad=0.2)


def main():
    common.setup()
    (OUT / 'volcano').mkdir(parents=True, exist_ok=True)
    counts, panel = [], {}
    for lb, c in C.items():
        if not (DE / f'{lb}.csv').exists():
            continue
        d = load(lb)
        for pa, lf in THRESHOLDS:
            s = (d.padj < pa) & (d.log2FoldChange.abs() >= lf)
            counts.append({'label': lb, 'padj': pa, 'lfc': lf, 'up': int((s & (d.log2FoldChange > 0)).sum()),
                           'down': int((s & (d.log2FoldChange < 0)).sum()), 'tested': len(d)})
        di = d.set_index('gene')
        panel[(lb, 'log2FC')] = di.log2FoldChange.reindex(ALL_GENES)
        panel[(lb, 'padj')] = di.padj.reindex(ALL_GENES)
        fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(70)))
        plot(ax, d, f"{lb}: {c['title']}")
        common.save(fig, f'volcano/{lb}', formats=('png', 'pdf'))
        plt.close(fig)
    cnt = pd.DataFrame(counts)
    cnt.to_csv(OUT / 'deg_counts.csv', index=False)
    pt = pd.DataFrame(panel)
    pt.index.name = 'gene'
    pt.insert(0, ('panel', ''), ['inflam+apop' if g == 'Il6' else ('inflam' if g in INFLAM else 'apop')
                                 for g in pt.index])
    pt.columns = [f'{a}|{b}' if b else a for a, b in pt.columns]
    pt.to_csv(OUT / 'panel_table.csv')

    fig, axes = plt.subplots(1, 3, figsize=(common.mm(170), common.mm(62)))
    for ax, (k, lb), pl in zip(axes, MAIN.items(), 'ABC'):
        plot(ax, load(lb), C[lb]['title'])
        common.panel_label(ax, pl, dx=-0.22)
    fig.tight_layout()
    common.save(fig, 'volcano/ABC', formats=('png', 'pdf'))
    pd.set_option('display.width', 200)
    print(cnt.pivot_table(index='label', columns=['padj', 'lfc'], values=['up', 'down']).to_string())


if __name__ == '__main__':
    main()

"""ヒト contrast の volcano・DEG 数・注目パネル表。prelim/volcano.py の plot() と閾値をそのまま使う
（パネル遺伝子リストだけヒトシンボルに差し替える）。

  uv run python -I prelim/human/volcano_h.py
出力: out/prelim/human/volcano/<label>.png|pdf, main3.png|pdf, deg_counts.csv, panel_table.csv
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
import volcano  # noqa: E402
import hgenes  # noqa: E402
from hcontrasts import HC, OUT, MAIN3  # noqa: E402

DE = OUT / 'de'
# volcano.plot はモジュール変数の INFLAM / APOP を参照する。APOP 側からは IL6 を除く（inflam 色で描く。元と同じ扱い）
volcano.INFLAM = hgenes.INFLAM
volcano.APOP = [g for g in hgenes.APOP if g != 'IL6']


def load(lb):
    d = pd.read_csv(DE / f'{lb}.csv')
    return d.dropna(subset=['log2FoldChange', 'pvalue'])


def main():
    common.setup()
    (OUT / 'volcano').mkdir(parents=True, exist_ok=True)
    counts, panel = [], {}
    for lb, c in HC.items():
        d = load(lb)
        for pa, lf in volcano.THRESHOLDS:
            s = (d.padj < pa) & (d.log2FoldChange.abs() >= lf)
            counts.append({'label': lb, 'padj': pa, 'lfc': lf,
                           'up': int((s & (d.log2FoldChange > 0)).sum()),
                           'down': int((s & (d.log2FoldChange < 0)).sum()), 'tested': len(d)})
        di = d.drop_duplicates('gene').set_index('gene')
        panel[(lb, 'log2FC')] = di.log2FoldChange.reindex(hgenes.ALL_GENES)
        panel[(lb, 'padj')] = di.padj.reindex(hgenes.ALL_GENES)
        fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(70)))
        volcano.plot(ax, d, f'{lb}: {c["title"]}')
        ax.title.set_fontsize(7)
        common.save(fig, f'human/volcano/{lb}', formats=('png', 'pdf'))
        plt.close(fig)
    cnt = pd.DataFrame(counts)
    cnt.to_csv(OUT / 'deg_counts.csv', index=False)
    pt = pd.DataFrame(panel)
    pt.index.name = 'gene'
    pt.insert(0, ('panel', ''), ['inflam+apop' if g == 'IL6' else ('inflam' if g in hgenes.INFLAM else 'apop')
                                 for g in pt.index])
    pt.columns = [f'{a}|{b}' if b else a for a, b in pt.columns]
    pt.to_csv(OUT / 'panel_table.csv')

    fig, axes = plt.subplots(1, 3, figsize=(common.mm(170), common.mm(62)))
    for ax, lb, pl in zip(axes, MAIN3, 'ABC'):
        volcano.plot(ax, load(lb), HC[lb]['title'])
        ax.title.set_fontsize(7)
        common.panel_label(ax, pl, dx=-0.22)
    fig.tight_layout()
    common.save(fig, 'human/volcano/main3', formats=('png', 'pdf'))
    pd.set_option('display.width', 220)
    print(cnt.pivot_table(index='label', columns=['padj', 'lfc'], values=['up', 'down']).to_string())


if __name__ == '__main__':
    main()

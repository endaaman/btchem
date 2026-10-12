"""volcano の共通描画（fig4.py / figS1.py が使う）。

入力は既存の DE 表（prelim の出力 csv。再計算しない）: gene, log2FoldChange, pvalue, padj。
- 縦軸 −log10(padj)。padj が NaN の遺伝子（DESeq2 の独立フィルタ / Cook 外れ値）は描かない。
- DEG = padj < 0.05 かつ |log2FC| ≥ 1。上昇 = 赤 (COLOR_HRO)、低下 = 青 (COLOR_LRO)、それ以外 = 薄灰。
- 注目遺伝子パネル（prelim/genes.py、事前定義）は padj < 0.05 のものだけ縁取り付きで強調しラベル（FC は問わない）。
  記号: 炎症・破骨 = 丸、アポトーシス・老化 = 菱形。塗りは変化の向き（赤 = 上昇, 青 = 低下）。
  ラベルが太字 = 既報の向き（炎症・破骨の促進 / アポトーシス・老化の促進）と一致、細字 = 逆向き。
"""
import os
import sys

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'prelim'))
import genes as G  # noqa: E402  (prelim/genes.py: 読むだけ)

PADJ = 0.05
LFC = 1.0
C_UP = common.COLOR_HRO
C_DOWN = common.COLOR_LRO
C_NS = '#D4D4D4'
C_LINE = '#8C8C8C'

# 既報で「炎症・破骨促進 / アポトーシス・老化促進」となる向き（+1 = 上昇, −1 = 低下）。
# OPG（Tnfrsf11b）は RANKL のデコイ受容体なので低下が破骨促進。Bcl2 / Bcl2l1 は抗アポトーシスなので低下がアポトーシス方向。
# Sost / Dkk1 は Wnt 阻害（骨形成抑制・RANKL/OPG 比上昇）として上昇を「破骨促進」側に置く。
EXPECTED = {g: +1 for g in G.ALL_GENES}
EXPECTED.update({'Tnfrsf11b': -1, 'Bcl2': -1, 'Bcl2l1': -1})

# マウス → ヒトのシンボル（大文字化で合わないもの）
HUMAN_ALIAS = {'Fasl': 'FASLG', 'Trp53': 'TP53'}


def human_symbol(g):
    return HUMAN_ALIAS.get(g, g.upper())


def panel_map(species):
    """{表示シンボル: (mouse シンボル, カテゴリ)}。Il6 は両方に属するので inflam として扱う。"""
    out = {}
    for g in G.ALL_GENES:
        cat = 'inflam' if g in G.INFLAM else 'apop'
        sym = g if species == 'mouse' else human_symbol(g)
        out[sym] = (g, cat)
    return out


def load(path):
    d = pd.read_csv(path)
    assert d.gene.is_unique, path
    return d


def deg_counts(d):
    s = (d.padj < PADJ) & (d.log2FoldChange.abs() >= LFC)
    return int((s & (d.log2FoldChange > 0)).sum()), int((s & (d.log2FoldChange < 0)).sum())


def sig_panel(d, species):
    """padj < 0.05 のパネル遺伝子の表（ラベル対象）。"""
    pm = panel_map(species)
    x = d[d.gene.isin(pm) & (d.padj < PADJ)].copy()
    x['mouse'] = x.gene.map(lambda s: pm[s][0])
    x['category'] = x.gene.map(lambda s: pm[s][1])
    x['consistent'] = np.sign(x.log2FoldChange) == x.mouse.map(EXPECTED)
    x['deg'] = x.log2FoldChange.abs() >= LFC
    return x.sort_values('padj')


def plot(ax, d, title, species='mouse', label_fs=7, xlim=None, ylim=None, label_kw=None):
    """1 パネル描画。戻り値は (up, down, ラベルした遺伝子の表)。"""
    d = d.dropna(subset=['padj', 'log2FoldChange']).copy()
    d['y'] = -np.log10(d.padj.clip(lower=1e-300))
    up = (d.padj < PADJ) & (d.log2FoldChange >= LFC)
    dn = (d.padj < PADJ) & (d.log2FoldChange <= -LFC)
    ns = ~(up | dn)

    kw = dict(s=2.5, linewidths=0, rasterized=True)
    ax.scatter(d.log2FoldChange[ns], d.y[ns], c=C_NS, zorder=1, **kw)
    ax.scatter(d.log2FoldChange[dn], d.y[dn], c=C_DOWN, alpha=0.75, zorder=2, **kw)
    ax.scatter(d.log2FoldChange[up], d.y[up], c=C_UP, alpha=0.75, zorder=2, **kw)

    # 軸範囲: x は対称、y は 0 から
    if xlim is None:
        m = max(d.log2FoldChange.abs().max(), 1.5)
        xlim = np.ceil(m * 1.08)
    ax.set_xlim(-xlim, xlim)
    ymax = max(d.y.max(), -np.log10(PADJ)) if ylim is None else ylim
    ax.set_ylim(0, ymax * 1.10 + 0.2)

    lw = 0.5
    for xv in (-LFC, LFC):
        ax.axvline(xv, c=C_LINE, lw=lw, ls=(0, (3, 2)), zorder=0)
    ax.axhline(-np.log10(PADJ), c=C_LINE, lw=lw, ls=(0, (3, 2)), zorder=0)

    # 注目遺伝子（padj < 0.05 のみ）。ラベルは上側の空き領域に縦に並べ、引き出し線で点とつなぐ
    sp = sig_panel(d, species)
    for _, r in sp.iterrows():
        mk = 'o' if r.category == 'inflam' else 'D'
        ax.scatter(r.log2FoldChange, r.y, s=16 if mk == 'o' else 12, marker=mk,
                   facecolor=C_UP if r.log2FoldChange > 0 else C_DOWN,
                   edgecolor='black', linewidth=0.6, zorder=4)
    callout(ax, sp, label_fs, **(label_kw or {}))
    if sp.empty:
        ax.text(0.5, 0.80, 'No panel gene with\nadj. $P$ < 0.05', transform=ax.transAxes,
                ha='center', va='top', fontsize=7, color='#555555',
                bbox=dict(boxstyle='square,pad=0.2', fc='white', ec='none'))

    n_up, n_dn = int(up.sum()), int(dn.sum())
    ax.text(0.98, 0.98, f'Up {n_up:,}', transform=ax.transAxes, ha='right', va='top', fontsize=7)
    ax.text(0.02, 0.98, f'Down {n_dn:,}', transform=ax.transAxes, ha='left', va='top', fontsize=7)
    # 色見本（数の左に小さい四角）
    ax.set_title(title, fontsize=8.5, linespacing=1.15, pad=4)
    ax.set_xlabel('log$_2$ fold change')
    ax.set_ylabel('$-$log$_{10}$ adjusted $P$')
    ax.tick_params(length=2.5, width=0.6)
    for s in ('left', 'bottom'):
        ax.spines[s].set_linewidth(0.6)
    return n_up, n_dn, sp


def callout(ax, sp, fs, col_x=0.38, top=0.70, step=None, inner=None):
    """上昇側は右、低下側は左に、ラベルを点の y の高い順に縦一列で並べる（重なり回避を手で決める）。

    col_x: 列の x 位置（x 軸半幅に対する比）。top: 列の先頭の y（軸高さに対する比）。
    step: 行送り（軸高さに対する比）。None ならフォントサイズから決める。
    inner: {+1 or -1: x 位置（半幅比）}。その側のラベル列を点の内側（中央寄り）に置く（外側に余白がないとき）。
    """
    inner = inner or {}
    if sp.empty:
        return
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    h_pt = ax.get_position().height * ax.figure.get_figheight() * 72
    step = step or fs * 1.25 / h_pt
    for sign in (+1, -1):
        g = sp[np.sign(sp.log2FoldChange) == sign].assign(a=lambda t: t.log2FoldChange.abs())
        # y の高い順。同じ高さなら、内側列では外側の点を上に（引き出し線が交差しないように）
        g = g.sort_values(['y', 'a'], ascending=[False, sign not in inner])
        if g.empty:
            continue
        # 列は点より外側に置く（ラベルが自分の点に乗らないように）
        tx = sign * max(col_x * x1, g.log2FoldChange.abs().max() + 0.15 * x1)
        ha = 'left' if sign > 0 else 'right'
        if sign in inner:
            tx = sign * inner[sign] * x1
            ha = 'right' if sign > 0 else 'left'
        for k, (_, r) in enumerate(g.iterrows()):
            ty = y0 + (top - k * step) * (y1 - y0)
            ax.annotate(r.gene, xy=(r.log2FoldChange, r.y), xytext=(tx, ty), textcoords='data',
                        fontsize=fs, style='italic', fontweight='bold' if r.consistent else 'normal',
                        color='black' if r.consistent else '#555555',
                        ha=ha, va='center', zorder=5,
                        bbox=dict(boxstyle='square,pad=0.05', fc='white', ec='none', alpha=0.7),
                        arrowprops=dict(arrowstyle='-', color='#444444', lw=0.4, shrinkA=1, shrinkB=2.5,
                                        relpos=(0, 0.5) if ha == 'left' else (1, 0.5)))


def legend_handles():
    m = dict(linestyle='none', markersize=4)
    return [
        Line2D([], [], marker='o', color=C_UP, markeredgewidth=0, label=f'Up (adj. $P$ < {PADJ}, log$_2$FC ≥ {LFC:g})', **m),
        Line2D([], [], marker='o', color=C_DOWN, markeredgewidth=0, label=f'Down (adj. $P$ < {PADJ}, log$_2$FC ≤ −{LFC:g})', **m),
        Line2D([], [], marker='o', color=C_NS, markeredgewidth=0, label='Not significant', **m),
        Line2D([], [], marker='o', mfc='white', mec='black', mew=0.6, label='Inflammation / osteoclastogenesis gene*', **m),
        Line2D([], [], marker='D', mfc='white', mec='black', mew=0.6, label='Apoptosis / senescence gene*', markersize=3.4,
               linestyle='none'),
    ]


LEGEND_NOTE = ('*Pre-defined panel genes with adj. $P$ < 0.05 (any fold change); fill = direction of change.\n'
               'Bold label = direction consistent with a pro-inflammatory / pro-osteoclastogenic '
               '(or pro-apoptotic / pro-senescent) shift; grey label = opposite direction.')

"""in vitro contrast の volcano・注目遺伝子ヒートマップ・DEG 数・QC（PCA, 群別 CPM）。

  uv run python -I prelim/invitro/figs_iv.py      # repo root から（de_iv.py の後）

閾値は prelim/volcano.py と同じ（|log2FC| >= 1 かつ padj < 0.05、THRESHOLDS も同じ 6 段）。
volcano: 注目遺伝子（ivcontrasts.FOCUS）をカテゴリ色の白抜き丸で全点表示し、padj < 0.05 かつ |log2FC| >= 1 のものにラベル。
出力（out/prelim/invitro/）: volcano/<label>.png|pdf, volcano_all.png|pdf|tiff, focus_heatmap.png|pdf|tiff,
  focus_table.csv, deg_counts.csv, qc/pca.png|pdf, qc/group_cpm.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/invitro'))
import common  # noqa: E402
from volcano import THRESHOLDS, COL_UP, COL_DOWN, COL_NS  # noqa: E402
from ivcontrasts import IC, OUT, FOCUS, EXTRA, CAT  # noqa: E402

OUT = ROOT / OUT
DE = OUT / 'de'
LFC, PADJ = 1.0, 0.05
CATCOL = {'apoptosis': '#1F4E79', 'senescence': '#7D3C98', 'SASP/inflam.': '#B03A2E',
          'OC support': '#B9770E', 'bone formation': '#1E8449', 'ER': '#566573'}


def load(lb):
    return pd.read_csv(DE / f'{lb}.csv').dropna(subset=['log2FoldChange', 'pvalue'])


def save(fig, name, formats=('png', 'pdf', 'tiff')):
    for ext in formats:
        p = OUT / f'{name}.{ext}'
        p.parent.mkdir(parents=True, exist_ok=True)
        kw = {'pil_kwargs': {'compression': 'tiff_lzw'}} if ext == 'tiff' else {}
        fig.savefig(p, dpi=common.DPI, bbox_inches='tight', pad_inches=0.02, **kw)
    print('saved', OUT / f'{name}.png')


def plot(ax, d, title):
    p = d.pvalue.clip(lower=1e-300)
    y = -np.log10(p)
    sig = (d.padj < PADJ) & (d.log2FoldChange.abs() >= LFC)
    up, dn = sig & (d.log2FoldChange > 0), sig & (d.log2FoldChange < 0)
    ax.scatter(d.log2FoldChange[~sig], y[~sig], s=1.2, c=COL_NS, lw=0, rasterized=True)
    ax.scatter(d.log2FoldChange[up], y[up], s=1.5, c=COL_UP, lw=0, rasterized=True, label=f'Up ({up.sum()})')
    ax.scatter(d.log2FoldChange[dn], y[dn], s=1.5, c=COL_DOWN, lw=0, rasterized=True, label=f'Down ({dn.sum()})')
    ymax = np.nanpercentile(y, 99.95) * 1.15
    for cat, genes in FOCUS.items():
        m = d.gene.isin(genes)
        col = CATCOL[cat]
        yy = y[m].clip(upper=ymax)
        ax.scatter(d.log2FoldChange[m], yy, s=10, facecolor='none', edgecolor=col, lw=0.7)
        lab = m & (d.padj < PADJ) & (d.log2FoldChange.abs() >= LFC)
        for gx, gy, gn in zip(d.log2FoldChange[lab], y[lab].clip(upper=ymax), d.gene[lab]):
            ax.annotate(gn, (gx, gy), xytext=(2, 1.5), textcoords='offset points', fontsize=5.5,
                        color=col, fontstyle='italic')
    for v in (-LFC, LFC):
        ax.axvline(v, ls='--', lw=0.5, c='k')
    if (d.padj < PADJ).any():
        ax.axhline(-np.log10(d.pvalue[d.padj < PADJ].max()), ls='--', lw=0.5, c='k')
    xm = max(np.nanpercentile(d.log2FoldChange.abs(), 99.9) * 1.1, 3)
    ax.set_xlim(-xm, xm)
    ax.set_ylim(0, ymax * 1.05)
    ax.set_xlabel('log$_2$ fold change')
    ax.set_ylabel('$-$log$_{10}$ P')
    ax.set_title(title, fontsize=8)
    ax.legend(fontsize=5.5, frameon=False, loc='upper left', markerscale=2, handletextpad=0.2)


def volcanoes():
    common.setup()
    for lb, c in IC.items():
        fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(70)))
        plot(ax, load(lb), c['title'])
        save(fig, f'volcano/{lb}', ('png', 'pdf'))
        plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(common.mm(170), common.mm(120)))
    for ax, (lb, c), pl in zip(axes.flat, IC.items(), 'ABCDE'):
        plot(ax, load(lb), c['short'])
        common.panel_label(ax, pl, dx=-0.25)
    ax = axes.flat[5]
    ax.axis('off')
    for i, (cat, col) in enumerate(CATCOL.items()):
        ax.scatter([0.05], [0.9 - i * 0.12], s=18, facecolor='none', edgecolor=col, lw=0.8,
                   transform=ax.transAxes)
        ax.text(0.12, 0.9 - i * 0.12, f'{cat}: ' + ', '.join(FOCUS[cat]), transform=ax.transAxes,
                fontsize=6, va='center', color=col, wrap=True)
    ax.text(0.0, 0.12, 'Labels: focus genes with padj < 0.05 & |log$_2$FC| ≥ 1.\nDashed: |log$_2$FC| = 1, padj = 0.05.\n'
            'E: TempO-Seq S1500+ panel (1,337 genes tested).', transform=ax.transAxes, fontsize=6, va='top')
    fig.tight_layout()
    save(fig, 'volcano_all')
    plt.close(fig)


def tables():
    rows, L, P = [], {}, {}
    genes = list(CAT)
    for lb in IC:
        d = load(lb)
        for pa, lf in THRESHOLDS:
            s = (d.padj < pa) & (d.log2FoldChange.abs() >= lf)
            rows.append(dict(label=lb, padj=pa, lfc=lf, up=int((s & (d.log2FoldChange > 0)).sum()),
                             down=int((s & (d.log2FoldChange < 0)).sum()), tested=len(d)))
        di = d.drop_duplicates('gene').set_index('gene')
        L[lb] = di.log2FoldChange.reindex(genes)
        P[lb] = di.padj.reindex(genes)
    cnt = pd.DataFrame(rows)
    cnt.to_csv(OUT / 'deg_counts.csv', index=False)
    print(cnt[(cnt.padj == 0.05) & (cnt.lfc == 1.0)].to_string(index=False))
    L, P = pd.DataFrame(L), pd.DataFrame(P)
    cpm = pd.read_csv(OUT / 'qc/group_cpm.csv', index_col=0) if (OUT / 'qc/group_cpm.csv').exists() else None
    t = pd.concat({'log2FC': L, 'padj': P}, axis=1)
    if cpm is not None:
        cpm.columns = pd.MultiIndex.from_tuples([('groupCPM', c) for c in cpm.columns])
        t = t.join(cpm)
    t.insert(0, 'category', [CAT[g] for g in genes])
    t.index.name = 'gene'
    t.to_csv(OUT / 'focus_table.csv')
    return L, P


def focus_heatmap(L, P):
    common.setup()
    genes = list(L.index)
    fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(150)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(L.values), cmap=cmap, vmin=-4, vmax=4, aspect='auto')
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            v, p = L.iloc[i, j], P.iloc[i, j]
            if np.isnan(v):
                continue
            star = '**' if p < 0.05 else ('*' if p < 0.1 else '')
            ax.text(j, i, f'{v:.1f}{star}', ha='center', va='center', fontsize=5.5,
                    color='white' if abs(v) > 3 else 'black')
    ax.set_yticks(range(len(genes)), [f'$\\it{{{g}}}$' for g in genes], fontsize=6.5)
    ax.set_xticks(range(L.shape[1]), [IC[c]['short'] for c in L.columns], rotation=45, ha='right', fontsize=6.5)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    prev = None
    for i, g in enumerate(genes):
        if prev is not None and CAT[g] != prev:
            ax.axhline(i - 0.5, c='k', lw=0.5)
            ax.text(L.shape[1] - 0.4, i - 0.5, '', fontsize=6)
        prev = CAT[g]
    # カテゴリ名を右端に
    cats = pd.Series([CAT[g] for g in genes])
    for cat, idx in cats.groupby(cats).groups.items():
        ax.text(L.shape[1] - 0.35, np.mean(list(idx)), cat, fontsize=6, va='center', ha='left',
                color=CATCOL.get(cat, '#333333'))
    ax.axvline(3.5, c='k', lw=1.0)
    fig.tight_layout()
    cax = fig.add_axes([0.80, 0.03, 0.04, 0.10])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('log$_2$FC (treated / control)', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title('Focus genes (** padj<0.05, * <0.1;\ngrey: filtered / not on TempO-Seq panel)', fontsize=7, loc='left')
    save(fig, 'focus_heatmap')
    plt.close(fig)


def qc():
    common.setup()
    (OUT / 'qc').mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(common.mm(170), common.mm(55)))
    gc = []
    panels = [('GSE284666', 'NHO', 'NHO'), ('GSE284666', 'BMMSC', 'BM-MSC'), ('GSE135323', None, 'TE-Bone')]
    for ax, (gse, cell, title) in zip(axes, panels):
        c = pd.read_csv(ROOT / f'data/seq/{gse}/counts.tsv.gz', sep='\t', index_col=0)
        s = pd.read_csv(ROOT / f'data/seq/{gse}/samples.csv')
        if cell:
            s = s[s.cell == cell]
        else:
            s = s[s.group.isin(['Dox', 'DoxCtrl'])]
        x = c[s['sample']]
        cpm = x / x.sum() * 1e6
        g = cpm.T.groupby(s.group.values).mean().T
        g.columns = [f'{gse}:{k}' for k in g.columns]
        gc.append(g.reindex(list(CAT)))
        lg = np.log2(cpm[(cpm >= 10).sum(axis=1) >= 3] + 1)
        lg = lg.loc[lg.var(axis=1).sort_values(ascending=False).index[:2000]]
        z = (lg.T - lg.T.mean()).values
        u, sv, _ = np.linalg.svd(z, full_matrices=False)
        pc = u * sv
        ve = sv ** 2 / (sv ** 2).sum()
        for grp, col in zip(sorted(s.group.unique()), ['#4A7BA7', '#C44E52', '#55A868']):
            m = (s.group == grp).values
            ax.scatter(pc[m, 0], pc[m, 1], s=14, c=col, label=grp, lw=0)
        ax.set_xlabel(f'PC1 ({ve[0]:.0%})')
        ax.set_ylabel(f'PC2 ({ve[1]:.0%})')
        ax.set_title(f'{title} (top 2,000 var. genes)', fontsize=8)
        ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'qc/pca.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    pd.concat(gc, axis=1).round(1).to_csv(OUT / 'qc/group_cpm.csv')


def main():
    qc()
    volcanoes()
    L, P = tables()
    focus_heatmap(L, P)


if __name__ == '__main__':
    main()

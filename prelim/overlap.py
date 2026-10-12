"""主 3 contrast（A: GSE222752 cortical OVX vs Sham / B: GSE205237 Dox / C: GSE276391 DEX）の重なり。

- 遺伝子レベル: 閾値 3 段階で up/down 別の集合を作り、Venn 相当の領域ごとの数と遺伝子リスト
    strict : padj < 0.05 & |log2FC| >= 1
    lenient: padj < 0.1  & |log2FC| >= 0.585
    nominal: pvalue < 0.01 & |log2FC| >= 0.585（探索用。多重比較補正なし）
- 各ペアの log2FC 散布と相関（共通に検定された遺伝子、Spearman / Pearson）
- パスウェイレベル: Hallmark prerank GSEA の NES を全 contrast 横並び（heatmap）

出力: out/prelim/overlap/
"""
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from contrasts import C, MAIN  # noqa: E402
from genes import INFLAM, APOP  # noqa: E402

DE = Path(common.OUT_DIR) / 'de'
EN = Path(common.OUT_DIR) / 'enrich'
OUT = Path(common.OUT_DIR) / 'overlap'
THR = {
    'strict': lambda d: (d.padj < 0.05) & (d.log2FoldChange.abs() >= 1),
    'lenient': lambda d: (d.padj < 0.1) & (d.log2FoldChange.abs() >= 0.585),
    'nominal': lambda d: (d.pvalue < 0.01) & (d.log2FoldChange.abs() >= 0.585),
}


def load(lb):
    d = pd.read_csv(DE / f'{lb}.csv').dropna(subset=['gene', 'log2FoldChange'])
    return d.drop_duplicates('gene').set_index('gene')


def main():
    common.setup()
    OUT.mkdir(parents=True, exist_ok=True)
    D = {k: load(lb) for k, lb in MAIN.items()}
    common_genes = sorted(set.intersection(*[set(d.index) for d in D.values()]))
    print('genes tested in all 3:', len(common_genes))

    venn, lists = [], []
    for tn, f in THR.items():
        for direction, sign in [('up', 1), ('down', -1)]:
            S = {k: set(d.index[f(d).fillna(False) & (np.sign(d.log2FoldChange) == sign)])
                 for k, d in D.items()}
            # Venn の 7 領域（共通に検定された遺伝子に限らず、各 contrast の DEG 集合そのもの）
            for r in range(1, 4):
                for combo in combinations('ABC', r):
                    inside = set.intersection(*[S[k] for k in combo])
                    outside = set().union(*[S[k] for k in 'ABC' if k not in combo])
                    region = sorted(inside - outside)
                    venn.append({'threshold': tn, 'direction': direction,
                                 'region': '&'.join(combo) + '_only' if r < 3 else 'A&B&C',
                                 'n': len(region)})
                    if r >= 2:
                        for g in region:
                            lists.append({'threshold': tn, 'direction': direction,
                                          'region': '&'.join(combo), 'gene': g,
                                          **{f'lfc_{k}': D[k].log2FoldChange.get(g) for k in 'ABC'},
                                          **{f'padj_{k}': D[k].padj.get(g) for k in 'ABC'}})
            for k in 'ABC':
                venn.append({'threshold': tn, 'direction': direction, 'region': f'{k}_total',
                             'n': len(S[k])})
    v = pd.DataFrame(venn)
    v.to_csv(OUT / 'venn_counts.csv', index=False)
    pd.DataFrame(lists).to_csv(OUT / 'shared_genes.csv', index=False)
    pd.set_option('display.width', 200)
    print(v.pivot_table(index='region', columns=['threshold', 'direction'], values='n',
                        aggfunc='sum').fillna(0).astype(int).to_string())

    # ペア散布
    fig, axes = plt.subplots(1, 3, figsize=(common.mm(170), common.mm(58)))
    cors = []
    panel = set(INFLAM + APOP)
    for ax, (a, b) in zip(axes, combinations('ABC', 2)):
        x = D[a].log2FoldChange.reindex(common_genes)
        y = D[b].log2FoldChange.reindex(common_genes)
        ok = x.notna() & y.notna()
        x, y = x[ok], y[ok]
        rs, ps = stats.spearmanr(x, y)
        rp, pp = stats.pearsonr(x, y)
        cors.append({'pair': f'{a}-{b}', 'n': len(x), 'spearman': rs, 'spearman_p': ps,
                     'pearson': rp, 'pearson_p': pp})
        ax.scatter(x, y, s=1, c='#AAAAAA', lw=0, rasterized=True)
        m = x.index.isin(panel)
        ax.scatter(x[m], y[m], s=8, facecolor='none', edgecolor='#B03A2E', lw=0.6)
        for g in x.index[m]:
            ax.annotate(g, (x[g], y[g]), fontsize=5, fontstyle='italic', color='#B03A2E')
        lim = np.nanpercentile(np.abs(np.r_[x, y]), 99.9)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.axhline(0, lw=0.4, c='k')
        ax.axvline(0, lw=0.4, c='k')
        ax.set_xlabel(f'{a}: {MAIN[a]} log$_2$FC')
        ax.set_ylabel(f'{b}: {MAIN[b]} log$_2$FC')
        ax.set_title(f'Spearman ρ = {rs:.2f} (n = {len(x)})', fontsize=8)
    fig.tight_layout()
    common.save(fig, 'overlap/lfc_scatter', formats=('png', 'pdf'))
    plt.close(fig)
    cr = pd.DataFrame(cors)
    cr.to_csv(OUT / 'lfc_correlation.csv', index=False)
    print(cr.round(4).to_string(index=False))

    # Hallmark NES heatmap（全 contrast）
    nes, fdr = {}, {}
    for lb in C:
        p = EN / f'gsea_{lb}.csv'
        if p.exists():
            g = pd.read_csv(p)
            nes[lb] = g.set_index('Term').NES.astype(float)
            fdr[lb] = g.set_index('Term')['FDR q-val'].astype(float)
    if nes:
        N = pd.DataFrame(nes)
        F = pd.DataFrame(fdr)
        N.to_csv(OUT / 'hallmark_nes.csv')
        F.to_csv(OUT / 'hallmark_fdr.csv')
        order = N[list(MAIN.values())].mean(axis=1).sort_values().index
        N, F = N.loc[order], F.loc[order]
        fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(200)))
        vmax = np.nanmax(np.abs(N.values))
        im = ax.imshow(N.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
        for i in range(N.shape[0]):
            for j in range(N.shape[1]):
                f = F.values[i, j]
                if f < 0.05:
                    ax.text(j, i, '**', ha='center', va='center', fontsize=6)
                elif f < 0.25:
                    ax.text(j, i, '*', ha='center', va='center', fontsize=6)
        ax.set_xticks(range(N.shape[1]))
        ax.set_xticklabels(N.columns, rotation=45, ha='right', fontsize=6)
        ax.set_yticks(range(N.shape[0]))
        ax.set_yticklabels(N.index, fontsize=6)
        cb = fig.colorbar(im, ax=ax, shrink=0.3)
        cb.set_label('NES')
        ax.set_title('Hallmark prerank GSEA (* FDR<0.25, ** FDR<0.05)')
        fig.tight_layout()
        common.save(fig, 'overlap/hallmark_nes', formats=('png', 'pdf'))


if __name__ == '__main__':
    main()

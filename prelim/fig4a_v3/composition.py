"""組成（細胞種マーカー）: デノスマブ・閉経後で何が「組成の変化」か。

  uv run python -I prelim/fig4a_v3/composition.py

対象:
- GSE141595 post30（消化後、閉経後 placebo 15 / denosumab 15、全遺伝子）
- GSE141610（遠心のみ、同じ 30 人の別の生検片。pairing.py 参照）
- GSE141595 共通遺伝子版（閉経前 15 を含む 45 例。閉経前ファイルの欠落で chr12 後半〜22 のマーカーは無い）
- GSE72815（若年 19 / 高齢 19、バッチの揃った独立コホート。比較の基準）
スコア = カテゴリのマーカーの log2(CPM+1) 平均（common4a.MARKERS_H）。群差は中央値差と Mann-Whitney、Hedges g。
骨髄スコア（HBB/HBA1/MPO/ELANE の z）は DE の共変量に使う（dmab_de.py）。

出力: out/prelim/fig4a_v3/composition/{scores_long.csv, composition_tests.csv, marker_genes.csv, composition.png|pdf}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common4a import OUT, counts, scores, MARKERS_H, cpm  # noqa: E402
import common  # noqa: E402  prelim/common.py

D = OUT / 'composition'
SETS = [  # (counts dir, 表示名, [(test, ref)])
    ('GSE141595/post30', 'GSE141595 digested (post30)', [('Post_Dmab', 'Post')]),
    ('GSE141610', 'GSE141610 centrifuged', [('Post_Dmab', 'Post')]),
    ('GSE141595', 'GSE141595 digested (common genes, 45)', [('Post', 'Pre'), ('Post_Dmab', 'Post')]),
    ('GSE72815', 'GSE72815 whole biopsy', [('Old', 'Young'), ('Old_E', 'Old')]),
]


def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    g = (a.mean() - b.mean()) / sp if sp > 0 else np.nan
    return g * (1 - 3 / (4 * (na + nb) - 9))


def main():
    D.mkdir(parents=True, exist_ok=True)
    longs, tests, genes = [], [], []
    for sub, name, comps in SETS:
        x, s = counts(sub)
        sc, used = scores(x)
        c = cpm(x)
        for k, gs in used.items():
            for g in gs:
                for grp in s.group.unique():
                    v = c.loc[g, s.index[s.group == grp]]
                    genes.append(dict(set=sub, cat=k, gene=g, group=grp, cpm_median=v.median()))
        L = sc.stack().rename('score').reset_index().rename(columns={'level_0': 'sample', 'level_1': 'cat'})
        L['group'] = s.group.reindex(L['sample']).values
        L['set'] = sub
        longs.append(L)
        for t, r in comps:
            for k in sc.columns:
                a, b = sc.loc[s.index[s.group == t], k], sc.loc[s.index[s.group == r], k]
                tests.append(dict(set=sub, name=name, test=t, ref=r, cat=k, genes=','.join(used[k]),
                                  n_genes=len(used[k]), n_markers_defined=len(MARKERS_H[k]),
                                  median_test=a.median(), median_ref=b.median(),
                                  diff_log2=a.median() - b.median(), hedges_g=hedges_g(a, b),
                                  mw_p=mannwhitneyu(a, b).pvalue))
    T = pd.DataFrame(tests)
    from statsmodels.stats.multitest import multipletests  # noqa: E402
    T['mw_fdr'] = np.nan
    for key, idx in T.groupby(['set', 'test', 'ref']).groups.items():
        T.loc[idx, 'mw_fdr'] = multipletests(T.loc[idx, 'mw_p'], method='fdr_bh')[1]
    T.to_csv(D / 'composition_tests.csv', index=False)
    pd.concat(longs).to_csv(D / 'scores_long.csv', index=False)
    pd.DataFrame(genes).to_csv(D / 'marker_genes.csv', index=False)
    show = T[['set', 'test', 'ref', 'cat', 'n_genes', 'diff_log2', 'hedges_g', 'mw_p', 'mw_fdr']]
    with pd.option_context('display.width', 200, 'display.max_rows', 200):
        print(show.round(3).to_string(index=False))
    plot(T)


def plot(T):
    common.setup()
    cats = list(MARKERS_H)
    cols = [('GSE141595/post30', 'Post_Dmab', 'Post', 'Dmab vs\nplacebo\n(digested)'),
            ('GSE141610', 'Post_Dmab', 'Post', 'Dmab vs\nplacebo\n(centrifuged)'),
            ('GSE141595', 'Post', 'Pre', 'Postmeno. vs\npremeno.†\n(digested)'),
            ('GSE72815', 'Old', 'Young', 'Old vs\nyoung\n(GSE72815)')]
    fig, ax = plt.subplots(figsize=(common.mm(120), common.mm(115)))
    M = np.full((len(cats), len(cols)), np.nan)
    P = np.full_like(M, np.nan)
    for j, (sub, t, r, _) in enumerate(cols):
        for i, k in enumerate(cats):
            q = T[(T.set == sub) & (T.test == t) & (T.ref == r) & (T.cat == k)]
            if len(q):
                M[i, j], P[i, j] = q.diff_log2.iloc[0], q.mw_fdr.iloc[0]
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-1.5, vmax=1.5, aspect='auto')
    for i in range(len(cats)):
        for j in range(len(cols)):
            if np.isnan(M[i, j]):
                ax.text(j, i, 'n/a', ha='center', va='center', fontsize=6, color='#888888')
                continue
            star = '**' if P[i, j] < 0.01 else '*' if P[i, j] < 0.05 else ''
            ax.text(j, i, f'{M[i, j]:+.2f}{star}'.replace('-', '−'), ha='center', va='center', fontsize=6,
                    color='white' if abs(M[i, j]) > 1.0 else '#222222')
    ax.set_yticks(range(len(cats)), cats, fontsize=7)
    ax.set_xticks(range(len(cols)), [c[3] for c in cols], fontsize=6.5)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.02)
    cb.set_label('median difference of\nmarker score (log2)', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title('Cell-type marker scores (Mann-Whitney FDR * < 0.05, ** < 0.01)', fontsize=7, loc='left', pad=40)
    fig.text(0.02, 0.005, '† batch-confounded; only markers on genes present in the premenopausal file (see composition_tests.csv n_genes)',
             fontsize=5.5, color='#444444')
    from common4a import save
    save(fig, 'composition', 'composition')
    plt.close(fig)


if __name__ == '__main__':
    main()

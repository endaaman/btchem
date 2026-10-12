"""ヒトデータの QC: PCA（log2 CPM, 上位 500 変動遺伝子）と組成マーカー（骨細胞・骨芽・破骨・骨髄/血液・筋・脂肪・性・ER/GC 応答）。

  uv run python -I prelim/human/qc_h.py
出力: out/prelim/human/qc/
  pca.png|pdf          データセットごとの PCA
  markers.csv          サンプル x マーカーの CPM（long）
  composition.csv      カテゴリ別スコア（マーカーの log2(CPM+1) 平均）の群中央値と Mann-Whitney P（各比較）
  libsize.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common  # noqa: E402

SEQ = Path('data/seq')
OUT = Path('out/prelim/human/qc')
SETS = {  # QC 単位: (counts dir, 表示名)
    'GSE141595': 'GSE141595 (common genes, 45)',
    'GSE141595/post30': 'GSE141595 post30 (all genes)',
    'GSE141610': 'GSE141610 (centrifuged)',
    'GSE72815': 'GSE72815',
    'GSE276529': 'GSE276529',
}
MARKERS = {
    'osteocyte': ['SOST', 'DMP1', 'PHEX', 'MEPE', 'FGF23'],
    'osteoblast': ['BGLAP', 'COL1A1', 'SP7', 'RUNX2', 'IBSP'],
    'osteoclast': ['ACP5', 'CTSK', 'MMP9'],
    'marrow/blood': ['PTPRC', 'HBB', 'HBA1', 'CD19', 'MPO', 'CD79A', 'ELANE'],
    'muscle': ['ACTA1', 'MYH1', 'MYH2', 'MYH7', 'CKM', 'DES', 'MB'],
    'adipocyte': ['ADIPOQ', 'PLIN1', 'LEP', 'CIDEA'],
    'sex': ['XIST', 'RPS4Y1', 'DDX3Y'],
    'ER target': ['PGR', 'GREB1'],
    'GC response': ['FKBP5', 'TSC22D3', 'PER1'],
}
COMPARE = {  # QC 単位ごとの比較（test, ref）
    'GSE141595': [('Post', 'Pre'), ('Post_Dmab', 'Post')],
    'GSE141595/post30': [('Post_Dmab', 'Post')],
    'GSE141610': [('Post_Dmab', 'Post')],
    'GSE72815': [('Old', 'Young'), ('Old_E', 'Old')],
    'GSE276529': [('GIOP', 'NonGIOP')],
}


def main():
    common.setup()
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(common.mm(170), common.mm(110)))
    rows, libs, comp = [], [], []
    for ax, (g, name) in zip(axes.flat, SETS.items()):
        x = pd.read_csv(SEQ / g / 'counts.tsv.gz', sep='\t', index_col=0)
        s = pd.read_csv(SEQ / g / 'samples.csv', dtype=str).set_index('sample')
        lab = s.group
        cpm = x / x.sum() * 1e6
        for c in x.columns:
            libs.append({'set': g, 'sample': c, 'group': lab[c], 'libsize': int(x[c].sum())})
        score = {}
        for cat, gs in MARKERS.items():
            have = [gn for gn in gs if gn in cpm.index]
            for gene in gs:
                for c in x.columns:
                    rows.append({'set': g, 'sample': c, 'group': lab[c], 'cat': cat, 'gene': gene,
                                 'cpm': cpm.loc[gene, c] if gene in cpm.index else np.nan})
            if have:
                score[cat] = np.log2(cpm.loc[have] + 1).mean()
        score = pd.DataFrame(score)
        for t, r in COMPARE[g]:
            for cat in score:
                a, b = score.loc[lab == t, cat], score.loc[lab == r, cat]
                comp.append({'set': g, 'test': t, 'ref': r, 'cat': cat,
                             'genes': ','.join(gn for gn in MARKERS[cat] if gn in cpm.index),
                             'median_test': a.median(), 'median_ref': b.median(),
                             'diff_log2': a.median() - b.median(),
                             'mw_p': mannwhitneyu(a, b).pvalue})
        lc = np.log2(cpm[(x.sum(axis=1) >= 10)] + 1)
        top = lc.loc[lc.var(axis=1).sort_values(ascending=False).index[:500]]
        z = (top.T - top.T.mean()).values
        u, sv, _ = np.linalg.svd(z, full_matrices=False)
        pc = u * sv
        ve = sv ** 2 / (sv ** 2).sum()
        for i, grp in enumerate(sorted(lab.unique())):
            m = (lab == grp).values
            ax.scatter(pc[m, 0], pc[m, 1], s=10, label=grp, color=f'C{i}')
        if len(x.columns) <= 12:
            for i, c in enumerate(x.columns):
                ax.annotate(c, (pc[i, 0], pc[i, 1]), fontsize=5, alpha=.7)
        ax.set_title(name, fontsize=8)
        ax.set_xlabel(f'PC1 ({ve[0]*100:.0f}%)')
        ax.set_ylabel(f'PC2 ({ve[1]*100:.0f}%)')
        ax.legend(fontsize=5.5, frameon=False)
    axes.flat[-1].axis('off')
    fig.tight_layout()
    common.save(fig, 'human/qc/pca', formats=('png', 'pdf'))
    m = pd.DataFrame(rows)
    m.to_csv(OUT / 'markers.csv', index=False)
    pd.DataFrame(libs).to_csv(OUT / 'libsize.csv', index=False)
    c = pd.DataFrame(comp)
    c.to_csv(OUT / 'composition.csv', index=False)
    pd.set_option('display.width', 250)
    print(c.round(3).to_string())
    print(m.groupby(['set', 'group', 'gene']).cpm.median().unstack(['set', 'group']).round(1).to_string())


if __name__ == '__main__':
    main()

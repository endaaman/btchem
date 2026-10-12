"""データセットごとの QC: PCA（log2 CPM, 上位 500 変動遺伝子）と細胞組成マーカー・性別マーカー。

出力:
  out/prelim/qc/pca.png|pdf          データセットごとの PCA
  out/prelim/qc/markers.csv          サンプル x マーカーの CPM
  out/prelim/qc/libsize.csv          ライブラリサイズ
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402

SEQ = Path('data/seq')
OUT = Path(common.OUT_DIR) / 'qc'
GSES = ['GSE222752', 'GSE238244', 'GSE205237', 'GSE296757', 'GSE276391']
MARKERS = {
    'osteocyte': ['Sost', 'Dmp1', 'Phex', 'Mepe'],
    'osteoblast': ['Bglap', 'Col1a1', 'Sp7'],
    'osteoclast': ['Acp5', 'Ctsk'],
    'marrow/blood': ['Ptprc', 'Hbb-bs', 'Hba-a1', 'Cd19', 'Mpo'],
    'muscle': ['Acta1', 'Myh4', 'Ckm'],
    'sex': ['Xist', 'Ddx3y', 'Kdm5d', 'Uty'],
    'GC response': ['Fkbp5', 'Tsc22d3', 'Per1'],
    'ER target': ['Pgr', 'Greb1'],
}


def main():
    common.setup()
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(common.mm(170), common.mm(110)))
    rows, libs = [], []
    for ax, g in zip(axes.flat, GSES):
        x = pd.read_csv(SEQ / g / 'counts.tsv.gz', sep='\t', index_col=0)
        s = pd.read_csv(SEQ / g / 'samples.csv', dtype=str).set_index('sample')
        lab = s.group if g != 'GSE222752' else s.compartment.str[:4] + ':' + s.group
        if g == 'GSE296757':
            lab = s.group + ' b' + s.batch
        cpm = x / x.sum() * 1e6
        for c in x.columns:
            libs.append({'gse': g, 'sample': c, 'group': lab[c], 'libsize': int(x[c].sum())})
        for cat, gs in MARKERS.items():
            for gene in gs:
                if gene in cpm.index:
                    for c in x.columns:
                        rows.append({'gse': g, 'sample': c, 'group': lab[c], 'cat': cat,
                                     'gene': gene, 'cpm': cpm.loc[gene, c]})
        lc = np.log2(cpm[(x.sum(axis=1) >= 10)] + 1)
        top = lc.loc[lc.var(axis=1).sort_values(ascending=False).index[:500]]
        z = (top.T - top.T.mean()).values
        u, sv, vt = np.linalg.svd(z, full_matrices=False)
        pc = u * sv
        ve = sv ** 2 / (sv ** 2).sum()
        for i, grp in enumerate(sorted(lab.unique())):
            m = (lab == grp).values
            ax.scatter(pc[m, 0], pc[m, 1], s=14, label=grp, color=f'C{i}')
        for i, c in enumerate(x.columns):
            ax.annotate(c, (pc[i, 0], pc[i, 1]), fontsize=5, alpha=.7)
        ax.set_title(g)
        ax.set_xlabel(f'PC1 ({ve[0]*100:.0f}%)')
        ax.set_ylabel(f'PC2 ({ve[1]*100:.0f}%)')
        ax.legend(fontsize=5, frameon=False)
    axes.flat[-1].axis('off')
    fig.tight_layout()
    common.save(fig, 'qc/pca', formats=('png', 'pdf'))
    m = pd.DataFrame(rows)
    m.to_csv(OUT / 'markers.csv', index=False)
    pd.DataFrame(libs).to_csv(OUT / 'libsize.csv', index=False)
    w = m.pivot_table(index=['gse', 'group', 'sample'], columns='gene', values='cpm')
    cols = [g for gs in MARKERS.values() for g in gs if g in w.columns]
    pd.set_option('display.width', 250)
    print(w[cols].round(1).to_string())


if __name__ == '__main__':
    main()

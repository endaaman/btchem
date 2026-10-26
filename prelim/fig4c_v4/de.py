"""4C 案1: GSE304562（MC3T3-E1 × コルチコステロン 100 nM 24 h vs 溶媒、各 3 well）の DESeq2 と QC。

  uv run python -I prelim/fig4c_v4/de.py      # repo root から（prepare.py の後）

DE は既存 prelim/de.py の run() をそのまま呼ぶ（pydeseq2、~ cond、MLE、edgeR filterByExpr 相当の低発現フィルタ）。
contrast は prelim/contrasts.py に登録せず、ここで dict を渡す（既存の登録簿は触らない）。

出力（out/prelim/fig4c_v4/）:
  de/Cort_MC3T3.csv        gene, baseMean, log2FoldChange, lfcSE, pvalue, padj（Cort / Ctrl）
  qc/pca.png|pdf           log2 CPM 上位 500 変動遺伝子の PCA
  qc/markers.csv           GC 応答・骨芽/骨細胞・FGF23 などのマーカー CPM（サンプル別）と DE 結果
  deg_counts.csv           閾値ごとの DEG 数（prelim/volcano.py の THRESHOLDS）
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
import common  # noqa: E402
import de as DE  # noqa: E402   prelim/de.py（run() だけ使う）
from volcano import THRESHOLDS  # noqa: E402

OUT = ROOT / 'out/prelim/fig4c_v4'
GSE = 'GSE304562'
CONTRAST = dict(label='Cort_MC3T3', gse=GSE, test='Cort', ref='Ctrl', subset={}, covariates=[],
                title='Corticosterone vs vehicle (MC3T3-E1, 24 h)')

# 事前に決めたマーカー（結果を見て足し引きしない）
MARKERS = {
    'GC response (positive control)': ['Fkbp5', 'Tsc22d3', 'Per1', 'Zbtb16', 'Sgk1', 'Dusp1', 'Klf15', 'Glul'],
    'paper target': ['Fgf23', 'Nr3c1', 'Hsd11b1', 'Vdr', 'Cyp24a1'],
    'osteoblast': ['Runx2', 'Sp7', 'Alpl', 'Col1a1', 'Col1a2', 'Bglap', 'Bglap2', 'Ibsp', 'Spp1', 'Sparc'],
    'osteocyte': ['Dmp1', 'Phex', 'Mepe', 'Sost', 'Pdpn'],
    'RANKL/OPG': ['Tnfsf11', 'Tnfrsf11b', 'Csf1'],
    'Wnt inhibitor': ['Dkk1', 'Sost', 'Sfrp1', 'Sfrp4'],
    'adipogenic': ['Pparg', 'Cebpa', 'Adipoq', 'Lpl', 'Fabp4'],
    'proliferation': ['Mki67', 'Top2a', 'Ccnb1'],
}


def pca(x, s):
    cpm = x / x.sum() * 1e6
    lc = np.log2(cpm[x.sum(axis=1) >= 10] + 1)
    top = lc.loc[lc.var(axis=1).sort_values(ascending=False).index[:500]]
    z = (top.T - top.T.mean()).values
    u, sv, _ = np.linalg.svd(z, full_matrices=False)
    pc = u * sv
    ve = sv ** 2 / (sv ** 2).sum()
    common.setup()
    fig, ax = plt.subplots(figsize=(common.mm(70), common.mm(60)))
    for grp, col in [('Ctrl', common.COLOR_LRO), ('Cort', common.COLOR_HRO)]:
        m = (s.group == grp).values
        ax.scatter(pc[m, 0], pc[m, 1], s=18, color=col, label=grp)
    for i, c in enumerate(x.columns):
        ax.annotate(c, (pc[i, 0], pc[i, 1]), fontsize=6, xytext=(2, 2), textcoords='offset points')
    ax.set_xlabel(f'PC1 ({ve[0]*100:.0f}%)')
    ax.set_ylabel(f'PC2 ({ve[1]*100:.0f}%)')
    ax.set_title('GSE304562 MC3T3-E1', fontsize=8)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    (OUT / 'qc').mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'qc/pca.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    return ve


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'de').mkdir(parents=True, exist_ok=True)
    r = DE.run(CONTRAST)
    r.to_csv(OUT / 'de/Cort_MC3T3.csv', index=False)

    x = pd.read_csv(ROOT / 'data/seq' / GSE / 'counts.tsv.gz', sep='\t', index_col=0)
    s = pd.read_csv(ROOT / 'data/seq' / GSE / 'samples.csv', dtype=str).set_index('sample').loc[x.columns]
    ve = pca(x, s)
    print(f'PCA: PC1 {ve[0]*100:.1f}%, PC2 {ve[1]*100:.1f}%')

    cpm = x / x.sum() * 1e6
    di = r.set_index('gene')
    rows = []
    for cat, gs in MARKERS.items():
        for g in gs:
            row = dict(cat=cat, gene=g)
            for c in x.columns:
                row[f'cpm_{c}'] = cpm.loc[g, c] if g in cpm.index else np.nan
            row['log2FC'] = di.log2FoldChange.get(g, np.nan)
            row['padj'] = di.padj.get(g, np.nan)
            rows.append(row)
    m = pd.DataFrame(rows)
    m.to_csv(OUT / 'qc/markers.csv', index=False)
    pd.set_option('display.width', 250)
    print(m.round(3).to_string(index=False))

    d = r.dropna(subset=['log2FoldChange', 'pvalue'])
    cnt = []
    for pa, lf in THRESHOLDS:
        sg = (d.padj < pa) & (d.log2FoldChange.abs() >= lf)
        cnt.append(dict(label=CONTRAST['label'], padj=pa, lfc=lf, up=int((sg & (d.log2FoldChange > 0)).sum()),
                        down=int((sg & (d.log2FoldChange < 0)).sum()), tested=len(d)))
    cnt = pd.DataFrame(cnt)
    cnt.to_csv(OUT / 'deg_counts.csv', index=False)
    print(cnt.to_string(index=False))


if __name__ == '__main__':
    main()

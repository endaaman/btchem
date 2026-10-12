"""化療データを並べたパネル遺伝子 log2FC ヒートマップ（4B 案 (c) の叩き台）。

bulk（DESeq2）は padj で印、scRNA（GSE289491, n=1）は擬似バルク log2FC のみ（印なし・斜線枠で区別）。
  uv run python -I prelim/chemo/panel_heatmap.py
出力: out/prelim/chemo/panel_heatmap.png|pdf, panel_heatmap.csv
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

OUT = ROOT / 'out/prelim/chemo'
GENES = {
    'Apoptosis / p53': ['Fas', 'Casp3', 'Bax', 'Bbc3', 'Pmaip1', 'Cdkn1a', 'Mdm2', 'Phlda3', 'Eda2r'],
    'Senescence / SASP': ['Cdkn2a', 'Serpine1', 'Il6', 'Mki67'],
    'Inflam. / OC': ['Tnfsf11', 'Tnfrsf11b', 'Il1b', 'Ccl2', 'Cxcl12'],
}
BULK = [  # (列名, csv, 説明)
    ('GSE205237 bone ♂ 30 d', ROOT / 'out/prelim/de/B_Dox.csv'),
    ('GSE296757 OB ♀ 8 wk (4x)', ROOT / 'out/prelim/de/B2_Dox4x_OB.csv'),
    ('GSE296758 MSPC ♀ 8 wk (4x)', OUT / 'de/B3_Dox4x_MSPC.csv'),
    ('GSE296758 MSPC ♀ 8 wk (1x)', OUT / 'de/B3_Dox1x_MSPC.csv'),
]
SC = ['CAR (LepR+)', 'Osteoblast', 'Osteoprogenitor', 'Chondrocyte']
EXTRA = [  # 追加 bulk（ヒト in vitro, day 0-7 連続曝露, IC50）
    ('GSE324491 hESC-OB 5-FU 7 d', OUT / 'de/H_FU5_hESC-OB.csv'),
    ('GSE324491 hESC-OB CPA 7 d', OUT / 'de/H_CPA_hESC-OB.csv'),
]


def main():
    common.setup()
    genes = sum(GENES.values(), [])
    lfc, sig = {}, {}
    sc = pd.read_csv(OUT / 'sc/pseudobulk_all_genes.csv.gz')
    for ct in SC:
        d = sc[sc.celltype == ct].set_index('gene')
        v = d.log2FC.reindex(genes)
        v[d.meanCPM.reindex(genes) < 10] = np.nan  # 低発現は空欄
        lfc[f'GSE289491 {ct} ♀ 96 h (n=1)'] = v
        sig[f'GSE289491 {ct} ♀ 96 h (n=1)'] = pd.Series('', index=genes)
    for name, f in BULK + EXTRA:
        if not Path(f).exists():
            continue
        d = pd.read_csv(f).set_index('gene')
        lfc[name] = d.log2FoldChange.reindex(genes)
        p = d.padj.reindex(genes)
        sig[name] = np.where(p < 0.05, '**', np.where(p < 0.1, '*', ''))
    L = pd.DataFrame(lfc)
    S = pd.DataFrame(sig, index=genes)
    L.to_csv(OUT / 'panel_heatmap.csv')
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(125)))
    m = np.ma.masked_invalid(L.values)
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(m, cmap=cmap, vmin=-2, vmax=2, aspect='auto')
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            if S.iloc[i, j]:
                ax.text(j, i, S.iloc[i, j], ha='center', va='center', fontsize=7)
    ax.set_yticks(range(len(genes)), [f'$\\it{{{g}}}$' for g in genes])
    ax.set_xticks(range(L.shape[1]), L.columns, rotation=40, ha='right')
    y = -0.5
    for k, v in GENES.items():
        y += len(v)
        ax.axhline(y, c='k', lw=0.6)
    ax.axvline(len(SC) - 0.5, c='k', lw=1.0)
    cb = fig.colorbar(im, ax=ax, shrink=0.6)
    cb.set_label('log$_2$ FC (drug / control)')
    ax.set_title('Chemotherapy and bone-resident cells: panel genes (Dox unless noted)\n'
                 '(left: scRNA pseudobulk, n=1 per arm, no test; right: bulk DESeq2, ** padj<0.05, * <0.1; grey = low expression)',
                 fontsize=8)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'panel_heatmap.{ext}', dpi=common.DPI, bbox_inches='tight')
    print(L.round(2).to_string())


if __name__ == '__main__':
    main()

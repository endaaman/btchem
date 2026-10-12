"""GSE289491: 細胞種ラベル付け → 細胞種 x 群の擬似バルクで Dox vs Veh を記述する。

各群 1 匹（n=1）なので統計検定はしない。出すのは
  - 擬似バルク CPM（細胞種ごとに UMI を合算）と log2FC = log2((CPM_Dox+1)/(CPM_Veh+1))
  - 発現細胞割合（UMI>0 の細胞の割合）
  - 細胞ブートストラップ（細胞を群内で復元抽出、500 回）の log2FC 2.5–97.5%:
    「同じ 1 匹から細胞をどれだけ拾ったかの揺らぎ」だけで、個体差（生物学的反復）は含まない
  - 細胞種内での全遺伝子の log2FC 分布におけるパーセンタイル（その遺伝子の変化がどれくらい外れているか）
  - 事前定義の p53 標的スコア・SASP スコア（細胞ごとの平均 log1p 発現）

  uv run python -I prelim/chemo/sc_pseudobulk.py
出力: out/prelim/chemo/sc/ 以下（panel_pseudobulk.csv, celltype_counts.csv, umap.png, dotplot_celltype.png,
      ma_<celltype>.png, scores.csv, partek_panel.csv）
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
import common  # noqa: E402
from genes import ALL_GENES  # noqa: E402

OUT = ROOT / 'out/prelim/chemo/sc'
G = ROOT / 'data/seq/GSE289491'

# クラスタ → 細胞種（sc_cluster.py の Leiden 0.8 とマーカー平均 cluster_marker_means.csv から付与）
ANNOT = {
    '0': 'CAR (LepR+)', '8': 'CAR (LepR+)',
    '10': 'Osteoblast',
    '4': 'Osteoprogenitor',
    '7': 'Hypertrophic chondro/osteo',
    '6': 'Chondrocyte', '11': 'Chondrocyte',
    '1': 'Fibroblast (Prg4+)', '2': 'Fibroblast (Prg4+)', '3': 'Fibroblast (Prg4+)', '12': 'Fibroblast (Prg4+)',
    '5': 'Fibroblast (Ly6a+)',
    '13': 'Tendon/ligament fibroblast',
    '9': 'Endothelial',
}
ORDER = ['CAR (LepR+)', 'Osteoblast', 'Osteoprogenitor', 'Hypertrophic chondro/osteo', 'Chondrocyte',
         'Fibroblast (Prg4+)', 'Fibroblast (Ly6a+)', 'Tendon/ligament fibroblast', 'Endothelial']

# 事前定義パネル（タスク指定 + prelim/genes.py + 元論文の SASP/老化マーカー + p53 標的）
TASK = {
    'apoptosis': ['Fas', 'Casp3', 'Bax', 'Bbc3', 'Pmaip1', 'Cdkn1a'],
    'senescence': ['Cdkn2a', 'Serpine1', 'Il6'],
    'inflam/OC': ['Tnfsf11', 'Il1b', 'Ccl2', 'Cxcl12'],
}
PAPER = ['Cdkn2a', 'Cdkn1a', 'Mki67', 'Il6', 'Spp1', 'Cstb', 'Serpine1', 'Mmp3']  # GSE289491 summary 記載
P53 = ['Cdkn1a', 'Mdm2', 'Ccng1', 'Eda2r', 'Phlda3', 'Aen', 'Zmat3', 'Bax', 'Ei24', 'Trp53inp1',
       'Sesn2', 'Gdf15', 'Plk2', 'Polk', 'Bbc3', 'Pmaip1', 'Fas', 'Tnfrsf10b', 'Apaf1']
SASP = ['Il6', 'Cxcl1', 'Ccl2', 'Ccl7', 'Il1a', 'Il1b', 'Serpine1', 'Mmp3', 'Mmp13', 'Timp1', 'Igfbp3',
        'Igfbp7', 'Spp1', 'Cstb', 'Gdf15', 'Tnfrsf11b', 'Csf1', 'Vegfa']
OTHER = ['Tnfrsf11b', 'Mki67', 'Top2a', 'Spp1', 'Cstb', 'Mmp3', 'Mdm2', 'Eda2r', 'Phlda3', 'Ccng1', 'Gdf15',
         'Trp53inp1', 'Sost', 'Bglap', 'Runx2', 'Sp7']
PANEL = list(dict.fromkeys(sum(TASK.values(), []) + ALL_GENES + PAPER + OTHER))
NBOOT = 500
MINCELLS = 20


def cpm(v):
    return v / v.sum() * 1e6


def main():
    common.setup()
    a = sc.read_h5ad(OUT / 'adata_unlabeled.h5ad')
    a.obs['celltype'] = a.obs.leiden.astype(str).map(ANNOT)
    assert a.obs.celltype.notna().all()
    a.obs['celltype'] = pd.Categorical(a.obs.celltype, ORDER)
    ct = pd.crosstab(a.obs.celltype, a.obs['sample'])[['Veh', 'Doxo']]
    ct['Veh_%'] = (ct.Veh / ct.Veh.sum() * 100).round(1)
    ct['Doxo_%'] = (ct.Doxo / ct.Doxo.sum() * 100).round(1)
    ct.to_csv(OUT / 'celltype_counts.csv')
    print(ct.to_string())

    # UMAP
    fig, axes = plt.subplots(1, 2, figsize=(common.mm(170), common.mm(70)))
    sc.pl.umap(a, color='celltype', ax=axes[0], show=False, frameon=False, legend_fontsize=6, s=4,
               title='GSE289491 cell types')
    sc.pl.umap(a, color='sample', ax=axes[1], show=False, frameon=False, s=4, title='Sample (n=1 each)',
               palette={'Veh': '#7F7F7F', 'Doxo': '#C44E52'})
    fig.tight_layout()
    fig.savefig(OUT / 'umap.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    mk = ['Lepr', 'Cxcl12', 'Adipoq', 'Bglap', 'Ifitm5', 'Alpl', 'Ibsp', 'Col10a1', 'Mmp13', 'Dmp1', 'Sost',
          'Col2a1', 'Acan', 'Prg4', 'Cilp', 'Ly6a', 'Tnmd', 'Cdh5', 'Ptprc', 'Acp5', 'Ctsk']
    sc.pl.dotplot(a, mk, groupby='celltype', show=False, save=None)
    plt.savefig(OUT / 'dotplot_celltype.png', dpi=300, bbox_inches='tight')
    plt.close('all')

    X = a.layers['counts'].tocsr()
    genes = np.array(a.var_names)
    gi = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(0)
    rows, ma = [], {}
    for cty in ORDER:
        idx = {s: np.where((a.obs.celltype == cty).values & (a.obs['sample'] == s).values)[0]
               for s in ['Veh', 'Doxo']}
        if min(len(v) for v in idx.values()) < MINCELLS:
            print('skip (cells < %d)' % MINCELLS, cty, {k: len(v) for k, v in idx.items()})
            continue
        sums = {s: np.asarray(X[v].sum(axis=0)).ravel() for s, v in idx.items()}
        c = {s: cpm(v) for s, v in sums.items()}
        lfc_all = np.log2((c['Doxo'] + 1) / (c['Veh'] + 1))
        expressed = (c['Doxo'] + c['Veh']) / 2 >= 10  # 分布の母集団（平均 CPM >= 10）
        ref = np.sort(lfc_all[expressed])
        ma[cty] = pd.DataFrame({'gene': genes, 'meanCPM': (c['Doxo'] + c['Veh']) / 2, 'log2FC': lfc_all})
        frac = {s: np.asarray((X[v] > 0).mean(axis=0)).ravel() for s, v in idx.items()}
        pi = [gi[g] for g in PANEL if g in gi]
        boots = []
        Xp = {s: X[v][:, pi].toarray() for s, v in idx.items()}
        lib = {s: np.asarray(X[v].sum(axis=1)).ravel() for s, v in idx.items()}
        for _ in range(NBOOT):
            cc = {}
            for s in ['Veh', 'Doxo']:
                r = rng.integers(0, len(idx[s]), len(idx[s]))
                cc[s] = Xp[s][r].sum(axis=0) / lib[s][r].sum() * 1e6
            boots.append(np.log2((cc['Doxo'] + 1) / (cc['Veh'] + 1)))
        boots = np.array(boots)
        for j, g in enumerate([g for g in PANEL if g in gi]):
            i = gi[g]
            rows.append(dict(celltype=cty, gene=g, n_veh=len(idx['Veh']), n_dox=len(idx['Doxo']),
                             cpm_veh=c['Veh'][i], cpm_dox=c['Doxo'][i], log2FC=lfc_all[i],
                             boot_lo=np.percentile(boots[:, j], 2.5), boot_hi=np.percentile(boots[:, j], 97.5),
                             frac_veh=frac['Veh'][i], frac_dox=frac['Doxo'][i],
                             pctile=(np.searchsorted(ref, lfc_all[i]) / len(ref) * 100
                                     if expressed[i] else np.nan)))
    p = pd.DataFrame(rows)
    p.to_csv(OUT / 'panel_pseudobulk.csv', index=False)

    # MA プロット（細胞種ごと）
    hl = sum(TASK.values(), [])
    for cty, d in ma.items():
        d = d[d.meanCPM >= 1]
        fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(65)))
        ax.scatter(np.log10(d.meanCPM), d.log2FC, s=1.5, c='#CFCFCF', lw=0, rasterized=True)
        m = d.gene.isin(hl)
        ax.scatter(np.log10(d.meanCPM[m]), d.log2FC[m], s=10, facecolor='none', edgecolor='#1F4E79', lw=0.6)
        for x, y, g in zip(np.log10(d.meanCPM[m]), d.log2FC[m], d.gene[m]):
            ax.annotate(g, (x, y), xytext=(2, 2), textcoords='offset points', fontsize=5.5, fontstyle='italic')
        ax.axhline(0, lw=0.5, c='k')
        ax.set_xlabel('log$_{10}$ mean pseudobulk CPM')
        ax.set_ylabel('log$_2$ (Dox / Veh), n=1 each')
        ax.set_title(f'GSE289491 {cty}')
        fig.savefig(OUT / f"ma_{cty.split(' ')[0].replace('/', '_')}.png", dpi=300, bbox_inches='tight')
        plt.close(fig)
    pd.concat({k: v for k, v in ma.items()}, names=['celltype']).reset_index(level=0).to_csv(
        OUT / 'pseudobulk_all_genes.csv.gz', index=False)

    # スコア（細胞ごとの平均 log1p 正規化発現）
    sc_rows = []
    for name, gl in [('p53_targets', P53), ('SASP', SASP)]:
        gl = [g for g in gl if g in gi]
        sc.tl.score_genes(a, gl, score_name=name, random_state=0)
    s = a.obs.groupby(['celltype', 'sample'], observed=True)[['p53_targets', 'SASP']].mean().unstack()
    s.to_csv(OUT / 'scores.csv')
    print(s.round(3).to_string())

    # 著者（Partek）の全細胞比較表から同じパネルを抜く（細胞を反復として扱った検定 = 擬似反復）
    x = pd.read_excel(G / 'GSE289491_Processed_data_GSA.xlsx', header=0)
    x.columns = ['gene_id', 'ensembl', 'gene', 'p', 'fdr', 'ratio', 'fc', 'lsm_veh', 'lsm_dox']
    x = x[x.gene.isin(PANEL)].copy()
    x['log2_dox_over_veh'] = np.log2(x.lsm_dox.astype(float) / x.lsm_veh.astype(float))
    x.to_csv(OUT / 'partek_panel.csv', index=False)

    pd.set_option('display.width', 250)
    show = p[p.gene.isin(hl + ['Mdm2', 'Eda2r', 'Phlda3', 'Gdf15', 'Mki67', 'Tnfrsf11b', 'Spp1', 'Mmp3'])]
    print(show.round(2).to_string())


if __name__ == '__main__':
    main()

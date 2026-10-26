"""GSE296752（GSE296761 の scRNA。雌 Lepr-Cre;tdTomato、8–10 週で PBS×4 / Dox 5 mg/kg×4 i.v.、最終投与 8 週後、
骨髄 + 骨片を消化 → 系統（B220/CD3/CD4/CD8/CD11b/Gr1/Ter119/CD71）と CD45/cKit を MACS で 2 回除去した非血球細胞、
10x 3' v3.1。各群 2 匹をプールした 1 ライブラリずつ）の QC → クラスタリング。

手順（事前に決めたもの。prelim/aging_f/sc_cluster.py に合わせた）:
  QC: 検出遺伝子 >= 300、UMI >= 500、ミトコンドリア率 < 20%。scrublet（sample ごと）で doublet 判定された細胞を除く。
  正規化 1e4 → log1p、HVG 3000（sample を batch_key）、scale → PCA 40、近傍 15、UMAP、Leiden（res 1.0）。
  バッチ補正はしない（PBS / Dox が 1 ライブラリずつで、補正すると処置効果ごと消える）。
細胞種の割り当ては、ここで出すクラスタごとのマーカー平均（log1p 正規化）を見て sc_annotate.py の ASSIGN に書く。

  uv run --with scikit-image python -I prelim/gse296761/sc_cluster.py      # repo root から（scrublet の自動閾値に scikit-image が要る。pyproject には足していない）
出力: out/prelim/gse296761/sc/{adata_clustered.h5ad, cluster_markers.csv, cluster_top_genes.csv, qc.csv, umap_clusters.png}
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
G = ROOT / 'data/seq/GSE296752'
OUT = ROOT / 'out/prelim/gse296761/sc'
FILES = {'PBS': 'GSM8975917_PBS_filtered_feature_bc_matrix.h5', 'DOX': 'GSM8975918_DOX_filtered_feature_bc_matrix.h5'}

MARK = [
    # 造血の残り
    'Ptprc', 'Hbb-bs', 'Lyz2', 'Cd79a', 'S100a8', 'Kit', 'Pf4',
    # T / NK（CD8+ TEM が居るかの確認用）
    'Cd3e', 'Cd3g', 'Cd8a', 'Cd8b1', 'Cd4', 'Cd44', 'Ccl5', 'Il7r', 'Nkg7', 'Gzmk', 'Klrk1', 'Ifng',
    # LepR+ MSC（Adipo- / Osteo-primed）
    'Lepr', 'Cxcl12', 'Kitl', 'Adipoq', 'Lpl', 'Apoe', 'Vcam1', 'Ebf3', 'Grem1', 'Wif1', 'Spp1', 'Alpl', 'Sp7', 'Runx2',
    # 骨芽・骨細胞
    'Bglap', 'Bglap2', 'Col1a1', 'Ibsp', 'Ifitm5', 'Dmp1', 'Sost', 'Phex',
    # 線維芽・軟骨
    'Pdgfra', 'Ly6a', 'Dcn', 'Pi16', 'Clec3b', 'Col2a1', 'Acan', 'Col10a1', 'Tnmd', 'Prg4',
    # 内皮（動脈 / 洞様）
    'Pecam1', 'Cdh5', 'Emcn', 'Gja5', 'Bmx', 'Sox17', 'Ly6c1', 'Stab2', 'Flt4', 'Vwf',
    # 周細胞・平滑筋
    'Acta2', 'Myh11', 'Cnn1', 'Rgs5', 'Kcnj8', 'Higd1b', 'Notch3', 'Pdgfrb',
    # 神経・その他
    'Plp1', 'Mpz', 'Acp5', 'Ctsk', 'Mki67', 'Top2a', 'Car3', 'Fabp4',
]


def load():
    obj = []
    for s, f in FILES.items():
        a = sc.read_10x_h5(G / f)
        a.var_names_make_unique()
        a.obs_names = [f'{s}_{x}' for x in a.obs_names]
        a.obs['sample'] = s
        obj.append(a)
    a = ad.concat(obj, join='inner')
    a.obs['sample'] = pd.Categorical(a.obs['sample'], ['PBS', 'DOX'])
    return a


def process(a, res=1.0, key='leiden'):
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=3000, batch_key='sample', flavor='seurat')
    b = a[:, a.var.highly_variable].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=40, random_state=0)
    a.obsm['X_pca'] = b.obsm['X_pca']
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=40, random_state=0)
    sc.tl.umap(a, random_state=0)
    sc.tl.leiden(a, resolution=res, key_added=key, random_state=0, flavor='igraph', n_iterations=2)
    return a


def cluster_table(a, key):
    mk = [g for g in dict.fromkeys(MARK) if g in a.var_names]
    X = pd.DataFrame(a[:, mk].X.toarray(), columns=mk, index=a.obs_names)
    m = X.groupby(a.obs[key].values).mean()
    m.insert(0, 'n', a.obs[key].value_counts().reindex(m.index).values)
    ct = pd.crosstab(a.obs[key], a.obs['sample'])
    m.insert(1, 'n_PBS', ct['PBS'].reindex(m.index).values)
    m.insert(2, 'n_DOX', ct['DOX'].reindex(m.index).values)
    return m


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = load()
    a.var['mt'] = a.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], inplace=True, percent_top=None, log1p=False)
    q0 = a.obs.groupby('sample', observed=True).size()
    med0 = a.obs.groupby('sample', observed=True)[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']].median()
    a = a[(a.obs.n_genes_by_counts >= 300) & (a.obs.total_counts >= 500) & (a.obs.pct_counts_mt < 20)].copy()
    q1 = a.obs.groupby('sample', observed=True).size()
    sc.pp.scrublet(a, batch_key='sample', random_state=0)
    dbl = a.obs.groupby('sample', observed=True).predicted_doublet.sum()
    a = a[~a.obs.predicted_doublet.astype(bool)].copy()
    sc.pp.filter_genes(a, min_cells=10)
    med = a.obs.groupby('sample', observed=True)[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']].median()
    q = pd.DataFrame({'cells_raw': q0, 'cells_qc': q1, 'doublets_removed': dbl,
                      'cells_final': a.obs.groupby('sample', observed=True).size(),
                      'median_UMI_raw': med0.total_counts, 'median_UMI_final': med.total_counts,
                      'median_genes_final': med.n_genes_by_counts, 'median_mt_final': med.pct_counts_mt})
    q.to_csv(OUT / 'qc.csv')
    print(q.to_string())
    a.layers['counts'] = a.X.copy()
    a = process(a)
    m = cluster_table(a, 'leiden')
    m.round(3).to_csv(OUT / 'cluster_markers.csv')
    pd.set_option('display.width', 400)
    pd.set_option('display.max_columns', 200)
    print(m.round(2).T.to_string())
    sc.tl.rank_genes_groups(a, 'leiden', method='wilcoxon', pts=True)
    rg = sc.get.rank_genes_groups_df(a, None)
    rg = rg[(rg.logfoldchanges > 1) & (rg.pct_nz_group > 0.3)]
    rg.groupby('group').head(30).to_csv(OUT / 'cluster_top_genes.csv', index=False)
    top = rg.groupby('group').head(15).groupby('group').names.apply(lambda s: ', '.join(s))
    for k, v in top.items():
        print(k, v)
    sc.pl.umap(a, color=['leiden', 'sample'], show=False, ncols=2, s=3, legend_loc='on data', legend_fontsize=6)
    plt.savefig(OUT / 'umap_clusters.png', dpi=200, bbox_inches='tight')
    plt.close('all')
    del a.uns['rank_genes_groups']
    a.write_h5ad(OUT / 'adata_clustered.h5ad', compression='gzip')


if __name__ == '__main__':
    main()

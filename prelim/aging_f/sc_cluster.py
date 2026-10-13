"""GSE210584（雌 C57BL/6、若齢 2–4 か月 n=5 / 中年 12–14 か月 n=9、個体別 10x）の QC → クラスタリング → 細胞種付け。
GSE299773（2 / 12 / 24 か月、各 1 ライブラリ）も同じ手順で処理する（引数で切り替え）。

手順（事前に決めたもの）:
  QC: 検出遺伝子 >= 300、UMI >= 500、ミトコンドリア率 < 20%（間質は mt が高めなので 10% にはしない）。
  正規化 1e4 → log1p、HVG 3000（sample を batch_key）、PCA 40、近傍 15、UMAP、Leiden（res 1.0）。
  バッチ補正はしない（個体 = 群の入れ子で、補正すると年齢効果も削る。sc_cluster.py（chemo）と同じ方針）。
  全細胞で大分類 → 非造血（Ptprc 低）を抜き出して同じ手順で再クラスタ → マーカー平均で細胞種を付ける。
細胞種の割り当ては、ここで出したクラスタごとのマーカー平均（log1p 正規化）を見て sc_pseudobulk.py の ASSIGN に書いた。

  uv run python -I prelim/aging_f/sc_cluster.py GSE210584
出力: out/prelim/aging_f/sc/<GSE>_{all,stroma}.h5ad, <GSE>_cluster_markers_{all,stroma}.csv, <GSE>_umap_*.png,
      <GSE>_celltype_counts.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'out/prelim/aging_f/sc'

MARK = ['Ptprc', 'Hbb-bs', 'Lyz2', 'Cd79a', 'Cd3e', 'S100a8', 'Kit', 'Pf4', 'Mki67',
        'Lepr', 'Cxcl12', 'Kitl', 'Adipoq', 'Vcam1', 'Ebf3', 'Pdgfra', 'Ly6a', 'Dcn', 'Pi16',
        'Runx2', 'Sp7', 'Bglap', 'Col1a1', 'Ibsp', 'Alpl', 'Spp1', 'Dmp1', 'Sost', 'Igf1',
        'Col2a1', 'Acan', 'Cdh5', 'Pecam1', 'Emcn', 'Ly6a', 'Acta2', 'Myh11', 'Rgs5', 'Kcnj8',
        'Plp1', 'Mpz', 'Acp5', 'Ctsk', 'Tnmd', 'Prg4', 'Ifitm5']

# 非造血の大分類: クラスタの Ptprc 平均 < 0.5 かつ Hbb-bs 平均 < 2（赤芽球除外）。実際の値は *_cluster_markers_all.csv


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
    return m


def main(gse):
    a = sc.read_h5ad(OUT / f'{gse}_raw.h5ad')
    a.var['mt'] = a.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], inplace=True, percent_top=None, log1p=False)
    q0 = a.obs.groupby('sample', observed=True).size()
    a = a[(a.obs.n_genes_by_counts >= 300) & (a.obs.total_counts >= 500) & (a.obs.pct_counts_mt < 20)].copy()
    sc.pp.filter_genes(a, min_cells=10)
    q = pd.DataFrame({'cells_raw': q0, 'cells_qc': a.obs.groupby('sample', observed=True).size()})
    print(q.to_string())
    a.layers['counts'] = a.X.copy()
    a = process(a)
    m = cluster_table(a, 'leiden')
    m.round(2).to_csv(OUT / f'{gse}_cluster_markers_all.csv')
    pd.set_option('display.width', 300)
    print(m[['n', 'Ptprc', 'Hbb-bs', 'Lepr', 'Cxcl12', 'Cdh5', 'Bglap', 'Col1a1', 'Acta2', 'Col2a1']].round(2).to_string())
    stro = m.index[(m.Ptprc < 0.5) & (m['Hbb-bs'] < 2)]
    print('non-hematopoietic clusters:', list(stro))
    a.obs['compartment'] = np.where(a.obs.leiden.isin(stro), 'non-hem', 'hem')
    sc.pl.umap(a, color=['leiden', 'sample', 'compartment'], show=False, ncols=3, s=2)
    plt.savefig(OUT / f'{gse}_umap_all.png', dpi=200, bbox_inches='tight')
    plt.close('all')
    a.write_h5ad(OUT / f'{gse}_all.h5ad', compression='gzip')

    s = a[a.obs.compartment == 'non-hem'].copy()
    s.X = s.layers['counts'].copy()
    s = process(s, res=1.0, key='sub')
    m2 = cluster_table(s, 'sub')
    ct = pd.crosstab(s.obs['sub'], s.obs['sample'])
    m2.round(2).to_csv(OUT / f'{gse}_cluster_markers_stroma.csv')
    ct.to_csv(OUT / f'{gse}_stroma_cluster_by_sample.csv')
    print(m2.round(2).to_string())
    print(ct.to_string())
    sc.pl.umap(s, color=['sub', 'sample', 'group'], show=False, ncols=2, s=4,
               legend_loc='on data', legend_fontsize=6)
    plt.savefig(OUT / f'{gse}_umap_stroma.png', dpi=200, bbox_inches='tight')
    plt.close('all')
    s.write_h5ad(OUT / f'{gse}_stroma.h5ad', compression='gzip')


if __name__ == '__main__':
    main(sys.argv[1])

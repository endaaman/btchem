"""GSE289491（雌 12 週, Dox 10 mg/kg 96 h, 骨常在細胞 CD45-/CD71-/Ter119-/CD31- の scRNA, 各 n=1）。

sc_counts.py で戻した UMI counts を統合 → QC → 正規化 → HVG/PCA → Leiden → マーカーで細胞種を付ける。
バッチ補正はしない（Veh / Dox が 1 ライブラリずつで、補正すると処置効果を消すため）。

  uv run python -I prelim/chemo/sc_cluster.py
出力: out/prelim/chemo/sc/{adata.h5ad, umap.png, markers_dotplot.png, cluster_markers.csv, cluster_counts.csv}
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
G = ROOT / 'data/seq/GSE289491'
OUT = ROOT / 'out/prelim/chemo/sc'

# 細胞種マーカー（事前定義）
MARKERS = {
    'CAR/LepR+': ['Lepr', 'Cxcl12', 'Kitl', 'Adipoq', 'Vcam1', 'Ebf3'],
    'Osteo-CAR': ['Spp1', 'Ibsp', 'Sp7', 'Alpl'],
    'Osteoblast': ['Bglap', 'Bglap2', 'Col1a1', 'Sp7', 'Ifitm5', 'Dmp1'],
    'Osteocyte': ['Sost', 'Dmp1', 'Phex', 'Mepe', 'Fgf23'],
    'Chondrocyte': ['Col2a1', 'Acan', 'Col10a1', 'Sox9', 'Comp'],
    'Fibroblast': ['Dcn', 'Pdgfra', 'Ly6a', 'Pi16', 'Prg4', 'Cilp', 'Tnmd'],
    'Pericyte/SMC': ['Acta2', 'Myh11', 'Rgs5', 'Kcnj8'],
    'Endothelial': ['Pecam1', 'Cdh5', 'Emcn'],
    'Hematopoietic': ['Ptprc', 'Lyz2', 'Cd74', 'Hbb-bs', 'Cd79a', 'Ms4a1'],
    'Osteoclast': ['Acp5', 'Ctsk', 'Mmp9', 'Oscar', 'Dcstamp', 'Calcr'],
    'Myeloid/Neutrophil': ['S100a8', 'S100a9', 'Ngp', 'Camp'],
    'Muscle': ['Acta1', 'Ckm', 'Myh4'],
    'Schwann/neural': ['Plp1', 'Mpz', 'S100b'],
    'Cycling': ['Mki67', 'Top2a'],
}


def load():
    genes = pd.read_csv(G / 'genes.txt', header=None)[0].astype(str).values
    obj = []
    for lab in ['Veh', 'Doxo']:
        m = sp.load_npz(G / f'counts_{lab}.npz').T.tocsr()  # cells x genes
        m.data = np.round(m.data * 2) / 2  # 0.5 刻み（重複シンボルの平均）は残す
        cells = pd.read_csv(G / f'cells_{lab}.txt', header=None)[0].astype(str)
        a = ad.AnnData(X=m.astype(np.float32), obs=pd.DataFrame(index=lab + '_' + cells),
                       var=pd.DataFrame(index=genes))
        a.obs['sample'] = lab
        obj.append(a)
    a = ad.concat(obj)
    a.var_names_make_unique()
    return a


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(OUT)
    a = load()
    a.var['mt'] = a.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], inplace=True, percent_top=None)
    print(a.obs.groupby('sample')[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']].median())
    a.layers['counts'] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    sc.tl.pca(a, n_comps=30, mask_var='highly_variable')
    sc.pp.neighbors(a, n_pcs=30)
    sc.tl.umap(a, random_state=0)
    sc.tl.leiden(a, resolution=0.8, random_state=0, flavor='igraph', n_iterations=2)
    sc.tl.rank_genes_groups(a, 'leiden', method='wilcoxon')
    rg = sc.get.rank_genes_groups_df(a, None)
    rg.groupby('group').head(25).to_csv(OUT / 'cluster_markers.csv', index=False)
    top = rg.groupby('group').head(12).groupby('group').names.apply(lambda s: ', '.join(s))
    for k, v in top.items():
        print(k, v)
    # クラスタごとのマーカー平均（log1p）
    mk = [g for v in MARKERS.values() for g in v if g in a.var_names]
    mk = list(dict.fromkeys(mk))
    df = pd.DataFrame(a[:, mk].X.toarray(), columns=mk, index=a.obs_names)
    df['leiden'] = a.obs.leiden.values
    mean = df.groupby('leiden').mean().T
    mean.to_csv(OUT / 'cluster_marker_means.csv')
    pd.set_option('display.width', 250)
    print(mean.round(2).to_string())
    ct = pd.crosstab(a.obs.leiden, a.obs['sample'])
    print(ct.to_string())
    a.obs.index.name = None
    a.var.index.name = None
    a.write_h5ad(OUT / 'adata_unlabeled.h5ad')


if __name__ == '__main__':
    main()

"""GSE255653（ソートした LepR+ MSC の bulk）の感度分析用: 「CAR 細胞より血球・内皮で高い遺伝子」の参照表。

GSE255653 は老齢 3 例すべてで Ptprc / Pecam1 / Cd74 が若齢より高く、ソート純度の差（血球・内皮の混入）が群と完全に
重なる（共変量では外せない）。そこで、雌マウス骨髄の単一細胞 GSE210584（既存の h5ad、prelim/aging_f/ が作ったもの）から
細胞種ごとの擬似バルク CPM を作り、混入で動きうる遺伝子を事前に定義して、GSEA から抜いた版も回せるようにする。

  定義: max(血球 CPM, 内皮 CPM) >= 10 かつ  (max + 1) / (CAR CPM + 1) >= 4
  血球 = 全細胞 h5ad の compartment == 'hem'、CAR / 内皮 = stroma h5ad の sub を sc_pseudobulk.py の ASSIGN で割り当てたもの。
  10x 3' と SMART-seq v4 の違いはあるが、細胞種間の比で使うだけなので許容とする。

  uv run python -I prelim/aging2/contam_ref.py
出力: out/prelim/aging2/ref/celltype_cpm_GSE210584.csv.gz, contam_genes.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

ROOT = Path(__file__).resolve().parents[2]
SC = ROOT / 'out/prelim/aging_f/sc'
OUT = ROOT / 'out/prelim/aging2/ref'
sys.path.insert(0, str(ROOT / 'prelim/aging_f'))
ASSIGN = {**{k: 'CAR' for k in ['2', '3', '4', '5', '6', '7', '8', '9', '10', '15']},
          '13': 'Osteo', '14': 'Osteo', '0': 'Endothelial', '1': 'Endothelial', '18': 'Endothelial'}
FOLD, MINCPM = 4, 10


def colsum(a, mask, chunk=5000):
    idx = np.where(mask)[0]
    tot = np.zeros(a.n_vars)
    X = a.layers['counts']
    for i in range(0, len(idx), chunk):
        x = X[idx[i:i + chunk]]
        tot += np.asarray(x.sum(axis=0)).ravel()
    return tot


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    allc = ad.read_h5ad(SC / 'GSE210584_all.h5ad', backed='r')
    hem = colsum(allc, (allc.obs.compartment == 'hem').values)
    genes = allc.var_names
    st = ad.read_h5ad(SC / 'GSE210584_stroma.h5ad', backed='r')
    assert (st.var_names == genes).all()
    ct = st.obs['sub'].astype(str).map(ASSIGN)
    res = {'hem': hem}
    for k in ['CAR', 'Endothelial', 'Osteo']:
        res[k] = colsum(st, (ct == k).values)
    d = pd.DataFrame(res, index=genes)
    cpm = d / d.sum() * 1e6
    cpm.round(3).to_csv(OUT / 'celltype_cpm_GSE210584.csv.gz')
    mx = cpm[['hem', 'Endothelial']].max(axis=1)
    ratio = (mx + 1) / (cpm.CAR + 1)
    flag = (mx >= MINCPM) & (ratio >= FOLD)
    out = pd.DataFrame({'CAR_cpm': cpm.CAR, 'hem_cpm': cpm.hem, 'endo_cpm': cpm.Endothelial,
                        'ratio_vs_CAR': ratio, 'source': np.where(cpm.hem >= cpm.Endothelial, 'hem', 'endo')})[flag]
    out.sort_values('ratio_vs_CAR', ascending=False).round(2).to_csv(OUT / 'contam_genes.csv')
    print('contamination-prone genes:', int(flag.sum()), out.source.value_counts().to_dict())
    print(out.loc[[g for g in ['Ptprc', 'Pecam1', 'Cd74', 'Ly6a', 'Il1b', 'Tnf', 'Il6', 'S100a8', 'Lyz2', 'Cdh5']
                   if g in out.index]].round(1).to_string())


if __name__ == '__main__':
    main()

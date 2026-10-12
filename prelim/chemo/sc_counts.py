"""GSE289491 の filtered matrix（Partek 出力）を UMI counts に逆変換する。

GEO の *filtered_feature_bc_matrix.txt は生 counts ではなく、各細胞で
  u = log(k * a + 1)（a = 1e6 / UMI 総数、k = UMI 数。つまり log(CPM+1)）
をもう一度 CPM 化して log2(x+1) したもの（v = log2(CPM(u) + 1)）になっている
（値の段階と頻度から逆算。docs/prelim_chemo.md 参照）。L = 2^v - 1 は u に比例するので、
最頻の段 L1（k=1）と次に多い段 L2（k=2）から a を解き、k = ((a+1)^(L/L1) - 1)/a で戻す。
同一シンボルに複数 Ensembl がある遺伝子は Partek が平均しているため k が 0.5 刻みになる → そのまま残す。

  uv run python -I prelim/chemo/sc_counts.py
出力: data/seq/GSE289491/counts_<Veh|Doxo>.npz（疎行列 genes x cells）, genes.txt, qc_inverse.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
G = ROOT / 'data/seq/GSE289491'
FILES = {'Veh': 'GSM8792668_Veh-filtered_feature_bc_matrix.txt.gz',
         'Doxo': 'GSM8792669_Doxo-filtered_feature_bc_matrix.txt.gz'}


def invert(v):
    L = np.power(2.0, v) - 1
    lv = pd.Series(np.round(L, 3)).value_counts()
    L1 = lv.index[0]
    above = lv[lv.index > L1 * 1.001]
    L2 = above.index[0]  # 次に多い段（value_counts は頻度降順）
    r = L2 / L1
    a = brentq(lambda a: np.log(2 * a + 1) / np.log(a + 1) - r, 1e-4, 1e10)
    k = (np.power(a + 1, L / L1) - 1) / a
    return k, a


def main():
    qc = []
    genes = None
    for lab, fn in FILES.items():
        d = pd.read_csv(G / 'raw' / fn, sep='\t', index_col=0)
        d = d.loc[d.index.notna()]
        if genes is None:
            genes = list(d.index)
        assert list(d.index) == genes
        x = d.values
        cols = []
        for j in range(x.shape[1]):
            v = x[:, j]
            nz = np.nonzero(v)[0]
            k, a = invert(v[nz])
            res = np.abs(k * 2 - np.round(k * 2))  # 0.5 刻みまで許容
            qc.append(dict(sample=lab, cell=d.columns[j], umi_from_a=1e6 / a, umi_sum=k.sum(),
                           n_genes=len(nz), frac_nonhalf=float((res > 0.05).mean()), max_k=k.max()))
            col = np.zeros(len(v))
            col[nz] = k
            cols.append(sp.csc_matrix(col[:, None]))
        m = sp.hstack(cols).tocsc()
        sp.save_npz(G / f'counts_{lab}.npz', m)
        pd.Series(d.columns).to_csv(G / f'cells_{lab}.txt', index=False, header=False)
        print(lab, m.shape)
    pd.Series(genes).to_csv(G / 'genes.txt', index=False, header=False)
    q = pd.DataFrame(qc)
    q.to_csv(G / 'qc_inverse.csv', index=False)
    q['ok'] = (q.frac_nonhalf < 0.01) & ((q.umi_sum / q.umi_from_a - 1).abs() < 0.01)
    print(q.groupby('sample').agg(n=('cell', 'size'), ok=('ok', 'mean'),
                                  umi_med=('umi_sum', 'median'), genes_med=('n_genes', 'median')))


if __name__ == '__main__':
    main()

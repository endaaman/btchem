"""GSE210584（若齢 5 / 中年 9、個体別 10x）と GSE299773（2 / 12 / 24 か月、各 1）の CellRanger filtered matrix を
AnnData にまとめ、性別マーカー（Xist と Y 遺伝子）をサンプルごとに集計する。

  uv run python -I prelim/aging_f/sc_load.py      # repo root から
出力: out/prelim/aging_f/sc/<GSE>_raw.h5ad, sex_markers_sc.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad

ROOT = Path(__file__).resolve().parents[2]
SEQ = ROOT / 'data/seq'
OUT = ROOT / 'out/prelim/aging_f/sc'
SEX = ['Xist', 'Tsix', 'Ddx3y', 'Eif2s3y', 'Kdm5d', 'Uty']

SAMPLES = {
    'GSE210584': [(f'GSM64328{57 + i}_MA{i + 1}', f'MA{i + 1}', 'middle') for i in range(9)]
    + [(f'GSM64328{66 + i}_Young{i + 1}', f'Y{i + 1}', 'young') for i in range(5)],
    'GSE299773': [('GSM9045893_2M', '2M', '2m'), ('GSM9045894_12M', '12M', '12m'), ('GSM9045895_24M', '24M', '24m')],
}


def read10x(prefix):
    m = scipy.io.mmread(f'{prefix}_matrix.mtx.gz').T.tocsr().astype(np.float32)  # cells x genes
    f = pd.read_csv(f'{prefix}_features.tsv.gz', sep='\t', header=None)
    b = pd.read_csv(f'{prefix}_barcodes.tsv.gz', header=None)[0]
    a = ad.AnnData(X=m, obs=pd.DataFrame(index=b.values), var=pd.DataFrame({'gene_id': f[0].values}, index=f[1].values))
    a.var_names_make_unique()
    return a


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for gse, ss in SAMPLES.items():
        objs = []
        for pre, name, grp in ss:
            a = read10x(SEQ / gse / 'raw' / pre)
            a.obs_names = name + '_' + a.obs_names
            a.obs['sample'] = name
            a.obs['group'] = grp
            objs.append(a)
            lib = np.asarray(a.X.sum(axis=1)).ravel()
            tot = lib.sum()
            r = dict(gse=gse, sample=name, group=grp, cells=a.n_obs, umi_median=np.median(lib))
            for g in SEX:
                if g in a.var_names:
                    r[g + '_CPM'] = a[:, g].X.sum() / tot * 1e6
            rows.append(r)
            print(gse, name, a.n_obs)
        a = ad.concat(objs, merge='first')
        a.obs['sample'] = a.obs['sample'].astype('category')
        a.write_h5ad(OUT / f'{gse}_raw.h5ad', compression='gzip')
    s = pd.DataFrame(rows)
    s.to_csv(OUT / 'sex_markers_sc.csv', index=False)
    print(s.round(1).to_string())


if __name__ == '__main__':
    main()

"""GSE289491 の filtered matrix の値の形式を調べる（使い捨て確認）。"""
import sys
import numpy as np
import pandas as pd
p = sys.argv[1]
d = pd.read_csv(p, sep='\t', index_col=0, nrows=None)
print(d.shape)
x = d.values
for j in range(5):
    v = x[:, j]; nz = v[v > 0]
    mn = nz.min()
    lib = 1e6 / (2 ** mn - 1)
    c = (2 ** nz - 1) * lib / 1e6
    print(j, 'min', mn, 'max', nz.max(), 'nnz', len(nz), 'lib_est', round(lib), 'resid', np.abs(c - c.round()).max(), 'sum_c', c.sum().round())
v = x[:, 0]; nz = np.sort(v[v > 0])
u = np.unique(np.round(nz, 4))
print('unique vals', len(u), u[:15], u[-5:])
for base in [2, np.e, 10]:
    cp = (base ** u - 1)
    print(base, (cp / cp[0])[:15].round(3))

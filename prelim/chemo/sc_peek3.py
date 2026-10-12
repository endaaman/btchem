"""GSE289491: 二重正規化（log(CPM+1) → CPM → log2(+1)）の逆変換で整数 counts に戻せるか確認（使い捨て）。"""
import sys
import numpy as np
import pandas as pd
from scipy.optimize import brentq
d = pd.read_csv(sys.argv[1], sep='\t', index_col=0, usecols=list(range(0, 41)))
for j in range(40):
    v = d.iloc[:, j]; v = v[v > 0]
    L = 2 ** v - 1
    lv = L.round(3).value_counts()
    L1 = lv.index[0]                      # 最頻 = k=1
    L2 = L[L > L1 * 1.0001].min()         # 次の段 = k=2
    r = L2 / L1
    a = brentq(lambda a: np.log(2 * a + 1) / np.log(a + 1) - r, 1e-3, 1e9)
    k = ((a + 1) ** (L / L1) - 1) / a
    res = np.abs(k - k.round())
    print(j, d.columns[j], 'lib', round(1e6 / a), 'sumk', round(k.sum()), 'maxres', res.max().round(3),
          'n_bad', int((res > 0.05).sum()), 'below1', [(g, round(x, 3)) for g, x in k[k < 0.95].items()][:3],
          'top', k.idxmax(), round(k.max(), 2))

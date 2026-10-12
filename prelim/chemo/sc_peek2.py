"""GSE289491 の値の段階を調べる（使い捨て確認）。"""
import sys
import numpy as np
import pandas as pd
d = pd.read_csv(sys.argv[1], sep='\t', index_col=0, usecols=[0, 1, 2])
for j in range(2):
    v = d.iloc[:, j]; v = v[v > 0].round(4)
    vc = v.value_counts().sort_index()
    print(vc.to_string())
    print(v.sort_values().tail(8))

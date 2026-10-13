"""prelim/de.py と同じ低発現フィルタ・pydeseq2 設定で DE を回す共通関数（de.py は触らずに中身を写した）。

フィルタ: CPM >= 10 / (中央値ライブラリ[M]) を満たすサンプルが小さい方の群の n 以上、かつ counts 合計 >= 15。
pydeseq2 MLE、shrink なし、refit_cooks=True。連続変数の contrast（年齢の傾き）にも対応する。
"""
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

COLS = ['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']


def expr_filter(x, n_min):
    lib = x.sum(axis=1)
    cpm_cut = 10 / (lib.median() / 1e6)
    cpm = x.div(lib, axis=0) * 1e6
    return ((cpm >= cpm_cut).sum(axis=0) >= n_min) & (x.sum(axis=0) >= 15)


def deseq(x, meta, design, contrast, n_min):
    """x: samples x genes の整数 counts。contrast: ['cond','test','ref'] または連続変数名（文字列）。"""
    keep = expr_filter(x, n_min)
    x = x.loc[:, keep].round().astype(int)
    inf = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(counts=x, metadata=meta, design=design, refit_cooks=True, inference=inf, quiet=True)
    dds.deseq2()
    if isinstance(contrast, str):
        # 連続変数: 係数をそのまま検定（LFC は 1 単位あたり）
        cols = list(dds.varm['LFC'].columns)
        c = np.array([1.0 if k == contrast else 0.0 for k in cols])
        st = DeseqStats(dds, contrast=c, inference=inf, quiet=True)
    else:
        st = DeseqStats(dds, contrast=contrast, inference=inf, quiet=True)
    st.summary()
    r = st.results_df.reset_index()
    r = r.rename(columns={r.columns[0]: 'gene'})[COLS]
    return r, int(keep.sum())


def ndeg(r, lfc=1.0, p=0.05):
    up = ((r.padj < p) & (r.log2FoldChange >= lfc)).sum()
    dn = ((r.padj < p) & (r.log2FoldChange <= -lfc)).sum()
    return int(up), int(dn)

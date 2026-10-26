"""prerank GSEA（Hallmark+SenMayo / GO BP / Reactome）と、要約モジュール図と同じ形式の列。

  uv run python -I prelim/fig4a_v3/gsea4a.py [label ...]

順位: Wald 統計量（log2FC / lfcSE）。gseapy prerank, 1000 perm, seed 0, size 15–500。
FDR はコレクションごと（Hallmark 51 / GO BP / Reactome）に計算。既存の out/prelim/human/enrich/ は Enrichr 版
Hallmark 2020 + SenMayo なので値が少し違う（例: H_Dmab の TNF-α は −1.82 → 本解析の値を参照）。
要約モジュール形式（summary.py が使う）は prelim/invitro/gsea_iv.py の collection() と同じ set 集合・設定
（1000 perm, seed 0, min 4, max 1000）で別に回す（out/prelim/fig4a_v3/gsea/summaryset_<label>.csv）。

出力: out/prelim/fig4a_v3/gsea/<collection>_<label>.csv, summaryset_<label>.csv
"""
import sys
from pathlib import Path

import pandas as pd
import gseapy as gp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common4a import OUT, HDE, collections, wald_rank  # noqa: E402

D = OUT / 'gsea'
# label -> DE 表
SRC = {
    'H_Dmab': HDE / 'H_Dmab.csv',
    'H_Dmab_common': HDE / 'H_Dmab_common.csv',
    'H_Dmab_CB': HDE / 'H_Dmab_CB.csv',
    'H_PostMeno': HDE / 'H_PostMeno.csv',
    'H_Old': HDE / 'H_Old.csv',
    'H_OldE': HDE / 'H_OldE.csv',
    'H_Dmab_marrowadj': OUT / 'de/H_Dmab_marrowadj.csv',
    'H_Dmab_CB_marrowadj': OUT / 'de/H_Dmab_CB_marrowadj.csv',
    'H_Dmab_ocadj': OUT / 'de/H_Dmab_ocadj.csv',
    'H_Dmab_CB_ocadj': OUT / 'de/H_Dmab_CB_ocadj.csv',
    'H_PostMeno_procadj': OUT / 'de/H_PostMeno_procadj.csv',
}
SUMMARY_LABELS = ['H_Dmab', 'H_Dmab_CB', 'H_Dmab_marrowadj', 'H_Dmab_CB_marrowadj', 'H_PostMeno_procadj']


def prerank(rnk, sets, lo, hi):
    r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=lo, max_size=hi,
                   seed=0, threads=16, outdir=None, no_plot=True, verbose=False).res2d
    for k in ['ES', 'NES', 'NOM p-val', 'FDR q-val', 'FWER p-val']:
        r[k] = r[k].astype(float)
    return r


def main(labels):
    D.mkdir(parents=True, exist_ok=True)
    cols = collections()
    import gsea_iv  # noqa: E402  prelim/invitro（collection() を借りるだけ）
    sumsets = None
    for lb in labels:
        de = pd.read_csv(SRC[lb])
        rnk = wald_rank(de)
        for cn, sets in cols.items():
            f = D / f'{cn}_{lb}.csv'
            if f.exists():
                continue
            r = prerank(rnk, sets, 15, 500).assign(label=lb, collection=cn, n_ranked=len(rnk))
            r.to_csv(f, index=False)
            print(lb, cn, len(r), 'sets;', int((r['FDR q-val'] < 0.05).sum()), 'FDR<0.05')
        if lb in SUMMARY_LABELS:
            f = D / f'summaryset_{lb}.csv'
            if not f.exists():
                if sumsets is None:
                    sumsets, _ = gsea_iv.collection()
                r = prerank(rnk, sumsets, 4, 1000).assign(label=lb, n_ranked=len(rnk))
                r.to_csv(f, index=False)
                print(lb, 'summary set', len(r))


if __name__ == '__main__':
    main(sys.argv[1:] or list(SRC))

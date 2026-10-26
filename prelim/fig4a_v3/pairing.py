"""GSE141595（消化後）と GSE141610（遠心のみ）は同じ 30 人か？ — 遺伝型で決まる発現から照合する。

  uv run python -I prelim/fig4a_v3/pairing.py

論文（Weivoda 2020 Nat Commun, PMID 31911667）では、1 人から 4 本の針生検を取り、
#1 を遠心のみ（GSE141610）、#2–4 をまとめて liberase 消化（GSE141595）に回している。
サンプル番号（SB / CB）は対応していない。同じ人なら、遺伝型で発現が決まる遺伝子（欠失多型・強い cis-eQTL・HLA の
アリル特異的なマッピング）の個人差が 2 つのデータで一致するはず。

方法:
 (a) 既知の遺伝型依存遺伝子（RHD 欠失、ERAP2 null アリル、SIGLEC14 欠失、HLA-DRB5 ハプロタイプ、HLA-DQ、CHIT1 null など）
     の log2 CPM を各データ内で z 化し、SB × CB の相関行列から最適 1 対 1 対応（Hungarian）と相互ベスト対を出す。
 (b) 感度: 各データで上位 5 主成分（組成）を除いた残差のうち、両データで残差分散の大きい 100 遺伝子で同じことをする。
判定: 相互ベスト対の数と、その対で群（placebo / denosumab）が一致する割合（同じ人なら 100% になるはず）。
帰無: CB の列ラベルをランダムに入れ替えて 1,000 回（相互ベスト対の数・群一致の分布）。

出力: out/prelim/fig4a_v3/pairing/pairing_{curated,residual}.csv, pairing_summary.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common4a import OUT, counts  # noqa: E402

GENO = ['RHD', 'ERAP2', 'SIGLEC14', 'HLA-DRB5', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DQA2', 'HLA-DQB2', 'HLA-DRB1', 'HLA-H',
        'CHIT1', 'FAM118A', 'PSMD5-AS1', 'GABPB1-AS1', 'PILRB', 'STAG3L5P', 'PVRIG2P', 'THNSL2', 'CHI3L2', 'XRRA1',
        'PEX6', 'POMZP3', 'LRRC37A4P', 'KANSL1-AS1', 'NBPF3', 'APOBEC3B', 'FCGR3B', 'CCL3L3', 'GSTM1', 'UGT2B17',
        'ZNF890P', 'CYP4F11', 'CRYBB2', 'NQO2', 'MICA', 'MICB', 'HLA-C', 'HLA-B', 'HLA-A', 'ZNF37BP']


def match(ra, rb):
    C = np.corrcoef(ra.T.values, rb.T.values)[:ra.shape[1], ra.shape[1]:]
    r, c = linear_sum_assignment(-C)
    rowmax, colmax = C.argmax(1), C.argmax(0)
    mutual = [i for i in range(C.shape[0]) if colmax[rowmax[i]] == i]
    return C, r, c, rowmax, mutual


def run(ra, rb, ga, gb, name, rng):
    C, r, c, rowmax, mutual = match(ra, rb)
    pairs = pd.DataFrame({'SB': ra.columns[r], 'CB': rb.columns[c], 'r': C[r, c],
                          'group_SB': ga[ra.columns[r]].values, 'group_CB': gb[rb.columns[c]].values})
    pairs['mutual_best'] = [i in mutual and rowmax[i] == c[k] for k, i in enumerate(r)]
    pairs['same_group'] = pairs.group_SB == pairs.group_CB
    obs_mut = len(mutual)
    obs_agree = float(np.mean([ga.iloc[i] == gb.iloc[rowmax[i]] for i in mutual])) if mutual else np.nan
    null_mut, null_agree = [], []
    for _ in range(1000):
        perm = rng.permutation(rb.shape[1])
        # 遺伝子方向の対応を壊す（各 CB 列の中で遺伝子をシャッフル）
        rbp = pd.DataFrame(np.column_stack([rng.permutation(rb.iloc[:, j].values) for j in range(rb.shape[1])]),
                           index=rb.index, columns=rb.columns[perm])
        _, _, _, rm, mu = match(ra, rbp)
        null_mut.append(len(mu))
        gbp = gb[rbp.columns]
        null_agree.append(np.mean([ga.iloc[i] == gbp.iloc[rm[i]] for i in mu]) if mu else np.nan)
    null_mut = np.array(null_mut)
    s = dict(method=name, n_genes=ra.shape[0], mutual_best_pairs=obs_mut,
             mutual_best_same_group=obs_agree,
             hungarian_same_group=pairs.same_group.mean(),
             matched_r_median=pairs.r.median(), offdiag_r_median=float(np.median(C[~np.eye(30, dtype=bool)[:C.shape[0], :C.shape[1]]])),
             null_mutual_mean=null_mut.mean(), null_mutual_p=(1 + (null_mut >= obs_mut).sum()) / (1 + len(null_mut)),
             null_same_group_mean=np.nanmean(null_agree))
    return pairs, s


def main():
    d = OUT / 'pairing'
    d.mkdir(parents=True, exist_ok=True)
    a, sa = counts('GSE141595/post30')
    b, sb = counts('GSE141610')
    la = np.log2(a / a.sum() * 1e6 + 1)
    lb = np.log2(b / b.sum() * 1e6 + 1)
    rng = np.random.default_rng(0)
    out = []
    # (a) 遺伝型依存の既知遺伝子
    g = [x for x in GENO if x in la.index and x in lb.index and la.loc[x].mean() > 0.5 and lb.loc[x].mean() > 0.5]
    za = la.loc[g].sub(la.loc[g].mean(1), axis=0).div(la.loc[g].std(1), axis=0)
    zb = lb.loc[g].sub(lb.loc[g].mean(1), axis=0).div(lb.loc[g].std(1), axis=0)
    p, s = run(za, zb, sa.group, sb.group, 'curated genotype genes', rng)
    s['genes'] = ','.join(g)
    p.to_csv(d / 'pairing_curated.csv', index=False)
    out.append(s)
    # (b) 組成（上位 5 PC）を除いた残差
    common = la.index.intersection(lb.index)
    keep = common[(la.loc[common].mean(1) > 4) & (lb.loc[common].mean(1) > 4)]

    def resid(l, k=5):
        X = l.loc[keep].T.values
        X = X - X.mean(0)
        U, S, Vt = np.linalg.svd(X, full_matrices=False)
        return pd.DataFrame((X - (U[:, :k] * S[:k]) @ Vt[:k]).T, index=keep, columns=l.columns)
    ra, rb = resid(la), resid(lb)
    v = ra.var(1).rank(pct=True) + rb.var(1).rank(pct=True)
    top = v.sort_values(ascending=False).index[:100]
    p, s = run(ra.loc[top], rb.loc[top], sa.group, sb.group, 'residual (5 PCs removed), top-100 genes', rng)
    s['genes'] = ','.join(top[:40]) + ',...'
    p.to_csv(d / 'pairing_residual.csv', index=False)
    out.append(s)
    S = pd.DataFrame(out)
    S.to_csv(d / 'pairing_summary.csv', index=False)
    print(S.drop(columns='genes').T.to_string())


if __name__ == '__main__':
    main()

"""デノスマブ vs placebo（閉経後）の中身: DEG・GSEA 上位（上昇・低下とも）・2 つの画分での再現・組成（破骨細胞の減少）との関係。

  uv run python -I prelim/fig4a_v3/dmab_sets.py      # composition.py, dmab_de.py, gsea4a.py の後

- DEG: 既存 out/prelim/human/de/H_Dmab.csv（消化後）と H_Dmab_CB.csv（遠心のみ）。閾値は padj < 0.05（|lfc| ≥ 1 と ≥ 0 の両方を数える）。
- GSEA: gsea4a.py の Hallmark / Reactome / GO BP。FDR < 0.05 の set を上昇・低下とも全部書き出す。
- 2 画分の一致: 同じ 30 人の別の生検片（pairing.py）なので、**独立コホートでの再現ではなく、画分（処理法）を変えた再現**。
- 組成との関係: 消化後 30 例で、各遺伝子と破骨マーカースコア・骨髄スコアの群内相関（placebo / denosumab それぞれで中心化して
  プールした Pearson）を計算し、有意 set の leading edge が「破骨細胞と一緒に動く遺伝子」でどれだけ占められるかを見る。
  共変量補正版（dmab_de.py の marrowadj / ocadj）の NES も並べる。ocadj は破骨スコアが群とほぼ分離している（下に点双列相関）ので
  不安定になりうる。

出力: out/prelim/fig4a_v3/dmab/{deg_*.csv, sets_fdr05.csv, fraction_concordance.csv, le_osteoclast.csv, covariate_nes.csv,
      dmab.png|pdf}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common4a import OUT, HDE, counts, logcpm, keep_genes, scores, zscore, MARROW_GENES, cpm, save  # noqa: E402
import common  # noqa: E402

D = OUT / 'dmab'
G = OUT / 'gsea'
COLLS = ['hallmark', 'reactome', 'gobp']


def within_corr(L, grp, score):
    """群ごとに中心化してプールした Pearson r（遺伝子 x スコア）。"""
    Lc = L.copy()
    sc = score.copy()
    for g in grp.unique():
        idx = grp.index[grp == g]
        Lc[idx] = Lc[idx].sub(Lc[idx].mean(1), axis=0)
        sc[idx] = sc[idx] - sc[idx].mean()
    X = Lc.values
    y = sc[Lc.columns].values
    num = X @ y
    den = np.sqrt((X ** 2).sum(1) * (y ** 2).sum())
    return pd.Series(num / den, index=L.index)


def main():
    D.mkdir(parents=True, exist_ok=True)
    # ---- DEG ----
    deg = {}
    for lb in ['H_Dmab', 'H_Dmab_CB', 'H_Dmab_marrowadj', 'H_Dmab_CB_marrowadj']:
        f = HDE / f'{lb}.csv' if lb in ('H_Dmab', 'H_Dmab_CB') else OUT / 'de' / f'{lb}.csv'
        d = pd.read_csv(f)
        deg[lb] = d.set_index('gene')
        q = d[d.padj < 0.05].sort_values('log2FoldChange')
        q.to_csv(D / f'deg_{lb}.csv', index=False)
        print(lb, 'padj<0.05:', len(q), f'(up {int((q.log2FoldChange > 0).sum())} / down {int((q.log2FoldChange < 0).sum())});',
              '|lfc|>=1:', int((q.log2FoldChange.abs() >= 1).sum()))
    a, b = deg['H_Dmab'], deg['H_Dmab_CB']
    sa, sb = set(a.index[a.padj < 0.05]), set(b.index[b.padj < 0.05])
    both = sorted(sa & sb)
    same = [g for g in both if np.sign(a.loc[g, 'log2FoldChange']) == np.sign(b.loc[g, 'log2FoldChange'])]
    print('DEG padj<0.05 digested', len(sa), 'centrifuged', len(sb), 'both', len(both), 'same direction', len(same))
    print('  both:', ','.join(both))

    # ---- GSEA の有意 set（上昇・低下とも） ----
    rows, conc = [], []
    for cn in COLLS:
        R = {lb: pd.read_csv(G / f'{cn}_{lb}.csv').set_index('Term') for lb in
             ['H_Dmab', 'H_Dmab_CB', 'H_Dmab_marrowadj', 'H_Dmab_CB_marrowadj', 'H_Dmab_ocadj', 'H_Dmab_CB_ocadj']}
        A, B = R['H_Dmab'], R['H_Dmab_CB']
        terms = A.index.intersection(B.index)
        rho = stats.spearmanr(A.loc[terms, 'NES'], B.loc[terms, 'NES'])[0]
        sigA = A.loc[terms, 'FDR q-val'] < 0.05
        sigB = B.loc[terms, 'FDR q-val'] < 0.05
        sameS = sigA & sigB & (np.sign(A.loc[terms, 'NES']) == np.sign(B.loc[terms, 'NES']))
        oppS = sigA & sigB & (np.sign(A.loc[terms, 'NES']) != np.sign(B.loc[terms, 'NES']))
        conc.append(dict(collection=cn, n_sets=len(terms), spearman_NES=rho,
                         sig_digested_up=int((sigA & (A.loc[terms, 'NES'] > 0)).sum()),
                         sig_digested_down=int((sigA & (A.loc[terms, 'NES'] < 0)).sum()),
                         sig_centrifuged_up=int((sigB & (B.loc[terms, 'NES'] > 0)).sum()),
                         sig_centrifuged_down=int((sigB & (B.loc[terms, 'NES'] < 0)).sum()),
                         sig_both_same=int(sameS.sum()), sig_both_opposite=int(oppS.sum()),
                         both_same_down=','.join(t for t in terms[sameS & (A.loc[terms, 'NES'] < 0)][:40]),
                         both_same_up=','.join(t for t in terms[sameS & (A.loc[terms, 'NES'] > 0)]),
                         both_opposite=','.join(t for t in terms[oppS][:40])))
        for t in terms:
            if not (sigA[t] or sigB[t]):
                continue
            row = dict(collection=cn, term=t, tag_digested=A.loc[t, 'Tag %'], lead_digested=A.loc[t, 'Lead_genes'])
            for lb, r in R.items():
                row[f'NES_{lb}'] = r.loc[t, 'NES'] if t in r.index else np.nan
                row[f'FDR_{lb}'] = r.loc[t, 'FDR q-val'] if t in r.index else np.nan
            rows.append(row)
    S = pd.DataFrame(rows)
    C = pd.DataFrame(conc)
    C.to_csv(D / 'fraction_concordance.csv', index=False)
    print(C.drop(columns=['both_same_down', 'both_same_up', 'both_opposite']).to_string(index=False))
    for _, r in C.iterrows():
        print(r.collection, 'both same up:', r.both_same_up)
        print(r.collection, 'both opposite:', r.both_opposite)

    # ---- 組成との関係（消化後） ----
    x, s = counts('GSE141595/post30')
    sc, _ = scores(x)
    c = np.log2(cpm(x) + 1)
    marrow = zscore(c.loc[MARROW_GENES].mean())
    oc = zscore(sc['osteoclast'])
    ob = zscore(sc['osteoblast'])
    grp = s.group
    pb = {k: stats.pointbiserialr((grp == 'Post_Dmab').astype(int), v)[0] for k, v in
          [('osteoclast', oc), ('osteoblast', ob), ('marrow', marrow)]}
    print('point-biserial r with denosumab:', {k: round(v, 3) for k, v in pb.items()})
    kk = keep_genes(x, grp)
    L = logcpm(x.loc[kk])
    L.index = L.index.str.upper()
    L = L[~L.index.duplicated()]
    r_oc = within_corr(L, grp, oc)
    r_ob = within_corr(L, grp, ob)
    r_mw = within_corr(L, grp, marrow)
    pd.DataFrame({'r_osteoclast': r_oc, 'r_osteoblast': r_ob, 'r_marrow': r_mw}).to_csv(D / 'gene_within_corr.csv')
    bg = dict(r_oc=r_oc.mean(), r_ob=r_ob.mean(), r_mw=r_mw.mean(), oc_hi=(r_oc > 0.5).mean())
    le = []
    for _, row in S[S.FDR_H_Dmab < 0.05].iterrows():
        genes = [g for g in str(row.lead_digested).split(';') if g in r_oc.index]
        if not genes:
            continue
        le.append(dict(collection=row.collection, term=row.term, NES=row.NES_H_Dmab, FDR=row.FDR_H_Dmab, n_lead=len(genes),
                       mean_r_osteoclast=r_oc[genes].mean(), frac_r_osteoclast_gt05=(r_oc[genes] > 0.5).mean(),
                       mean_r_osteoblast=r_ob[genes].mean(), mean_r_marrow=r_mw[genes].mean(),
                       NES_marrowadj=row.NES_H_Dmab_marrowadj, FDR_marrowadj=row.FDR_H_Dmab_marrowadj,
                       NES_ocadj=row.NES_H_Dmab_ocadj, FDR_ocadj=row.FDR_H_Dmab_ocadj,
                       NES_centrifuged=row.NES_H_Dmab_CB, FDR_centrifuged=row.FDR_H_Dmab_CB))
    LE = pd.DataFrame(le)
    LE.to_csv(D / 'le_osteoclast.csv', index=False)
    print('background: mean r_oc %.3f, r_ob %.3f, r_marrow %.3f, frac r_oc>0.5 %.3f' % (bg['r_oc'], bg['r_ob'], bg['r_mw'], bg['oc_hi']))
    with pd.option_context('display.width', 250, 'display.max_colwidth', 60):
        print(LE[LE.collection == 'hallmark'].round(3).to_string(index=False))
        for cn in ['reactome', 'gobp']:
            q = LE[LE.collection == cn]
            print(cn, 'n sig', len(q), 'median mean_r_oc', round(q.mean_r_osteoclast.median(), 3),
                  'frac sets mean_r_oc>0.3', round((q.mean_r_osteoclast > 0.3).mean(), 3),
                  'median mean_r_marrow', round(q.mean_r_marrow.median(), 3))
    S.to_csv(D / 'sets_fdr05.csv', index=False)
    # 共変量補正の NES（Hallmark 全体）
    cov = []
    for lb in ['H_Dmab', 'H_Dmab_marrowadj', 'H_Dmab_ocadj', 'H_Dmab_CB', 'H_Dmab_CB_marrowadj', 'H_Dmab_CB_ocadj']:
        r = pd.read_csv(G / f'hallmark_{lb}.csv')
        cov.append(r[['Term', 'NES', 'FDR q-val']].assign(contrast=lb))
    CV = pd.concat(cov)
    W = CV.pivot(index='Term', columns='contrast', values='NES')
    W.to_csv(D / 'covariate_nes.csv')
    plot(deg, R_h=pd.read_csv(G / 'hallmark_H_Dmab.csv').set_index('Term'),
         R_hc=pd.read_csv(G / 'hallmark_H_Dmab_CB.csv').set_index('Term'))


def spread(ys, gap):
    """ラベルの y を、上から順に最低 gap だけ離す（単純な重なり回避）。"""
    order = np.argsort(ys)[::-1]
    out = np.array(ys, dtype=float)
    last = None
    for i in order:
        if last is not None and out[i] > last - gap:
            out[i] = last - gap
        last = out[i]
    return out


LABEL_GENES = ['SIGLEC15', 'ACP5', 'CTSK', 'ATP6V0D2', 'MMP13', 'COL1A1', 'COL1A2', 'SPARC', 'TNN', 'PRSS35', 'IFITM5',
               'PILRB', 'STAG3L5P']


def volcano(ax, d, title):
    d = d.dropna(subset=['padj'])
    y = -np.log10(d.padj.clip(lower=1e-30))
    sig = d.padj < 0.05
    col = np.where(~sig, '#CCCCCC', np.where(d.log2FoldChange > 0, common.COLOR_HRO, common.COLOR_LRO))
    ax.scatter(d.log2FoldChange, y, s=2, c=col, lw=0, rasterized=True)
    m = max(3.5, d.log2FoldChange.abs().quantile(0.999))
    ax.set_xlim(-m, m)
    ymax = y.max()
    lab = d[sig & d.index.isin(LABEL_GENES)]
    for side in (-1, 1):
        q = lab[np.sign(lab.log2FoldChange) == side]
        if not len(q):
            continue
        yy = -np.log10(q.padj.clip(lower=1e-30)).values
        ty = spread(yy, ymax * 0.07)
        tx = -m * 0.97 if side < 0 else m * 0.97
        for (g, r), y0, y1 in zip(q.iterrows(), yy, ty):
            ax.annotate(g, (r.log2FoldChange, y0), xytext=(tx, y1), fontsize=5.5, style='italic',
                        ha='left' if side < 0 else 'right', va='center',
                        arrowprops=dict(arrowstyle='-', lw=0.3, color='#555555'))
    ax.axhline(-np.log10(0.05), c='k', lw=0.4, ls='--')
    ax.set_xlabel('log2FC denosumab vs placebo')
    ax.set_ylabel('−log10 padj')
    ax.set_title(f'{title}\npadj<0.05: up {int((sig & (d.log2FoldChange > 0)).sum())} / down {int((sig & (d.log2FoldChange < 0)).sum())}',
                 fontsize=7)


def plot(deg, R_h, R_hc):
    common.setup()
    fig, axs = plt.subplots(1, 3, figsize=(common.mm(170), common.mm(70)))
    volcano(axs[0], deg['H_Dmab'], 'GSE141595 digested (n = 15 vs 15)')
    common.panel_label(axs[0], 'A', dx=-0.28)
    volcano(axs[1], deg['H_Dmab_CB'], 'GSE141610 centrifuged (same women)')
    common.panel_label(axs[1], 'B', dx=-0.28)
    ax = axs[2]
    t = R_h.index.intersection(R_hc.index)
    xa, xb = R_h.loc[t, 'NES'], R_hc.loc[t, 'NES']
    fa, fb = R_h.loc[t, 'FDR q-val'], R_hc.loc[t, 'FDR q-val']
    both = (fa < 0.05) & (fb < 0.05)
    one = ((fa < 0.05) | (fb < 0.05)) & ~both
    ax.scatter(xa[~(both | one)], xb[~(both | one)], s=6, c='#CCCCCC', lw=0)
    ax.scatter(xa[one], xb[one], s=8, c='#888888', lw=0)
    ax.scatter(xa[both], xb[both], s=10, c='k', lw=0)
    short = lambda k: k.replace('HALLMARK_', '').replace('_', ' ').lower().replace('saul sen mayo', 'SenMayo')
    lab = [k for k in t if both[k]] + ['SAUL_SEN_MAYO']
    for side in (1, -1):
        ks = [k for k in lab if np.sign(xb[k]) == side]
        yy = np.array([xb[k] for k in ks])
        ty = spread(yy, 0.36)
        for k, y0, y1 in zip(ks, yy, ty):
            ax.annotate(short(k), (xa[k], y0), xytext=(0.3, y1), fontsize=5.3, va='center',
                        arrowprops=dict(arrowstyle='-', lw=0.3, color='#777777'))
    ax.axhline(0, c='k', lw=0.4)
    ax.axvline(0, c='k', lw=0.4)
    ax.set_xlim(-3, 3.6)
    ax.set_ylim(-3.4, 3.4)
    ax.set_xlabel('NES, digested')
    ax.set_ylabel('NES, centrifuged')
    ax.set_title(f'Hallmark: black = FDR<0.05 in both\nSpearman ρ = {stats.spearmanr(xa, xb)[0]:.2f}', fontsize=7)
    common.panel_label(ax, 'C', dx=-0.28)
    fig.tight_layout()
    save(fig, 'dmab', 'dmab')
    plt.close(fig)


if __name__ == '__main__':
    main()

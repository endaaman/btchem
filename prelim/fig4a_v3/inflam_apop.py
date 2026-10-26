"""Inflammatory response・Apoptosis の上昇は、どの比較で、どのくらいの強さで言えるか。

  uv run python -I prelim/fig4a_v3/inflam_apop.py      # gsea4a.py の後

3 つの物差しで並べる:
 (1) prerank GSEA の NES / FDR（gsea4a.py の Hallmark コレクション。遺伝子ラベルの置換なので、遺伝子間相関があると楽観的）
 (2) サンプル単位のセットスコア: 各遺伝子を比較内で z 化し、セット内の平均をサンプルのスコアにする。
     群差を Hedges g（ブートストラップ 95% CI, 2,000 回）と**サンプルラベルの置換 P**（10,000 回）で出す。
     遺伝子間相関があっても正しい（自己完結型の検定）。効果の「強さ」はこちらで読む。
 (3) 中身: GSEA の leading edge に占める即時初期遺伝子（FOS/JUN/EGR1/IER/DUSP1/NR4A1/KLF… = 検体処理で上がる）と
     熱ショック遺伝子の割合。閉経後 vs 閉経前はバッチ（処理）と完全に交絡しているので、ここが多ければ処理の痕跡を疑う。
     加えて、即時初期・熱ショック遺伝子をセットから除いた版で NES を出し直す（noIEG）。
比較: GSE141595 閉経後 vs 閉経前（共通遺伝子）、同 デノスマブ vs placebo（消化後 全遺伝子）、GSE141610 デノスマブ vs placebo（遠心のみ）、
      GSE72815 高齢 vs 若年・高齢+E vs 高齢。注目パネル（prelim/genes.py の炎症・破骨 15 / アポトーシス・老化 17 をヒト化）も同じ形で出す。

出力: out/prelim/fig4a_v3/inflam_apop/{nes_table.csv, setscore_effects.csv, leading_edge.csv, nes_noIEG.csv, inflam_apop.png|pdf}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common4a import (OUT, HDE, counts, logcpm, keep_genes, collections, KEY_SETS, IEG, HSP, wald_rank,  # noqa: E402
                      save, MARKERS_H)
import common  # noqa: E402
import hgenes  # noqa: E402  prelim/human/hgenes.py

D = OUT / 'inflam_apop'
CONTRASTS = [  # (label, counts dir, test, ref, 表示名)
    ('H_PostMeno', 'GSE141595', 'Post', 'Pre', 'Postmeno. vs premeno.\n(GSE141595; batch-confounded)'),
    ('H_Old', 'GSE72815', 'Old', 'Young', 'Old vs young\n(GSE72815)'),
    ('H_OldE', 'GSE72815', 'Old_E', 'Old', 'Old+E2 3 wk vs old\n(GSE72815)'),
    ('H_Dmab', 'GSE141595/post30', 'Post_Dmab', 'Post', 'Denosumab vs placebo\n(GSE141595 digested)'),
    ('H_Dmab_CB', 'GSE141610', 'Post_Dmab', 'Post', 'Denosumab vs placebo\n(GSE141610 centrifuged)'),
]
EXTRA = ['H_PostMeno_procadj', 'H_Dmab_marrowadj', 'H_Dmab_CB_marrowadj', 'H_Dmab_ocadj', 'H_Dmab_CB_ocadj']


def hedges(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    return (np.mean(a) - np.mean(b)) / sp * (1 - 3 / (4 * (na + nb) - 9))


def setsets():
    h = collections()['hallmark']
    S = {k: h[k] for k in KEY_SETS}
    S['PANEL_INFLAM_OC'] = hgenes.INFLAM
    S['PANEL_APOP_SEN'] = hgenes.APOP
    return S


LAB = dict(KEY_SETS, PANEL_INFLAM_OC='Panel: inflam./osteoclast (15)', PANEL_APOP_SEN='Panel: apoptosis/senescence (17)')


def main():
    D.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    S = setsets()
    # ---- (1) NES 表 ----
    rows = []
    for lb in [c[0] for c in CONTRASTS] + EXTRA:
        f = OUT / 'gsea' / f'hallmark_{lb}.csv'
        if not f.exists():
            print('missing', f)
            continue
        r = pd.read_csv(f).set_index('Term')
        for t in KEY_SETS:
            if t in r.index:
                rows.append(dict(contrast=lb, term=t, label=KEY_SETS[t], NES=r.loc[t, 'NES'], FDR=r.loc[t, 'FDR q-val'],
                                 lead=r.loc[t, 'Lead_genes'], tag=r.loc[t, 'Tag %']))
    NES = pd.DataFrame(rows)
    NES.drop(columns='lead').to_csv(D / 'nes_table.csv', index=False)
    print(NES.pivot(index='label', columns='contrast', values='NES').round(2).to_string())
    print(NES.pivot(index='label', columns='contrast', values='FDR').round(3).to_string())

    # ---- (3) leading edge の中身 ----
    comp_markers = {g for k, v in MARKERS_H.items() if 'processing' not in k for g in v}
    le = []
    for _, r in NES.iterrows():
        genes = [g for g in str(r.lead).split(';') if g and g != 'nan']
        n = len(genes)
        le.append(dict(contrast=r.contrast, label=r.label, NES=r.NES, FDR=r.FDR, n_lead=n,
                       n_ieg=sum(g in IEG for g in genes), n_hsp=sum(g in HSP for g in genes),
                       frac_ieg_hsp=(sum(g in IEG or g in HSP for g in genes) / n) if n else np.nan,
                       n_celltype_marker=sum(g in comp_markers for g in genes),
                       ieg_hsp_genes=','.join(g for g in genes if g in IEG or g in HSP),
                       top_lead=','.join(genes[:25])))
    LE = pd.DataFrame(le)
    LE.to_csv(D / 'leading_edge.csv', index=False)
    print(LE[['contrast', 'label', 'NES', 'FDR', 'n_lead', 'n_ieg', 'n_hsp', 'frac_ieg_hsp']].round(3).to_string(index=False))

    # ---- noIEG 版 NES（即時初期・熱ショックを除いたセットで prerank を回し直す） ----
    no = []
    sets_full = {k: v for k, v in S.items()}
    sets_no = {k + '__noIEG': [g for g in v if g not in IEG and g not in HSP] for k, v in S.items()}
    for lb, *_ in CONTRASTS:
        src = HDE / f'{lb}.csv'
        rnk = wald_rank(pd.read_csv(src))
        r = gp.prerank(rnk=rnk, gene_sets={**sets_full, **sets_no}, permutation_num=1000, min_size=4, max_size=1000,
                       seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
        for _, q in r.iterrows():
            base = q.Term.replace('__noIEG', '')
            no.append(dict(contrast=lb, term=base, label=LAB.get(base, base), version='noIEG' if 'noIEG' in q.Term else 'full',
                           NES=float(q.NES), NOM_p=float(q['NOM p-val'])))
    NO = pd.DataFrame(no)
    NO.to_csv(D / 'nes_noIEG.csv', index=False)
    print(NO.pivot_table(index=['label'], columns=['contrast', 'version'], values='NES').round(2).to_string())

    # ---- (2) サンプル単位のセットスコア ----
    eff = []
    for lb, sub, t, r_, name in CONTRASTS:
        x, s = counts(sub)
        s = s[s.group.isin([t, r_])]
        x = x[s.index]
        kk = keep_genes(x, s.group)
        L = logcpm(x.loc[kk])
        L.index = L.index.str.upper()
        L = L[~L.index.duplicated()]
        Z = L.sub(L.mean(1), axis=0).div(L.std(1).replace(0, np.nan), axis=0).dropna()
        a_idx = (s.group == t).values
        for k, genes in S.items():
            gs = [g for g in genes if g in Z.index]
            if len(gs) < 5:
                continue
            sc = Z.loc[gs].mean().values
            a, b = sc[a_idx], sc[~a_idx]
            g = hedges(a, b)
            null = np.empty(10000)
            for i in range(len(null)):
                p = rng.permutation(a_idx)
                null[i] = sc[p].mean() - sc[~p].mean()
            obs = a.mean() - b.mean()
            pperm = (1 + (np.abs(null) >= abs(obs)).sum()) / (1 + len(null))
            boots = []
            for _ in range(2000):
                boots.append(hedges(rng.choice(a, len(a)), rng.choice(b, len(b))))
            gnoi = [x_ for x_ in gs if x_ not in IEG and x_ not in HSP]
            sc2 = Z.loc[gnoi].mean().values
            eff.append(dict(contrast=lb, name=name.replace('\n', ' '), term=k, label=LAB.get(k, k), n_genes=len(gs),
                            n_test=int(a_idx.sum()), n_ref=int((~a_idx).sum()), mean_diff_z=obs, hedges_g=g,
                            g_ci_lo=np.quantile(boots, 0.025), g_ci_hi=np.quantile(boots, 0.975), perm_p=pperm,
                            hedges_g_noIEG=hedges(sc2[a_idx], sc2[~a_idx]), n_genes_noIEG=len(gnoi)))
    E = pd.DataFrame(eff)
    from statsmodels.stats.multitest import multipletests  # noqa: E402
    E['perm_fdr'] = np.nan
    for c, idx in E.groupby('contrast').groups.items():
        E.loc[idx, 'perm_fdr'] = multipletests(E.loc[idx, 'perm_p'], method='fdr_bh')[1]
    E.to_csv(D / 'setscore_effects.csv', index=False)
    print(E[['contrast', 'label', 'n_genes', 'hedges_g', 'g_ci_lo', 'g_ci_hi', 'perm_p', 'perm_fdr', 'hedges_g_noIEG']]
          .round(3).to_string(index=False))
    plot(E, NES)


def plot(E, NES):
    common.setup()
    order = ['HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_IL6_JAK_STAT3_SIGNALING',
             'PANEL_INFLAM_OC', 'HALLMARK_APOPTOSIS', 'HALLMARK_P53_PATHWAY', 'SAUL_SEN_MAYO', 'PANEL_APOP_SEN',
             'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE']
    cons = [c for c in CONTRASTS]
    fig, axs = plt.subplots(1, len(cons), figsize=(common.mm(170), common.mm(80)), sharey=True)
    for j, (lb, _, _, _, name) in enumerate(cons):
        ax = axs[j]
        q = E[E.contrast == lb].set_index('term')
        for i, t in enumerate(order):
            if t not in q.index:
                ax.text(0, i, 'n/a', ha='center', va='center', fontsize=6, color='#999999')
                continue
            r = q.loc[t]
            col = (common.COLOR_HRO if r.hedges_g > 0 else common.COLOR_LRO) if r.perm_fdr < 0.05 else '#888888'
            ax.plot([r.g_ci_lo, r.g_ci_hi], [i, i], c=col, lw=1.0)
            ax.plot(r.hedges_g, i, 'o', ms=3.5, c=col, mfc=col if r.perm_fdr < 0.05 else 'white')
            nes = NES[(NES.contrast == lb) & (NES.term == t)]
            if len(nes):
                n = nes.iloc[0]
                ax.text(6.3, i, f'{n.NES:+.1f}{"*" if n.FDR < 0.05 else ""}'.replace('-', '−'), fontsize=5.5,
                        va='center', ha='right', color='#444444')
        ax.axvline(0, c='k', lw=0.4)
        ax.set_xlim(-3, 6.4)
        ax.axvspan(4.0, 6.4, color='#F4F4F4', lw=0)
        ax.set_xticks([-2, 0, 2])
        ax.set_title(name, fontsize=6.5)
        ax.set_xlabel('Hedges g (95% CI)', fontsize=7)
        ax.tick_params(labelsize=6.5)
        if j == 0:
            ax.set_yticks(range(len(order)), [LAB[t] for t in order], fontsize=6.5)
        ax.set_ylim(len(order) - 0.5, -0.8)
        ax.text(6.3, -0.75, 'NES', fontsize=5.5, ha='right', color='#444444')
    fig.suptitle('Sample-level set scores (filled: permutation FDR < 0.05). Right: prerank NES (* FDR < 0.05)',
                 fontsize=7, x=0.01, ha='left')
    fig.tight_layout()
    save(fig, 'inflam_apop', 'inflam_apop')
    plt.close(fig)


if __name__ == '__main__':
    main()

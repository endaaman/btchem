"""「閉経後に起こる骨の変化の一部をデノスマブが抑える」を直接検証する。

  uv run python -I prelim/fig4a_v4/reversal.py

問い: 閉経後 vs 閉経前（GSE141595）で動いた遺伝子・遺伝子セットのうち、デノスマブ vs placebo（閉経後）で
逆向きに戻るもの / 戻らないものはどれだけあるか。

注意（設計上の落とし穴）:
1. **共有参照のバイアス**: 「閉経後 placebo − 閉経前」と「デノスマブ − 閉経後 placebo」は同じ placebo 15 例を
   逆符号で使う。placebo 群の平均の偶然の揺れが両方に逆符号で入るので、何もしなくても「戻る」ように見える
   （平均への回帰）。そこで 3 通りで見る:
   (a) naive: 既存 DESeq2（H_PostMeno と H_Dmab_common）をそのまま比べる（バイアス込みの見かけの値）
   (b) split-half: placebo 15 例をランダムに 8 / 7 に分け、8 例を「閉経後」側（vs 閉経前）、7 例を「デノスマブの対照」側に使う。
       2 つの比較が別のサンプルになるのでバイアスが消える。200 回の分割の平均を報告する。
       速度のため DESeq2 ではなく log2 正規化カウント（median-of-ratios）の平均差と Welch t を使う。
   (c) 置換の帰無分布: 閉経後 30 例の placebo / denosumab ラベルを入れ替えて (a)(b) と同じ統計量を計算（各 500 回）。
       (a) の観察値はこの帰無（バイアス込み）と比べて初めて解釈できる。
2. **バッチ交絡**: 閉経前 15 例は別ファイル（2022 年追加・ライブラリ深度が浅い）で、閉経後 vs 閉経前の差には
   バッチの差が混ざる。デノスマブはバッチの差を戻せないので、「戻る割合」は過小評価の方向に偏る。
   そこで、(d) バッチの揃った独立コホート GSE72815（高齢 vs 若年女性）の変化を「加齢・閉経軸」として、
   デノスマブ（GSE141595 消化後・GSE141610 遠心のみ）がそれを戻すかも見る。サンプルの重なりが無いのでバイアスは無い。
   帰無は GSE141595 の placebo / denosumab ラベル置換（500 回）。
3. 閉経前ファイルの欠落（chr12 後半〜22, Y, M）で、閉経後 vs 閉経前は 16,823 遺伝子しか比べられない。

統計量（「閉経で動いた遺伝子」= DEG 集合ごと、上昇・低下別）:
  frac_opposite   デノスマブの log2FC が閉経の変化と逆符号の割合（帰無で 0.5 付近、バイアスがあるとそれ以上）
  slope           DEG 集合での回帰 lfc_dmab ~ lfc_meno の傾き（−1 = 完全に戻る、0 = 戻らない）。restoration = −slope
  frac_rev_sig    逆符号かつデノスマブ側で名目 p < 0.05 の割合
遺伝子セット（Hallmark 50 + SenMayo）: セット内遺伝子の lfc 平均 − 全遺伝子の lfc 平均（競合的な平均差）で同じことをする。

出力: out/prelim/fig4a_v4/reversal/{summary.csv, genes.csv, sets_splithalf.csv, sets_nes.csv, null_*.csv, reversal.png|pdf}
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
from common4a import OUT, HDE, counts, logcpm, keep_genes, collections, KEY_SETS, save  # noqa: E402
import common  # noqa: E402

D = OUT / 'reversal'
NSPLIT, NPERM = 200, 500
LFC_MIN = 0.585


def welch(L, a, b):
    A, B = L[a].values, L[b].values
    na, nb = A.shape[1], B.shape[1]
    ma, mb = A.mean(1), B.mean(1)
    va, vb = A.var(1, ddof=1) / na, B.var(1, ddof=1) / nb
    se = np.sqrt(va + vb)
    with np.errstate(divide='ignore', invalid='ignore'):
        t = (ma - mb) / se
        df = (va + vb) ** 2 / (va ** 2 / (na - 1) + vb ** 2 / (nb - 1))
    p = 2 * stats.t.sf(np.abs(t), df)
    return pd.DataFrame({'lfc': ma - mb, 't': t, 'p': p}, index=L.index)


def metrics(l1, l2, p2, sel_up, sel_dn):
    """l1: 閉経（など）の lfc、l2: デノスマブの lfc、p2: デノスマブの p。sel_*: bool Series。"""
    out = {}
    for nm, sel in [('up', sel_up), ('down', sel_dn), ('all', sel_up | sel_dn)]:
        a, b, pb = l1[sel], l2[sel], p2[sel]
        n = int(sel.sum())
        if n < 3:
            out.update({f'{nm}_n': n})
            continue
        opp = np.sign(a) != np.sign(b)
        out[f'{nm}_n'] = n
        out[f'{nm}_frac_opposite'] = float(opp.mean())
        out[f'{nm}_frac_rev_sig'] = float((opp & (pb < 0.05)).mean())
        out[f'{nm}_frac_same_sig'] = float((~opp & (pb < 0.05)).mean())
        out[f'{nm}_median_lfc_dmab'] = float(np.median(b))
        if nm == 'all':
            out['all_slope'] = float(np.polyfit(a, b, 1)[0])
    return out


def setlevel(l, sets):
    """競合的なセット平均差: セット内遺伝子の lfc 平均 − 全体平均。"""
    g = l.index
    mu = l.mean()
    return pd.Series({k: l[g.isin(v)].mean() - mu for k, v in sets.items() if g.isin(v).sum() >= 10})


def main():
    D.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    sets = collections()['hallmark']

    # ---------- データ ----------
    x, s = counts('GSE141595')               # 45 例・共通遺伝子
    pre = s.index[s.group == 'Pre']
    pla = s.index[s.group == 'Post']
    dmb = s.index[s.group == 'Post_Dmab']
    keep = keep_genes(x, pd.Series(np.where(s.group == 'Pre', 'a', 'b'), index=s.index))
    L = logcpm(x.loc[keep])
    L.index = L.index.str.upper()
    post_all = np.array(list(pla) + list(dmb))

    rows, null_rows, gene_rows, set_rows = [], [], [], []

    # ---------- (a) naive: DESeq2 の既存表 ----------
    m = pd.read_csv(HDE / 'H_PostMeno.csv').set_index('gene')
    dm = pd.read_csv(HDE / 'H_Dmab_common.csv').set_index('gene')
    g = m.index.intersection(dm.index)
    m, dm = m.loc[g], dm.loc[g]
    ok = m.padj.notna() & dm.pvalue.notna()
    m, dm = m[ok], dm[ok]
    for thr_name, lfc_min in [('strict', 1.0), ('lenient', LFC_MIN)]:
        up = (m.padj < 0.05) & (m.log2FoldChange >= lfc_min)
        dn = (m.padj < 0.05) & (m.log2FoldChange <= -lfc_min)
        r = metrics(m.log2FoldChange, dm.log2FoldChange, dm.pvalue, up, dn)
        rows.append(dict(method='naive DESeq2 (H_PostMeno vs H_Dmab_common)', axis='PostMeno', dmab='digested',
                         threshold=thr_name, **r))
    r_all = stats.spearmanr(m.log2FoldChange, dm.log2FoldChange)[0]
    rows.append(dict(method='naive DESeq2 (H_PostMeno vs H_Dmab_common)', axis='PostMeno', dmab='digested',
                     threshold='genome-wide', spearman=r_all, n=len(m)))
    gene_tab = pd.DataFrame({'lfc_postmeno': m.log2FoldChange, 'padj_postmeno': m.padj,
                             'lfc_dmab_naive': dm.log2FoldChange, 'p_dmab_naive': dm.pvalue, 'padj_dmab_naive': dm.padj})

    # ---------- (a') naive の fast 版と、その置換帰無（バイアス込みの帰無） ----------
    def naive_fast(pl, dx):
        a = welch(L, list(pl), list(pre))
        b = welch(L, list(dx), list(pl))
        sel = (a.p < 0.001) & (a.lfc.abs() >= LFC_MIN)
        r = metrics(a.lfc, b.lfc, b.p, sel & (a.lfc > 0), sel & (a.lfc < 0))
        r['spearman'] = stats.spearmanr(a.lfc, b.lfc)[0]
        return r, a, b
    obs_naive, a0, b0 = naive_fast(pla, dmb)
    rows.append(dict(method='naive fast (logCPM, Welch p<0.001 & |lfc|>=0.585)', axis='PostMeno', dmab='digested',
                     threshold='fast', **obs_naive))
    for i in range(NPERM):
        p = rng.permutation(post_all)
        r, _, _ = naive_fast(p[:15], p[15:])
        null_rows.append(dict(kind='naive_fast', i=i, **r))

    # ---------- (b) split-half ----------
    def split_half(pl, dx, rng, nrep):
        acc, sacc, lfc1s, lfc2s = [], [], [], []
        for _ in range(nrep):
            pp = rng.permutation(np.array(pl))
            A, B = pp[:8], pp[8:]
            a = welch(L, list(A), list(pre))
            b = welch(L, list(dx), list(B))
            sel = (a.p < 0.001) & (a.lfc.abs() >= LFC_MIN)
            r = metrics(a.lfc, b.lfc, b.p, sel & (a.lfc > 0), sel & (a.lfc < 0))
            r['spearman'] = stats.spearmanr(a.lfc, b.lfc)[0]
            acc.append(r)
            sa, sb = setlevel(a.lfc, sets), setlevel(b.lfc, sets)
            sacc.append(pd.DataFrame({'meno': sa, 'dmab': sb}))
            lfc1s.append(a.lfc)
            lfc2s.append(b.lfc)
        R = pd.DataFrame(acc).mean().to_dict()
        S = pd.concat(sacc, keys=range(nrep)).groupby(level=1).mean()
        return R, S, pd.concat(lfc1s, axis=1).mean(1), pd.concat(lfc2s, axis=1).mean(1)
    obs_split, S_obs, l1_split, l2_split = split_half(pla, dmb, rng, NSPLIT)
    rows.append(dict(method=f'split-half (placebo 8 vs pre / dmab vs placebo 7; {NSPLIT} splits)', axis='PostMeno',
                     dmab='digested', threshold='fast', **obs_split))
    gene_tab = gene_tab.join(pd.DataFrame({'lfc_postmeno_split': l1_split, 'lfc_dmab_split': l2_split}), how='outer')
    S_null = []
    for i in range(NPERM):
        p = rng.permutation(post_all)
        r, S, _, _ = split_half(p[:15], p[15:], rng, 1)
        null_rows.append(dict(kind='split_half', i=i, **r))
        S_null.append(S['dmab'].rename(i))
    S_null = pd.concat(S_null, axis=1)

    # set-level split-half
    for k in S_obs.index:
        nd = S_null.loc[k].values
        set_rows.append(dict(term=k, label=KEY_SETS.get(k, k), meno_setdiff=S_obs.loc[k, 'meno'],
                             dmab_setdiff=S_obs.loc[k, 'dmab'],
                             dmab_perm_p=(1 + (np.abs(nd) >= abs(S_obs.loc[k, 'dmab'])).sum()) / (1 + len(nd)),
                             reverses=np.sign(S_obs.loc[k, 'meno']) != np.sign(S_obs.loc[k, 'dmab'])))
    SET = pd.DataFrame(set_rows)
    SET['dmab_perm_fdr'] = stats.false_discovery_control(SET.dmab_perm_p)
    SET.sort_values('meno_setdiff', ascending=False).to_csv(D / 'sets_splithalf.csv', index=False)

    # ---------- (d) 独立コホート GSE72815（高齢 vs 若年）を軸にする ----------
    old = pd.read_csv(HDE / 'H_Old.csv').dropna(subset=['padj']).set_index('gene')
    old.index = old.index.str.upper()
    for sub, nm in [('GSE141595/post30', 'digested'), ('GSE141610', 'centrifuged')]:
        xx, ss = counts(sub)
        pl2, dx2 = ss.index[ss.group == 'Post'], ss.index[ss.group == 'Post_Dmab']
        kk = keep_genes(xx, ss.group)
        LL = logcpm(xx.loc[kk])
        LL.index = LL.index.str.upper()
        LL = LL[~LL.index.duplicated()]
        gg = old.index.intersection(LL.index)
        o = old.loc[gg]
        for thr_name, lfc_min in [('strict', 1.0), ('lenient', LFC_MIN)]:
            up = (o.padj < 0.05) & (o.log2FoldChange >= lfc_min)
            dn = (o.padj < 0.05) & (o.log2FoldChange <= -lfc_min)

            def stat(pl, dx):
                b = welch(LL.loc[gg], list(dx), list(pl))
                r = metrics(o.log2FoldChange, b.lfc, b.p, up, dn)
                r['spearman'] = stats.spearmanr(o.log2FoldChange, b.lfc)[0]
                return r
            obs = stat(pl2, dx2)
            nul = []
            pa = np.array(list(pl2) + list(dx2))
            for i in range(NPERM):
                p = rng.permutation(pa)
                nul.append(stat(p[:15], p[15:]))
            nul = pd.DataFrame(nul)
            nul.assign(kind=f'old_axis_{nm}_{thr_name}').to_csv(D / f'null_old_axis_{nm}_{thr_name}.csv', index=False)
            r = dict(method=f'independent axis: GSE72815 old vs young (DESeq2) x Dmab {nm} (fast)', axis='Old (GSE72815)',
                     dmab=nm, threshold=thr_name, **obs)
            for k in ['all_frac_opposite', 'all_slope', 'spearman', 'up_frac_opposite', 'down_frac_opposite']:
                if k in obs:
                    nv = nul[k].dropna().values
                    r[f'{k}_null_mean'] = nv.mean()
                    r[f'{k}_null_q025'] = np.quantile(nv, 0.025)
                    r[f'{k}_null_q975'] = np.quantile(nv, 0.975)
                    # 片側: 「戻る」方向（frac_opposite 大、slope・spearman 小）
                    r[f'{k}_perm_p'] = ((1 + (nv >= obs[k]).sum()) if 'opposite' in k else (1 + (nv <= obs[k]).sum())) / (1 + len(nv))
            rows.append(r)

    # ---------- 帰無との比較（(a') と (b)） ----------
    N = pd.DataFrame(null_rows)
    N.to_csv(D / 'null_postmeno_axis.csv', index=False)
    for r in rows:
        kind = 'naive_fast' if r['method'].startswith('naive fast') else 'split_half' if r['method'].startswith('split-half') else None
        if kind is None:
            continue
        nn = N[N.kind == kind]
        for k in ['all_frac_opposite', 'all_slope', 'spearman', 'up_frac_opposite', 'down_frac_opposite',
                  'all_frac_rev_sig']:
            if k in r and k in nn:
                nv = nn[k].dropna().values
                r[f'{k}_null_mean'] = nv.mean()
                r[f'{k}_null_q025'] = np.quantile(nv, 0.025)
                r[f'{k}_null_q975'] = np.quantile(nv, 0.975)
                r[f'{k}_perm_p'] = ((1 + (nv >= r[k]).sum()) if ('opposite' in k or 'rev_sig' in k)
                                    else (1 + (nv <= r[k]).sum())) / (1 + len(nv))
    R = pd.DataFrame(rows)
    R.to_csv(D / 'summary.csv', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 60):
        print(R.T.to_string())

    # ---------- 遺伝子表（閉経 DEG の内訳） ----------
    gene_tab['meno_deg'] = np.where((gene_tab.padj_postmeno < 0.05) & (gene_tab.lfc_postmeno >= LFC_MIN), 'up',
                                    np.where((gene_tab.padj_postmeno < 0.05) & (gene_tab.lfc_postmeno <= -LFC_MIN), 'down', ''))
    o2 = old[['log2FoldChange', 'padj']].rename(columns={'log2FoldChange': 'lfc_old', 'padj': 'padj_old'})
    gene_tab.index = gene_tab.index.str.upper()
    gene_tab = gene_tab[~gene_tab.index.duplicated()].join(o2[~o2.index.duplicated()], how='left')
    dmf = pd.read_csv(HDE / 'H_Dmab.csv').set_index('gene')[['log2FoldChange', 'pvalue', 'padj']]
    dmf.columns = ['lfc_dmab_full', 'p_dmab_full', 'padj_dmab_full']
    dmf.index = dmf.index.str.upper()
    gene_tab = gene_tab.join(dmf[~dmf.index.duplicated()], how='left')
    gene_tab['call'] = ''
    sel = gene_tab.meno_deg != ''
    opp = np.sign(gene_tab.lfc_postmeno) != np.sign(gene_tab.lfc_dmab_split)
    gene_tab.loc[sel & opp & (gene_tab.p_dmab_naive < 0.05), 'call'] = 'reversed (Dmab p<0.05)'
    gene_tab.loc[sel & opp & ~(gene_tab.p_dmab_naive < 0.05), 'call'] = 'opposite, n.s.'
    gene_tab.loc[sel & ~opp, 'call'] = 'same direction'
    gene_tab.index.name = 'gene'
    gene_tab.sort_values('lfc_postmeno', ascending=False).to_csv(D / 'genes.csv')

    nes_classify()
    plot(gene_tab, R, N, SET)


def nes_classify():
    """セット単位: 閉経後 vs 閉経前（naive、共有参照のバイアス込み）または高齢 vs 若年（独立）で FDR < 0.05 の set が、
    デノスマブ（消化後 / 遠心のみ）でどちら向きに動くか。デノスマブは多くの set を一様に下げる傾向があるので、
    「動かなかった set」での同じ割合を基準にして Fisher 検定で比べる（baseline 列）。"""
    from scipy.stats import fisher_exact
    nes = []
    for cn in ['hallmark', 'reactome', 'gobp']:
        fs = {lb: OUT / 'gsea' / f'{cn}_{lb}.csv' for lb in ['H_PostMeno', 'H_Old', 'H_Dmab', 'H_Dmab_CB']}
        if not all(f.exists() for f in fs.values()):
            continue
        T = {lb: pd.read_csv(f).set_index('Term')[['NES', 'FDR q-val']] for lb, f in fs.items()}
        terms = T['H_PostMeno'].index.union(T['H_Old'].index)
        for t in terms:
            row = dict(collection=cn, term=t)
            for lb in T:
                if t in T[lb].index:
                    row[f'NES_{lb}'] = T[lb].loc[t, 'NES']
                    row[f'FDR_{lb}'] = T[lb].loc[t, 'FDR q-val']
            nes.append(row)
    if not nes:
        return
    NE = pd.DataFrame(nes)

    def one(r, ax, dm):
        n0, f0 = r.get(f'NES_{ax}', np.nan), r.get(f'FDR_{ax}', np.nan)
        f, n = r.get(f'FDR_{dm}', np.nan), r.get(f'NES_{dm}', np.nan)
        if np.isnan(n) or np.isnan(n0):
            return 'n/a'
        opp = np.sign(n) != np.sign(n0)
        if opp and f < 0.05:
            return 'reversed'
        if opp and f < 0.25:
            return 'reversed (trend)'
        if not opp and f < 0.05:
            return 'same direction'
        if not opp and f < 0.25:
            return 'same (trend)'
        return 'not reversed'
    summ = []
    for ax, axn in [('H_PostMeno', 'postmeno'), ('H_Old', 'old')]:
        for dm, dmn in [('H_Dmab', 'digested'), ('H_Dmab_CB', 'centrifuged')]:
            NE[f'call_{axn}_{dmn}'] = [one(r, ax, dm) if r.get(f'FDR_{ax}', 1) < 0.05 else '' for _, r in NE.iterrows()]
            for cn in NE.collection.unique():
                q = NE[NE.collection == cn].dropna(subset=[f'NES_{ax}', f'NES_{dm}'])
                for sgn, nm in [(1, 'up'), (-1, 'down')]:
                    ch = q[(q[f'FDR_{ax}'] < 0.05) & (np.sign(q[f'NES_{ax}']) == sgn)]
                    un = q[q[f'FDR_{ax}'] >= 0.25]
                    # 「戻る」= デノスマブで逆向きに FDR < 0.05
                    rev_ch = ((np.sign(ch[f'NES_{dm}']) == -sgn) & (ch[f'FDR_{dm}'] < 0.05)).sum()
                    rev_un = ((np.sign(un[f'NES_{dm}']) == -sgn) & (un[f'FDR_{dm}'] < 0.05)).sum()
                    same_ch = ((np.sign(ch[f'NES_{dm}']) == sgn) & (ch[f'FDR_{dm}'] < 0.05)).sum()
                    p = fisher_exact([[rev_ch, len(ch) - rev_ch], [rev_un, len(un) - rev_un]])[1] if len(ch) and len(un) else np.nan
                    summ.append(dict(axis=axn, dmab=dmn, collection=cn, direction=nm, n_changed=len(ch),
                                     reversed_fdr05=int(rev_ch), frac_reversed=rev_ch / len(ch) if len(ch) else np.nan,
                                     same_dir_fdr05=int(same_ch), n_unchanged=len(un),
                                     baseline_frac_unchanged=rev_un / len(un) if len(un) else np.nan, fisher_p=p))
    NE.to_csv(D / 'sets_nes.csv', index=False)
    # 「戻った」set の中身: 加齢・閉経側の leading edge 遺伝子が、デノスマブでどう動いたか
    stat = {}
    for lb in ['H_Dmab', 'H_Dmab_CB']:
        d = pd.read_csv(HDE / f'{lb}.csv').dropna(subset=['lfcSE'])
        d = d[d.lfcSE > 0]
        stat[lb] = (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean()
    le_rows = []
    for ax, axn in [('H_PostMeno', 'postmeno'), ('H_Old', 'old')]:
        for dm, dmn in [('H_Dmab', 'digested'), ('H_Dmab_CB', 'centrifuged')]:
            q = NE[NE[f'call_{axn}_{dmn}'] == 'reversed']
            for _, r in q.iterrows():
                A = pd.read_csv(OUT / 'gsea' / f'{r.collection}_{ax}.csv').set_index('Term')
                B = pd.read_csv(OUT / 'gsea' / f'{r.collection}_{dm}.csv').set_index('Term')
                la = [g for g in str(A.loc[r.term, 'Lead_genes']).split(';') if g]
                lb_ = [g for g in str(B.loc[r.term, 'Lead_genes']).split(';') if g]
                st = stat[dm].reindex(la).dropna()
                sgn = np.sign(r[f'NES_{ax}'])
                le_rows.append(dict(axis=axn, dmab=dmn, collection=r.collection, term=r.term, NES_axis=r[f'NES_{ax}'],
                                    NES_dmab=r[f'NES_{dm}'], FDR_dmab=r[f'FDR_{dm}'], n_le_axis=len(la), n_le_dmab=len(lb_),
                                    n_le_overlap=len(set(la) & set(lb_)),
                                    jaccard=len(set(la) & set(lb_)) / max(1, len(set(la) | set(lb_))),
                                    axisLE_mean_dmab_stat=st.mean(),
                                    axisLE_frac_reversed=float((np.sign(st) == -sgn).mean()) if len(st) else np.nan,
                                    axisLE_frac_reversed_p05=float(((np.sign(st) == -sgn) & (st.abs() > 1.96)).mean()) if len(st) else np.nan,
                                    le_axis_top=','.join(la[:12]), le_dmab_top=','.join(lb_[:12])))
    LR = pd.DataFrame(le_rows)
    LR.to_csv(D / 'sets_nes_reversed_le.csv', index=False)
    if len(LR):
        with pd.option_context('display.width', 250, 'display.max_colwidth', 50):
            print(LR.drop(columns=['le_axis_top', 'le_dmab_top']).round(3).to_string(index=False))
    SM = pd.DataFrame(summ)
    SM.to_csv(D / 'sets_nes_summary.csv', index=False)
    with pd.option_context('display.width', 200):
        print(SM.round(3).to_string(index=False))


def plot(G, R, N, SET):
    common.setup()
    fig, axs = plt.subplots(1, 3, figsize=(common.mm(170), common.mm(66)),
                            gridspec_kw=dict(width_ratios=[1, 1.05, 1], wspace=0.75))
    # A: split-half の遺伝子散布（閉経 DEG を色付け）
    ax = axs[0]
    g = G.dropna(subset=['lfc_postmeno_split', 'lfc_dmab_split'])
    deg = g.meno_deg.fillna('') != ''
    ax.scatter(g.lfc_postmeno_split[~deg], g.lfc_dmab_split[~deg], s=1.5, c='#CCCCCC', lw=0, rasterized=True)
    ax.scatter(g.lfc_postmeno_split[deg], g.lfc_dmab_split[deg], s=3, lw=0, rasterized=True,
               c=np.where(g.lfc_postmeno[deg] > 0, common.COLOR_HRO, common.COLOR_LRO))
    for gene, dx, dy in [('CDKN1A', 6, 14), ('IL6', -4, -16), ('HSPA1A', 4, 6)]:
        if gene in g.index:
            ax.annotate(gene, (g.loc[gene, 'lfc_postmeno_split'], g.loc[gene, 'lfc_dmab_split']), fontsize=6,
                        xytext=(dx, dy), textcoords='offset points', style='italic',
                        arrowprops=dict(arrowstyle='-', lw=0.4, color='k'))
    ax.axhline(0, c='k', lw=0.4)
    ax.axvline(0, c='k', lw=0.4)
    rs = R[R.method.str.startswith('split-half')].iloc[0]
    ax.set_title(f'slope {rs.all_slope:+.2f} (95% null {rs.all_slope_null_q025:+.2f} to {rs.all_slope_null_q975:+.2f})'
                 .replace('-', '−'), fontsize=6.5)
    ax.set_xlabel('log2FC postmeno. vs premeno.\n(placebo half A)', fontsize=7.5)
    ax.set_ylabel('log2FC denosumab vs placebo\n(placebo half B)', fontsize=7.5)
    common.panel_label(ax, 'A', dx=-0.42, dy=1.08)
    # B: frac_opposite の観察値と帰無
    ax = axs[1]
    items = [('naive fast', 'naive (shared\nplacebo)'), ('split-half', 'split-half'),
             ('independent axis: GSE72815 old vs young (DESeq2) x Dmab digested', 'old-vs-young axis\n× Dmab, digested'),
             ('independent axis: GSE72815 old vs young (DESeq2) x Dmab centrifuged', 'old-vs-young axis\n× Dmab, centrifuged')]
    for i, (key, lab) in enumerate(items):
        q = R[R.method.str.startswith(key) & R.threshold.isin(['fast', 'lenient'])]
        if not len(q):
            continue
        q = q.iloc[0]
        ax.plot([q.all_frac_opposite_null_q025, q.all_frac_opposite_null_q975], [i, i], c='#BBBBBB', lw=4, solid_capstyle='butt')
        ax.plot(q.all_frac_opposite, i, 'o', c='k', ms=4)
        ax.text(q.all_frac_opposite, i - 0.3, f'{q.all_frac_opposite:.2f} (P = {q.all_frac_opposite_perm_p:.2f})',
                fontsize=5.5, ha='center')
    ax.axvline(0.5, c='k', lw=0.4, ls=':')
    ax.set_yticks(range(len(items)), [i[1] for i in items], fontsize=6.5)
    ax.set_ylim(len(items) - 0.5, -0.7)
    ax.set_xlim(0, 1)
    ax.set_xlabel('fraction of age/menopause DEGs\nmoving the opposite way with denosumab', fontsize=7.5)
    ax.set_title('grey bar: 95% permutation null', fontsize=6.5)
    common.panel_label(ax, 'B', dx=-0.85, dy=1.08)
    # C: セット単位（split-half のセット平均差）
    ax = axs[2]
    S = SET.copy()
    ax.scatter(S.meno_setdiff, S.dmab_setdiff, s=6, c=np.where(S.dmab_perm_fdr < 0.05, 'k', '#999999'), lw=0)
    lab = {'HALLMARK_INFLAMMATORY_RESPONSE': ('Inflam. resp.', 4, 8), 'HALLMARK_TNFA_SIGNALING_VIA_NFKB': ('TNF-α', -10, -12),
           'HALLMARK_APOPTOSIS': ('Apoptosis', -30, 10), 'SAUL_SEN_MAYO': ('SenMayo', 6, -10),
           'HALLMARK_P53_PATHWAY': ('p53', -12, 10), 'HALLMARK_MYOGENESIS': ('Myogenesis', -20, 10)}
    for _, r in S.iterrows():
        if r.term in lab:
            t, dx, dy = lab[r.term]
            ax.annotate(t, (r.meno_setdiff, r.dmab_setdiff), fontsize=5.5, xytext=(dx, dy), textcoords='offset points',
                        arrowprops=dict(arrowstyle='-', lw=0.3, color='#444444'))
    ax.axhline(0, c='k', lw=0.4)
    ax.axvline(0, c='k', lw=0.4)
    ax.set_xlabel('set mean log2FC − all genes,\npostmeno. vs premeno. (split-half)', fontsize=7.5)
    ax.set_ylabel('same, denosumab vs placebo', fontsize=7.5)
    ax.set_title(f'Hallmark + SenMayo: perm. FDR < 0.05 in {int((S.dmab_perm_fdr < 0.05).sum())}/{len(S)}', fontsize=6.5)
    common.panel_label(ax, 'C', dx=-0.42, dy=1.08)
    save(fig, 'reversal', 'reversal')
    plt.close(fig)


def replot():
    G = pd.read_csv(D / 'genes.csv').set_index('gene')
    R = pd.read_csv(D / 'summary.csv')
    SET = pd.read_csv(D / 'sets_splithalf.csv')
    plot(G, R, None, SET)


if __name__ == '__main__':
    if '--plot' in sys.argv:
        replot()
    elif '--nes' in sys.argv:
        nes_classify()
    else:
        main()

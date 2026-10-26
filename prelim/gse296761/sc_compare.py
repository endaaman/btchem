"""GSE296752: 細胞種ごとに DOX vs PBS を記述する（各群 2 匹プールの 1 ライブラリ = n=1 対 n=1）。

統計的な反復が無いので、検定の P 値は出さないか、出すときは「細胞を反復とした擬似反復」と明記する。出すもの:
  1. 擬似バルク: 細胞種ごとに UMI を合算 → CPM → log2((CPM_DOX+1)/(CPM_PBS+1))（prelim/chemo/sc_pseudobulk.py と同じ）。
     パネル遺伝子には細胞ブートストラップ（群内で細胞を復元抽出、500 回）の 2.5–97.5% を付ける
     = 細胞の拾い方の揺らぎだけで、個体差は含まない。
  2. 擬似バルク log2FC（平均 CPM >= 10）を順位にした prerank GSEA（prelim/chemo2/gsea.py の sc_CAR と同じ作法）。
     72 set（chemo2 コレクション、要約図と同じ FDR の母集団）と、追加 set + 手製 set（panels.extra_collection）を別々に回す。
     FDR は参考値（順位が統計量ではない・n=1）。
  3. 細胞ごとの gene set スコア（scanpy score_genes）の DOX − PBS 差と AUC（DOX の細胞が PBS の細胞より高い確率）。
     Wilcoxon の P も出すが、細胞を反復とした擬似反復で、個体差を含まない（探索的）。
     DOX ライブラリは細胞あたり UMI が 1.3–1.8 倍（骨芽細胞は 3.3 倍）深いので、スコアと 4 の前に
     細胞種ごとに DOX 細胞の UMI を二項間引き（p = PBS の UMI 中央値 / DOX の UMI 中央値）して深さを揃える（thin_factors.csv）。
  4. 細胞レベルの DE（Wilcoxon, DOX vs PBS, 細胞種ごと、間引き後）: 探索的。上位だけ保存。
  5. bulk との整合: Lepr+ MSC（MSC-Adipo + MSC-Osteo）の擬似バルク log2FC と GSE296758 Dox×4（LepR+ MSPC bulk, n=4/4）、
     Osteoblast 擬似バルクと GSE296757（Col2.3+ 骨芽細胞 bulk）の log2FC の Spearman。

  uv run python -I prelim/gse296761/sc_compare.py      # repo root から（sc_annotate.py の後）
出力: out/prelim/gse296761/sc/{pseudobulk_all_genes.csv.gz, panel_pseudobulk.csv, gsea_<ct>.csv, gsea_extra_<ct>.csv,
      sc_nes_long.csv, scores_by_celltype.csv, thin_factors.csv, cellde_<ct>.csv, bulk_concordance.csv},
      out/prelim/gse296761/sc_scores.png|pdf, sc_panel_heatmap.png|pdf
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/gse296761'))
import common  # noqa: E402
import panels as PN  # noqa: E402

OUT = ROOT / 'out/prelim/gse296761'
SC = OUT / 'sc'
NBOOT = 500
MINCELLS = 50
# 擬似バルクと GSEA を回す細胞種（依頼の細胞種 + 束ねた Lepr+ MSC と EC 全体）
GROUPS = {
    'Lepr+ MSC (all)': ['MSC-Adipo', 'MSC-Osteo'],
    'MSC-Adipo': ['MSC-Adipo'],
    'MSC-Osteo': ['MSC-Osteo'],
    'Osteoprogenitor': ['Osteoprogenitor'],
    'Osteoblast': ['Osteoblast'],
    'EC (all)': ['EC arteriolar', 'EC capillary', 'EC venular', 'EC sinusoidal'],
    'EC arteriolar': ['EC arteriolar'],
    'EC capillary': ['EC capillary'],
    'EC sinusoidal': ['EC sinusoidal'],
    'Pericyte': ['Pericyte'],
    'SMC': ['SMC'],
    'Fibroblast': ['Fibroblast'],
}
KEY = {k: 'sc296752_' + k.replace('+', '').replace(' (all)', '_all').replace(' ', '_').replace('-', '_') for k in GROUPS}
SCORE_SETS = ['HAND_ADIPO', 'HALLMARK_ADIPOGENESIS', 'HAND_ECM', 'NABA_CORE_MATRISOME', 'NABA_COLLAGENS', 'HAND_OSTEO',
              'GOBP_OSSIFICATION', 'HAND_IFNG', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_INTERFERON_ALPHA_RESPONSE',
              'HALLMARK_INFLAMMATORY_RESPONSE', 'SAUL_SEN_MAYO', 'HALLMARK_P53_PATHWAY', 'HALLMARK_ANGIOGENESIS',
              'HAND_HSC_NICHE']
PANEL_GENES = list(dict.fromkeys(sum(PN.PANELS.values(), []) + ['Cdkn1a', 'Cdkn2a', 'Bax', 'Mki67', 'Gzmk', 'Il6',
                                                               'Ccl2', 'Pgr', 'Greb1', 'Esr1']))


def cpm(v):
    return v / v.sum() * 1e6


def pseudobulk(a):
    X = a.layers['counts'].tocsr()
    genes = np.array(a.var_names)
    gi = {g: i for i, g in enumerate(genes)}
    pi = [gi[g] for g in PANEL_GENES if g in gi]
    rng = np.random.default_rng(0)
    rows, allg = [], {}
    for name, cts in GROUPS.items():
        idx = {s: np.where(a.obs.celltype.isin(cts).values & (a.obs['sample'] == s).values)[0] for s in ['PBS', 'DOX']}
        if min(len(v) for v in idx.values()) < MINCELLS:
            print('skip', name, {k: len(v) for k, v in idx.items()})
            continue
        sums = {s: np.asarray(X[v].sum(axis=0)).ravel() for s, v in idx.items()}
        c = {s: cpm(v) for s, v in sums.items()}
        lfc = np.log2((c['DOX'] + 1) / (c['PBS'] + 1))
        mean = (c['DOX'] + c['PBS']) / 2
        frac = {s: np.asarray((X[v] > 0).mean(axis=0)).ravel() for s, v in idx.items()}
        allg[name] = pd.DataFrame({'gene': genes, 'cpm_PBS': c['PBS'], 'cpm_DOX': c['DOX'], 'meanCPM': mean,
                                   'log2FC': lfc, 'frac_PBS': frac['PBS'], 'frac_DOX': frac['DOX']})
        ref = np.sort(lfc[mean >= 10])
        Xp = {s: X[v][:, pi].toarray() for s, v in idx.items()}
        lib = {s: np.asarray(X[v].sum(axis=1)).ravel() for s, v in idx.items()}
        boots = []
        for _ in range(NBOOT):
            cc = {}
            for s in ['PBS', 'DOX']:
                r = rng.integers(0, len(idx[s]), len(idx[s]))
                cc[s] = Xp[s][r].sum(axis=0) / lib[s][r].sum() * 1e6
            boots.append(np.log2((cc['DOX'] + 1) / (cc['PBS'] + 1)))
        boots = np.array(boots)
        for j, i in enumerate(pi):
            rows.append(dict(group=name, gene=genes[i], n_PBS=len(idx['PBS']), n_DOX=len(idx['DOX']),
                             cpm_PBS=c['PBS'][i], cpm_DOX=c['DOX'][i], log2FC=lfc[i],
                             boot_lo=np.percentile(boots[:, j], 2.5), boot_hi=np.percentile(boots[:, j], 97.5),
                             frac_PBS=frac['PBS'][i], frac_DOX=frac['DOX'][i],
                             pctile=(np.searchsorted(ref, lfc[i]) / len(ref) * 100 if mean[i] >= 10 else np.nan)))
    P = pd.DataFrame(rows)
    P.to_csv(SC / 'panel_pseudobulk.csv', index=False)
    A = pd.concat(allg, names=['group']).reset_index(level=0)
    A.to_csv(SC / 'pseudobulk_all_genes.csv.gz', index=False)
    return P, A


def gsea(A):
    base, _ = PN.base_collection()
    extra, _ = PN.extra_collection()
    rows = []
    for name in A.group.unique():
        d = A[(A.group == name) & (A.meanCPM >= 10)]
        rnk = d.groupby(d.gene.str.upper()).log2FC.mean().sort_values(ascending=False)
        for tag, sets in [('base', base), ('extra', extra)]:
            r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=4, max_size=1000,
                           seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
            r = r.assign(group=name, contrast=KEY[name], collection=tag, n_ranked=len(rnk),
                         tie_frac=1 - rnk.round(6).nunique() / len(rnk))
            fn = SC / f"gsea{'_extra' if tag == 'extra' else ''}_{KEY[name]}.csv"
            r.to_csv(fn, index=False)
            rows.append(r)
        print(name, len(rnk), 'genes ranked')
    L = pd.concat(rows, ignore_index=True).rename(columns={'Term': 'term', 'FDR q-val': 'fdr', 'NOM p-val': 'pval',
                                                           'Lead_genes': 'lead_genes'})
    for k in ['NES', 'ES', 'fdr', 'pval']:
        L[k] = L[k].astype(float)
    L = L[['group', 'contrast', 'collection', 'term', 'NES', 'fdr', 'pval', 'lead_genes', 'n_ranked', 'tie_frac']]
    L.to_csv(SC / 'sc_nes_long.csv', index=False)
    return L


def thin(a):
    """細胞種ごとに DOX 細胞の counts を二項間引きして UMI 中央値を PBS に揃え、正規化（1e4, log1p）し直したコピーを返す。"""
    import scipy.sparse as sp
    rng = np.random.default_rng(0)
    X = a.layers['counts'].tocsr().astype(np.float64).copy()
    rows = []
    for ct in a.obs.celltype.cat.categories:
        m = (a.obs.celltype == ct).values
        p_, d_ = m & (a.obs['sample'] == 'PBS').values, m & (a.obs['sample'] == 'DOX').values
        if p_.sum() == 0 or d_.sum() == 0:
            continue
        mp, md = np.median(a.obs.total_counts[p_]), np.median(a.obs.total_counts[d_])
        f = min(1.0, mp / md)
        rows.append(dict(celltype=ct, median_UMI_PBS=mp, median_UMI_DOX=md, thin_p=f))
        if f < 1:
            for i in np.where(d_)[0]:
                s0, s1 = X.indptr[i], X.indptr[i + 1]
                X.data[s0:s1] = rng.binomial(X.data[s0:s1].astype(np.int64), f)
    X.eliminate_zeros()
    pd.DataFrame(rows).to_csv(SC / 'thin_factors.csv', index=False)
    t = a.copy()
    t.X = sp.csr_matrix(X, dtype=np.float32)
    t.obs['total_counts_thin'] = np.asarray(X.sum(axis=1)).ravel()
    sc.pp.normalize_total(t, target_sum=1e4)
    sc.pp.log1p(t)
    return t


def scores(a):
    base, _ = PN.base_collection()
    extra, _ = PN.extra_collection()
    allsets = {**base, **extra}
    up = {g.upper(): g for g in a.var_names}
    rows = []
    for s in SCORE_SETS:
        gl = [up[g] for g in allsets[s] if g in up]
        sc.tl.score_genes(a, gl, score_name=s, random_state=0, use_raw=False)
    for ct in [c for c in a.obs.celltype.cat.categories if c not in ('Low quality', 'Doublet/unclear')]:
        m = (a.obs.celltype == ct).values
        p, d = m & (a.obs['sample'] == 'PBS').values, m & (a.obs['sample'] == 'DOX').values
        if min(p.sum(), d.sum()) < MINCELLS:
            continue
        for s in SCORE_SETS:
            x, y = a.obs.loc[d, s].values, a.obs.loc[p, s].values
            u = mannwhitneyu(x, y, alternative='two-sided')
            rows.append(dict(celltype=ct, set=s, n_PBS=p.sum(), n_DOX=d.sum(), mean_PBS=y.mean(), mean_DOX=x.mean(),
                             diff=x.mean() - y.mean(), sd_pooled=np.sqrt((x.var() + y.var()) / 2),
                             auc=u.statistic / (len(x) * len(y)), p_cells=u.pvalue))
    S = pd.DataFrame(rows)
    S['d'] = S['diff'] / S.sd_pooled
    S.to_csv(SC / 'scores_by_celltype.csv', index=False)
    return S


def cell_de(a):
    for name in ['MSC-Adipo', 'MSC-Osteo', 'Osteoblast', 'EC arteriolar', 'EC sinusoidal', 'Pericyte']:
        b = a[a.obs.celltype == name].copy()
        b.obs['sample'] = b.obs['sample'].astype(str)
        sc.tl.rank_genes_groups(b, 'sample', groups=['DOX'], reference='PBS', method='wilcoxon', pts=True)
        r = sc.get.rank_genes_groups_df(b, 'DOX')
        for s_ in ['PBS', 'DOX']:
            f = pd.Series(np.asarray((b[b.obs['sample'] == s_].X > 0).mean(axis=0)).ravel(), index=b.var_names)
            r[f'frac_{s_}'] = r.names.map(f).values
        r = r[(r.frac_PBS > 0.1) | (r.frac_DOX > 0.1)]
        r = pd.concat([r.head(100), r.tail(100)])
        r.to_csv(SC / f"cellde_{KEY.get(name, name).replace(' ', '_')}.csv", index=False)


def bulk_concordance(A):
    out = []
    pairs = [('Lepr+ MSC (all)', 'out/prelim/chemo/de/B3_Dox4x_MSPC.csv', 'GSE296758 LepR+ MSPC Dox×4 (bulk, n=4/4)'),
             ('MSC-Adipo', 'out/prelim/chemo/de/B3_Dox4x_MSPC.csv', 'GSE296758 LepR+ MSPC Dox×4 (bulk, n=4/4)'),
             ('MSC-Osteo', 'out/prelim/chemo/de/B3_Dox4x_MSPC.csv', 'GSE296758 LepR+ MSPC Dox×4 (bulk, n=4/4)'),
             ('Osteoblast', 'out/prelim/de/B2_Dox4x_OB.csv', 'GSE296757 Col2.3+ OB Dox×4 (bulk, n=5/3)'),
             ('Lepr+ MSC (all)', 'out/prelim/de/B2_Dox4x_OB.csv', 'GSE296757 Col2.3+ OB Dox×4 (bulk, n=5/3)'),
             ('Osteoblast', 'out/prelim/chemo/de/B3_Dox4x_MSPC.csv', 'GSE296758 LepR+ MSPC Dox×4 (bulk, n=4/4)')]
    for g, f, lab in pairs:
        b = pd.read_csv(ROOT / f).dropna(subset=['log2FoldChange']).drop_duplicates('gene').set_index('gene')
        s = A[(A.group == g) & (A.meanCPM >= 10)].set_index('gene')
        j = s.join(b, how='inner')
        r_all = spearmanr(j.log2FC, j.log2FoldChange).statistic
        sig = j[(j.padj < 0.05)]
        r_sig = spearmanr(sig.log2FC, sig.log2FoldChange).statistic if len(sig) > 5 else np.nan
        agree = (np.sign(sig.log2FC) == np.sign(sig.log2FoldChange)).mean() if len(sig) else np.nan
        out.append(dict(sc_group=g, bulk=lab, n_genes=len(j), spearman_all=r_all, n_bulk_padj05=len(sig),
                        spearman_bulk_sig=r_sig, sign_agreement_bulk_sig=agree))
    C = pd.DataFrame(out)
    C.to_csv(SC / 'bulk_concordance.csv', index=False)
    print(C.round(3).to_string())


def plot_scores(S):
    common.setup()
    cts = ['MSC-Adipo', 'MSC-Osteo', 'Osteoprogenitor', 'Osteoblast', 'EC arteriolar', 'EC capillary', 'EC venular',
           'EC sinusoidal', 'Pericyte', 'SMC', 'Fibroblast']
    lab = PN.all_labels()
    M = S.pivot(index='set', columns='celltype', values='d').reindex(index=SCORE_SETS, columns=cts)
    U = S.pivot(index='set', columns='celltype', values='auc').reindex(index=SCORE_SETS, columns=cts)
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(95)))
    im = ax.imshow(M.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.iloc[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, f'{v:+.2f}'.replace('-', '−'), ha='center', va='center', fontsize=5.5,
                    color='white' if abs(v) > 0.6 else 'black',
                    fontweight='bold' if abs(U.iloc[i, j] - 0.5) >= 0.1 else 'normal')
    ax.set_yticks(range(len(SCORE_SETS)), [lab.get(s, s) for s in SCORE_SETS], fontsize=6.5)
    n = S.drop_duplicates('celltype').set_index('celltype')
    ax.set_xticks(range(len(cts)), [f'{c}\n({n.loc[c, "n_PBS"]}/{n.loc[c, "n_DOX"]})' if c in n.index else c
                                    for c in cts], rotation=50, ha='right', fontsize=6)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.01)
    cb.set_label("Cohen's d (DOX − PBS cells)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    ax.set_title('GSE296752 per-cell gene-set scores, DOX vs PBS (n = 1 pooled library each; cells as units = '
                 'exploratory).\nBold: |AUC − 0.5| ≥ 0.1. Column labels: cells PBS/DOX.', fontsize=6.5, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'sc_scores.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)


def plot_panel(P):
    common.setup()
    groups = ['MSC-Adipo', 'MSC-Osteo', 'Osteoprogenitor', 'Osteoblast', 'EC arteriolar', 'EC capillary',
              'EC sinusoidal', 'Pericyte']
    genes, sep = [], []
    for k in ['adipo', 'ecm', 'osteo', 'oc_support', 'hsc_niche', 'ifng']:
        gl = [g for g in PN.PANELS[k] if g in set(P.gene)]
        genes += gl
        sep.append(len(genes))
    M = P.pivot(index='gene', columns='group', values='log2FC').reindex(index=genes, columns=groups)
    lo = P.pivot(index='gene', columns='group', values='cpm_PBS').reindex(index=genes, columns=groups)
    hi = P.pivot(index='gene', columns='group', values='cpm_DOX').reindex(index=genes, columns=groups)
    low = ((lo + hi) / 2) < 10
    M = M.mask(low)
    fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(200)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(M.values), cmap=cmap, vmin=-2, vmax=2, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.iloc[i, j]
            if not np.isnan(v):
                ax.text(j, i, f'{v:+.1f}'.replace('-', '−'), ha='center', va='center', fontsize=4.8,
                        color='white' if abs(v) > 1.4 else 'black')
    for s in sep[:-1]:
        ax.axhline(s - 0.5, c='k', lw=0.6)
    ax.set_yticks(range(len(genes)), genes, fontsize=5.5, fontstyle='italic')
    ax.set_xticks(range(len(groups)), groups, rotation=50, ha='right', fontsize=6)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, shrink=0.3, pad=0.02)
    cb.set_label('log$_2$ DOX/PBS (pseudobulk CPM)', fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    ax.set_title('GSE296752 pseudobulk, n = 1 vs 1\n(grey: mean CPM < 10)', fontsize=6.5, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'sc_panel_heatmap.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)


def main():
    a = sc.read_h5ad(SC / 'adata_annot.h5ad')
    a.obs['sample'] = a.obs['sample'].astype(str)
    P, A = pseudobulk(a)
    plot_panel(P)
    bulk_concordance(A)
    L = gsea(A)
    pd.set_option('display.width', 250)
    show = ['HALLMARK_ADIPOGENESIS', 'GOBP_FAT_CELL_DIFFERENTIATION', 'NABA_CORE_MATRISOME', 'GOBP_OSSIFICATION',
            'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_ANGIOGENESIS', 'HAND_IFNG',
            'HAND_ADIPO', 'HAND_OSTEO', 'NABA_COLLAGENS']
    print(L[L.term.isin(show)].pivot(index='term', columns='group', values='NES').round(2).to_string())
    t = thin(a)
    S = scores(t)
    plot_scores(S)
    print(S.pivot(index='set', columns='celltype', values='d').round(2).to_string())
    cell_de(t)


if __name__ == '__main__':
    main()

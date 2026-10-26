"""4C 案1（GSE304562 コルチコステロン MC3T3）と案2（GSE276391 DEX 脛骨）を、互いに・4A/4B と比べる。

  uv run python -I prelim/fig4c_v4/compare.py      # repo root から（de.py, gsea.py の後）

既存の DE / GSEA 表は読むだけ（再計算しない）。出力はすべて out/prelim/fig4c_v4/ の下:
  gc_targets.csv            GR（グルココルチコイド受容体）標的遺伝子の log2FC/padj（各ステロイド列）＝ステロイドが効いているかの陽性対照
  gene_concordance.csv      案1 / 案2 と各 contrast の遺伝子単位の一致（Wald 統計量の Spearman、共通 DEG の同方向/逆方向と超幾何 P）
  module_concordance.csv    同じ 72 set の NES ベクトルの Spearman と、両方 FDR<0.05 の同方向/逆方向の set 数
  gsea_noGR.csv             感度分析: GR 標的（下の GR_TARGETS、事前定義）を順位から抜いて要約図の 12 行を prerank し直した NES/FDR
  dex_environment.csv       案2 の脛骨（骨髄込み）を区画マーカーで読む: サンプル別スコア（平均 log2 CPM）の DEX − 対照と Welch t
  dex_environment_genes.csv 同じマーカーの遺伝子単位 log2FC/padj（C_DEX, C_DEX_bgadj）と、H_GIOP（大文字化）
  dex_inflam_vs_cells.csv   脛骨 10 サンプルで、炎症パネル遺伝子（S100a8/9, Mmp9, Tnf, Il1b）と区画スコアの Spearman
  panel_genes.csv           注目遺伝子パネル（prelim/genes.py の INFLAM + APOP、chemo2 の手製 niche パネル）の log2FC/padj
  summary_preview.png|pdf   要約図の 8 行 × [4A 4 列 | 4B 3 列 | 4C 案1・案2] のプレビュー（要約図本体には触らない）
  dex_environment.png|pdf   区画スコアの図（案2）
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[2]
for p in ('prelim', 'prelim/chemo2', 'prelim/invitro', 'prelim/summary'):
    sys.path.insert(0, str(ROOT / p))
import common  # noqa: E402
import genesets as GS  # noqa: E402
import summary_modules as SM  # noqa: E402   ROWS / DROPPED / SIG / TREND / 列ラベルだけ借りる

P = ROOT / 'out/prelim'
OUT = P / 'fig4c_v4'

# ---------- DE 表（読むだけ） ----------
DE = {
    'Cort_MC3T3': OUT / 'de/Cort_MC3T3.csv',
    'C_DEX': P / 'de/C_DEX.csv',
    'C_DEX_bgadj': P / 'de/C_DEX_bgadj.csv',
    'H_GIOP': P / 'human/de/H_GIOP.csv',
    'H_GIOP_marrowadj': P / 'human/de/H_GIOP_marrowadj.csv',
    'B2_Dox4x_OB': P / 'de/B2_Dox4x_OB.csv',
    'B3_Dox4x_MSPC': P / 'chemo/de/B3_Dox4x_MSPC.csv',
    'B_Dox': P / 'de/B_Dox.csv',
    'B_Dox_compadj': P / 'chemo2/de/B_Dox_compadj.csv',
    'NHO_Eto': P / 'invitro/de/NHO_Eto.csv',
    'Age_2.5y': P / 'de/Age_2.5y.csv',
    'sc210584_CAR': P / 'aging_f/de/sc210584_CAR.csv',
    'H_PostMeno': P / 'human/de/H_PostMeno.csv',
    'H_Old': P / 'human/de/H_Old.csv',
}
NAME = {
    'Cort_MC3T3': '4C-1 Cort MC3T3 (in vitro)', 'C_DEX': '4C-2 DEX tibia ♂', 'C_DEX_bgadj': '4C-2 DEX tibia (bg adj.)',
    'H_GIOP': 'Hs GIOP cortical ♀', 'H_GIOP_marrowadj': 'Hs GIOP (marrow adj.)',
    'B2_Dox4x_OB': '4B Dox×4 Col2.3+ OB ♀', 'B3_Dox4x_MSPC': '4B Dox×4 LepR+ MSPC ♀', 'B_Dox': '4B Dox×1 bone ♂',
    'B_Dox_compadj': '4B Dox×1 bone (comp. adj.)', 'NHO_Eto': '4B Hs NHO etoposide',
    'Age_2.5y': '4A aging bone ♂', 'sc210584_CAR': '4A mid-age CAR ♀', 'H_PostMeno': '4A Hs postmeno.', 'H_Old': '4A Hs old',
}
HUMAN = {'H_GIOP', 'H_GIOP_marrowadj', 'NHO_Eto', 'H_PostMeno', 'H_Old'}

# GR 標的（事前定義）: GC 応答の陽性対照として文献で繰り返し使われる直接標的。解析結果を見て足し引きしない。
GR_TARGETS = ['Fkbp5', 'Tsc22d3', 'Per1', 'Sgk1', 'Dusp1', 'Klf9', 'Klf15', 'Zbtb16', 'Nfkbia', 'Zfp36',
              'Errfi1', 'Ddit4', 'Pdk4', 'Mt1', 'Mt2', 'Glul', 'Cdkn1c', 'Sult1a1', 'Angptl4', 'Bcl2l11']

# 骨（骨髄込み）の区画マーカー（事前定義）。炎症パネルの遺伝子（S100a8/9, Mmp9, Tnf, Il1b）は区画マーカーに入れない（循環を避ける）
COMPART = {
    'GC response': ['Fkbp5', 'Tsc22d3', 'Per1', 'Sgk1', 'Zbtb16', 'Klf9'],
    'Osteoblast': ['Runx2', 'Sp7', 'Bglap', 'Bglap2', 'Ibsp', 'Col1a1', 'Alpl'],
    'Osteocyte': ['Sost', 'Dmp1', 'Phex', 'Mepe'],
    'Osteoclast': ['Acp5', 'Ctsk', 'Dcstamp', 'Oscar', 'Calcr'],
    'LepR+ stroma': ['Lepr', 'Cxcl12', 'Kitl', 'Ebf3'],
    'Adipocyte': ['Adipoq', 'Cfd', 'Lpl', 'Plin1', 'Cidec'],
    'Endothelial': ['Cdh5', 'Emcn', 'Pecam1', 'Kdr'],
    'Neutrophil': ['Ngp', 'Camp', 'Ltf', 'Lcn2', 'Chil3', 'Retnlg', 'Mpo', 'Elane'],
    'Mono/macrophage': ['Csf1r', 'Adgre1', 'Cd68', 'Lyz2', 'C1qa'],
    'B cell': ['Cd79a', 'Cd79b', 'Cd19', 'Ms4a1', 'Vpreb3', 'Rag1', 'Dntt'],
    'T cell': ['Cd3e', 'Cd3g', 'Cd8a', 'Cd4'],
    'Erythroid': ['Hba-a1', 'Hbb-bs', 'Gypa', 'Klf1', 'Alas2'],
    'Megakaryocyte': ['Pf4', 'Ppbp', 'Itga2b'],
    'Muscle': ['Acta1', 'Ckm', 'Myh4'],
}
INFLAM_Q = ['S100a8', 'S100a9', 'Mmp9', 'Tnf', 'Il1b']


def load_de(k):
    d = pd.read_csv(DE[k]).dropna(subset=['gene'])
    d['gene'] = d.gene.astype(str)
    d = d.drop_duplicates('gene')
    d['key'] = d.gene.str.upper()
    d = d.drop_duplicates('key').set_index('key')
    d['stat'] = d.log2FoldChange / d.lfcSE
    return d


# ---------- 1. GR 標的 ----------
def gc_targets(D):
    cols = ['Cort_MC3T3', 'C_DEX', 'C_DEX_bgadj', 'H_GIOP']
    rows = []
    for g in GR_TARGETS:
        r = dict(gene=g)
        for c in cols:
            d = D[c]
            k = g.upper()
            r[f'{c}|log2FC'] = d.log2FoldChange.get(k, np.nan)
            r[f'{c}|padj'] = d.padj.get(k, np.nan)
        rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(OUT / 'gc_targets.csv', index=False)
    s = {}
    for c in cols:
        x = T[[f'{c}|log2FC', f'{c}|padj']].dropna()
        s[c] = dict(tested=len(x), up_sig=int(((x.iloc[:, 0] > 0) & (x.iloc[:, 1] < 0.05)).sum()),
                    down_sig=int(((x.iloc[:, 0] < 0) & (x.iloc[:, 1] < 0.05)).sum()),
                    median_lfc=x.iloc[:, 0].median())
    print('GR targets:\n', pd.DataFrame(s).T.round(2).to_string())
    return T


def panel_genes(D):
    sys.path.insert(0, str(ROOT / 'prelim'))
    from genes import INFLAM, APOP
    order = list(dict.fromkeys(INFLAM + APOP + [g for v in GS.NICHE_PANEL.values() for g in v]))
    cols = ['Cort_MC3T3', 'C_DEX', 'C_DEX_bgadj', 'H_GIOP', 'B2_Dox4x_OB', 'B3_Dox4x_MSPC']
    tr = {'Fasl': 'FASLG', 'Trp53': 'TP53'}
    rows = []
    for g in order:
        panel = 'inflam' if g in INFLAM else 'apop' if g in APOP else 'niche'
        r = dict(gene=g, panel=panel)
        for c in cols:
            k = tr.get(g, g).upper() if c.startswith('H_') else g.upper()
            r[f'{c}|log2FC'] = D[c].log2FoldChange.get(k, np.nan)
            r[f'{c}|padj'] = D[c].padj.get(k, np.nan)
        rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(OUT / 'panel_genes.csv', index=False)
    v = T[['gene', 'panel']].copy()
    for c in cols:
        v[c] = [('' if np.isnan(a) else f'{a:+.2f}' + ('**' if p < 0.05 else '*' if p < 0.1 else ''))
                for a, p in zip(T[f'{c}|log2FC'], T[f'{c}|padj'].fillna(1))]
    print(v.to_string(index=False))
    return T


# ---------- 2. 遺伝子単位 ----------
def hyper(a, b, n_all, n_ov):
    return stats.hypergeom.sf(n_ov - 1, n_all, a, b) if n_ov > 0 else 1.0


def gene_concordance(D):
    rows = []
    for ref in ['Cort_MC3T3', 'C_DEX']:
        for k in DE:
            if k == ref:
                continue
            a, b = D[ref], D[k]
            com = a.index.intersection(b.index)
            com = [g for g in com if np.isfinite(a.stat.get(g, np.nan)) and np.isfinite(b.stat.get(g, np.nan))]
            x, y = a.loc[com], b.loc[com]
            rho = stats.spearmanr(x.stat, y.stat).statistic
            rho_lfc = stats.spearmanr(x.log2FoldChange, y.log2FoldChange).statistic
            r = dict(ref=ref, other=k, other_name=NAME[k], species='human' if k in HUMAN else 'mouse',
                     n_common=len(com), spearman_stat=rho, spearman_lfc=rho_lfc)
            for tag, pa, lf in [('strict', 0.05, 1.0), ('lenient', 0.1, 0.585)]:
                sx = (x.padj < pa) & (x.log2FoldChange.abs() >= lf)
                sy = (y.padj < pa) & (y.log2FoldChange.abs() >= lf)
                for dname, dx, dy in [('up', 1, 1), ('down', -1, -1), ('opp_xup_ydown', 1, -1), ('opp_xdown_yup', -1, 1)]:
                    mx = sx & (np.sign(x.log2FoldChange) == dx)
                    my = sy & (np.sign(y.log2FoldChange) == dy)
                    ov = int((mx & my).sum())
                    r[f'{tag}_{dname}'] = ov
                    if dname in ('up', 'down'):
                        r[f'{tag}_{dname}_P'] = hyper(int(mx.sum()), int(my.sum()), len(com), ov)
                r[f'{tag}_same'] = r[f'{tag}_up'] + r[f'{tag}_down']
                r[f'{tag}_opposite'] = r[f'{tag}_opp_xup_ydown'] + r[f'{tag}_opp_xdown_yup']
                if tag == 'lenient':
                    both = sx & sy
                    r['lenient_shared_genes_same'] = ','.join(
                        sorted(x.index[both & (np.sign(x.log2FoldChange) == np.sign(y.log2FoldChange))],
                               key=lambda g: -abs(x.stat[g]))[:25])
            rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(OUT / 'gene_concordance.csv', index=False)
    pd.set_option('display.width', 250)
    print(T[['ref', 'other', 'n_common', 'spearman_stat', 'spearman_lfc', 'strict_same', 'strict_opposite',
             'lenient_same', 'lenient_opposite', 'lenient_up_P', 'lenient_down_P']].round(3).to_string(index=False))
    return T


# ---------- 3. モジュール単位 ----------
def nes_table():
    cols = {}
    w = pd.read_csv(P / 'chemo2/nes_wide.csv', header=[0, 1], index_col=0)
    for k in ['Age_2.5y', 'B_Dox', 'B_Dox_compadj', 'B3_Dox4x_MSPC', 'B2_Dox4x_OB', 'C_DEX', 'C_DEX_bgadj']:
        cols[k] = (w[('NES', k)].astype(float), w[('FDR', k)].astype(float))
    w = pd.read_csv(P / 'aging_f/nes_wide_all.csv', header=[0, 1], index_col=0)
    cols['sc210584_CAR'] = (w[('NES', 'sc210584_CAR')].astype(float), w[('FDR', 'sc210584_CAR')].astype(float))
    w = pd.read_csv(P / 'invitro/nes_wide.csv', header=[0, 1], index_col=0)
    cols['NHO_Eto'] = (w[('NES', 'NHO_Eto')].astype(float), w[('FDR', 'NHO_Eto')].astype(float))
    for k, f in [('H_PostMeno', P / 'summary/gsea_human/H_PostMeno.csv'), ('H_Old', P / 'summary/gsea_human/H_Old.csv'),
                 ('H_GIOP', OUT / 'gsea_human/H_GIOP.csv'), ('H_GIOP_marrowadj', OUT / 'gsea_human/H_GIOP_marrowadj.csv'),
                 ('Cort_MC3T3', OUT / 'gsea/Cort_MC3T3.csv')]:
        r = pd.read_csv(f).set_index('Term')
        cols[k] = (r.NES.astype(float), r['FDR q-val'].astype(float))
    N = pd.DataFrame({k: v[0] for k, v in cols.items()})
    F = pd.DataFrame({k: v[1] for k, v in cols.items()})
    _, meta = GS.collection()
    terms = [t for t in meta.term if t in N.index]
    return N.loc[terms], F.loc[terms]


def module_concordance(N, F):
    rows = []
    for ref in ['Cort_MC3T3', 'C_DEX']:
        for k in N.columns:
            if k == ref:
                continue
            m = N[[ref, k]].dropna()
            sa, sb = F.loc[m.index, ref] < 0.05, F.loc[m.index, k] < 0.05
            same = (sa & sb & (np.sign(m[ref]) == np.sign(m[k])))
            opp = (sa & sb & (np.sign(m[ref]) != np.sign(m[k])))
            rows.append(dict(ref=ref, other=k, other_name=NAME.get(k, k), n_sets=len(m),
                             spearman_NES=stats.spearmanr(m[ref], m[k]).statistic,
                             both_sig_same=int(same.sum()), both_sig_opposite=int(opp.sum()),
                             ref_only=int((sa & ~sb).sum()), other_only=int((~sa & sb).sum()),
                             same_sets=';'.join(m.index[same]), opposite_sets=';'.join(m.index[opp])))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / 'module_concordance.csv', index=False)
    print(T[['ref', 'other', 'n_sets', 'spearman_NES', 'both_sig_same', 'both_sig_opposite', 'ref_only',
             'other_only']].round(2).to_string(index=False))
    return T


# ---------- 4. GR 標的を抜いた感度分析 ----------
def gsea_nogr():
    terms = [t for t, *_ in SM.ROWS + SM.DROPPED]
    sets, _ = GS.collection()
    sub = {t: sets[t] for t in terms}
    drop = {g.upper() for g in GR_TARGETS}
    rows = []
    for k in ['Cort_MC3T3', 'C_DEX']:
        d = pd.read_csv(DE[k]).dropna(subset=['lfcSE', 'gene'])
        d = d[d.lfcSE > 0]
        rnk = (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean()
        for tag, rr in [('all', rnk), ('noGR', rnk[~rnk.index.isin(drop)])]:
            # 12 set だけで回すと FDR の母集団が変わるので、FDR は全 collection で回して取り出す
            r = gp.prerank(rnk=rr.sort_values(ascending=False), gene_sets=sets, permutation_num=1000, min_size=4,
                           max_size=1000, seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
            r = r.set_index('Term')
            for t in terms:
                rows.append(dict(contrast=k, ranking=tag, term=t, NES=float(r.loc[t, 'NES']),
                                 FDR=float(r.loc[t, 'FDR q-val']), n_removed=int(rnk.index.isin(drop).sum())))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / 'gsea_noGR.csv', index=False)
    v = T.assign(cell=T.NES.map('{:+.2f}'.format) + ' (' + T.FDR.map('{:.3f}'.format) + ')')
    print(v.pivot_table(index='term', columns=['contrast', 'ranking'], values='cell', aggfunc='first').to_string())
    return T


# ---------- 5. 案2 の脛骨を区画で読む ----------
def dex_environment():
    seq = ROOT / 'data/seq/GSE276391'
    x = pd.read_csv(seq / 'counts.tsv.gz', sep='\t', index_col=0)
    s = pd.read_csv(seq / 'samples.csv', dtype=str).set_index('sample')
    s = s[s.group.isin(['DEX_Ctrl', 'Vehicle_Ctrl'])]
    x = x[s.index]
    lc = np.log2(x / x.sum() * 1e6 + 1)
    rows, score = [], {}
    for cat, gs in COMPART.items():
        g = [z for z in gs if z in lc.index]
        sc = lc.loc[g].mean()
        score[cat] = sc
        a, b = sc[s.group == 'DEX_Ctrl'], sc[s.group == 'Vehicle_Ctrl']
        t = stats.ttest_ind(a, b, equal_var=False)
        rows.append(dict(compartment=cat, genes=','.join(g), n_genes=len(g), mean_log2cpm_vehicle=b.mean(),
                         mean_log2cpm_DEX=a.mean(), diff_DEX_minus_vehicle=a.mean() - b.mean(), welch_t=t.statistic,
                         welch_p=t.pvalue))
    T = pd.DataFrame(rows)
    T['welch_p_BH'] = stats.false_discovery_control(T.welch_p)
    T.to_csv(OUT / 'dex_environment.csv', index=False)
    print(T.drop(columns='genes').round(3).to_string(index=False))
    S = pd.DataFrame(score)
    S['group'] = s.group
    S['bg'] = s.bg.astype(float)
    S.to_csv(OUT / 'dex_environment_scores.csv')

    # 遺伝子単位（既存 DE 表）
    D = {k: load_de(k) for k in ['C_DEX', 'C_DEX_bgadj', 'H_GIOP', 'Cort_MC3T3']}
    g_rows = []
    for cat, gs in list(COMPART.items()) + [('Inflam. panel (question)', INFLAM_Q), ('OPG/RANKL', ['Tnfrsf11b', 'Tnfsf11'])]:
        for g in gs:
            r = dict(compartment=cat, gene=g)
            for k, d in D.items():
                kk = g.upper() if k.startswith('H_') else g
                kk = kk.upper()
                r[f'{k}|log2FC'] = d.log2FoldChange.get(kk, np.nan)
                r[f'{k}|padj'] = d.padj.get(kk, np.nan)
            g_rows.append(r)
    G = pd.DataFrame(g_rows)
    G.to_csv(OUT / 'dex_environment_genes.csv', index=False)

    # 炎症パネル遺伝子と区画スコアの相関（10 サンプル）
    c_rows = []
    for g in INFLAM_Q + ['Tnfrsf11b']:
        if g not in lc.index:
            continue
        for cat in COMPART:
            rho, p = stats.spearmanr(lc.loc[g], S[cat])
            c_rows.append(dict(gene=g, compartment=cat, spearman=rho, p=p))
    C = pd.DataFrame(c_rows)
    C.to_csv(OUT / 'dex_inflam_vs_cells.csv', index=False)
    print(C.pivot(index='compartment', columns='gene', values='spearman').round(2).to_string())
    plot_env(T, S)
    return T, G, C


def plot_env(T, S):
    common.setup()
    cats = list(COMPART)
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(62)))
    for i, c in enumerate(cats):
        for grp, col, dx in [('Vehicle_Ctrl', common.COLOR_LRO, -0.17), ('DEX_Ctrl', common.COLOR_HRO, 0.17)]:
            v = S.loc[S.group == grp, c] - S.loc[S.group == 'Vehicle_Ctrl', c].mean()
            ax.bar(i + dx, v.mean(), width=0.3, facecolor='none', edgecolor=col, lw=0.8)
            ax.scatter(np.full(len(v), i + dx) + np.linspace(-0.06, 0.06, len(v)), v, s=6, color=col, zorder=3,
                       label=('Vehicle' if grp == 'Vehicle_Ctrl' else 'DEX') if i == 0 else None)
        r = T.set_index('compartment').loc[c]
        star = '**' if r.welch_p_BH < 0.01 else '*' if r.welch_p_BH < 0.05 else ''
        if star:
            top = max((S[c] - S.loc[S.group == 'Vehicle_Ctrl', c].mean()).max(), 0)
            ax.text(i, top + 0.08, star, ha='center', va='bottom', fontsize=8)
    ax.axhline(0, c='k', lw=0.5)
    ax.set_xticks(range(len(cats)), cats, rotation=35, ha='right', fontsize=7)
    ax.set_ylabel('Marker score\n(Δ mean log$_2$ CPM vs vehicle)', fontsize=8)
    ax.set_title('GSE276391 tibia (with marrow), DEX 5 d vs vehicle, n = 5 vs 5; * BH-adjusted Welch P < 0.05, ** < 0.01',
                 fontsize=7.5, loc='left')
    ax.legend(frameon=False, fontsize=7, loc='lower left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'dex_environment.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)


# ---------- 6. 要約図のプレビュー（4A | 4B | 4C） ----------
def summary_preview():
    base = pd.read_csv(P / 'summary/summary_modules_values.csv')
    new = pd.read_csv(OUT / 'summary_modules_cols.csv')
    keep_new = ['Cort_MC3T3', 'C_DEX']
    T = pd.concat([base[base.shown], new[new.shown & new.contrast.isin(keep_new)]], ignore_index=True)
    cols = [(k, lab, g) for k, _, _, lab, g in SM.COLS] + [
        ('Cort_MC3T3', 'Mouse\nMC3T3-E1 OB\ncorticosterone\n24 h (in vitro)', 'steroid'),
        ('C_DEX', 'Mouse ♂\ntibia\nDEX 5 d', 'steroid')]
    rows = [r[1] for r in SM.ROWS]
    common.setup()
    nr, nc = len(rows), len(cols)
    fig = plt.figure(figsize=(common.mm(170), common.mm(85)))
    L, R, B, Tp = 0.2, 0.83, 0.035, 0.62
    ax = fig.add_axes([L, B, R - L, Tp - B])
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(nr - 0.5, -0.5)
    ax.axis('off')
    cmap, norm, pad = plt.get_cmap('RdBu_r'), Normalize(-3, 3), 0.06
    for _, x in T.iterrows():
        i = rows.index(x.module)
        j = [c[0] for c in cols].index(x.contrast)
        box, s_ = (j - 0.5 + pad, i - 0.5 + pad), 1 - 2 * pad
        if x.call in ('up', 'down'):
            ax.add_patch(Rectangle(box, s_, s_, fc=cmap(norm(x.NES)), ec='none'))
            ax.text(j, i, f'{x.NES:+.1f}'.replace('-', '−'), ha='center', va='center', fontsize=6,
                    color='white' if abs(x.NES) >= 2.0 else '#222222')
        elif str(x.call).startswith('trend'):
            ax.add_patch(Rectangle(box, s_, s_, fc='none', ec=cmap(norm(1.4 if x.NES > 0 else -1.4)), lw=0.45))
    rg = [r[2] for r in SM.ROWS]
    for i in range(1, nr):
        if rg[i] != rg[i - 1]:
            ax.plot([-0.5, nc - 0.5], [i - 0.5, i - 0.5], c='#888888', lw=0.4, clip_on=False)
    grp = [c[2] for c in cols]
    bounds = [j for j in range(1, nc) if grp[j] != grp[j - 1]]
    for jb in bounds:
        ax.plot([jb - 0.5, jb - 0.5], [-0.5, nr - 0.5], c='black', lw=0.8, clip_on=False)
    for i, rl in enumerate(rows):
        ax.text(-0.62, i, rl, ha='right', va='center', fontsize=8)
    for j, c in enumerate(cols):
        ax.text(j, -0.62, c[1], ha='center', va='bottom', fontsize=5.3, linespacing=1.05)
    edges = [0] + bounds + [nc]
    labs = {'aging': '4A  Menopause / aging', 'chemo': '4B  Chemotherapy', 'steroid': '4C  Glucocorticoid'}
    for a, b in zip(edges[:-1], edges[1:]):
        x0, x1, y = a - 0.42, b - 1 + 0.42, -3.4
        ax.plot([x0, x0, x1, x1], [y + 0.18, y, y, y + 0.18], c='black', lw=0.7, clip_on=False)
        ax.text((x0 + x1) / 2, y - 0.12, labs[grp[a]], ha='center', va='bottom', fontsize=8, fontweight='bold')
    lx = fig.add_axes([0.86, 0.36, 0.02, 0.24])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=lx)
    cb.set_ticks([-3, 0, 3])
    cb.ax.tick_params(labelsize=6, length=2, width=0.5)
    cb.outline.set_linewidth(0.4)
    lx.set_title('NES\n(FDR < 0.05)', fontsize=6.5, loc='left', pad=3)
    fig.text(0.85, 0.12, 'Open box:\nFDR 0.05–0.25\nBlank: ≥ 0.25\n\nPreview only\n(not the summary\nfigure itself)', fontsize=6)
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'summary_preview.{ext}', dpi=common.DPI)
    plt.close(fig)


def main():
    D = {k: load_de(k) for k in DE}
    gc_targets(D)
    panel_genes(D)
    gene_concordance(D)
    N, F = nes_table()
    module_concordance(N, F)
    gsea_nogr()
    dex_environment()
    summary_preview()


if __name__ == '__main__':
    main()

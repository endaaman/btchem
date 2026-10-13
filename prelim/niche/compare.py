"""LepR+ / CAR 間質で「加齢で動く軸」と「化療で動く軸」をモジュール単位・遺伝子単位で並べる。

  uv run python -I prelim/niche/compare.py     # repo root から。先に design.py → de.py

contrast（すべて log2FC = test / ref）:
  加齢   Age_LepRPa（GSE130885 LepR+PDGFRa+, 24-25 mo vs 8-12 wk, 3/3）, Age_CAR（GSE133922 CAR, 2 y vs 2 mo, 4/4）
         Age_LepRPaN（GSE130885 LepR+PDGFRa-, 感度）
  炎症刺激（参考）  LPS_CAR, PIC_CAR（GSE133922, 2 mo）
  化療   B3_Dox4x_MSPC, B3_Dox1x_MSPC（GSE296758, 雌 LepR+ MSPC, 8 週後。既存 out/prelim/chemo/de/ を読むだけ）
gene set: prelim/chemo2/genesets.py の collection()（Hallmark 50 + niche 17 + 手製 5）をそのまま使う。
GSEA の設定も chemo2/gsea.py と同じ（Wald 統計量で順位、gseapy prerank 1000 perm, seed 0, min 4, max 1000、FDR は 72 set 全体）。

分類（モジュール）: 加齢の向き = 主 2 本（Age_LepRPa, Age_CAR）のどちらかで FDR<0.05、かつもう片方の NES が逆符号で FDR<0.05 でない。
  2 本とも FDR<0.05 なら 'both'。逆向きに有意なら 'discordant'。
  化療の向き = Dox×4 か Dox×1 で FDR<0.05、かつもう片方が逆向きに有意でない。
これは向きの並置で、交互作用の検定ではない（加齢 × 化療を同じ実験で組んだデータは無い）。
DEG: |log2FC|>=1 & padj<0.05（prelim/de.py と同じ）。重なりは両 contrast で検定された遺伝子を母集団にした超幾何検定。

出力（out/prelim/niche/）: gsea/<label>.csv, nes_long.csv, nes_wide.csv, heatmap_niche.png|pdf, heatmap_hallmark.png|pdf,
  concordance_sets.csv, concordance_counts.csv, spearman_nes.csv, spearman_gene_stat.csv, module_class.csv,
  panel_genes.csv, deg_overlap.csv, deg_overlap_genes.csv, fig_niche.png|pdf|tiff
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
import common  # noqa: E402
import genesets as GS  # noqa: E402

OUT = ROOT / 'out/prelim/niche'
DEN = OUT / 'de'
DEC = ROOT / 'out/prelim/chemo/de'

# (label, 図の名前, 軸, 入力)
CONTRASTS = [
    ('Age_LepRPa', 'Aging LepR+PDGFRa+\n24mo vs 2-3mo (GSE130885)', 'aging', DEN / 'Age_LepRPa.csv'),
    ('Age_CAR', 'Aging CAR\n2y vs 2mo (GSE133922)', 'aging', DEN / 'Age_CAR.csv'),
    ('Age_LepRPaN', 'Aging LepR+PDGFRa-\n(sensitivity)', 'aging-sens', DEN / 'Age_LepRPaN.csv'),
    ('LPS_CAR', 'LPS CAR\nvs 2mo (ref.)', 'inflam', DEN / 'LPS_CAR.csv'),
    ('PIC_CAR', 'poly I:C CAR\nvs 2mo (ref.)', 'inflam', DEN / 'PIC_CAR.csv'),
    ('B3_Dox4x_MSPC', 'Dox 4x LepR+ MSPC\n8wk (GSE296758)', 'chemo', DEC / 'B3_Dox4x_MSPC.csv'),
    ('B3_Dox1x_MSPC', 'Dox 1x LepR+ MSPC\n8wk (GSE296758)', 'chemo', DEC / 'B3_Dox1x_MSPC.csv'),
]
NAME = {a: b for a, b, *_ in CONTRASTS}
SHORT = {'Age_LepRPa': 'Aging LepR+', 'Age_CAR': 'Aging CAR', 'Age_LepRPaN': 'Aging LepR+Pa− (sens.)',
         'LPS_CAR': 'LPS (CAR)', 'PIC_CAR': 'poly I:C (CAR)', 'B3_Dox4x_MSPC': 'Dox ×4 (MSPC)', 'B3_Dox1x_MSPC': 'Dox ×1 (MSPC)'}
AXIS = {a: c for a, _, c, _ in CONTRASTS}
AGING = ['Age_LepRPa', 'Age_CAR']
CHEMO = ['B3_Dox4x_MSPC', 'B3_Dox1x_MSPC']

# 遺伝子パネル（事前定義。依頼の niche 因子 + ER 標的 + 骨芽 + 脂肪 + 炎症 + 老化）
GENE_PANEL = {
    'niche factor': ['Cxcl12', 'Kitl', 'Angpt1', 'Vcam1', 'Lepr', 'Spp1'],
    'osteoclast support': ['Tnfsf11', 'Tnfrsf11b', 'Csf1'],
    'inflammation': ['Il6', 'Ccl2', 'Il1b', 'Cxcl1', 'Tnf'],
    'adipocyte': ['Adipoq', 'Pparg', 'Cebpa', 'Lpl'],
    'osteoblast': ['Runx2', 'Sp7', 'Alpl', 'Col1a1', 'Ibsp', 'Bglap', 'Dmp1'],
    'ER target': ['Esr1', 'Pgr', 'Greb1'],
    'senescence': ['Cdkn1a', 'Cdkn2a', 'Serpine1', 'Igfbp3'],
    'IFN': ['Isg15', 'Ifit1', 'Irf7', 'Stat1'],
}


def read_de(src):
    d = pd.read_csv(src).dropna(subset=['gene'])
    d['G'] = d.gene.str.upper()
    return d.drop_duplicates('G').set_index('G')


def ranking(d):
    d = d.dropna(subset=['lfcSE'])
    d = d[d.lfcSE > 0]
    return (d.log2FoldChange / d.lfcSE).sort_values(ascending=False)


def run_gsea(sets, DE):
    (OUT / 'gsea').mkdir(parents=True, exist_ok=True)
    rows = []
    for lb, *_ in CONTRASTS:
        rnk = ranking(DE[lb])
        r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=4, max_size=1000,
                       seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
        r = r.assign(label=lb, n_ranked=len(rnk))
        r.to_csv(OUT / 'gsea' / f'{lb}.csv', index=False)
        rows.append(r)
        print(lb, len(rnk), 'genes;', (r['FDR q-val'].astype(float) < 0.05).sum(), 'sets FDR<0.05')
    L = pd.concat(rows, ignore_index=True).rename(columns={'Term': 'term', 'FDR q-val': 'fdr', 'NOM p-val': 'pval',
                                                           'Lead_genes': 'lead_genes', 'Tag %': 'tag_pct'})
    for k in ['NES', 'ES', 'fdr', 'pval']:
        L[k] = L[k].astype(float)
    L['matched'] = L.tag_pct.astype(str).str.split('/').str[1].astype(int)
    return L[['label', 'term', 'NES', 'ES', 'pval', 'fdr', 'matched', 'lead_genes', 'n_ranked']]


def heatmap_ax(ax, W, F, labels, cats, vmax=3, fs=6):
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(W.values), cmap=cmap, vmin=-vmax, vmax=vmax, aspect='auto')
    for i in range(W.shape[0]):
        for j in range(W.shape[1]):
            v, f = W.iloc[i, j], F.iloc[i, j]
            if np.isnan(v):
                continue
            if f < 0.05:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=fs - 0.5, fontweight='bold',
                        color='white' if abs(v) > 0.67 * vmax else 'black')
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, lw=0.7, ec='black'))
            elif f < 0.25:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=fs - 0.5, color='#333333')
    ax.set_yticks(range(len(labels)), labels, fontsize=fs)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    prev = None
    for i, c in enumerate(cats):
        if prev is not None and c != prev:
            ax.axhline(i - 0.5, c='k', lw=0.5)
        prev = c
    ax_ = [AXIS[c] for c in W.columns]
    for j in range(1, len(ax_)):
        if ax_[j] != ax_[j - 1]:
            ax.axvline(j - 0.5, c='k', lw=1.0)
    return im


def heatmap(W, F, labels, cats, fname, title, height_mm):
    common.setup()
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(height_mm)))
    im = heatmap_ax(ax, W, F, labels, cats)
    ax.set_xticks(range(W.shape[1]), [NAME[c].replace('\n', ', ') for c in W.columns], rotation=60, ha='right', fontsize=6)
    cb = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.01)
    cb.set_label('NES', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=7.5, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'{fname}.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    print('saved', OUT / f'{fname}.png')


def sgn(W, F, t, c):
    v, f = W.loc[t, c], F.loc[t, c]
    if np.isnan(v) or f >= 0.05:
        return 0
    return int(np.sign(v))


def combine(W, F, t, cols):
    """2 本の contrast をまとめた向き: (dir, support)。"""
    s = [sgn(W, F, t, c) for c in cols]
    nz = [v for v in s if v != 0]
    if not nz:
        return 0, 'none'
    if len(set(nz)) > 1:
        return 0, 'discordant'
    d = nz[0]
    other = [c for c, v in zip(cols, s) if v == 0]
    if not other:
        return d, 'both'
    # 片方のみ有意: もう片方の NES が同符号なら 'one (+other same sign)'、逆符号なら 'one (other opposite NS)'
    o = W.loc[t, other[0]]
    return d, 'one, other same sign' if np.sign(o) == d else 'one, other opposite (NS)'


def module_class(W, F, meta):
    arrow = {1: '↑', -1: '↓', 0: '·'}
    rows = []
    for t in W.index:
        a, asup = combine(W, F, t, AGING)
        c, csup = combine(W, F, t, CHEMO)
        if a and c:
            k = f'aging{arrow[a]} chemo{arrow[c]} (' + ('same' if a == c else 'OPPOSITE') + ')'
        elif a:
            k = f'aging{arrow[a]} only (chemo not sig.)'
        elif c:
            k = f'chemo{arrow[c]} only (aging not sig.)'
        else:
            k = 'neither' + (' [aging discordant]' if asup == 'discordant' else '')
        r = dict(term=t, cls=k, aging=arrow[a], aging_support=asup, chemo=arrow[c], chemo_support=csup,
                 LPS=arrow[sgn(W, F, t, 'LPS_CAR')], PIC=arrow[sgn(W, F, t, 'PIC_CAR')],
                 aging_sens=arrow[sgn(W, F, t, 'Age_LepRPaN')])
        for col in W.columns:
            r[f'NES:{col}'] = W.loc[t, col]
            r[f'FDR:{col}'] = F.loc[t, col]
        rows.append(r)
    d = pd.DataFrame(rows).merge(meta[['term', 'label', 'category', 'source']], on='term')
    d = d.sort_values(['cls', 'category', 'term'])
    d.to_csv(OUT / 'module_class.csv', index=False)
    print(d.cls.value_counts().to_string())
    return d


def concordance(W, F, meta):
    rows = []
    for a in AGING + ['Age_LepRPaN', 'LPS_CAR', 'PIC_CAR']:
        for c in CHEMO + ([x for x in AGING if x != a] if a in AGING else []):
            for t in W.index:
                x, y = sgn(W, F, t, a), sgn(W, F, t, c)
                k = ('same' if x == y else 'opposite') if x and y else ('A only' if x else ('B only' if y else 'neither'))
                rows.append(dict(A=a, B=c, term=t, NES_A=W.loc[t, a], FDR_A=F.loc[t, a], NES_B=W.loc[t, c],
                                 FDR_B=F.loc[t, c], call=k))
    d = pd.DataFrame(rows).merge(meta[['term', 'label', 'category']], on='term')
    d.to_csv(OUT / 'concordance_sets.csv', index=False)
    s = d.pivot_table(index=['A', 'B'], columns='call', values='term', aggfunc='count', fill_value=0)
    s['spearman_NES'] = [W[[a, b]].dropna().corr(method='spearman').iloc[0, 1] for a, b in s.index]
    s.to_csv(OUT / 'concordance_counts.csv')
    print(s.round(2).to_string())
    W.corr(method='spearman').to_csv(OUT / 'spearman_nes.csv')


def gene_level(DE):
    # 遺伝子単位の Wald 統計量の Spearman（両方で検定された遺伝子）
    S = pd.DataFrame({lb: ranking(DE[lb]) for lb in DE})
    rows = []
    for a in S.columns:
        for b in S.columns:
            if a == b:
                continue
            j = S[[a, b]].dropna()
            rows.append(dict(A=a, B=b, n=len(j), spearman=j[a].corr(j[b], method='spearman')))
    pd.DataFrame(rows).to_csv(OUT / 'spearman_gene_stat.csv', index=False)

    # DEG の重なり（向き付き）
    def deg(d):
        return d[(d.padj < 0.05) & (d.log2FoldChange.abs() >= 1)]
    ov, genes = [], []
    for a in AGING + ['Age_LepRPaN', 'LPS_CAR', 'PIC_CAR']:
        for c in CHEMO:
            A, C = DE[a], DE[c]
            U = A.dropna(subset=['padj']).index.intersection(C.dropna(subset=['padj']).index)
            da, dc = deg(A.loc[U]), deg(C.loc[U])
            both = da.index.intersection(dc.index)
            same = [g for g in both if np.sign(A.loc[g, 'log2FoldChange']) == np.sign(C.loc[g, 'log2FoldChange'])]
            p = hypergeom.sf(len(both) - 1, len(U), len(da), len(dc)) if len(both) else 1.0
            ov.append(dict(aging_or_ref=a, chemo=c, universe=len(U), deg_A=len(da), deg_chemo=len(dc), overlap=len(both),
                           expected=len(da) * len(dc) / len(U), hypergeom_p=p, same_dir=len(same), opposite_dir=len(both) - len(same)))
            for g in both:
                genes.append(dict(A=a, chemo=c, gene=A.loc[g, 'gene'], lfc_A=A.loc[g, 'log2FoldChange'], padj_A=A.loc[g, 'padj'],
                                  lfc_chemo=C.loc[g, 'log2FoldChange'], padj_chemo=C.loc[g, 'padj'],
                                  direction=('both ↑' if A.loc[g, 'log2FoldChange'] > 0 and C.loc[g, 'log2FoldChange'] > 0 else
                                             'both ↓' if A.loc[g, 'log2FoldChange'] < 0 and C.loc[g, 'log2FoldChange'] < 0 else
                                             'A ↑ chemo ↓' if A.loc[g, 'log2FoldChange'] > 0 else 'A ↓ chemo ↑')))
    ov = pd.DataFrame(ov)
    ov.to_csv(OUT / 'deg_overlap.csv', index=False)
    print(ov.round(4).to_string())
    pd.DataFrame(genes).sort_values(['A', 'chemo', 'direction', 'gene']).to_csv(OUT / 'deg_overlap_genes.csv', index=False)

    # パネル遺伝子
    order = [g for v in GENE_PANEL.values() for g in v]
    Lf = pd.DataFrame({lb: DE[lb].log2FoldChange.reindex([g.upper() for g in order]).values for lb in DE}, index=order)
    Pf = pd.DataFrame({lb: DE[lb].padj.reindex([g.upper() for g in order]).values for lb in DE}, index=order)
    out = pd.concat({'log2FC': Lf, 'padj': Pf}, axis=1)
    out.index.name = 'gene'
    out.to_csv(OUT / 'panel_genes.csv')
    return Lf, Pf


def gene_heatmap_ax(ax, Lf, Pf, fs=6):
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(Lf.values), cmap=cmap, vmin=-3, vmax=3, aspect='auto')
    for i in range(Lf.shape[0]):
        for j in range(Lf.shape[1]):
            p = Pf.iloc[i, j]
            if p < 0.05:
                ax.text(j, i, '**', ha='center', va='center', fontsize=fs)
            elif p < 0.1:
                ax.text(j, i, '*', ha='center', va='center', fontsize=fs)
    sub = {g: k for k, v in GENE_PANEL.items() for g in v}
    ax.set_yticks(range(len(Lf)), [f'$\\it{{{g}}}$' for g in Lf.index], fontsize=fs)
    prev = None
    for i, g in enumerate(Lf.index):
        if prev is not None and sub[g] != prev:
            ax.axhline(i - 0.5, c='k', lw=0.5)
        prev = sub[g]
    ax_ = [AXIS[c] for c in Lf.columns]
    for j in range(1, len(ax_)):
        if ax_[j] != ax_[j - 1]:
            ax.axvline(j - 0.5, c='k', lw=1.0)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # 右側にカテゴリ名
    starts = {}
    for i, g in enumerate(Lf.index):
        starts.setdefault(sub[g], []).append(i)
    for k, ii in starts.items():
        ax.text(Lf.shape[1] - 0.3, np.mean(ii), k, fontsize=fs - 0.5, va='center', ha='left', color='#444444')
    return im


# 主図に載せる module（事前に決めた niche 軸。chemo2 の niche 図から主要なものを抜粋。結果を見て選んでいない）
FIG_SETS = [
    'HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_IL6_JAK_STAT3_SIGNALING',
    'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'NICHE_INFLAM_CHEMOKINE',
    'SAUL_SEN_MAYO', 'REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP', 'HALLMARK_P53_PATHWAY',
    'HALLMARK_ADIPOGENESIS', 'GOBP_FAT_CELL_DIFFERENTIATION', 'NICHE_ADIPO',
    'GOBP_OSSIFICATION', 'GOBP_OSTEOBLAST_DIFFERENTIATION', 'GOBP_BONE_MINERALIZATION',
    'NICHE_HSC_RETENTION', 'GOBP_HEMATOPOIETIC_STEM_CELL_HOMEOSTASIS',
    'HALLMARK_ESTROGEN_RESPONSE_EARLY', 'HALLMARK_ESTROGEN_RESPONSE_LATE',
    'HALLMARK_ANGIOGENESIS', 'GOBP_VASCULATURE_DEVELOPMENT',
    'NABA_CORE_MATRISOME', 'GOBP_COLLAGEN_FIBRIL_ORGANIZATION', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
    'HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT',
]
FIG_CAT = ['inflammation'] * 6 + ['senescence/SASP'] * 3 + ['adipo'] * 3 + ['bone formation'] * 3 + ['HSC niche'] * 2 + \
    ['estrogen'] * 2 + ['vascular'] * 2 + ['ECM'] * 3 + ['proliferation'] * 2
FIG_COLS = ['Age_LepRPa', 'Age_CAR', 'LPS_CAR', 'PIC_CAR', 'B3_Dox4x_MSPC', 'B3_Dox1x_MSPC']


def main_figure(W, F, meta, Lf, Pf):
    common.setup()
    lab = dict(zip(meta.term, meta.label))
    fig = plt.figure(figsize=(common.mm(143), common.mm(150)))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.8], wspace=0.75)
    ax = fig.add_subplot(gs[0])
    Ws, Fs = W.loc[FIG_SETS, FIG_COLS], F.loc[FIG_SETS, FIG_COLS]
    im = heatmap_ax(ax, Ws, Fs, [lab[t] for t in FIG_SETS], FIG_CAT, fs=6)
    ax.set_xticks(range(len(FIG_COLS)), [SHORT[c] for c in FIG_COLS], rotation=55, ha='right', fontsize=6)
    cb = fig.colorbar(im, ax=ax, orientation='horizontal', shrink=0.5, pad=0.17, aspect=20)
    cb.set_label('NES (bold/box: FDR<0.05)', fontsize=6)
    cb.ax.tick_params(labelsize=6)
    common.panel_label(ax, 'A', dx=-0.95, dy=1.01)
    ax2 = fig.add_subplot(gs[1])
    im2 = gene_heatmap_ax(ax2, Lf[FIG_COLS], Pf[FIG_COLS], fs=5.5)
    ax2.set_xticks(range(len(FIG_COLS)), [SHORT[c] for c in FIG_COLS], rotation=55, ha='right', fontsize=6)
    cb2 = fig.colorbar(im2, ax=ax2, orientation='horizontal', shrink=0.5, pad=0.17, aspect=20)
    cb2.set_label('log$_2$FC (** padj<0.05, * <0.1)', fontsize=6)
    cb2.ax.tick_params(labelsize=6)
    common.panel_label(ax2, 'B', dx=-0.40, dy=1.01)
    for ext in ('png', 'pdf', 'tiff'):
        fig.savefig(OUT / f'fig_niche.{ext}', dpi=common.DPI, bbox_inches='tight', pad_inches=0.02)
    from PIL import Image
    w, h = Image.open(OUT / 'fig_niche.png').size
    print(f'fig_niche: {w / common.DPI * 25.4:.0f} x {h / common.DPI * 25.4:.0f} mm')
    plt.close(fig)
    print('saved', OUT / 'fig_niche.png')


def main():
    DE = {lb: read_de(src) for lb, _, _, src in CONTRASTS}
    sets, meta = GS.collection()
    L = run_gsea(sets, DE)
    L = L.merge(meta[['term', 'label', 'category', 'source', 'niche']].rename(columns={'label': 'set_label'}), on='term', how='left')
    L.to_csv(OUT / 'nes_long.csv', index=False)
    order = [c for c, *_ in CONTRASTS]
    W = L.pivot(index='term', columns='label', values='NES')[order]
    F = L.pivot(index='term', columns='label', values='fdr')[order]
    pd.concat({'NES': W, 'FDR': F}, axis=1).to_csv(OUT / 'nes_wide.csv')
    catorder = ['inflammation', 'myeloid', 'osteoclast', 'senescence/SASP', 'HSC niche', 'vascular', 'ECM',
                'bone formation', 'adipo', 'proliferation', 'estrogen', 'composition?', 'hand panel']
    nm = meta[meta.niche].copy()
    nm['co'] = nm.category.map({c: i for i, c in enumerate(catorder)})
    nm = nm.sort_values(['co', 'source', 'term'])
    nm = nm[nm.term.isin(W.index)]
    heatmap(W.loc[nm.term], F.loc[nm.term], [f'{l}  [{c}]' for l, c in zip(nm.label, nm.category)], list(nm.category),
            'heatmap_niche', 'LepR+/CAR stroma: prerank GSEA NES (bold+box: FDR<0.05; plain: FDR<0.25)', 205)
    hm = meta[meta.source == 'Hallmark'].sort_values('term')
    heatmap(W.loc[hm.term], F.loc[hm.term], list(hm.label), ['h'] * len(hm), 'heatmap_hallmark',
            'Hallmark 50 in LepR+/CAR stroma (bold+box: FDR<0.05; plain: FDR<0.25)', 205)
    concordance(W, F, meta)
    module_class(W, F, meta)
    Lf, Pf = gene_level(DE)
    main_figure(W, F, meta, Lf, Pf)


if __name__ == '__main__':
    main()

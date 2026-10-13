"""雌の加齢 contrast（GSE210584 擬似バルク、GSE132040 Bone）を prerank GSEA にかけ、
既存の加齢♂皮質骨・OVX・化療（chemo2）と加齢 LepR+/CAR（niche）の NES と並べる。niche 因子パネルの log2FC 表と volcano も出す。

GSEA の設定は prelim/chemo2/gsea.py と同じ（gene set = genesets.collection()、Wald 統計量の prerank、gseapy 1000 perm、
seed 0、min 4、max 1000、FDR は 72 set で計算）。既存の列は nes_wide.csv を読むだけ（再計算しない）。

  uv run python -I prelim/aging_f/compare.py
出力（out/prelim/aging_f/）: gsea/<label>.csv, nes_long.csv, nes_wide_new.csv, nes_wide_all.csv,
  heatmap_niche.png|pdf, heatmap_hallmark.png|pdf, concordance_nes_spearman.csv, concordance_sets.csv,
  gene_stat_spearman.csv, panel_genes.csv, heatmap_panel.png|pdf, volcano_4A.png|pdf, deg_counts.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
import common  # noqa: E402
import genesets as GS  # noqa: E402

OUT = ROOT / 'out/prelim/aging_f'
DEN = OUT / 'de'
C2 = ROOT / 'out/prelim/chemo2'
NI = ROOT / 'out/prelim/niche'

# 新しい contrast: (label, 表示名, 軸)
NEW = [
    ('sc210584_CAR', 'Middle-age vs young\nCAR (scRNA PB) ♀', 'aging♀'),
    ('sc210584_Osteo', 'Middle-age vs young\nOsteo-lineage (PB) ♀', 'aging♀'),
    ('sc210584_Endothelial', 'Middle-age vs young\nEndothelial (PB) ♀', 'aging♀'),
    ('TMS_F_old', 'TMS bone 12-21 vs 1-6m ♀', 'aging♀'),
    ('TMS_F_old_comp', 'TMS bone ♀\n+composition adj.', 'aging♀'),
    ('TMS_F_old_ex1', 'TMS bone ♀ 12-21 vs 3-6m', 'aging♀'),
    ('TMS_F_trend', 'TMS bone ♀ age trend', 'aging♀'),
    ('TMS_M_old', 'TMS bone 12-21 vs 1-6m ♂', 'aging♂'),
    ('TMS_M_old_comp', 'TMS bone ♂\n+composition adj.', 'aging♂'),
    ('TMS_M_trend', 'TMS bone ♂ age trend', 'aging♂'),
]
# 既存列（nes_wide から読む）
OLD = [
    (NI, 'Age_LepRPa', 'Aging LepR+PDGFRa+\n(GSE130885) ♀', 'aging (existing)'),
    (NI, 'Age_CAR', 'Aging CAR\n(GSE133922) ♂', 'aging (existing)'),
    (C2, 'Age_2.5y', 'Aging 2.5y cortical ♂\n(current 4A)', 'aging (existing)'),
    (C2, 'A_OVX_cort', 'OVX cortical ♀', 'OVX'),
    (C2, 'A_OVX_trab', 'OVX trabecular ♀', 'OVX'),
    (C2, 'B3_Dox4x_MSPC', 'Dox 4x LepR+ MSPC ♀', 'chemo'),
    (C2, 'B3_Dox1x_MSPC', 'Dox 1x LepR+ MSPC ♀', 'chemo'),
    (C2, 'B2_Dox4x_OB', 'Dox 4x Col2.3+ OB ♀', 'chemo'),
    (C2, 'B_Dox', 'Dox 1x bone ♂', 'chemo'),
]
OLD_DE = {'Age_LepRPa': NI / 'de/Age_LepRPa.csv', 'Age_CAR': NI / 'de/Age_CAR.csv',
          'Age_2.5y': ROOT / 'out/prelim/de/Age_2.5y.csv', 'A_OVX_cort': ROOT / 'out/prelim/de/A_OVX_cort.csv',
          'A_OVX_trab': ROOT / 'out/prelim/de/A_OVX_trab.csv',
          'B3_Dox4x_MSPC': ROOT / 'out/prelim/chemo/de/B3_Dox4x_MSPC.csv',
          'B3_Dox1x_MSPC': ROOT / 'out/prelim/chemo/de/B3_Dox1x_MSPC.csv',
          'B2_Dox4x_OB': ROOT / 'out/prelim/de/B2_Dox4x_OB.csv', 'B_Dox': ROOT / 'out/prelim/de/B_Dox.csv'}
NAME = {**{a: b for a, b, _ in NEW}, **{a: b for _, a, b, _ in OLD}}
AXIS = {**{a: c for a, _, c in NEW}, **{a: c for _, a, _, c in OLD}}
PANEL = ['Cxcl12', 'Kitl', 'Tnfsf11', 'Tnfrsf11b', 'Il6', 'Il1b', 'Ccl2', 'Cdkn1a', 'Cdkn2a', 'Serpine1', 'Spp1',
         'Dmp1', 'Runx2', 'Sp7', 'Pgr', 'Greb1', 'Esr1', 'Igf1', 'Lepr', 'Adipoq', 'Bglap', 'Col1a1']
REF_AGING = 'Age_2.5y'
CHEMO_REF = ['B3_Dox4x_MSPC', 'B3_Dox1x_MSPC']


def de_path(lb):
    return OLD_DE[lb] if lb in OLD_DE else DEN / f'{lb}.csv'


def ranking(lb):
    d = pd.read_csv(de_path(lb)).dropna(subset=['lfcSE', 'gene'])
    d = d[d.lfcSE > 0]
    return (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean().sort_values(ascending=False)


def run_gsea(sets):
    (OUT / 'gsea').mkdir(parents=True, exist_ok=True)
    rows = []
    for lb, _, _ in NEW:
        rnk = ranking(lb)
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
    return L[['label', 'term', 'NES', 'ES', 'pval', 'fdr', 'lead_genes', 'n_ranked']]


def heatmap(W, F, labels, cats, fname, title, height_mm):
    common.setup()
    cols = list(W.columns)
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(height_mm)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(W.values), cmap=cmap, vmin=-3, vmax=3, aspect='auto')
    for i in range(W.shape[0]):
        for j in range(W.shape[1]):
            v, f = W.iloc[i, j], F.iloc[i, j]
            if np.isnan(v):
                continue
            if f < 0.05:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=5.5, fontweight='bold',
                        color='white' if abs(v) > 2 else 'black')
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, lw=0.7, ec='black'))
            elif f < 0.25:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=5.5, color='#333333')
    ax.set_yticks(range(len(labels)), labels, fontsize=6)
    ax.set_xticks(range(len(cols)), [NAME[c].replace('\n', ', ') for c in cols], rotation=60, ha='right', fontsize=5.5)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    prev = None
    for i, c in enumerate(cats):
        if prev is not None and c != prev:
            ax.axhline(i - 0.5, c='k', lw=0.5)
        prev = c
    ax_ = [AXIS[c] for c in cols]
    for j in range(1, len(cols)):
        if ax_[j] != ax_[j - 1]:
            ax.axvline(j - 0.5, c='k', lw=1.0)
    cb = fig.colorbar(im, ax=ax, shrink=0.3, pad=0.01)
    cb.set_label('NES', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=7.5, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'{fname}.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    print('saved', OUT / f'{fname}.png')


def lfc_padj(lb):
    d = pd.read_csv(de_path(lb)).dropna(subset=['gene']).drop_duplicates('gene').set_index('gene')
    return d.log2FoldChange, d.padj


def panel(cols):
    L, P = {}, {}
    for lb in cols:
        l, p = lfc_padj(lb)
        L[lb], P[lb] = l.reindex(PANEL), p.reindex(PANEL)
    L, P = pd.DataFrame(L), pd.DataFrame(P)
    out = pd.concat({'log2FC': L, 'padj': P}, axis=1)
    out.index.name = 'gene'
    out.to_csv(OUT / 'panel_genes.csv')
    common.setup()
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(120)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(L.values), cmap=cmap, vmin=-2, vmax=2, aspect='auto')
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            v, p = L.iloc[i, j], P.iloc[i, j]
            if np.isnan(v):
                continue
            star = '**' if p < 0.05 else ('*' if p < 0.1 else '')
            ax.text(j, i, f'{v:.1f}{star}', ha='center', va='center', fontsize=5,
                    color='white' if abs(v) > 1.5 else 'black', fontweight='bold' if star else 'normal')
    ax.set_yticks(range(len(PANEL)), [f'$\\it{{{g}}}$' for g in PANEL], fontsize=6.5)
    ax.set_xticks(range(L.shape[1]), [NAME[c].replace('\n', ', ') for c in L.columns], rotation=60, ha='right', fontsize=5.5)
    ax_ = [AXIS[c] for c in L.columns]
    for j in range(1, len(ax_)):
        if ax_[j] != ax_[j - 1]:
            ax.axvline(j - 0.5, c='k', lw=1.0)
    cb = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.01)
    cb.set_label('log$_2$FC (trend: per 10 months)', fontsize=6.5)
    ax.set_title('Niche factor panel (** padj<0.05, * <0.1; grey = filtered)', fontsize=7.5, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'heatmap_panel.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)


def volcano(labels):
    """4A 候補の見栄え確認。prelim/volcano.py と同じ閾値（|log2FC|>=1, padj<0.05）・配色。"""
    common.setup()
    fig, axes = plt.subplots(1, len(labels), figsize=(common.mm(170), common.mm(62)))
    hl = PANEL
    for ax, lb, pl in zip(axes, labels, 'ABCD'):
        d = pd.read_csv(de_path(lb)).dropna(subset=['log2FoldChange', 'pvalue'])
        y = -np.log10(d.pvalue.clip(lower=1e-300))
        sig = (d.padj < 0.05) & (d.log2FoldChange.abs() >= 1)
        up, dn = sig & (d.log2FoldChange > 0), sig & (d.log2FoldChange < 0)
        ax.scatter(d.log2FoldChange[~sig], y[~sig], s=1.5, c='#CFCFCF', lw=0, rasterized=True)
        ax.scatter(d.log2FoldChange[up], y[up], s=2, c='#D98880', lw=0, rasterized=True, label=f'Up ({up.sum()})')
        ax.scatter(d.log2FoldChange[dn], y[dn], s=2, c='#7FB3D5', lw=0, rasterized=True, label=f'Down ({dn.sum()})')
        m = d.gene.isin(hl) & (d.padj < 0.1)
        ax.scatter(d.log2FoldChange[m], y[m], s=9, facecolor='none', edgecolor='#1F4E79', lw=0.6)
        for gx, gy, gn in zip(d.log2FoldChange[m], y[m], d.gene[m]):
            ax.annotate(gn, (gx, gy), xytext=(2, 2), textcoords='offset points', fontsize=5.5, fontstyle='italic')
        if (d.padj < 0.05).any():
            ax.axhline(-np.log10(d.pvalue[d.padj < 0.05].max()), ls='--', lw=0.5, c='k')
        for v in (-1, 1):
            ax.axvline(v, ls='--', lw=0.5, c='k')
        xm = max(np.nanpercentile(d.log2FoldChange.abs(), 99.9) * 1.1, 3)
        ax.set_xlim(-xm, xm)
        inx = d.log2FoldChange.abs() <= xm  # 範囲外の点で y 軸が伸びないように（範囲外は描かれない）
        ax.set_ylim(0, y[inx].max() * 1.08)
        ax.set_xlabel('log$_2$ fold change')
        ax.set_ylabel('$-$log$_{10}$ P')
        ax.set_title(NAME[lb].replace('\n', ' '), fontsize=7)
        ax.legend(fontsize=5.5, frameon=False, loc='upper left', handletextpad=0.2)
        common.panel_label(ax, pl)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'volcano_4A.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)


def main():
    sets, meta = GS.collection()
    L = run_gsea(sets)
    L = L.merge(meta[['term', 'label', 'category', 'source', 'niche']].rename(columns={'label': 'set_label'}),
                on='term', how='left')
    L.to_csv(OUT / 'nes_long.csv', index=False)
    order = [a for a, *_ in NEW]
    Wn = L.pivot(index='term', columns='label', values='NES')[order]
    Fn = L.pivot(index='term', columns='label', values='fdr')[order]
    pd.concat({'NES': Wn, 'FDR': Fn}, axis=1).to_csv(OUT / 'nes_wide_new.csv')
    Wo, Fo = {}, {}
    for d, lb, *_ in OLD:
        w = pd.read_csv(d / 'nes_wide.csv', header=[0, 1], index_col=0)
        Wo[lb], Fo[lb] = w[('NES', lb)], w[('FDR', lb)]
    Wo, Fo = pd.DataFrame(Wo), pd.DataFrame(Fo)
    W = Wn.join(Wo)
    F = Fn.join(Fo)
    # 列順: 新（雌加齢 sc → TMS 雌 → TMS 雄）→ 既存加齢 → OVX → 化療
    W, F = W[order + list(Wo.columns)], F[order + list(Wo.columns)]
    pd.concat({'NES': W, 'FDR': F}, axis=1).to_csv(OUT / 'nes_wide_all.csv')

    catorder = ['inflammation', 'myeloid', 'osteoclast', 'senescence/SASP', 'HSC niche', 'vascular', 'ECM',
                'bone formation', 'adipo', 'proliferation', 'estrogen', 'composition?', 'hand panel']
    nm = meta[meta.niche].copy()
    nm['co'] = nm.category.map({c: i for i, c in enumerate(catorder)})
    nm = nm.sort_values(['co', 'source', 'term'])
    nm = nm[nm.term.isin(W.index)]
    heatmap(W.loc[nm.term], F.loc[nm.term], [f'{l}  [{c}]' for l, c in zip(nm.label, nm.category)],
            list(nm.category), 'heatmap_niche',
            'Female aging vs existing aging / OVX / Dox: prerank GSEA NES (bold+box FDR<0.05; plain FDR<0.25)', 205)
    hm = meta[meta.source == 'Hallmark'].sort_values('term')
    heatmap(W.loc[hm.term], F.loc[hm.term], list(hm.label), ['h'] * len(hm), 'heatmap_hallmark',
            'Hallmark 50, prerank GSEA NES (bold+box FDR<0.05; plain FDR<0.25)', 205)

    # NES ベクトルの一致（全 72 set、set は重なるので独立ではない。記述用）
    refs = ['Age_LepRPa', 'Age_CAR', REF_AGING, 'A_OVX_cort'] + CHEMO_REF + ['B_Dox']
    sp = pd.DataFrame({r: {c: W[[c, r]].dropna().corr(method='spearman').iloc[0, 1] for c in order} for r in refs})
    sp.to_csv(OUT / 'concordance_nes_spearman.csv')
    print(sp.round(2).to_string())
    # set ごとの向き一致（両方 FDR<0.05）
    rows = []
    for c in order:
        for r in refs:
            a, b = W[c], W[r]
            sa, sb = F[c] < 0.05, F[r] < 0.05
            both = sa & sb
            rows.append(dict(new=c, ref=r, n_sig_new=int(sa.sum()), n_sig_ref=int(sb.sum()), both_sig=int(both.sum()),
                             same=int((both & (np.sign(a) == np.sign(b))).sum()),
                             opposite=int((both & (np.sign(a) != np.sign(b))).sum())))
    cs = pd.DataFrame(rows)
    cs.to_csv(OUT / 'concordance_sets.csv', index=False)

    # 遺伝子単位: Wald 統計量の Spearman（共通遺伝子）と DEG の重なり
    R = {lb: ranking(lb) for lb in order + refs}
    rows = []
    for c in order:
        for r in refs:
            j = pd.concat([R[c], R[r]], axis=1, join='inner').dropna()
            dc, dr = lfc_padj(c), lfc_padj(r)
            degc = set(dc[0][(dc[1] < 0.05) & (dc[0].abs() >= 1)].index)
            degr = set(dr[0][(dr[1] < 0.05) & (dr[0].abs() >= 1)].index)
            ov = degc & degr
            same = sum(np.sign(dc[0][g]) == np.sign(dr[0][g]) for g in ov)
            rows.append(dict(new=c, ref=r, n_genes=len(j), spearman_stat=j.corr(method='spearman').iloc[0, 1],
                             deg_new=len(degc), deg_ref=len(degr), deg_overlap=len(ov), overlap_same_sign=int(same),
                             overlap_genes=','.join(sorted(ov))[:500]))
    gs = pd.DataFrame(rows)
    gs.to_csv(OUT / 'gene_stat_spearman.csv', index=False)
    print(gs[['new', 'ref', 'spearman_stat', 'deg_new', 'deg_ref', 'deg_overlap', 'overlap_same_sign']].round(2).to_string())

    # DEG 数（volcano.py と同じ閾値の組）
    rows = []
    for lb in order:
        d = pd.read_csv(de_path(lb))
        for p, l in [(0.05, 1.0), (0.05, 0.585), (0.05, 0.0), (0.1, 1.0)]:
            s = (d.padj < p) & (d.log2FoldChange.abs() >= l)
            rows.append(dict(label=lb, padj=p, lfc=l, up=int((s & (d.log2FoldChange > 0)).sum()),
                             down=int((s & (d.log2FoldChange < 0)).sum()), tested=int(d.padj.notna().sum())))
    pd.DataFrame(rows).to_csv(OUT / 'deg_counts.csv', index=False)

    panel(['sc210584_CAR', 'sc210584_Osteo', 'sc210584_Endothelial', 'TMS_F_old', 'TMS_F_old_comp', 'TMS_F_trend',
           'TMS_M_old', 'Age_LepRPa', 'Age_CAR', 'Age_2.5y', 'A_OVX_cort', 'B3_Dox4x_MSPC', 'B3_Dox1x_MSPC'])
    volcano(['sc210584_CAR', 'TMS_F_old', 'Age_2.5y'])


if __name__ == '__main__':
    main()

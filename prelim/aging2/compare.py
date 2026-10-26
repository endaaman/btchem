"""Figure 4A 案2・案3 を既存の加齢・閉経・化療の contrast と並べる（log2FC の相関、遺伝子セットの向き、図）。

  uv run python -I prelim/aging2/compare.py          # repo root から（先に de.py, gsea.py）

既存の値は読むだけ（再計算しない）:
  DE:  out/prelim/de/{Age_2.5y,A_OVX_cort,B_Dox,B2_Dox4x_OB}.csv, out/prelim/chemo2/de/B_Dox_compadj.csv,
       out/prelim/chemo/de/B3_Dox{4x,1x}_MSPC.csv, out/prelim/aging_f/de/{sc210584_CAR,TMS_F_old,TMS_F_old_comp}.csv,
       out/prelim/niche/de/{Age_LepRPa,Age_CAR}.csv, out/prelim/human/de/{H_PostMeno,H_Old}.csv
  NES: out/prelim/chemo2/nes_wide.csv, out/prelim/aging_f/nes_wide_all.csv（マウス、同じ 72 set・同じ設定）、
       out/prelim/summary/gsea_human/{H_PostMeno,H_Old}.csv（ヒト、MSigDB Hs の同名 set。手製パネルはヒト側に 3 つだけ）
マウス → ヒトは既存 compare_h.py と同じく大文字化（Fasl→FASLG, Trp53→TP53）、大文字化後の重複は落とす。

出力（out/prelim/aging2/）:
  corr_log2fc.csv          新 contrast × 既存 contrast の log2FC Spearman（全共通遺伝子 / どちらかで padj<0.1 の遺伝子）
  corr_nes.csv             同じ組の NES Spearman と、両方 FDR<0.05 で同方向 / 逆方向の set 数
  module_grid.csv|png|pdf  要約モジュール 8+4 行 × 加齢・新規・化療の列（FDR<0.05 は色、0.05–0.25 は枠）
  heatmap_niche.png|pdf    niche 関連 set の NES（新規 + 既存）
  corr_matrix.png|pdf      log2FC Spearman の行列（新規 × 既存）
  volcano_aging2.png|pdf   主 3 contrast の volcano（prelim/volcano.py の plot をそのまま使う）
  panel_genes.csv, heatmap_panel.png|pdf  注目遺伝子パネルの log2FC
  deg_counts.csv           閾値ごとの DEG 数
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[2]
for p in ['prelim', 'prelim/chemo2', 'prelim/summary', 'prelim/human']:
    sys.path.insert(0, str(ROOT / p))
import common  # noqa: E402
import genesets as GS  # noqa: E402
import summary_modules as SM  # noqa: E402
import hgenes  # noqa: E402
from genes import INFLAM, APOP  # noqa: E402
import volcano as V  # noqa: E402  （plot() を借りる。volcano.py の main は呼ばない）

OUT = ROOT / 'out/prelim/aging2'
P = ROOT / 'out/prelim'

# 新しい contrast（行）: label -> (表示名, DE csv)
NEW = {
    'G235_F_age': ('♀ tibia 22 vs 6 mo (GSE235942)', OUT / 'de/G235_F_age.csv'),
    'G235_F_age_comp': ('♀ tibia 22 vs 6 mo, comp. adj.', OUT / 'de/G235_F_age_comp.csv'),
    'G235_M_age': ('♂ tibia 22 vs 6 mo (GSE235942)', OUT / 'de/G235_M_age.csv'),
    'G255_MSC_age': ('♀ LepR+ MSC 18 vs 3 mo (GSE255653)', OUT / 'de/G255_MSC_age.csv'),
}
MAIN_NEW = ['G235_F_age', 'G235_M_age', 'G255_MSC_age']
# 既存の contrast（列）: label -> (表示名, 軸, DE csv, NES の出どころ)
REF = {
    'H_PostMeno': ('Human ♀ bone, postmeno. vs premeno. (GSE141595)', 'human aging', P / 'human/de/H_PostMeno.csv', 'human'),
    'H_Old': ('Human ♀ bone, old vs young (GSE72815)', 'human aging', P / 'human/de/H_Old.csv', 'human'),
    'Age_2.5y': ('♂ cortical bone 30 vs 2 mo (Kaya 2022; current 4A)', 'mouse aging', P / 'de/Age_2.5y.csv', 'chemo2'),
    'TMS_F_old': ('♀ bone 12–21 vs 1–6 mo (TMS GSE132040)', 'mouse aging', P / 'aging_f/de/TMS_F_old.csv', 'aging_f'),
    'TMS_F_old_comp': ('♀ bone TMS, comp. adj.', 'mouse aging', P / 'aging_f/de/TMS_F_old_comp.csv', 'aging_f'),
    'sc210584_CAR': ('♀ CAR 12–14 vs 2–4 mo (GSE210584)', 'mouse aging', P / 'aging_f/de/sc210584_CAR.csv', 'aging_f'),
    'Age_LepRPa': ('LepR+PDGFRa+ 24 vs 2–3 mo (GSE130885)', 'mouse aging', P / 'niche/de/Age_LepRPa.csv', 'aging_f'),
    'Age_CAR': ('♂ CAR 24 vs 2 mo (GSE133922)', 'mouse aging', P / 'niche/de/Age_CAR.csv', 'aging_f'),
    'A_OVX_cort': ('♀ cortical bone OVX vs sham (GSE222752)', 'mouse OVX', P / 'de/A_OVX_cort.csv', 'chemo2'),
    'B_Dox': ('♂ bone Dox ×1, 30 d (GSE205237)', 'Dox', P / 'de/B_Dox.csv', 'chemo2'),
    'B_Dox_compadj': ('♂ bone Dox ×1, comp. adj.', 'Dox', P / 'chemo2/de/B_Dox_compadj.csv', 'chemo2'),
    'B3_Dox4x_MSPC': ('♀ LepR+ MSPC Dox ×4 (GSE296758)', 'Dox', P / 'chemo/de/B3_Dox4x_MSPC.csv', 'chemo2'),
    'B3_Dox1x_MSPC': ('♀ LepR+ MSPC Dox ×1 (GSE296758)', 'Dox', P / 'chemo/de/B3_Dox1x_MSPC.csv', 'chemo2'),
    'B2_Dox4x_OB': ('♀ Col2.3+ OB Dox ×4 (GSE296757)', 'Dox', P / 'de/B2_Dox4x_OB.csv', 'chemo2'),
}
WIDE = {'chemo2': P / 'chemo2/nes_wide.csv', 'aging_f': P / 'aging_f/nes_wide_all.csv'}
SHORT = {
    'H_PostMeno': 'Human ♀ bone\npostmeno.', 'H_Old': 'Human ♀ bone\nold', 'Age_2.5y': '♂ cortical\n30 mo (4A)',
    'TMS_F_old': '♀ bone TMS\n12–21 mo', 'TMS_F_old_comp': '♀ bone TMS\ncomp. adj.', 'sc210584_CAR': '♀ CAR\n12–14 mo',
    'Age_LepRPa': 'LepR+ (♀?)\n24 mo', 'Age_CAR': '♂ CAR\n24 mo', 'A_OVX_cort': '♀ cortical\nOVX',
    'G235_F_age': '♀ tibia 22 mo\nGSE235942', 'G235_F_age_comp': '♀ tibia 22 mo\ncomp. adj.',
    'G235_F_age_excl': '♀ tibia 22 mo\nblood/endo out', 'G235_M_age': '♂ tibia 22 mo\nGSE235942',
    'G235_M_age_comp': '♂ tibia 22 mo\ncomp. adj.', 'G255_MSC_age': '♀ LepR+ MSC\n18 mo GSE255653',
    'G255_MSC_age_excl': '♀ LepR+ MSC\nblood/endo out',
    'B_Dox': '♂ bone Dox ×1', 'B_Dox_compadj': '♂ bone Dox ×1\ncomp. adj.', 'B3_Dox4x_MSPC': '♀ LepR+ MSPC\nDox ×4',
    'B3_Dox1x_MSPC': '♀ LepR+ MSPC\nDox ×1', 'B2_Dox4x_OB': '♀ Col2.3+ OB\nDox ×4',
}


# ---------- 読み込み ----------
def load_de(path, human=False):
    d = pd.read_csv(path).dropna(subset=['gene', 'log2FoldChange'])
    d['gene'] = d.gene.astype(str)
    if human:
        d['gene'] = d.gene.map(lambda g: hgenes.RENAME.get(g, g.upper()))
    d = d[~d.gene.duplicated(keep=False)]
    return d.set_index('gene')


def nes_tables():
    """NES / FDR の wide（term × contrast）。新規 + 既存。"""
    new = pd.read_csv(OUT / 'nes_wide.csv', header=[0, 1], index_col=0)
    W, F = {k: new[('NES', k)] for k in new['NES'].columns}, {k: new[('FDR', k)] for k in new['FDR'].columns}
    wide = {k: pd.read_csv(v, header=[0, 1], index_col=0) for k, v in WIDE.items()}
    for k, (_, _, _, src) in REF.items():
        if src == 'human':
            h = pd.read_csv(P / f'summary/gsea_human/{k}.csv').set_index('Term')
            W[k], F[k] = h.NES.astype(float), h['FDR q-val'].astype(float)
        else:
            W[k], F[k] = wide[src][('NES', k)], wide[src][('FDR', k)]
    return pd.DataFrame(W), pd.DataFrame(F)


# ---------- 相関 ----------
def corr_log2fc():
    rows = []
    for a, (an, ap) in NEW.items():
        for b, (bn, ax, bp, src) in REF.items():
            hum = src == 'human'
            x, y = load_de(ap, human=hum), load_de(bp)
            g = x.index.intersection(y.index)
            r = dict(new=a, ref=b, ref_axis=ax, n_common=len(g),
                     spearman=stats.spearmanr(x.loc[g, 'log2FoldChange'], y.loc[g, 'log2FoldChange'])[0])
            sel = ((x.loc[g, 'padj'] < 0.1) | (y.loc[g, 'padj'] < 0.1)).fillna(False)
            r['n_sig_either'] = int(sel.sum())
            if sel.sum() >= 10:
                r['spearman_sig_either'] = stats.spearmanr(x.loc[g][sel].log2FoldChange, y.loc[g][sel].log2FoldChange)[0]
            both = ((x.loc[g, 'padj'] < 0.1) & (y.loc[g, 'padj'] < 0.1)).fillna(False)
            sx, sy = np.sign(x.loc[g][both].log2FoldChange), np.sign(y.loc[g][both].log2FoldChange)
            r['deg_both_same'] = int((sx == sy).sum())
            r['deg_both_opposite'] = int((sx != sy).sum())
            rows.append(r)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / 'corr_log2fc.csv', index=False)
    print(d.pivot(index='ref', columns='new', values='spearman').reindex(list(REF)).round(2).to_string())
    return d


def corr_nes(W, F):
    rows = []
    for a in list(NEW) + ['G235_F_age_excl', 'G255_MSC_age_excl', 'G235_M_age_comp']:
        for b, (_, ax, _, _) in REF.items():
            t = W[[a, b]].dropna().index
            s = stats.spearmanr(W.loc[t, a], W.loc[t, b])[0]
            sa, sb = F.loc[t, a] < 0.05, F.loc[t, b] < 0.05
            same = (sa & sb & (np.sign(W.loc[t, a]) == np.sign(W.loc[t, b]))).sum()
            opp = (sa & sb & (np.sign(W.loc[t, a]) != np.sign(W.loc[t, b]))).sum()
            rows.append(dict(new=a, ref=b, ref_axis=ax, n_sets=len(t), spearman_NES=s, both_sig_same=int(same),
                             both_sig_opposite=int(opp), sig_new=int(sa.sum()), sig_ref=int(sb.sum())))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / 'corr_nes.csv', index=False)
    print(d.pivot(index='ref', columns='new', values='spearman_NES').reindex(list(REF)).round(2).to_string())
    return d


# ---------- 図 ----------
def module_grid(W, F):
    """要約モジュール 8 行（+ 外した 4 行）を、加齢の既存列 | 新規列 | 化療列 で並べる。表現は summary_modules.py と同じ。"""
    groups = [
        ('Human aging', ['H_PostMeno', 'H_Old']),
        ('Mouse aging (existing)', ['Age_2.5y', 'TMS_F_old', 'sc210584_CAR', 'Age_LepRPa', 'Age_CAR']),
        ('New (Fig 4A plan 2 / 3)', ['G235_F_age', 'G235_F_age_comp', 'G235_F_age_excl', 'G235_M_age',
                                     'G255_MSC_age', 'G255_MSC_age_excl']),
        ('Doxorubicin', ['B_Dox_compadj', 'B3_Dox4x_MSPC', 'B3_Dox1x_MSPC', 'B2_Dox4x_OB']),
    ]
    cols = [c for _, cs in groups for c in cs]
    rows = SM.ROWS + SM.DROPPED
    recs = []
    for term, rl, rg in rows:
        for c in cols:
            v, f = (W.loc[term, c], F.loc[term, c]) if term in W.index else (np.nan, np.nan)
            call = ('n/a' if np.isnan(v) else 'up' if (f < SM.SIG and v > 0) else 'down' if f < SM.SIG
                    else 'trend_up' if (f < SM.TREND and v > 0) else 'trend_down' if f < SM.TREND else 'ns')
            recs.append(dict(module=rl, term=term, shown_in_summary=(term, rl, rg) in SM.ROWS, contrast=c, NES=v, FDR=f, call=call))
    T = pd.DataFrame(recs)
    T.to_csv(OUT / 'module_grid.csv', index=False)

    common.setup()
    nr, nc = len(rows), len(cols)
    fig = plt.figure(figsize=(common.mm(170), common.mm(105)))
    L, R, B, Tp = 0.17, 0.90, 0.03, 0.66
    ax = fig.add_axes([L, B, R - L, Tp - B])
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(nr - 0.5, -0.5)
    ax.axis('off')
    cmap, norm, pad = plt.get_cmap('RdBu_r'), Normalize(-3, 3), 0.06
    for _, x in T.iterrows():
        i = [r[1] for r in rows].index(x.module)
        j = cols.index(x.contrast)
        box, s = (j - 0.5 + pad, i - 0.5 + pad), 1 - 2 * pad
        if x.call in ('up', 'down'):
            ax.add_patch(Rectangle(box, s, s, fc=cmap(norm(x.NES)), ec='none'))
            ax.text(j, i, f'{x.NES:+.1f}'.replace('-', '−'), ha='center', va='center', fontsize=5.5,
                    color='white' if abs(x.NES) >= 2.0 else '#222222')
        elif x.call.startswith('trend'):
            ax.add_patch(Rectangle(box, s, s, fc='none', ec=cmap(norm(1.4 if x.NES > 0 else -1.4)), lw=0.45))
        elif x.call == 'n/a':
            ax.text(j, i, '–', ha='center', va='center', fontsize=6, color='#999999')
    ax.plot([-0.5, nc - 0.5], [len(SM.ROWS) - 0.5] * 2, c='black', lw=0.6, clip_on=False)
    for i, (_, rl, _) in enumerate(rows):
        ax.text(-0.62, i, rl + ('' if i < len(SM.ROWS) else ' †'), ha='right', va='center', fontsize=7)
    j0 = 0
    for k, (gname, cs) in enumerate(groups):
        if k:
            ax.plot([j0 - 0.5] * 2, [-0.5, nr - 0.5], c='black', lw=0.8, clip_on=False)
        x0, x1 = j0 - 0.42, j0 + len(cs) - 1 + 0.42
        ax.plot([x0, x0, x1, x1], [-4.25, -4.4, -4.4, -4.25], c='black', lw=0.6, clip_on=False)
        ax.text((x0 + x1) / 2, -4.5, gname, ha='center', va='bottom', fontsize=6.5, fontweight='bold')
        j0 += len(cs)
    for j, c in enumerate(cols):
        ax.text(j, -0.62, SHORT[c], ha='left', va='bottom', fontsize=5.5, rotation=60, rotation_mode='anchor',
                linespacing=1.0)
    lx = fig.add_axes([0.92, 0.38, 0.015, 0.22])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=lx)
    cb.set_ticks([-3, 0, 3])
    cb.ax.tick_params(labelsize=5.5, length=2, width=0.5)
    lx.set_title('NES\nFDR<0.05', fontsize=5.5, loc='left')
    fig.text(0.92, 0.30, 'box: FDR\n0.05–0.25\nblank: ≥0.25\n† not in the\nsummary fig.', fontsize=5.5, va='top')
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'module_grid.{ext}', dpi=common.DPI)
    plt.close(fig)
    print('saved module_grid')
    return T


def heatmap_niche(W, F):
    _, meta = GS.collection()
    catorder = ['inflammation', 'myeloid', 'osteoclast', 'senescence/SASP', 'HSC niche', 'vascular', 'ECM',
                'bone formation', 'adipo', 'proliferation', 'estrogen', 'composition?', 'hand panel']
    nm = meta[meta.niche].copy()
    nm['co'] = nm.category.map({c: i for i, c in enumerate(catorder)})
    nm = nm.sort_values(['co', 'source', 'term'])
    cols = ['H_PostMeno', 'H_Old', 'Age_2.5y', 'TMS_F_old', 'sc210584_CAR', 'Age_LepRPa', 'Age_CAR', 'A_OVX_cort',
            'G235_F_age', 'G235_F_age_comp', 'G235_F_age_excl', 'G235_M_age', 'G235_M_age_comp', 'G255_MSC_age',
            'G255_MSC_age_excl', 'B_Dox_compadj', 'B3_Dox4x_MSPC', 'B3_Dox1x_MSPC', 'B2_Dox4x_OB']
    w, f = W.reindex(index=nm.term, columns=cols), F.reindex(index=nm.term, columns=cols)
    common.setup()
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(200)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(w.values), cmap=cmap, vmin=-3, vmax=3, aspect='auto')
    for i in range(w.shape[0]):
        for j in range(w.shape[1]):
            v, q = w.iloc[i, j], f.iloc[i, j]
            if np.isnan(v):
                continue
            if q < 0.05:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=5, fontweight='bold',
                        color='white' if abs(v) > 2 else 'black')
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, lw=0.7, ec='black'))
            elif q < 0.25:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=5, color='#333333')
    ax.set_yticks(range(len(nm)), [f'{l}  [{c}]' for l, c in zip(nm.label, nm.category)], fontsize=6)
    ax.set_xticks(range(len(cols)), [SHORT[c].replace('\n', ' ') for c in cols], rotation=60, ha='right', fontsize=5.5)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    prev = None
    for i, c in enumerate(nm.category):
        if prev is not None and c != prev:
            ax.axhline(i - 0.5, c='k', lw=0.5)
        prev = c
    for j in (2, 8, 15):
        ax.axvline(j - 0.5, c='k', lw=1.0)
    cb = fig.colorbar(im, ax=ax, shrink=0.3, pad=0.01)
    cb.set_label('NES', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title('Prerank GSEA NES (bold+box: FDR<0.05; plain: FDR<0.25). Human columns: MSigDB Hs, same set names\n'
                 '| human | mouse aging / OVX (existing) | new (GSE235942, GSE255653) | doxorubicin |', fontsize=7, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'heatmap_niche.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    print('saved heatmap_niche')


def corr_matrix(C, N):
    common.setup()
    fig, axes = plt.subplots(1, 2, figsize=(common.mm(170), common.mm(85)))
    for ax, d, val, title in [(axes[0], C, 'spearman', 'log$_2$FC Spearman (all shared genes)'),
                              (axes[1], N, 'spearman_NES', 'NES Spearman (72 gene sets)')]:
        m = d[d.new.isin(MAIN_NEW)].pivot(index='ref', columns='new', values=val).reindex(index=list(REF), columns=MAIN_NEW)
        im = ax.imshow(m.values, cmap='RdBu_r', vmin=-0.8, vmax=0.8, aspect='auto')
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                v = m.iloc[i, j]
                ax.text(j, i, f'{v:+.2f}'.replace('-', '−'), ha='center', va='center', fontsize=6,
                        color='white' if abs(v) > 0.55 else 'black')
        ax.set_xticks(range(len(MAIN_NEW)), [SHORT[c] for c in MAIN_NEW], fontsize=6)
        ax.set_yticks(range(len(REF)), [REF[k][0] if ax is axes[0] else '' for k in REF], fontsize=6)
        for y in (1.5, 7.5, 8.5):
            ax.axhline(y, c='k', lw=0.6)
        ax.set_title(title, fontsize=7.5)
        ax.tick_params(length=0)
    fig.colorbar(im, ax=axes, shrink=0.5, pad=0.01).ax.tick_params(labelsize=6)
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'corr_matrix.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    print('saved corr_matrix')


def volcanoes():
    common.setup()
    titles = {'G235_F_age': '♀ tibia 22 vs 6 mo (GSE235942)', 'G235_M_age': '♂ tibia 22 vs 6 mo (GSE235942)',
              'G255_MSC_age': '♀ LepR+ MSC 18 vs 3 mo (GSE255653)'}
    fig, axes = plt.subplots(1, 3, figsize=(common.mm(170), common.mm(62)))
    for ax, lb, pl in zip(axes, MAIN_NEW, 'ABC'):
        d = pd.read_csv(OUT / f'de/{lb}.csv').dropna(subset=['log2FoldChange', 'pvalue'])
        V.plot(ax, d, titles[lb])
        ax.title.set_fontsize(7.5)
        common.panel_label(ax, pl, dx=-0.22)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'volcano_aging2.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    rows = []
    for lb in ['G235_F_age', 'G235_F_age_comp', 'G235_M_age', 'G235_M_age_comp', 'G235_M6_HLU', 'G235_M22_HLU',
               'G255_MSC_age', 'G255_MSC_KO']:
        d = pd.read_csv(OUT / f'de/{lb}.csv').dropna(subset=['log2FoldChange', 'pvalue'])
        for pa, lf in V.THRESHOLDS:
            s = (d.padj < pa) & (d.log2FoldChange.abs() >= lf)
            rows.append(dict(label=lb, padj=pa, lfc=lf, up=int((s & (d.log2FoldChange > 0)).sum()),
                             down=int((s & (d.log2FoldChange < 0)).sum()), tested=len(d)))
    pd.DataFrame(rows).to_csv(OUT / 'deg_counts.csv', index=False)
    print('saved volcano_aging2, deg_counts')


PANEL = list(dict.fromkeys(INFLAM + APOP + ['Kitl', 'Lepr', 'Vcam1', 'Igf1', 'Spp1', 'Adipoq', 'Pparg', 'Bglap', 'Sp7',
                                            'Dmp1', 'Col1a1', 'Esr1', 'Greb1', 'Ly6a', 'Cd74']))


def panel():
    cols = ['H_PostMeno', 'H_Old', 'Age_2.5y', 'TMS_F_old', 'sc210584_CAR', 'Age_LepRPa', 'Age_CAR',
            'G235_F_age', 'G235_M_age', 'G255_MSC_age', 'B_Dox_compadj', 'B3_Dox4x_MSPC', 'B2_Dox4x_OB']
    L, Pv = {}, {}
    for c in cols:
        path = NEW[c][1] if c in NEW else REF[c][2]
        d = load_de(path)
        idx = [hgenes.RENAME.get(g, g.upper()) for g in PANEL] if c.startswith('H_') else PANEL
        L[c] = pd.Series(d.log2FoldChange.reindex(idx).values, index=PANEL)
        Pv[c] = pd.Series(d.padj.reindex(idx).values, index=PANEL)
    L, Pv = pd.DataFrame(L), pd.DataFrame(Pv)
    out = pd.concat({'log2FC': L, 'padj': Pv}, axis=1)
    out.index.name = 'gene'
    out.to_csv(OUT / 'panel_genes.csv')
    common.setup()
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(170)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(L.values), cmap=cmap, vmin=-2, vmax=2, aspect='auto')
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            v, p = L.iloc[i, j], Pv.iloc[i, j]
            if np.isnan(v):
                continue
            star = '**' if p < 0.05 else ('*' if p < 0.1 else '')
            ax.text(j, i, f'{v:.1f}{star}', ha='center', va='center', fontsize=4.8,
                    color='white' if abs(v) > 1.5 else 'black', fontweight='bold' if star else 'normal')
    ax.set_yticks(range(len(PANEL)), [f'$\\it{{{g}}}$' for g in PANEL], fontsize=6)
    ax.set_xticks(range(len(cols)), [SHORT[c].replace('\n', ' ') for c in cols], rotation=60, ha='right', fontsize=5.5)
    for j in (2, 7, 10):
        ax.axvline(j - 0.5, c='k', lw=1.0)
    cb = fig.colorbar(im, ax=ax, shrink=0.3, pad=0.01)
    cb.set_label('log$_2$FC', fontsize=6.5)
    ax.set_title('Panel genes (** padj<0.05, * <0.1; grey = filtered / not detected)', fontsize=7.5, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'heatmap_panel.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    print('saved heatmap_panel')


CONTAM_MARKERS = {
    'plasma / B': ['Igkc', 'Ighm', 'Jchain', 'Cd79a', 'Cd79b', 'Ms4a1'],
    'mast / basophil': ['Cpa3', 'Fcer1a'],
    'leukocyte': ['Ptprc', 'Mpo', 'Lyz2'],
    'endothelial': ['Pecam1', 'Cdh5', 'Emcn'],
    'cytokine': ['Il1b', 'Tnf', 'Il6', 'Cxcl2'],
    'stromal': ['Lepr', 'Cxcl12', 'Kitl', 'Adipoq', 'Bglap', 'Sp7', 'Ly6a', 'Cd74'],
}
COMP_MARKERS = {'muscle': ['Acta1', 'Ckm', 'Myh4'], 'blood': ['Hbb-bs', 'Hba-a1'], 'leukocyte': ['Ptprc', 'S100a8', 'Lyz2'],
                'osteocyte': ['Sost', 'Dmp1'], 'cytokine': ['Il1b', 'Tnf', 'Il6'], 'p16 / p21': ['Cdkn2a', 'Cdkn1a']}


def sample_scores():
    """サンプルごとのモジュールスコア（群平均の差が何例で支えられているかを見る）と、混入・組成マーカーの log2 CPM。
    スコア = set の遺伝子（平均 CPM >= 5、混入しうる 1,904 遺伝子は除く）の log2(正規化 count + 1)（median-of-ratios）を
    データセット内で z 化して平均。マーカーの図は log2(CPM+1)。"""
    sets, _ = GS.collection()
    drop = set(pd.read_csv(OUT / 'ref/contam_genes.csv', index_col=0).index.str.upper())
    rows = SM.ROWS + SM.DROPPED
    specs = [('GSE235942', ['F6_ctrl', 'F22_ctrl', 'M6_ctrl', 'M22_ctrl'], COMP_MARKERS),
             ('GSE255653', ['Young_Control', 'Aged_Control', 'Aged_TSPKO'], CONTAM_MARKERS)]
    common.setup()
    fig, axes = plt.subplots(2, 2, figsize=(common.mm(170), common.mm(185)),
                             gridspec_kw=dict(width_ratios=[1, 1], height_ratios=[1, 1.6]))
    recs = []
    for k, (gse, groups, markers) in enumerate(specs):
        c = pd.read_csv(ROOT / f'data/seq/{gse}/counts.tsv.gz', sep='\t', index_col=0)
        s = pd.read_csv(ROOT / f'data/seq/{gse}/samples.csv').set_index('sample')
        s = s[s.group.isin(groups)]
        s['o'] = s.group.map({g: i for i, g in enumerate(groups)})
        s = s.sort_values(['o']).copy()
        c = c[s.index]
        cpm = c / c.sum() * 1e6
        # サンプル間の比較は DESeq2 と同じ median-of-ratios で正規化（Ig など少数の巨大な転写物で CPM 全体が薄まるのを避ける）
        pos = c[(c > 0).all(axis=1)]
        lgm = np.log(pos).mean(axis=1)
        sf = np.exp(np.log(pos).sub(lgm, axis=0).median())
        lg = np.log2(c / sf + 1)
        lg.index = lg.index.str.upper()
        lg = lg[~lg.index.duplicated()]
        keep = (cpm.mean(axis=1) >= 5).values[~c.index.str.upper().duplicated()]
        z = lg[keep]
        z = z[~z.index.isin(drop)]
        z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1) + 1e-9, axis=0)
        S = pd.DataFrame({rl: z.reindex([g for g in sets[t] if g in z.index]).mean() for t, rl, _ in rows}).T
        S.columns = s.index
        for rl in S.index:
            for smp in S.columns:
                recs.append(dict(dataset=gse, module=rl, sample=smp, group=s.loc[smp, 'group'], score=S.loc[rl, smp],
                                 n_genes=int(z.index.isin(sets[[t for t, r, _ in rows if r == rl][0]]).sum())))
        ax = axes[0, k]
        im = ax.imshow(S.values, cmap='RdBu_r', vmin=-1.2, vmax=1.2, aspect='auto')
        ax.set_yticks(range(len(S)), [r + ('' if i < len(SM.ROWS) else ' †') for i, r in enumerate(S.index)], fontsize=6)
        xl = [x if g in x else f'{g}·{x}' for x, g in zip(S.columns, s.group)]
        ax.set_xticks(range(S.shape[1]), xl, rotation=90, fontsize=4.8)
        ax.axhline(len(SM.ROWS) - 0.5, c='k', lw=0.6)
        for j in range(1, len(s)):
            if s.group.iloc[j] != s.group.iloc[j - 1]:
                ax.axvline(j - 0.5, c='k', lw=0.8)
        ax.set_title(f'{gse}: module score per sample\n(mean z; contamination-prone genes removed)', fontsize=6.5, loc='left')
        fig.colorbar(im, ax=ax, shrink=0.6, pad=0.01).ax.tick_params(labelsize=5)
        common.panel_label(ax, 'AB'[k], dx=-0.5 if k == 0 else -0.45)
        # マーカー
        mk = [(cat, g) for cat, gs in markers.items() for g in gs if g in cpm.index]
        M = np.log2(cpm.loc[[g for _, g in mk]] + 1)
        M.to_csv(OUT / f'markers_{gse}.csv')
        ax = axes[1, k]
        im = ax.imshow(M.values, cmap='viridis', aspect='auto')
        ax.set_yticks(range(len(mk)), [f'{g}  ({cat})' for cat, g in mk], fontsize=5)
        ax.set_xticks(range(S.shape[1]), xl, rotation=90, fontsize=4.8)
        for j in range(1, len(s)):
            if s.group.iloc[j] != s.group.iloc[j - 1]:
                ax.axvline(j - 0.5, c='w', lw=0.8)
        ax.set_title(f'{gse}: marker genes, log$_2$(CPM+1)', fontsize=6.5, loc='left')
        fig.colorbar(im, ax=ax, shrink=0.6, pad=0.01).ax.tick_params(labelsize=5)
        common.panel_label(ax, 'CD'[k], dx=-0.5 if k == 0 else -0.45)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'sample_scores.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    pd.DataFrame(recs).to_csv(OUT / 'sample_scores.csv', index=False)
    print('saved sample_scores')


def main():
    sample_scores()
    W, F = nes_tables()
    C = corr_log2fc()
    N = corr_nes(W, F)
    module_grid(W, F)
    heatmap_niche(W, F)
    corr_matrix(C, N)
    volcanoes()
    panel()


if __name__ == '__main__':
    main()

"""Figure 4B 候補の横並び図と、要約モジュール図（prelim/summary）に足すための csv。

1. fig4b_modules.png|pdf: 依頼の軸（脂肪形成・ECM・骨形成・破骨支持・IFN/炎症・血管）× データ
   （bulk: GSE296757 OB, GSE296758 MSPC ×4 / ×1 | scRNA 擬似バルク GSE296752（n=1）| ヒト ALL 培養間質 GSE318370）。
   表現は要約図と同じ: FDR < 0.05 は塗り（濃さ = NES）、0.05–0.25 は枠線、それ以上は空白。
   FDR はコレクションごと（72 set の term は 72 set の中、追加 set・手製 set は追加コレクションの中）で計算した値。
   scRNA 列は順位が擬似バルクの log2FC（統計量ではない）・n=1 なので FDR は参考。
2. summary_modules_values_gse296761.csv: prelim/summary/summary_modules.py の summary_modules_values.csv と同じ列・同じ 8+4 行で、
   ここで新しく作った比較（GSE296752 の細胞種別擬似バルク、GSE318370）の NES/FDR。
   FDR は 72 set のコレクション（マウス = chemo2.genesets.collection()、ヒト = invitro.gsea_iv.collection()）で計算したもの。
   summary_modules_values_gse296761_existing.csv: 要約 csv に未収載の既存比較（B2_Dox4x_OB, B3_Dox1x_MSPC; chemo2 の値そのまま）。

  uv run python -I prelim/gse296761/modules.py      # repo root から（bulk.py, sc_compare.py, human.py の後）
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/gse296761'))
import common  # noqa: E402
import panels as PN  # noqa: E402

OUT = ROOT / 'out/prelim/gse296761'
SIG, TREND = 0.05, 0.25

ROWS = [
    ('Adipogenesis', 'HALLMARK_ADIPOGENESIS'), ('Adipogenesis', 'GOBP_FAT_CELL_DIFFERENTIATION'), ('Adipogenesis', 'HAND_ADIPO'),
    ('ECM', 'NABA_CORE_MATRISOME'), ('ECM', 'NABA_COLLAGENS'), ('ECM', 'GOBP_COLLAGEN_FIBRIL_ORGANIZATION'),
    ('ECM', 'HAND_ECM_REMODEL'),
    ('Bone formation', 'GOBP_OSSIFICATION'), ('Bone formation', 'GOBP_OSTEOBLAST_DIFFERENTIATION'),
    ('Bone formation', 'GOBP_BONE_MINERALIZATION'), ('Bone formation', 'HAND_OSTEO'),
    ('Osteoclast', 'GOBP_POSITIVE_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION'), ('Osteoclast', 'GOBP_BONE_RESORPTION'),
    ('IFN / inflammation', 'HALLMARK_INTERFERON_GAMMA_RESPONSE'), ('IFN / inflammation', 'HALLMARK_INTERFERON_ALPHA_RESPONSE'),
    ('IFN / inflammation', 'HAND_IFNG'), ('IFN / inflammation', 'HALLMARK_INFLAMMATORY_RESPONSE'),
    ('IFN / inflammation', 'SAUL_SEN_MAYO'), ('IFN / inflammation', 'HALLMARK_P53_PATHWAY'),
    ('Vascular / niche', 'HALLMARK_ANGIOGENESIS'), ('Vascular / niche', 'GOBP_ARTERY_DEVELOPMENT'),
    ('Vascular / niche', 'NICHE_HSC_RETENTION'),
]
# (key, 群, 図のラベル)
COLS = [
    ('B2_Dox4x_OB', 'bulk', 'Col2.3+ OB\nDOX×4 (5/3)'),
    ('B3_Dox4x_MSPC', 'bulk', 'LepR+ MSPC\nDOX×4 (4/4)'),
    ('B3_Dox1x_MSPC', 'bulk', 'LepR+ MSPC\nDOX×1 (3/4)'),
    ('Lepr+ MSC (all)', 'sc', 'Lepr+ MSC\n(all)'),
    ('MSC-Adipo', 'sc', 'MSC-Adipo'),
    ('MSC-Osteo', 'sc', 'MSC-Osteo'),
    ('Osteoblast', 'sc', 'Osteoblast'),
    ('EC arteriolar', 'sc', 'EC\narteriolar'),
    ('EC capillary', 'sc', 'EC\ncapillary'),
    ('EC sinusoidal', 'sc', 'EC\nsinusoidal'),
    ('Pericyte', 'sc', 'Pericyte'),
    ('H318370_Remission', 'human', 'Human ALL\nstroma, post-\ninduction (5 pairs)'),
]
GROUP_TITLE = {'bulk': 'Mouse ♀ bulk, 8 wk after DOX\n(GSE296757 / GSE296758)',
               'sc': 'Mouse ♀ scRNA pseudobulk, DOX×4, 8 wk\n(GSE296752; n = 1 pooled library each)',
               'human': 'Human (GSE318370)'}


def load_long():
    rows = []
    b = pd.read_csv(OUT / 'bulk/axes_nes.csv')
    for _, r in b.iterrows():
        rows.append(dict(col=r.label, term=r.term, NES=r.NES, fdr=r.fdr, collection=r.collection))
    s = pd.read_csv(OUT / 'sc/sc_nes_long.csv')
    for _, r in s.iterrows():
        rows.append(dict(col=r.group, term=r.term, NES=r.NES, fdr=r.fdr, collection=r.collection))
    for lb in ['H318370_Remission', 'H318370_PostInduction']:
        for tag, coll in [('', 'base'), ('extra_', 'extra')]:
            h = pd.read_csv(OUT / f'human/gsea_{tag}{lb}.csv')
            for _, r in h.iterrows():
                rows.append(dict(col=lb, term=r.Term, NES=float(r.NES), fdr=float(r['FDR q-val']), collection=coll))
    L = pd.DataFrame(rows)
    # 同じ term が 2 つのコレクションに出ることは無い（base と extra は重ならない）が、念のため base を優先
    L['pri'] = L.collection.str.startswith('base').map({True: 0, False: 1})
    return L.sort_values('pri').drop_duplicates(['col', 'term'])


def call(nes, fdr):
    if np.isnan(nes):
        return 'n/a'
    if fdr < SIG:
        return 'up' if nes > 0 else 'down'
    if fdr < TREND:
        return 'trend_up' if nes > 0 else 'trend_down'
    return 'ns'


def figure(L):
    common.setup()
    lab = PN.all_labels()
    nr, nc = len(ROWS), len(COLS)
    fig = plt.figure(figsize=(common.mm(170), common.mm(150)))
    Lf, R, Bt, Tp = 0.25, 0.86, 0.02, 0.80
    ax = fig.add_axes([Lf, Bt, R - Lf, Tp - Bt])
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(nr - 0.5, -0.5)
    ax.axis('off')
    cmap = plt.get_cmap('RdBu_r')
    norm = Normalize(-3, 3)
    pad = 0.06
    M = L.set_index(['col', 'term'])
    for i, (axis, term) in enumerate(ROWS):
        for j, (key, grp, _) in enumerate(COLS):
            if (key, term) in M.index:
                nes, fdr = M.loc[(key, term), 'NES'], M.loc[(key, term), 'fdr']
            else:
                nes, fdr = np.nan, np.nan
            c = call(nes, fdr)
            box = (j - 0.5 + pad, i - 0.5 + pad)
            s = 1 - 2 * pad
            if c in ('up', 'down'):
                ax.add_patch(Rectangle(box, s, s, fc=cmap(norm(nes)), ec='none'))
                ax.text(j, i, f'{nes:+.1f}'.replace('-', '−'), ha='center', va='center', fontsize=5.5,
                        color='white' if abs(nes) >= 2.0 else '#222222')
            elif c.startswith('trend'):
                ax.add_patch(Rectangle(box, s, s, fc='none', ec=cmap(norm(1.4 if nes > 0 else -1.4)), lw=0.5))
            elif c == 'n/a':
                ax.text(j, i, '–', ha='center', va='center', fontsize=6, color='#999999')
    axes = [r[0] for r in ROWS]
    for i in range(1, nr):
        if axes[i] != axes[i - 1]:
            ax.plot([-0.5, nc - 0.5], [i - 0.5, i - 0.5], c='#888888', lw=0.4, clip_on=False)
    prev = None
    for i, (a, term) in enumerate(ROWS):
        ax.text(-0.62, i, lab.get(term, term), ha='right', va='center', fontsize=6.3)
        if a != prev:
            ax.text(-5.6, i - 0.45, a, ha='left', va='top', fontsize=6.5, fontweight='bold')
            prev = a
    grps = [c[1] for c in COLS]
    for j in range(1, nc):
        if grps[j] != grps[j - 1]:
            ax.plot([j - 0.5, j - 0.5], [-0.5, nr - 0.5], c='black', lw=0.8, clip_on=False)
    for j, c in enumerate(COLS):
        ax.text(j, -0.62, c[2], ha='center', va='bottom', fontsize=5.5, linespacing=1.05, rotation=90)
    ytop = -4.6
    for g in ['bulk', 'sc', 'human']:
        js = [j for j, c in enumerate(COLS) if c[1] == g]
        x0, x1 = js[0] - 0.42, js[-1] + 0.42
        ax.plot([x0, x0, x1, x1], [ytop + 0.2, ytop, ytop, ytop + 0.2], c='black', lw=0.6, clip_on=False)
        ax.text((x0 + x1) / 2, ytop - 0.15, GROUP_TITLE[g], ha='center', va='bottom', fontsize=6, fontweight='bold')
    cx = fig.add_axes([0.885, 0.45, 0.018, 0.2])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cx)
    cb.set_ticks([-3, 0, 3])
    cb.set_ticklabels(['−3', '0', '+3'])
    cb.ax.tick_params(labelsize=5.5, length=2, width=0.5)
    cb.outline.set_linewidth(0.4)
    cx.set_title('NES\n(FDR<0.05)', fontsize=6, loc='left', pad=3)
    kx = fig.add_axes([0.875, 0.2, 0.12, 0.2])
    kx.set_xlim(0, 1)
    kx.set_ylim(0, 1)
    kx.axis('off')
    for k, (kind, t) in enumerate([('fill', 'FDR < 0.05'), ('edge', 'FDR 0.05–0.25'), ('blank', 'FDR ≥ 0.25'),
                                   ('na', 'not tested')]):
        y = 0.9 - k * 0.2
        if kind == 'fill':
            kx.add_patch(Rectangle((0.02, y - 0.06), 0.16, 0.12, fc=cmap(norm(2.2)), ec='none'))
        elif kind == 'edge':
            kx.add_patch(Rectangle((0.02, y - 0.06), 0.16, 0.12, fc='none', ec=cmap(norm(1.4)), lw=0.5))
        elif kind == 'blank':
            kx.add_patch(Rectangle((0.02, y - 0.06), 0.16, 0.12, fc='none', ec='#DDDDDD', lw=0.4, ls=':'))
        else:
            kx.text(0.1, y, '–', ha='center', va='center', fontsize=6, color='#999999')
        kx.text(0.24, y, t, va='center', fontsize=5.5)
    kx.text(0.0, 0.05, 'scRNA: rank = pseudobulk\nlog2FC, n = 1 → FDR\nfor reference only', fontsize=5, va='bottom',
            color='#444444')
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'fig4b_modules.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)


def summary_csv(L):
    from sc_compare import GROUPS, KEY
    rows, rows_old = [], []
    sc_cols = [(g, KEY[g], f'Mouse ♀ {g} (scRNA pseudobulk, n=1) Dox ×4 (8 wk)') for g in GROUPS]
    hum_cols = [('H318370_Remission', 'H318370_Remission', 'Human ALL BM stroma (cultured) post-induction remission vs diagnosis (5 pairs)'),
                ('H318370_PostInduction', 'H318370_PostInduction', 'Human ALL BM stroma (cultured) post-induction vs diagnosis (7 pairs)')]
    base_sc = pd.read_csv(OUT / 'sc/sc_nes_long.csv')
    base_sc = base_sc[base_sc.collection == 'base']
    for term, rl, rg in PN.SUMMARY_ROWS + PN.SUMMARY_DROPPED:
        shown = (term, rl, rg) in PN.SUMMARY_ROWS
        for g, key, cl in sc_cols:
            x = base_sc[(base_sc.group == g) & (base_sc.term == term)]
            nes, fdr = (float(x.NES.iloc[0]), float(x.fdr.iloc[0])) if len(x) else (np.nan, np.nan)
            rows.append(dict(shown=shown, module=rl, term=term, row_group=rg, contrast=key, dataset='GSE296752',
                             column_label=cl, group='chemo', NES=nes, FDR=fdr, call=call(nes, fdr),
                             source=f'out/prelim/gse296761/sc/gsea_{key}.csv'))
        for lb, key, cl in hum_cols:
            h = pd.read_csv(OUT / f'human/gsea_{lb}.csv').set_index('Term')
            nes, fdr = (float(h.loc[term, 'NES']), float(h.loc[term, 'FDR q-val'])) if term in h.index else (np.nan, np.nan)
            rows.append(dict(shown=shown, module=rl, term=term, row_group=rg, contrast=key, dataset='GSE318370',
                             column_label=cl, group='chemo', NES=nes, FDR=fdr, call=call(nes, fdr),
                             source=f'out/prelim/gse296761/human/gsea_{lb}.csv'))
        w = pd.read_csv(ROOT / 'out/prelim/chemo2/nes_wide.csv', header=[0, 1], index_col=0)
        for key, gse, cl in [('B2_Dox4x_OB', 'GSE296757', 'Mouse ♀ Col2.3+ osteoblasts Dox ×4 (8 wk)'),
                             ('B3_Dox1x_MSPC', 'GSE296758', 'Mouse ♀ LepR+ MSPC Dox ×1 (8 wk)')]:
            ok = term in w.index
            nes = float(w.loc[term, ('NES', key)]) if ok else np.nan
            fdr = float(w.loc[term, ('FDR', key)]) if ok else np.nan
            rows_old.append(dict(shown=shown, module=rl, term=term, row_group=rg, contrast=key, dataset=gse,
                                 column_label=cl, group='chemo', NES=nes, FDR=fdr, call=call(nes, fdr),
                                 source='out/prelim/chemo2/nes_wide.csv'))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / 'summary_modules_values_gse296761.csv', index=False)
    pd.DataFrame(rows_old).to_csv(OUT / 'summary_modules_values_gse296761_existing.csv', index=False)
    pd.set_option('display.width', 250)
    print(T[T.shown].pivot(index='module', columns='contrast', values='call').to_string())


def main():
    L = load_long()
    figure(L)
    summary_csv(L)


if __name__ == '__main__':
    main()

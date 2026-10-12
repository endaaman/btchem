"""マウス × ヒトの突き合わせ（log2FC の相関・共通 DEG・Hallmark NES）と、ヒト主 3 比較のパネル横並び図。

マウス → ヒトは遺伝子シンボルの大文字化で近似（パネル内の Fasl→FASLG, Trp53→TP53 だけ手当て）。
大文字化後に重複するマウス遺伝子は落とす。1:1 でない遺伝子（H2-*, Hbb-*, Ly6*, Olfr, Gm*, 多重遺伝子族）は
大文字化で一致しないか、一致しても対応がずれる（docs/prelim_human.md 参照）。

  uv run python -I prelim/human/compare_h.py
出力: out/prelim/human/compare/
  pairs.csv              ペアごとの n・Spearman/Pearson・共通 DEG 数（同方向 / 逆方向）
  shared_<pair>.csv      共通 DEG（lenient 以上）の遺伝子リスト
  scatter.png|pdf        log2FC 散布
  hallmark_nes.png|pdf, hallmark_nes.csv   マウス主要 + ヒトの NES
  panel_main3.png|pdf    ヒト主 3 比較のパネル log2FC ヒートマップ
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
import hgenes  # noqa: E402
from hcontrasts import HC, OUT, MAIN3  # noqa: E402

MDE = Path(common.OUT_DIR) / 'de'
MEN = Path(common.OUT_DIR) / 'enrich'
HDE = OUT / 'de'
HEN = OUT / 'enrich'
CO = OUT / 'compare'
THR = {
    'strict': lambda d: (d.padj < 0.05) & (d.log2FoldChange.abs() >= 1),
    'lenient': lambda d: (d.padj < 0.1) & (d.log2FoldChange.abs() >= 0.585),
    'nominal': lambda d: (d.pvalue < 0.01) & (d.log2FoldChange.abs() >= 0.585),
}
# (mouse, human) — 先頭 4 つが主
PAIRS = [
    ('Age_2.5y', 'H_PostMeno'), ('Age_2.5y', 'H_Old'),
    ('A_OVX_cort', 'H_PostMeno'), ('A_OVX_cort', 'H_Old'),
    ('Age_2y', 'H_PostMeno'), ('Age_2y', 'H_Old'),
    ('A_OVX_trab', 'H_PostMeno'), ('A_OVX_trab', 'H_Old'),
    ('A_ERT_cort', 'H_OldE'),
    ('C_DEX', 'H_GIOP'), ('C_DEX', 'H_GIOP_bgadj'), ('C_DEX', 'H_GIOP_marrowadj'),
    ('C_DEX_bgadj', 'H_GIOP_bgadj'),
    ('H_Old', 'H_PostMeno'),   # ヒト同士（独立データでの再現性）
    ('H_Dmab', 'H_Dmab_CB'),
]


def load(lb):
    if lb.startswith('H_'):
        d = pd.read_csv(HDE / f'{lb}.csv')
    else:
        d = pd.read_csv(MDE / f'{lb}.csv')
        d['gene'] = d.gene.astype(str).map(lambda g: hgenes.RENAME.get(g, g.upper()))
    d = d.dropna(subset=['gene', 'log2FoldChange'])
    d = d[~d.gene.duplicated(keep=False)]
    return d.set_index('gene')


def main():
    common.setup()
    CO.mkdir(parents=True, exist_ok=True)
    rows = []
    fig, axes = plt.subplots(3, 5, figsize=(common.mm(170), common.mm(110)))
    for ax, (m, h) in zip(axes.flat, PAIRS):
        a, b = load(m), load(h)
        g = a.index.intersection(b.index)
        x, y = a.loc[g, 'log2FoldChange'], b.loc[g, 'log2FoldChange']
        r = {'mouse': m, 'human': h, 'n_common': len(g),
             'mouse_mapped_frac': round(len(g) / len(a), 3),
             'spearman': stats.spearmanr(x, y)[0], 'pearson': stats.pearsonr(x, y)[0]}
        # どちらかで lenient 有意な遺伝子に限った相関
        sel = (THR['lenient'](a.loc[g]) | THR['lenient'](b.loc[g])).fillna(False)
        if sel.sum() >= 10:
            r['spearman_sig_either'] = stats.spearmanr(x[sel], y[sel])[0]
            r['n_sig_either'] = int(sel.sum())
        shared = []
        for tn, f in THR.items():
            sa, sb = f(a.loc[g]).fillna(False), f(b.loc[g]).fillna(False)
            both = sa & sb
            same = both & (np.sign(x) == np.sign(y))
            r[f'{tn}_mouse'] = int(sa.sum())
            r[f'{tn}_human'] = int(sb.sum())
            r[f'{tn}_both_same'] = int(same.sum())
            r[f'{tn}_both_opp'] = int((both & ~same).sum())
            # 超幾何（同方向一致の過剰）: up と down を別に
            for sgn, nm in [(1, 'up'), (-1, 'down')]:
                A = sa & (np.sign(x) == sgn)
                B = sb & (np.sign(y) == sgn)
                k = int((A & B).sum())
                r[f'{tn}_{nm}_shared'] = k
                r[f'{tn}_{nm}_p'] = stats.hypergeom.sf(k - 1, len(g), int(A.sum()), int(B.sum())) if k else 1.0
            if tn != 'nominal':
                for gn in g[both]:
                    shared.append({'threshold': tn, 'gene': gn, 'lfc_mouse': x[gn], 'padj_mouse': a.loc[gn, 'padj'],
                                   'lfc_human': y[gn], 'padj_human': b.loc[gn, 'padj'],
                                   'concordant': bool(np.sign(x[gn]) == np.sign(y[gn]))})
        rows.append(r)
        if shared:
            sh = pd.DataFrame(shared).drop_duplicates('gene')
            sh.sort_values(['concordant', 'lfc_human'], ascending=[False, False]).to_csv(
                CO / f'shared_{m}__{h}.csv', index=False)
        ax.scatter(x, y, s=1, c='#999999', lw=0, rasterized=True)
        ax.scatter(x[sel], y[sel], s=2, c='#B03A2E', lw=0, rasterized=True)
        pn = [p for p in hgenes.ALL_GENES if p in g]
        ax.scatter(x[pn], y[pn], s=8, facecolor='none', edgecolor='#1F4E79', lw=0.5)
        ax.axhline(0, lw=0.4, c='k')
        ax.axvline(0, lw=0.4, c='k')
        lim = 4
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_title(f'{m} vs {h}\nρ={r["spearman"]:.2f} (n={len(g)})', fontsize=6)
        ax.tick_params(labelsize=6)
    for ax in axes.flat[len(PAIRS):]:
        ax.axis('off')
    fig.supxlabel('log$_2$FC (first)', fontsize=8)
    fig.supylabel('log$_2$FC (second)', fontsize=8)
    fig.tight_layout()
    common.save(fig, 'human/compare/scatter', formats=('png', 'pdf'))
    plt.close(fig)
    pr = pd.DataFrame(rows)
    pr.to_csv(CO / 'pairs.csv', index=False)
    pd.set_option('display.width', 250)
    print(pr[['mouse', 'human', 'n_common', 'mouse_mapped_frac', 'spearman', 'spearman_sig_either', 'n_sig_either',
              'strict_mouse', 'strict_human', 'strict_both_same', 'strict_both_opp',
              'lenient_both_same', 'lenient_both_opp', 'lenient_up_p', 'lenient_down_p']].round(3).to_string())

    # Hallmark NES（マウス主要 + ヒト）
    labs = [('M', 'Age_2.5y'), ('M', 'Age_2y'), ('M', 'A_OVX_cort'), ('M', 'A_OVX_trab'), ('M', 'A_ERT_cort'),
            ('M', 'C_DEX'), ('M', 'C_DEX_bgadj')] + [('H', k) for k in HC]
    nes, fdr = {}, {}
    for sp, lb in labs:
        f = (MEN if sp == 'M' else HEN) / f'gsea_{lb}.csv'
        if not f.exists():
            continue
        d = pd.read_csv(f)
        nes[lb] = d.set_index('Term').NES.astype(float)
        fdr[lb] = d.set_index('Term')['FDR q-val'].astype(float)
    N, F = pd.DataFrame(nes), pd.DataFrame(fdr)
    pd.concat({'NES': N, 'FDR': F}, axis=1).to_csv(CO / 'hallmark_nes.csv')
    order = N.loc[N.index != 'SenMayo (SAUL_SEN_MAYO)'].abs().mean(axis=1).sort_values(ascending=False).index.tolist()
    if 'SenMayo (SAUL_SEN_MAYO)' in N.index:
        order = ['SenMayo (SAUL_SEN_MAYO)'] + order
    N, F = N.loc[order], F.loc[order]
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(200)))
    im = ax.imshow(N.values, cmap='RdBu_r', vmin=-3, vmax=3, aspect='auto')
    for i in range(N.shape[0]):
        for j in range(N.shape[1]):
            if F.iat[i, j] < 0.05:
                ax.text(j, i, '*', ha='center', va='center', fontsize=6)
    ax.set_xticks(range(N.shape[1]), N.columns, rotation=60, ha='right', fontsize=6)
    ax.set_yticks(range(N.shape[0]), N.index, fontsize=6)
    ax.axvline(6.5, c='k', lw=0.8)
    fig.colorbar(im, ax=ax, shrink=0.3, label='NES (* FDR<0.05)')
    fig.tight_layout()
    common.save(fig, 'human/compare/hallmark_nes', formats=('png', 'pdf'))
    plt.close(fig)

    # パネル（ヒト主 3 比較 + 感度分析）
    pt = pd.read_csv(OUT / 'panel_table.csv', index_col=0)
    cols = MAIN3 + ['H_GIOP_bgadj', 'H_GIOP_marrowadj', 'H_OldE', 'H_Dmab']
    L = pt[[f'{c}|log2FC' for c in cols]]
    P = pt[[f'{c}|padj' for c in cols]]
    fig, ax = plt.subplots(figsize=(common.mm(110), common.mm(150)))
    im = ax.imshow(L.values, cmap='RdBu_r', vmin=-1.5, vmax=1.5, aspect='auto')
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            v, p = L.iat[i, j], P.iat[i, j]
            t = '—' if pd.isna(v) else ('**' if p < 0.05 else ('*' if p < 0.1 else ''))
            ax.text(j, i, t, ha='center', va='center', fontsize=6)
    ax.set_xticks(range(len(cols)), cols, rotation=45, ha='right', fontsize=6)
    ax.set_yticks(range(L.shape[0]), [f'{g} ({pt.panel[g]})' for g in L.index], fontsize=6, fontstyle='italic')
    ax.axvline(2.5, c='k', lw=0.8)
    fig.colorbar(im, ax=ax, shrink=0.4, label='log$_2$FC (** padj<0.05, * <0.1, — not tested)')
    fig.tight_layout()
    common.save(fig, 'human/compare/panel_main3', formats=('png', 'pdf'))


if __name__ == '__main__':
    main()

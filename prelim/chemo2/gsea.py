"""化療 / 加齢 / OVX の contrast を同じ gene set（Hallmark 50 + niche 関連 MSigDB + 手製 niche パネル）で
prerank GSEA にかけ、NES/FDR を横並びにする。

  uv run python -I prelim/chemo2/gsea.py          # repo root から
順位: bulk は Wald 統計量（log2FC / lfcSE、既存 enrich.py と同じ）。scRNA（GSE289491, 各群 n=1）は
細胞種別擬似バルクの log2FC（平均 CPM >= 10 の遺伝子のみ）で、統計量ではない → FDR は参考値。
gseapy prerank, 1000 perm, seed 0, min_size 4, max_size 1000。FDR はコレクション全体（~72 set）で計算。

感度分析: B_Dox（GSE205237）を組成共変量（筋 = Acta1/Ckm/Myh4、血液 = Hbb-bs/Hba-a1 の log2 CPM 平均、z 化）
入りの design で再計算した B_Dox_compadj も回す（de.py の run と同じフィルタ・pydeseq2 設定。de.py は変えていない）。

出力（out/prelim/chemo2/）: gsea/<label>.csv, nes_long.csv, nes_wide.csv, heatmap_niche.png|pdf,
  heatmap_hallmark.png|pdf, concordance.csv, classify_aging_dex_chemo.csv, panel_genes.csv, heatmap_panel_genes.png|pdf, de/B_Dox_compadj.csv
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

OUT = ROOT / 'out/prelim/chemo2'
DE1 = ROOT / 'out/prelim/de'
DE2 = ROOT / 'out/prelim/chemo/de'
SCF = ROOT / 'out/prelim/chemo/sc/pseudobulk_all_genes.csv.gz'

# (label, 短縮名, 軸, 入力)
CONTRASTS = [
    ('Age_2.5y', 'Aging 2.5y vs 2m\ncortical ♂', 'aging', DE1 / 'Age_2.5y.csv'),
    ('Age_2y', 'Aging 2y vs 2m\ncortical ♂', 'aging', DE1 / 'Age_2y.csv'),
    ('A_OVX_cort', 'OVX vs Sham\ncortical ♀', 'OVX', DE1 / 'A_OVX_cort.csv'),
    ('A_OVX_trab', 'OVX vs Sham\ntrabecular ♀', 'OVX', DE1 / 'A_OVX_trab.csv'),
    ('B_Dox', 'Dox 1x 30d\nbone ♂ (GSE205237)', 'chemo', DE1 / 'B_Dox.csv'),
    ('B_Dox_compadj', 'Dox 1x 30d\n+composition adj.', 'chemo', OUT / 'de/B_Dox_compadj.csv'),
    ('B3_Dox4x_MSPC', 'Dox 4x 8wk\nLepR+ MSPC ♀', 'chemo', DE2 / 'B3_Dox4x_MSPC.csv'),
    ('B3_Dox1x_MSPC', 'Dox 1x 8wk\nLepR+ MSPC ♀', 'chemo', DE2 / 'B3_Dox1x_MSPC.csv'),
    ('B2_Dox4x_OB', 'Dox 4x 8wk\nCol2.3+ OB ♀', 'chemo', DE1 / 'B2_Dox4x_OB.csv'),
    ('sc_CAR', 'Dox 96h scRNA\nCAR ♀ (n=1)', 'chemo-sc', ('sc', 'CAR (LepR+)')),
    ('sc_Osteoprogenitor', 'Dox 96h scRNA\nOsteoprog. (n=1)', 'chemo-sc', ('sc', 'Osteoprogenitor')),
    ('sc_Osteoblast', 'Dox 96h scRNA\nOsteoblast (n=1)', 'chemo-sc', ('sc', 'Osteoblast')),
    ('C_DEX', 'DEX 5d\ntibia ♂', 'DEX', DE1 / 'C_DEX.csv'),
    ('C_DEX_bgadj', 'DEX 5d\n+bg adj.', 'DEX', DE1 / 'C_DEX_bgadj.csv'),
]
NAME = {a: b for a, b, *_ in CONTRASTS}
AXIS = {a: c for a, _, c, _ in CONTRASTS}
AGING_REF = 'Age_2.5y'
CHEMO = [a for a, _, c, _ in CONTRASTS if c.startswith('chemo')]
# 分類に使う化療 bulk（組成補正版・scRNA n=1 は分類に入れず、列として並べるだけ）
CHEMO_BULK = ['B_Dox', 'B3_Dox4x_MSPC', 'B3_Dox1x_MSPC', 'B2_Dox4x_OB']
DEX = 'C_DEX'
DEX_SENS = 'C_DEX_bgadj'


# ---------- 組成補正版 B_Dox ----------
def de_compadj():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    seq = ROOT / 'data/seq/GSE205237'
    counts = pd.read_csv(seq / 'counts.tsv.gz', sep='\t', index_col=0)
    s = pd.read_csv(seq / 'samples.csv', dtype=str)
    s = s[s.group.isin(['Dox', 'Vehicle'])].copy()
    s['cond'] = np.where(s.group == 'Dox', 'test', 'ref')
    s = s.set_index('sample')
    x = counts[s.index].T
    lib = x.sum(axis=1)
    cpm = x.div(lib, axis=0) * 1e6
    for k, g in [('muscle', ['Acta1', 'Ckm', 'Myh4']), ('blood', ['Hbb-bs', 'Hba-a1'])]:
        v = np.log2(cpm[g] + 1).mean(axis=1)
        s[k] = ((v - v.mean()) / v.std()).values
    print(s[['cond', 'muscle', 'blood']].round(2).to_string())
    # de.py と同じ低発現フィルタ
    cpm_cut = 10 / (lib.median() / 1e6)
    n_min = s.cond.value_counts().min()
    keep = ((cpm >= cpm_cut).sum(axis=0) >= n_min) & (x.sum(axis=0) >= 15)
    x = x.loc[:, keep]
    dds = DeseqDataSet(counts=x, metadata=s[['muscle', 'blood', 'cond']], design='~muscle + blood + cond',
                       refit_cooks=True, inference=DefaultInference(n_cpus=8), quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=['cond', 'test', 'ref'], inference=DefaultInference(n_cpus=8), quiet=True)
    st.summary()
    r = st.results_df.reset_index().rename(columns={'index': 'gene'})
    r = r.rename(columns={r.columns[0]: 'gene'})[['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj']]
    (OUT / 'de').mkdir(parents=True, exist_ok=True)
    r.to_csv(OUT / 'de/B_Dox_compadj.csv', index=False)
    print('B_Dox_compadj: padj<0.05&|lfc|>=1:', ((r.padj < 0.05) & (r.log2FoldChange.abs() >= 1)).sum())


# ---------- 順位 ----------
_sc = None


def ranking(src):
    global _sc
    if isinstance(src, tuple):
        if _sc is None:
            _sc = pd.read_csv(SCF)
        d = _sc[(_sc.celltype == src[1]) & (_sc.meanCPM >= 10)]
        return d.groupby(d.gene.str.upper()).log2FC.mean().sort_values(ascending=False)
    d = pd.read_csv(src).dropna(subset=['lfcSE', 'gene'])
    d = d[d.lfcSE > 0]
    return (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean().sort_values(ascending=False)


def lfc_table(src):
    """遺伝子単位の表用: (log2FC, padj)。sc は padj なし。"""
    if isinstance(src, tuple):
        d = _sc[(_sc.celltype == src[1]) & (_sc.meanCPM >= 10)].set_index('gene')
        return d.log2FC, pd.Series(np.nan, index=d.index)
    d = pd.read_csv(src).dropna(subset=['gene']).drop_duplicates('gene').set_index('gene')
    return d.log2FoldChange, d.padj


def run_gsea(sets):
    (OUT / 'gsea').mkdir(parents=True, exist_ok=True)
    rows = []
    for lb, _, _, src in CONTRASTS:
        rnk = ranking(src)
        r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=4, max_size=1000,
                       seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
        r = r.assign(label=lb, n_ranked=len(rnk))
        r.to_csv(OUT / 'gsea' / f'{lb}.csv', index=False)
        rows.append(r)
        print(lb, len(rnk), 'genes;', (r['FDR q-val'].astype(float) < 0.05).sum(), 'sets FDR<0.05')
    L = pd.concat(rows, ignore_index=True)
    L = L.rename(columns={'Term': 'term', 'FDR q-val': 'fdr', 'NOM p-val': 'pval', 'Lead_genes': 'lead_genes',
                          'Tag %': 'tag_pct'})
    for k in ['NES', 'ES', 'fdr', 'pval']:
        L[k] = L[k].astype(float)
    L['matched'] = L.tag_pct.astype(str).str.split('/').str[1].astype(int)
    return L[['label', 'term', 'NES', 'ES', 'pval', 'fdr', 'matched', 'lead_genes', 'n_ranked']]


# ---------- 図 ----------
def heatmap(W, F, labels, cats, fname, title, height_mm):
    common.setup()
    cols = list(W.columns)
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(height_mm)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    m = np.ma.masked_invalid(W.values)
    im = ax.imshow(m, cmap=cmap, vmin=-3, vmax=3, aspect='auto')
    for i in range(W.shape[0]):
        for j in range(W.shape[1]):
            v, f = W.iloc[i, j], F.iloc[i, j]
            if np.isnan(v):
                continue
            if f < 0.05:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=6, fontweight='bold',
                        color='white' if abs(v) > 2 else 'black')
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, lw=0.8, ec='black'))
            elif f < 0.25:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=6, color='#333333')
    ax.set_yticks(range(len(labels)), labels, fontsize=6.5)
    ax.set_xticks(range(len(cols)), [NAME[c].replace('\n', ', ') for c in cols], rotation=60, ha='right', fontsize=6)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # カテゴリ境界
    prev = None
    for i, c in enumerate(cats):
        if prev is not None and c != prev:
            ax.axhline(i - 0.5, c='k', lw=0.5)
        prev = c
    # 軸境界（aging | OVX | chemo bulk | chemo sc | DEX）
    ax_ = [AXIS[c] for c in cols]
    for j in range(1, len(cols)):
        if ax_[j] != ax_[j - 1]:
            ax.axvline(j - 0.5, c='k', lw=1.0)
    cb = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.01)
    cb.set_label('NES', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=7.5, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'{fname}.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)
    print('saved', OUT / f'{fname}.png')


def concordance(W, F, meta):
    """加齢（Age_2.5y）と各化療 contrast の向きの一致。交互作用の検定ではない（モジュールの向きの並置）。"""
    rows = []
    for t in W.index:
        a, fa = W.loc[t, AGING_REF], F.loc[t, AGING_REF]
        for c in CHEMO:
            b, fb = W.loc[t, c], F.loc[t, c]
            if np.isnan(a) or np.isnan(b):
                continue
            sa, sb = fa < 0.05, fb < 0.05
            if sa and sb:
                k = 'same (both FDR<0.05)' if np.sign(a) == np.sign(b) else 'OPPOSITE (both FDR<0.05)'
            elif sa:
                k = 'aging only'
            elif sb:
                k = 'chemo only'
            else:
                k = 'neither'
            rows.append(dict(term=t, chemo=c, NES_aging=a, FDR_aging=fa, NES_chemo=b, FDR_chemo=fb, call=k))
    d = pd.DataFrame(rows).merge(meta[['term', 'label', 'category']], on='term')
    d.to_csv(OUT / 'concordance.csv', index=False)
    s = d.pivot_table(index='chemo', columns='call', values='term', aggfunc='count', fill_value=0)
    s.to_csv(OUT / 'concordance_counts.csv')
    print(s.to_string())
    # NES ベクトルの Spearman（set は重なるので独立ではない。記述用）
    cor = {c: W[[AGING_REF, c]].dropna().corr(method='spearman').iloc[0, 1] for c in CHEMO}
    pd.Series(cor, name='spearman_NES_vs_Age2.5y').to_csv(OUT / 'concordance_nes_spearman.csv')
    print(pd.Series(cor).round(2).to_string())


def classify(W, F, meta):
    """各 gene set を「加齢 / DEX / 化療（Dox bulk 4 本）」で FDR<0.05 かつ向き付きで分類する。
    加齢 = Age_2.5y、DEX = C_DEX（C_DEX_bgadj で同じ判定になるかを dex_robust に残す）、
    化療 = CHEMO_BULK のどれかで FDR<0.05（向きは有意なものの多数決、割れたら 'mixed'）。
    交互作用・差分の検定ではない。モジュールの向きを並べた記述。
    """
    def sg(t, c):
        if t not in W.index or np.isnan(W.loc[t, c]) or F.loc[t, c] >= 0.05:
            return 0
        return int(np.sign(W.loc[t, c]))

    rows = []
    for t in W.index:
        a, d, d2 = sg(t, AGING_REF), sg(t, DEX), sg(t, DEX_SENS)
        cs = {c: sg(t, c) for c in CHEMO_BULK}
        sig = [v for v in cs.values() if v != 0]
        if not sig:
            c = 0
        elif all(v == sig[0] for v in sig):
            c = sig[0]
        else:
            c = 'mixed'
        arrow = {1: '↑', -1: '↓', 0: '·', 'mixed': '↕'}

        def rel(x, y):
            if 'mixed' in (x, y):
                return 'mixed'
            return 'same' if x == y else 'opposite'
        A, D, C = a != 0, d != 0, c != 0
        if A and D and C:
            k = f'aging∩DEX∩chemo (DEX {rel(a, d)}, chemo {rel(a, c)})'
        elif A and D:
            k = f'aging∩DEX ({rel(a, d)})'
        elif A and C:
            k = f'aging∩chemo ({rel(a, c)})'
        elif A:
            k = 'aging only'
        elif D and C:
            k = f'DEX∩chemo ({rel(d, c)})'
        elif D:
            k = 'DEX only'
        elif C:
            k = 'chemo only'
        else:
            k = 'none'
        r = dict(term=t, aging=arrow[a], DEX=arrow[d], DEX_bgadj=arrow[d2], dex_robust=(d == d2),
                 chemo=arrow[c], **{f'chemo:{x}': arrow[v] for x, v in cs.items()},
                 sc_CAR=arrow[sg(t, 'sc_CAR')], cls=k,
                 NES_aging=W.loc[t, AGING_REF], NES_DEX=W.loc[t, DEX], NES_DEX_bgadj=W.loc[t, DEX_SENS])
        rows.append(r)
    d = pd.DataFrame(rows).merge(meta[['term', 'label', 'category', 'source']], on='term')
    d = d.sort_values(['cls', 'category', 'term'])
    d.to_csv(OUT / 'classify_aging_dex_chemo.csv', index=False)
    print(d.cls.value_counts().to_string())
    return d


def panel_genes():
    order = [g for v in GS.NICHE_PANEL.values() for g in v]
    order = list(dict.fromkeys(order))
    L, P = {}, {}
    for lb, _, _, src in CONTRASTS:
        l, p = lfc_table(src)
        L[lb] = l.reindex(order)
        P[lb] = p.reindex(order)
    L, P = pd.DataFrame(L), pd.DataFrame(P)
    out = pd.concat({'log2FC': L, 'padj': P}, axis=1)
    out.index.name = 'gene'
    out.to_csv(OUT / 'panel_genes.csv')
    common.setup()
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(150)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(L.values), cmap=cmap, vmin=-2, vmax=2, aspect='auto')
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            p = P.iloc[i, j]
            if p < 0.05:
                ax.text(j, i, '**', ha='center', va='center', fontsize=7)
            elif p < 0.1:
                ax.text(j, i, '*', ha='center', va='center', fontsize=7)
    sub = {g: k for k, v in GS.NICHE_PANEL.items() for g in v}
    ax.set_yticks(range(len(order)), [f'$\\it{{{g}}}$  ' + sub[g].replace('NICHE_', '').lower() for g in order],
                  fontsize=6.5)
    ax.set_xticks(range(L.shape[1]), [NAME[c].replace('\n', ', ') for c in L.columns], rotation=60, ha='right', fontsize=6)
    prev = None
    for i, g in enumerate(order):
        if prev is not None and sub[g] != prev:
            ax.axhline(i - 0.5, c='k', lw=0.5)
        prev = sub[g]
    ax_ = [AXIS[c] for c in L.columns]
    for j in range(1, len(ax_)):
        if ax_[j] != ax_[j - 1]:
            ax.axvline(j - 0.5, c='k', lw=1.0)
    cb = fig.colorbar(im, ax=ax, shrink=0.35, pad=0.01)
    cb.set_label('log$_2$FC', fontsize=7)
    ax.set_title('Niche factor panel log2FC (bulk: ** padj<0.05, * <0.1)\nscRNA n=1: no test; grey = filtered/low',
                 fontsize=7.5, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'heatmap_panel_genes.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if not (OUT / 'de/B_Dox_compadj.csv').exists():
        de_compadj()
    sets, meta = GS.collection()
    meta.to_csv(OUT / 'genesets_meta.csv', index=False)
    L = run_gsea(sets)
    L = L.merge(meta[['term', 'label', 'category', 'source', 'niche']].rename(columns={'label': 'set_label'}), on='term', how='left')
    L.to_csv(OUT / 'nes_long.csv', index=False)
    order = [c for c, *_ in CONTRASTS]
    W = L.pivot(index='term', columns='label', values='NES')[order]
    F = L.pivot(index='term', columns='label', values='fdr')[order]
    wide = pd.concat({'NES': W, 'FDR': F}, axis=1)
    wide.to_csv(OUT / 'nes_wide.csv')

    # niche 図: Hallmark niche + curated + hand をカテゴリ順
    catorder = ['inflammation', 'myeloid', 'osteoclast', 'senescence/SASP', 'HSC niche', 'vascular', 'ECM',
                'bone formation', 'adipo', 'proliferation', 'estrogen', 'composition?', 'hand panel']
    nm = meta[meta.niche].copy()
    nm['co'] = nm.category.map({c: i for i, c in enumerate(catorder)})
    nm = nm.sort_values(['co', 'source', 'term'])
    nm = nm[nm.term.isin(W.index)]
    heatmap(W.loc[nm.term], F.loc[nm.term], [f'{l}  [{c}]' for l, c in zip(nm.label, nm.category)],
            list(nm.category), 'heatmap_niche',
            'Prerank GSEA NES (bold+box: FDR<0.05; plain: FDR<0.25)\nscRNA columns: ranked by n=1 log2FC, FDR indicative only',
            205)
    hm = meta[meta.source == 'Hallmark'].sort_values('term')
    heatmap(W.loc[hm.term], F.loc[hm.term], list(hm.label), ['h'] * len(hm), 'heatmap_hallmark',
            'Hallmark 50, prerank GSEA NES (bold+box: FDR<0.05; plain: FDR<0.25)', 205)
    concordance(W, F, meta)
    classify(W, F, meta)
    panel_genes()


if __name__ == '__main__':
    main()

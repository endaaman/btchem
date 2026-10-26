"""GSE296757（Col2.3-GFP+ 骨芽細胞 bulk）と GSE296758（Lepr-Cre;tdTomato+ 骨髄 MSPC bulk）で、依頼の
「Chemotherapy による osteoblast における骨基質環境の再編（脂肪形成偏向・ECM 再構築・骨量減少）」を軸ごとに確かめる。

DE は既存のものを読むだけ（再計算しない）:
  B2_Dox4x_OB   = out/prelim/de/B2_Dox4x_OB.csv（prelim/de.py, design ~ batch + cond, DOX×4 5 / PBS 3）
  B3_Dox4x_MSPC = out/prelim/chemo/de/B3_Dox4x_MSPC.csv（prelim/chemo/bulk.py, ~ cond, 4 / 4）
  B3_Dox1x_MSPC = out/prelim/chemo/de/B3_Dox1x_MSPC.csv（同, 3 / 4）
GSEA は 72 set（chemo2 コレクション）の結果を out/prelim/chemo2/gsea/<label>.csv から読むだけ。
追加 set + 手製 set（panels.extra_collection）だけここで回す（同じ順位 = Wald 統計量、同じ prerank 設定）。

RANKL/OPG: サンプルごとに log2((Tnfsf11 + 0.5) / (Tnfrsf11b + 0.5))（counts）を出し、群差を OLS（OB は ~ batch + group）で見る。
純度: 系統マーカーの CPM をサンプルごとに並べる（骨芽細胞 bulk に間質・血球がどれだけ混ざるか）。

  uv run python -I prelim/gse296761/bulk.py      # repo root から
出力: out/prelim/gse296761/bulk/{panel_genes.csv, gsea_extra_<label>.csv, axes_nes.csv, rankl_opg.csv,
      rankl_opg_test.csv, purity_cpm.csv}, out/prelim/gse296761/bulk_panel_heatmap.png|pdf
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/gse296761'))
import common  # noqa: E402
import panels as PN  # noqa: E402

OUT = ROOT / 'out/prelim/gse296761'
B = OUT / 'bulk'
CONTRASTS = [
    ('B2_Dox4x_OB', 'GSE296757', ROOT / 'out/prelim/de/B2_Dox4x_OB.csv', 'Col2.3+ OB, DOX×4 vs PBS (5/3)', 'DOX4X'),
    ('B3_Dox4x_MSPC', 'GSE296758', ROOT / 'out/prelim/chemo/de/B3_Dox4x_MSPC.csv', 'LepR+ MSPC, DOX×4 vs PBS (4/4)', 'DOX4X'),
    ('B3_Dox1x_MSPC', 'GSE296758', ROOT / 'out/prelim/chemo/de/B3_Dox1x_MSPC.csv', 'LepR+ MSPC, DOX×1 vs PBS (3/4)', 'DOX1X'),
]
PURITY = ['Bglap', 'Col1a1', 'Sp7', 'Ifitm5', 'Dmp1', 'Sost', 'Lepr', 'Cxcl12', 'Kitl', 'Adipoq', 'Lpl', 'Pparg', 'Fabp4',
          'Ptprc', 'Hbb-bs', 'Cd79a', 'Cdh5', 'Pecam1', 'Acta2', 'Rgs5', 'Col2a1', 'Acta1', 'Ckm']


def read_de(f):
    d = pd.read_csv(f).dropna(subset=['gene']).drop_duplicates('gene').set_index('gene')
    return d


def panel_table():
    rows = []
    for lb, gse, f, title, test in CONTRASTS:
        d = read_de(f)
        c = pd.read_csv(ROOT / f'data/seq/{gse}/counts.tsv.gz', sep='\t', index_col=0)
        s = pd.read_csv(ROOT / f'data/seq/{gse}/samples.csv', dtype=str)
        s = s[s.group.isin([test, 'PBS'])]
        cpm = c[s['sample']] / c[s['sample']].sum() * 1e6
        for k, gl in PN.PANELS.items():
            for g in gl:
                r = dict(contrast=lb, panel=k, gene=g)
                if g in cpm.index:
                    r['cpm_PBS'] = cpm.loc[g, s[s.group == 'PBS']['sample']].mean()
                    r['cpm_DOX'] = cpm.loc[g, s[s.group == test]['sample']].mean()
                if g in d.index:
                    r.update(log2FC=d.loc[g, 'log2FoldChange'], lfcSE=d.loc[g, 'lfcSE'], padj=d.loc[g, 'padj'],
                             baseMean=d.loc[g, 'baseMean'])
                else:
                    r['note'] = 'filtered (low expression)' if g in cpm.index else 'not in counts'
                rows.append(r)
    P = pd.DataFrame(rows)
    P.to_csv(B / 'panel_genes.csv', index=False)
    return P


def extra_gsea():
    sets, _ = PN.extra_collection()
    rows = []
    for lb, _, f, _, _ in CONTRASTS:
        d = pd.read_csv(f).dropna(subset=['lfcSE', 'gene'])
        d = d[d.lfcSE > 0]
        rnk = (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean().sort_values(ascending=False)
        r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=4, max_size=1000,
                       seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
        r = r.assign(label=lb, n_ranked=len(rnk))
        r.to_csv(B / f'gsea_extra_{lb}.csv', index=False)
        rows.append(r.assign(collection='extra'))
        base = pd.read_csv(ROOT / f'out/prelim/chemo2/gsea/{lb}.csv').assign(collection='base (chemo2 72 sets)')
        rows.append(base)
    L = pd.concat(rows, ignore_index=True).rename(columns={'Term': 'term', 'FDR q-val': 'fdr', 'Lead_genes': 'lead_genes'})
    L['NES'] = L.NES.astype(float)
    L['fdr'] = L.fdr.astype(float)
    terms = [t for v in PN.AXES.values() for t in v]
    ax = {t: k for k, v in PN.AXES.items() for t in v}
    lab = PN.all_labels()
    A = L[L.term.isin(terms)][['label', 'collection', 'term', 'NES', 'fdr', 'lead_genes']].copy()
    A['axis'] = A.term.map(ax)
    A['set_label'] = A.term.map(lab)
    A.to_csv(B / 'axes_nes.csv', index=False)
    return L


def rankl_opg():
    rows, tests = [], []
    for lb, gse, f, title, test in CONTRASTS:
        c = pd.read_csv(ROOT / f'data/seq/{gse}/counts.tsv.gz', sep='\t', index_col=0)
        s = pd.read_csv(ROOT / f'data/seq/{gse}/samples.csv', dtype=str)
        s = s[s.group.isin([test, 'PBS'])].copy()
        x = c[s['sample']]
        cpm = x / x.sum() * 1e6
        s['rankl_cpm'] = cpm.loc['Tnfsf11'].values
        s['opg_cpm'] = cpm.loc['Tnfrsf11b'].values
        s['log2_rankl_opg'] = np.log2((x.loc['Tnfsf11'].values + 0.5) / (x.loc['Tnfrsf11b'].values + 0.5))
        s['dox'] = (s.group == test).astype(int)
        s['contrast'] = lb
        rows.append(s)
        form = 'log2_rankl_opg ~ C(batch) + dox' if 'batch' in s and s.batch.notna().all() else 'log2_rankl_opg ~ dox'
        m = smf.ols(form, data=s).fit()
        tests.append(dict(contrast=lb, model=form, n_PBS=int((s.dox == 0).sum()), n_DOX=int(s.dox.sum()),
                          mean_PBS=s[s.dox == 0].log2_rankl_opg.mean(), mean_DOX=s[s.dox == 1].log2_rankl_opg.mean(),
                          effect=m.params['dox'], se=m.bse['dox'], p=m.pvalues['dox']))
    R = pd.concat(rows)
    R.to_csv(B / 'rankl_opg.csv', index=False)
    T = pd.DataFrame(tests)
    T.to_csv(B / 'rankl_opg_test.csv', index=False)
    print(T.round(3).to_string())


def purity():
    out = []
    for gse in ['GSE296757', 'GSE296758']:
        c = pd.read_csv(ROOT / f'data/seq/{gse}/counts.tsv.gz', sep='\t', index_col=0)
        s = pd.read_csv(ROOT / f'data/seq/{gse}/samples.csv', dtype=str)
        cpm = c[s['sample']] / c[s['sample']].sum() * 1e6
        t = cpm.reindex(PURITY).T
        t.insert(0, 'group', s.group.values)
        t.insert(0, 'gse', gse)
        if 'batch' in s:
            t.insert(2, 'batch', s.batch.values)
        out.append(t)
    P = pd.concat(out)
    P.round(1).to_csv(B / 'purity_cpm.csv')


def plot_panel(P):
    common.setup()
    labs = [c[0] for c in CONTRASTS]
    titles = {c[0]: c[3] for c in CONTRASTS}
    genes, sep = [], []
    for k in ['adipo', 'ecm', 'ecm_remodel', 'osteo', 'oc_support', 'hsc_niche', 'ifng']:
        genes += PN.PANELS[k]
        sep.append(len(genes))
    M = P.pivot(index='gene', columns='contrast', values='log2FC').reindex(index=genes, columns=labs)
    Q = P.pivot(index='gene', columns='contrast', values='padj').reindex(index=genes, columns=labs)
    fig, ax = plt.subplots(figsize=(common.mm(62), common.mm(210)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(M.values), cmap=cmap, vmin=-2, vmax=2, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v, q = M.iloc[i, j], Q.iloc[i, j]
            if np.isnan(v):
                continue
            star = '**' if q < 0.05 else ('*' if q < 0.1 else '')
            ax.text(j, i, f'{v:+.1f}{star}'.replace('-', '−'), ha='center', va='center', fontsize=4.8,
                    color='white' if abs(v) > 1.4 else 'black', fontweight='bold' if q < 0.05 else 'normal')
    for s in sep[:-1]:
        ax.axhline(s - 0.5, c='k', lw=0.6)
    ax.set_yticks(range(len(genes)), genes, fontsize=5.5, fontstyle='italic')
    ax.set_xticks(range(len(labs)), [titles[x].replace(', ', ',\n') for x in labs], rotation=50, ha='right', fontsize=5.5)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, shrink=0.25, pad=0.03)
    cb.set_label('log$_2$FC (DESeq2)', fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    ax.set_title('GSE296757 / GSE296758 bulk\n** padj<0.05, * <0.1; grey = filtered', fontsize=6.5, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'bulk_panel_heatmap.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)


def main():
    B.mkdir(parents=True, exist_ok=True)
    P = panel_table()
    plot_panel(P)
    pd.set_option('display.width', 250)
    print(P.pivot(index='gene', columns='contrast', values='log2FC').round(2).to_string())
    extra_gsea()
    rankl_opg()
    purity()


if __name__ == '__main__':
    main()

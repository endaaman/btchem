"""加齢側（LepR+ / CAR）の DE。prelim/de.py の run と同じ流儀（filterByExpr 相当フィルタ、pydeseq2、MLE、shrink なし）。

  uv run python -I prelim/niche/de.py        # repo root から。先に design.py

contrast（log2FC は test / ref）:
  Age_LepRPa   GSE130885 LEPR+_PDGFRA+_CBM  OLD (24-25 mo) vs YOUNG (8-12 wk), 3/3   ← 主（LepR+ の代表。design.py の根拠）
  Age_LepRPaN  GSE130885 LEPR+_PDGFRA-_CBM  OLD vs YOUNG, 3/3                       ← 感度
  Age_CAR      GSE133922 CAR  AGE (2 y) vs ADU (2 mo), 4/4                            ← 主
  LPS_CAR      GSE133922 CAR  LPS vs ADU, 4/4（参考: 炎症刺激）
  PIC_CAR      GSE133922 CAR  poly I:C vs ADU, 4/4（参考）
GSE133922 は normalized counts を size factor で整数に戻した counts（design.py）。DESeq2 で size factor を推定し直す。
出力: out/prelim/niche/de/<label>.csv, de_summary.csv, qc_pca.png|pdf, qc_pca.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
import common  # noqa: E402

OUT = ROOT / 'out/prelim/niche'
COLS = ['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj']
MIN_COUNT, MIN_TOTAL = 10, 15  # prelim/de.py と同じ

C = {
    'Age_LepRPa': dict(gse='GSE130885', pop='LEPR+_PDGFRA+_CBM', test='OLD', ref='YOUNG'),
    'Age_LepRPaN': dict(gse='GSE130885', pop='LEPR+_PDGFRA-_CBM', test='OLD', ref='YOUNG'),
    'Age_CAR': dict(gse='GSE133922', pop='CAR', test='AGE', ref='ADU'),
    'LPS_CAR': dict(gse='GSE133922', pop='CAR', test='LPS', ref='ADU'),
    'PIC_CAR': dict(gse='GSE133922', pop='CAR', test='PIC', ref='ADU'),
}


def load(gse):
    x = pd.read_csv(OUT / f'counts/{gse}_counts.tsv.gz', sep='\t', index_col=0)
    s = pd.read_csv(OUT / 'design_samples.csv', dtype=str)
    return x, s[s.dataset == gse]


def run(label, c):
    counts, s = load(c['gse'])
    s = s[(s.population == c['pop']) & s.group.isin([c['test'], c['ref']])].copy()
    s['cond'] = np.where(s.group == c['test'], 'test', 'ref')
    s = s.set_index('sample')
    x = counts[s.index].T
    lib = x.sum(axis=1)
    cpm_cut = MIN_COUNT / (lib.median() / 1e6)
    cpm = x.div(lib, axis=0) * 1e6
    n_min = s.cond.value_counts().min()
    keep = ((cpm >= cpm_cut).sum(axis=0) >= n_min) & (x.sum(axis=0) >= MIN_TOTAL)
    x = x.loc[:, keep]
    dds = DeseqDataSet(counts=x, metadata=s[['cond']], design='~cond', refit_cooks=True,
                       inference=DefaultInference(n_cpus=8), quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=['cond', 'test', 'ref'], inference=DefaultInference(n_cpus=8), quiet=True)
    st.summary()
    r = st.results_df.reset_index()
    r = r.rename(columns={r.columns[0]: 'gene'})[COLS]
    up = ((r.padj < 0.05) & (r.log2FoldChange >= 1)).sum()
    dn = ((r.padj < 0.05) & (r.log2FoldChange <= -1)).sum()
    n = s.cond.value_counts().to_dict()
    row = dict(label=label, gse=c['gse'], population=c['pop'], test=c['test'], n_test=n.get('test'), ref=c['ref'],
               n_ref=n.get('ref'), lib_median=int(lib.median()), genes_tested=int(keep.sum()), deg_up=int(up), deg_down=int(dn))
    print(row)
    return r, row


def pca_panel(ax, gse, pops, title):
    counts, s = load(gse)
    s = s[s.population.isin(pops)]
    x = counts[s['sample']]
    lc = np.log2(x.div(x.sum(axis=0), axis=1) * 1e6 + 1)
    lc = lc[(x > 0).sum(axis=1) >= 3]
    hv = lc.var(axis=1).sort_values().index[-2000:]
    z = lc.loc[hv].sub(lc.loc[hv].mean(axis=1), axis=0).T.values
    u, sv, _ = np.linalg.svd(z, full_matrices=False)
    pc = u[:, :2] * sv[:2]
    ve = sv ** 2 / (sv ** 2).sum()
    rows = []
    mk = {'OLD': 'o', 'YOUNG': 's', 'AGE': 'o', 'ADU': 's', 'YOU': 'v', 'LPS': '^', 'PIC': 'D'}
    col = {'OLD': '#C44E52', 'AGE': '#C44E52', 'YOUNG': '#4A7BA7', 'ADU': '#4A7BA7', 'YOU': '#8FB8DE', 'LPS': '#DD8452', 'PIC': '#937860'}
    for i, (_, r) in enumerate(s.reset_index(drop=True).iterrows()):
        fill = col[r.group] if r.population == pops[0] else 'none'
        ax.scatter(pc[i, 0], pc[i, 1], marker=mk[r.group], s=18, facecolor=fill, edgecolor=col[r.group], lw=0.8)
        rows.append(dict(dataset=gse, sample=r['sample'], population=r.population, group=r.group, PC1=pc[i, 0], PC2=pc[i, 1]))
    ax.set_xlabel(f'PC1 ({ve[0] * 100:.0f}%)')
    ax.set_ylabel(f'PC2 ({ve[1] * 100:.0f}%)')
    ax.set_title(title, fontsize=8, loc='left')
    from matplotlib.lines import Line2D
    gs = list(dict.fromkeys(s.group))
    h = [Line2D([], [], marker=mk[g], ls='', mec=col[g], mfc=col[g], ms=4, label=g) for g in gs]
    if len(pops) > 1:
        h.append(Line2D([], [], marker='o', ls='', mec='k', mfc='k', ms=4, label=pops[0]))
        h.append(Line2D([], [], marker='o', ls='', mec='k', mfc='none', ms=4, label=pops[1]))
    ax.legend(handles=h, fontsize=6, frameon=False, loc='best')
    return rows


def main():
    (OUT / 'de').mkdir(parents=True, exist_ok=True)
    rows = []
    for lb, c in C.items():
        r, row = run(lb, c)
        r.to_csv(OUT / 'de' / f'{lb}.csv', index=False)
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / 'de_summary.csv', index=False)
    common.setup()
    fig, axs = plt.subplots(1, 2, figsize=(common.mm(170), common.mm(70)))
    p = pca_panel(axs[0], 'GSE130885', ['LEPR+_PDGFRA+_CBM', 'LEPR+_PDGFRA-_CBM'], 'GSE130885 LepR+ (filled: PDGFRa+, open: PDGFRa-)')
    p += pca_panel(axs[1], 'GSE133922', ['CAR'], 'GSE133922 CAR cells')
    for a, l in zip(axs, 'AB'):
        common.panel_label(a, l)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'qc_pca.{ext}', dpi=common.DPI, bbox_inches='tight')
    pd.DataFrame(p).to_csv(OUT / 'qc_pca.csv', index=False)


if __name__ == '__main__':
    main()

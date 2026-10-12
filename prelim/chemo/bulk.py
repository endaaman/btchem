"""化療 × 骨の追加 bulk データを既存パイプライン（prelim/de.py の run, volcano.plot, enrich と同じ GSEA）で解析する。

  uv run python -I prelim/chemo/bulk.py prepare      # counts.tsv.gz / samples.csv を作る
  uv run python -I prelim/chemo/bulk.py de [label..] # DE + volcano + Hallmark GSEA + パネル表
（repo root から実行。de.py が相対パス data/seq を読むため）

出力: out/prelim/chemo/de/<label>.csv, volcano/<label>.png|pdf, gsea/<label>.csv,
      out/prelim/chemo/bulk_panel.csv, bulk_deg_counts.csv, qc/<GSE>_*.csv|png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
import common  # noqa: E402
import de  # noqa: E402
import prepare  # noqa: E402
import volcano  # noqa: E402
from genes import ALL_GENES  # noqa: E402

SEQ = ROOT / 'data/seq'
OUT = ROOT / 'out/prelim/chemo'
EXTRA = ['Mki67', 'Top2a', 'Mdm2', 'Eda2r', 'Phlda3', 'Ccng1', 'Gdf15', 'Trp53inp1', 'Spp1', 'Mmp3',
         'Lepr', 'Kitl', 'Adipoq', 'Sp7', 'Bglap', 'Runx2']
PANEL = list(dict.fromkeys(ALL_GENES + EXTRA))

C = {}


def reg(label, gse, test, ref, covariates=(), title='', tissue=''):
    C[label] = dict(label=label, gse=gse, test=test, ref=ref, subset={}, covariates=list(covariates),
                    title=title, tissue=tissue)


reg('B3_Dox4x_MSPC', 'GSE296758', 'DOX4X', 'PBS', title='Dox 5 mg/kg x4 vs PBS (LepR+ MSPC)',
    tissue='LepR+ BM-MSPC, female, 8 wk after last dose')
C['B3_Dox4x_MSPC_noC4'] = dict(label='B3_Dox4x_MSPC_noC4', gse='GSE296758', test='DOX4X', ref='PBS',
                               subset={'title_ok': 'yes'}, covariates=[],
                               title='Dox x4 vs PBS, C4 excluded (sensitivity)', tissue='')
reg('B3_Dox1x_MSPC', 'GSE296758', 'DOX1X', 'PBS', title='Dox 10 mg/kg x1 vs PBS (LepR+ MSPC)',
    tissue='LepR+ BM-MSPC, female, 8 wk after dose')

# GSE324491: ヒト H9 hESC の骨芽方向分化培養に day0–7 連続曝露（IC50）、day 7 回収、各 n=3。
# counts の対照列は UNT_I–III と UNT_IA–IIIA の 6 本あるが GSM は 3 本（対応は SOFT に記載なし）→
# 主解析は処理群と同じ命名の UNT_I–III、感度分析で UNT_IA–IIIA。
# CPA（シクロホスファミド）はプロドラッグで、肝代謝なしの in vitro では活性体がほぼできない点に注意。
for drug, name in [('CPA', 'cyclophosphamide'), ('FU5', '5-FU'), ('MTX', 'methotrexate')]:
    for st in ['I', 'A']:
        lb = f'H_{drug}_hESC-OB' + ('' if st == 'I' else '_altUNT')
        C[lb] = dict(label=lb, gse='GSE324491', test=drug, ref='UNT', subset={f'set_{st}': 'yes'},
                     covariates=[], title=f'{name} vs vehicle (hESC osteogenic d7)' + ('' if st == 'I' else ', alt UNT'),
                     tissue='human H9 hESC osteogenic culture, day 0-7 exposure')


def prep_gse324491():
    g = 'GSE324491'
    x = pd.read_excel(SEQ / g / 'GSE324491_Supplemental_metadatatable.Counts.DEGs_RNAseq.xlsx',
                      sheet_name='mirtox_fullcounts', header=0)
    cols = [c for c in x.columns if c not in ('Geneid', 'Length', 'hgnc_symbol', 'description')]
    raw = pd.read_csv(SEQ / g / 'GSE324491_allSample_RNAseq_count.txt.gz', sep='\t', index_col=0)
    sym = x.set_index('Geneid').hgnc_symbol.reindex(raw.index)
    # xlsx（シンボルのある 26,579 遺伝子のみ）の counts が txt と一致することを確認
    xi = x.set_index('Geneid')
    j = raw.index.intersection(xi.index)
    print('xlsx/txt counts identical on', len(j), 'genes:', bool((xi.loc[j, cols].values == raw.loc[j, cols].values).all()))
    sym = sym.where(sym.notna() & (sym.astype(str) != 'nan'), raw.index.to_series()).astype(str)
    # パネル（マウス表記）と突き合わせるため頭文字だけ大文字に（FAS -> Fas）。ENSG はそのまま
    sym = sym.map(lambda v: v if v.startswith('ENSG') else v[0] + v[1:].lower())
    counts = prepare.to_symbol_counts(raw[cols].astype(float), sym)
    counts.columns = [c.replace('5FU', 'FU5') for c in counts.columns]
    grp = [c.split('_')[0] for c in counts.columns]
    s = pd.DataFrame({'sample': counts.columns, 'group': grp})
    s['set_I'] = np.where((s.group != 'UNT') | s['sample'].isin(['UNT_I', 'UNT_II', 'UNT_III']), 'yes', 'no')
    s['set_A'] = np.where((s.group != 'UNT') | s['sample'].isin(['UNT_IA', 'UNT_IIA', 'UNT_IIIA']), 'yes', 'no')
    s['note'] = 'GSM 対応は SOFT になし（列名で群を決定）'
    prepare.write(g, counts, s)


def prep_gse296758():
    g = 'GSE296758'
    a = pd.read_csv(SEQ / g / 'GSE296758_raw_counts.csv.gz')
    c = a.drop(columns='gene_name').astype(float)
    c.columns = [x.split('_CKDL')[0] for x in c.columns]
    fc = {x.split('_CKDL')[0]: x.split('_')[-2] + '_' + x.split('_')[-1] for x in a.columns[1:]}
    counts = prepare.to_symbol_counts(c, a.gene_name.astype(str))
    s = prepare.read_soft(g)
    s['sample'] = s.description.str.replace('Library name: ', '', regex=False).str.replace('-', '_')
    # 群: library 名（C = control, DOX1X, DOX4X）と treatment 欄で決める。
    # GSM8975978 は title が "DOX bio rep 1" だが library 名 C4・treatment "4X PBS IV" → PBS とする
    # （他の Dox 試料の title は必ず DOX1X/DOX4X を含む。下の QC で発現からも確認する）
    s['group'] = np.select([s['sample'].str.startswith('C'), s['sample'].str.startswith('DOX1X'),
                            s['sample'].str.startswith('DOX4X')], ['PBS', 'DOX1X', 'DOX4X'], 'NA')
    s['flowcell'] = s['sample'].map(fc)
    s['note'] = np.where(s.title.str.startswith('DOX bio'),
                         'title は DOX だが library 名 C4・treatment PBS → PBS として扱う', '')
    s['title_ok'] = np.where(s.note == '', 'yes', 'no')  # 感度分析用（C4 を除く）
    s = s.set_index('sample').loc[list(counts.columns)].rename_axis('sample').reset_index()
    prepare.write(g, counts, s[['sample', 'gsm', 'title', 'group', 'treatment', 'flowcell', 'title_ok', 'note']])
    qc(g)


def qc(g):
    (OUT / 'qc').mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(SEQ / g / 'counts.tsv.gz', sep='\t', index_col=0)
    s = pd.read_csv(SEQ / g / 'samples.csv', dtype=str).set_index('sample')
    cpm = counts / counts.sum() * 1e6
    mk = ['Xist', 'Ddx3y', 'Lepr', 'Cxcl12', 'Kitl', 'Adipoq', 'Sp7', 'Bglap', 'Col1a1', 'Ptprc', 'Hbb-bs',
          'Pecam1', 'Acta1', 'Cdkn1a', 'Mdm2', 'Eda2r', 'Phlda3', 'Mki67']
    t = cpm.reindex(mk).T.round(1).join(s[['group']])
    t['libsize_M'] = (counts.sum() / 1e6).round(2)
    t.to_csv(OUT / 'qc' / f'{g}_markers.csv')
    print(t.to_string())
    lg = np.log2(cpm[(cpm > 5).sum(axis=1) >= 3] + 1)
    lg = lg.loc[lg.var(axis=1).sort_values().index[-2000:]]
    z = (lg.T - lg.T.mean())
    u, sv, _ = np.linalg.svd(z.values, full_matrices=False)
    pc = u[:, :2] * sv[:2]
    ve = sv ** 2 / (sv ** 2).sum()
    common.setup()
    fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(65)))
    cols = {'PBS': '#7F7F7F', 'DOX1X': '#E39C3A', 'DOX4X': '#C44E52'}
    for i, smp in enumerate(lg.columns):
        ax.scatter(pc[i, 0], pc[i, 1], c=cols.get(s.loc[smp, 'group'], 'k'), s=20)
        ax.annotate(smp, pc[i], fontsize=5)
    ax.set_xlabel(f'PC1 ({ve[0]*100:.0f}%)')
    ax.set_ylabel(f'PC2 ({ve[1]*100:.0f}%)')
    ax.set_title(f'{g} PCA (top 2000 var genes)')
    fig.savefig(OUT / 'qc' / f'{g}_pca.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def gsea(label, d):
    import gseapy as gp
    import enrich
    d = d.dropna(subset=['lfcSE'])
    d = d[d.lfcSE > 0]
    rnk = (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean().sort_values(ascending=False)
    r = gp.prerank(rnk=rnk, gene_sets=enrich.lib('MSigDB_Hallmark_2020'), permutation_num=1000,
                   min_size=15, max_size=500, seed=0, threads=8, outdir=None, no_plot=True,
                   verbose=False).res2d
    r = r.assign(label=label).sort_values('NES', ascending=False)
    (OUT / 'gsea').mkdir(parents=True, exist_ok=True)
    r.to_csv(OUT / 'gsea' / f'{label}.csv', index=False)
    return r


def run_de(labels):
    common.setup()
    for sub in ['de', 'volcano']:
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    for lb in labels:
        r = de.run(C[lb])
        r.to_csv(OUT / 'de' / f'{lb}.csv', index=False)
        d = r.dropna(subset=['log2FoldChange', 'pvalue'])
        fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(70)))
        volcano.plot(ax, d, f"{lb}: {C[lb]['title']}")
        for ext in ('png', 'pdf'):
            fig.savefig(OUT / 'volcano' / f'{lb}.{ext}', dpi=common.DPI, bbox_inches='tight', pad_inches=0)
        plt.close(fig)
        g = gsea(lb, d)
        print(g[['Term', 'NES', 'FDR q-val']].head(6).to_string())
        print(g[['Term', 'NES', 'FDR q-val']].tail(6).to_string())
    summarize()


def summarize():
    rows, panel = [], {}
    for lb in C:
        f = OUT / 'de' / f'{lb}.csv'
        if not f.exists():
            continue
        d = pd.read_csv(f).dropna(subset=['log2FoldChange', 'pvalue'])
        for pa, lf in volcano.THRESHOLDS:
            s = (d.padj < pa) & (d.log2FoldChange.abs() >= lf)
            rows.append(dict(label=lb, padj=pa, lfc=lf, up=int((s & (d.log2FoldChange > 0)).sum()),
                             down=int((s & (d.log2FoldChange < 0)).sum()), tested=len(d)))
        di = d.set_index('gene')
        panel[(lb, 'log2FC')] = di.log2FoldChange.reindex(PANEL)
        panel[(lb, 'padj')] = di.padj.reindex(PANEL)
    pd.DataFrame(rows).to_csv(OUT / 'bulk_deg_counts.csv', index=False)
    pt = pd.DataFrame(panel)
    pt.columns = [f'{a}|{b}' for a, b in pt.columns]
    pt.index.name = 'gene'
    pt.to_csv(OUT / 'bulk_panel.csv')
    print(pd.DataFrame(rows).query('padj == 0.05 and lfc == 1.0').to_string())


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'prepare':
        prep_gse296758()
        prep_gse324491()
    elif cmd == 'de':
        run_de(sys.argv[2:] or list(C))

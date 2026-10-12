"""Enrichment: DEG (up/down 別) の ORA（Enrichr ライブラリ、背景 = その contrast で検定した全遺伝子）と
Hallmark prerank GSEA（順位 = Wald 統計量 log2FC/lfcSE）。

DEG の定義: padj < 0.05 かつ |log2FC| >= 1。DEG が 10 未満なら padj < 0.1 かつ |log2FC| >= 0.585 に緩めて
実施し、使った閾値を出力に残す。

出力: out/prelim/enrich/ora_<label>.csv, gsea_<label>.csv, summary_top.csv
"""
import sys
from pathlib import Path

import pandas as pd
import gseapy as gp

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
from contrasts import C  # noqa: E402

DE = Path(common.OUT_DIR) / 'de'
OUT = Path(common.OUT_DIR) / 'enrich'
LIBS = ['GO_Biological_Process_2023', 'KEGG_2019_Mouse', 'MSigDB_Hallmark_2020']
_lib_cache = {}


def lib(name):
    if name not in _lib_cache:
        g = gp.get_library(name=name, organism='Mouse')
        _lib_cache[name] = {k: [x.upper() for x in v] for k, v in g.items()}
    return _lib_cache[name]


def degs(d):
    for pa, lf in [(0.05, 1.0), (0.1, 0.585)]:
        s = d[(d.padj < pa) & (d.log2FoldChange.abs() >= lf)]
        if len(s) >= 10:
            return s, f'padj<{pa},|lfc|>={lf}'
    return s, f'padj<{pa},|lfc|>={lf}'


def ora(label, d):
    bg = d.gene.str.upper().unique().tolist()
    s, thr = degs(d)
    rows = []
    for direction, sub in [('up', s[s.log2FoldChange > 0]), ('down', s[s.log2FoldChange < 0])]:
        genes = sub.gene.str.upper().unique().tolist()
        if len(genes) < 5:
            print(f'  {label} {direction}: {len(genes)} genes, ORA skipped')
            continue
        for ln in LIBS:
            try:
                r = gp.enrich(gene_list=genes, gene_sets=lib(ln), background=bg, outdir=None,
                              cutoff=1.0, no_plot=True).results
            except Exception as e:  # noqa: BLE001
                print(f'  {label} {direction} {ln}: failed {e}')
                continue
            if r is None or len(r) == 0:
                continue
            r = r.assign(label=label, direction=direction, library=ln, threshold=thr,
                         n_input=len(genes))
            rows.append(r)
    if rows:
        r = pd.concat(rows, ignore_index=True)
        r = r.sort_values(['direction', 'library', 'Adjusted P-value'])
        r.to_csv(OUT / f'ora_{label}.csv', index=False)
        return r


def gsea(label, d):
    d = d.dropna(subset=['lfcSE'])
    d = d[d.lfcSE > 0]
    rnk = (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean().sort_values(ascending=False)
    r = gp.prerank(rnk=rnk, gene_sets=lib('MSigDB_Hallmark_2020'), permutation_num=1000,
                   min_size=15, max_size=500, seed=0, threads=8, outdir=None, no_plot=True,
                   verbose=False).res2d
    r = r.assign(label=label).sort_values('NES', ascending=False)
    r.to_csv(OUT / f'gsea_{label}.csv', index=False)
    return r


def main(labels):
    OUT.mkdir(parents=True, exist_ok=True)
    tops = []
    for lb in labels:
        d = pd.read_csv(DE / f'{lb}.csv').dropna(subset=['padj', 'gene'])
        print(lb)
        r = ora(lb, d)
        if r is not None:
            sig = r[r['Adjusted P-value'] < 0.05]
            for (dr, ln), g in sig.groupby(['direction', 'library']):
                for _, x in g.head(5).iterrows():
                    tops.append({'label': lb, 'kind': 'ORA', 'direction': dr, 'library': ln,
                                 'term': x.Term, 'padj': x['Adjusted P-value'],
                                 'overlap': x.Overlap, 'threshold': x.threshold})
        g = gsea(lb, d)
        g['FDR q-val'] = g['FDR q-val'].astype(float)
        for _, x in g[g['FDR q-val'] < 0.25].iterrows():
            tops.append({'label': lb, 'kind': 'GSEA', 'direction': 'up' if x.NES > 0 else 'down',
                         'library': 'MSigDB_Hallmark_2020', 'term': x.Term, 'padj': x['FDR q-val'],
                         'overlap': f'NES={float(x.NES):.2f}', 'threshold': 'FDR<0.25'})
    pd.DataFrame(tops).to_csv(OUT / 'summary_top.csv', index=False)


if __name__ == '__main__':
    main(sys.argv[1:] or [lb for lb in C if (DE / f'{lb}.csv').exists()])

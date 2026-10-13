"""問い 1: GSE205237（Dox 単回 30 日, B_Dox）の DEG 9 個（padj<0.05, |log2FC|>=1）が、所属するモジュールで
一貫した向きを伴っているか（孤立ヒットか）を見る。

  uv run python -I prelim/chemo2/dox_deg.py      # gsea.py の後に repo root から
1. B_Dox の順位（Wald 統計量）で GO BP 全体（mouse M5, size 15–500）の prerank GSEA を追加で回す（1000 perm, seed 0）。
2. DEG ごとに、それを含む set（gsea.py のコレクション + GO BP）を列挙し、
   - その set の NES/FDR、DEG の向きと NES の向きが一致するか、DEG が leading edge に入るか
   - 「同じ set の他遺伝子」の Wald 統計量の平均（DEG 自身を除く）
   を表にする。
3. 組成（筋・血液）の補正版 B_Dox_compadj でその DEG が残るか、組成スコアとのサンプル間相関も付ける。
出力: out/prelim/chemo2/dox_deg/gsea_gobp_B_Dox.csv, deg_membership.csv, deg_summary.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
import genesets as GS  # noqa: E402

OUT = ROOT / 'out/prelim/chemo2'
O = OUT / 'dox_deg'


def main():
    O.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(ROOT / 'out/prelim/de/B_Dox.csv').dropna(subset=['lfcSE', 'gene'])
    d = d[d.lfcSE > 0]
    d['stat'] = d.log2FoldChange / d.lfcSE
    d['G'] = d.gene.str.upper()
    rnk = d.groupby('G').stat.mean().sort_values(ascending=False)
    deg = d[(d.padj < 0.05) & (d.log2FoldChange.abs() >= 1)].sort_values('log2FoldChange')
    print(deg[['gene', 'log2FoldChange', 'padj']].to_string())

    f = O / 'gsea_gobp_B_Dox.csv'
    if not f.exists():
        go = GS.gobp()
        r = gp.prerank(rnk=rnk, gene_sets=go, permutation_num=1000, min_size=15, max_size=500, seed=0,
                       threads=8, outdir=None, no_plot=True, verbose=False).res2d
        r.to_csv(f, index=False)
    go_r = pd.read_csv(f).rename(columns={'Term': 'term', 'FDR q-val': 'fdr', 'Lead_genes': 'lead_genes'})
    go_r['source'] = 'GO BP (all, 15-500)'
    col_r = pd.read_csv(OUT / 'gsea/B_Dox.csv').rename(columns={'Term': 'term', 'FDR q-val': 'fdr',
                                                                 'Lead_genes': 'lead_genes'})
    col_r['source'] = 'collection (Hallmark+niche)'
    res = pd.concat([col_r, go_r], ignore_index=True)
    sets, _ = GS.collection()
    sets.update(GS.gobp())
    print('GO BP sets tested:', len(go_r), '; FDR<0.05:', (go_r.fdr < 0.05).sum())

    # 組成スコア（gsea.py と同じ定義）
    seq = ROOT / 'data/seq/GSE205237'
    counts = pd.read_csv(seq / 'counts.tsv.gz', sep='\t', index_col=0)
    s = pd.read_csv(seq / 'samples.csv', dtype=str).set_index('sample')
    x = counts[s.index]
    cpm = x / x.sum() * 1e6
    comp = {k: np.log2(cpm.loc[g] + 1).mean() for k, g in
            [('muscle', ['Acta1', 'Ckm', 'Myh4']), ('blood', ['Hbb-bs', 'Hba-a1']),
             ('macrophage', ['Adgre1', 'Csf1r', 'Mrc1']), ('Bcell', ['Cd19', 'Cd79a', 'Ms4a1'])]}
    adj = pd.read_csv(OUT / 'de/B_Dox_compadj.csv').set_index('gene')

    mem, summ = [], []
    for _, g in deg.iterrows():
        G, sgn = g.G, np.sign(g.log2FoldChange)
        rows = []
        for _, r in res.iterrows():
            if G not in sets.get(r.term, []):
                continue
            others = [y for y in sets[r.term] if y != G and y in rnk.index]
            lead = str(r.lead_genes).split(';')
            rows.append(dict(gene=g.gene, log2FC=g.log2FoldChange, source=r.source, term=r.term,
                             NES=float(r.NES), FDR=float(r.fdr), same_dir=np.sign(float(r.NES)) == sgn,
                             in_leading_edge=G in lead, mean_stat_others=np.mean(rnk[others]) if others else np.nan,
                             n_others=len(others)))
        m = pd.DataFrame(rows)
        mem.append(m)
        lv = np.log2(cpm.loc[g.gene] + 1) if g.gene in cpm.index else None
        cr = {f'r_{k}': np.corrcoef(lv, v)[0, 1] for k, v in comp.items()} if lv is not None else {}
        sig = m[m.FDR < 0.05] if len(m) else m
        coh = sig[sig.same_dir & sig.in_leading_edge] if len(sig) else sig
        summ.append(dict(gene=g.gene, log2FC=round(g.log2FoldChange, 2), padj=g.padj, n_sets=len(m),
                         n_sets_FDR05=len(sig),
                         n_FDR05_same_dir_leading=len(coh),
                         n_FDR05_opposite=int((~sig.same_dir).sum()) if len(sig) else 0,
                         top_coherent_sets='; '.join(coh.sort_values('FDR').term.head(4)) if len(coh) else '',
                         compadj_log2FC=adj.log2FoldChange.get(g.gene), compadj_padj=adj.padj.get(g.gene),
                         **{k: round(v, 2) for k, v in cr.items()}))
    pd.concat(mem).to_csv(O / 'deg_membership.csv', index=False)
    S = pd.DataFrame(summ)
    S.to_csv(O / 'deg_summary.csv', index=False)
    pd.set_option('display.width', 300)
    pd.set_option('display.max_colwidth', 120)
    print(S.to_string())


if __name__ == '__main__':
    main()

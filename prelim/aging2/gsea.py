"""Figure 4A 案2・案3 の contrast を prerank GSEA にかけ、要約モジュール図と同じ形式の csv を出す。

設定は prelim/chemo2/gsea.py・prelim/aging_f/compare.py と同じ:
  gene set = prelim/chemo2/genesets.collection()（Hallmark 50 + MSigDB M2/M5 17 + 手製 5 = 72、MSigDB 2025.1.Mm）
  順位 = Wald 統計量（log2FC / lfcSE、大文字シンボルで重複は平均）、gseapy prerank 1000 perm, seed 0, min 4, max 1000、
  FDR は 72 set 全体で計算。
GO BP 全体（M5 GO:BP, size 15–500、chemo2/dox_deg.py と同じ）も主な contrast で回す（探索用。要約図には使わない）。

感度分析（順位から遺伝子を抜いた版、DE はやり直さない）:
  <label>_excl  = out/prelim/aging2/ref/contam_genes.csv（GSE210584 で CAR より血球・内皮で 4 倍以上高い 1,904 遺伝子）を抜いた順位。
  GSE255653 の老齢群の血球・内皮の混入、GSE235942 の雌 22 か月の骨髄残存の差に対する感度を見る。

  uv run python -I prelim/aging2/gsea.py          # repo root から（先に de.py, contam_ref.py）
出力（out/prelim/aging2/）: gsea/<label>.csv, gsea/gobp_<label>.csv, nes_long.csv, nes_wide.csv,
  summary_modules_values_aging2.csv（prelim/summary/summary_modules.py の values csv と同じ列）
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
sys.path.insert(0, str(ROOT / 'prelim/summary'))
import genesets as GS  # noqa: E402
import summary_modules as SM  # noqa: E402  （ROWS / DROPPED / 閾値を借りるだけ）

OUT = ROOT / 'out/prelim/aging2'
DE = OUT / 'de'
CONTAM = OUT / 'ref/contam_genes.csv'

# (label, DE の csv 名, 抜く遺伝子, 図・表の名前, 群, データ)
CONTRASTS = [
    ('G235_F_age', 'G235_F_age', None, 'Mouse ♀ tibia (marrow-flushed) 22 vs 6 mo', 'aging', 'GSE235942'),
    ('G235_F_age_comp', 'G235_F_age_comp', None, 'Mouse ♀ tibia 22 vs 6 mo, composition adj.', 'aging (sens.)', 'GSE235942'),
    ('G235_F_age_excl', 'G235_F_age', 'contam', 'Mouse ♀ tibia 22 vs 6 mo, blood/endo genes removed', 'aging (sens.)', 'GSE235942'),
    ('G235_M_age', 'G235_M_age', None, 'Mouse ♂ tibia (marrow-flushed) 22 vs 6 mo', 'aging', 'GSE235942'),
    ('G235_M_age_comp', 'G235_M_age_comp', None, 'Mouse ♂ tibia 22 vs 6 mo, composition adj.', 'aging (sens.)', 'GSE235942'),
    ('G235_M6_HLU', 'G235_M6_HLU', None, 'Mouse ♂ tibia 6 mo, hindlimb unloading 3 wk', 'disuse (ref.)', 'GSE235942'),
    ('G235_M22_HLU', 'G235_M22_HLU', None, 'Mouse ♂ tibia 22 mo, hindlimb unloading 3 wk', 'disuse (ref.)', 'GSE235942'),
    ('G255_MSC_age', 'G255_MSC_age', None, 'Mouse ♀ LepR+ MSC 18 vs 3 mo', 'aging', 'GSE255653'),
    ('G255_MSC_age_excl', 'G255_MSC_age', 'contam', 'Mouse ♀ LepR+ MSC 18 vs 3 mo, blood/endo genes removed', 'aging (sens.)', 'GSE255653'),
    ('G255_MSC_KO', 'G255_MSC_KO', None, 'Mouse ♀ LepR+ MSC aged Thbs1 KO vs aged ctrl', 'KO (ref.)', 'GSE255653'),
]
GOBP = ['G235_F_age', 'G235_F_age_comp', 'G235_M_age', 'G255_MSC_age', 'G255_MSC_age_excl']
# 要約モジュール図に列として足す候補（主 2 本 + 同データの雄 + 感度分析）
SUMMARY_COLS = ['G235_F_age', 'G235_F_age_comp', 'G235_F_age_excl', 'G235_M_age', 'G235_M_age_comp',
                'G255_MSC_age', 'G255_MSC_age_excl']


def ranking(de, excl):
    d = pd.read_csv(DE / f'{de}.csv').dropna(subset=['lfcSE', 'gene'])
    d = d[d.lfcSE > 0]
    if excl == 'contam':
        drop = set(pd.read_csv(CONTAM, index_col=0).index)
        d = d[~d.gene.isin(drop)]
    return (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean().sort_values(ascending=False)


def prerank(rnk, sets, lb, mn, mx):
    r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=mn, max_size=mx,
                   seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
    return r.assign(label=lb, n_ranked=len(rnk))


def main():
    (OUT / 'gsea').mkdir(parents=True, exist_ok=True)
    sets, meta = GS.collection()
    rows = []
    for lb, de, excl, *_ in CONTRASTS:
        rnk = ranking(de, excl)
        r = prerank(rnk, sets, lb, 4, 1000)
        r.to_csv(OUT / 'gsea' / f'{lb}.csv', index=False)
        rows.append(r)
        print(lb, len(rnk), 'genes;', (r['FDR q-val'].astype(float) < 0.05).sum(), 'sets FDR<0.05')
    L = pd.concat(rows, ignore_index=True).rename(columns={'Term': 'term', 'FDR q-val': 'fdr', 'NOM p-val': 'pval',
                                                           'Lead_genes': 'lead_genes', 'Tag %': 'tag_pct'})
    for k in ['NES', 'ES', 'fdr', 'pval']:
        L[k] = L[k].astype(float)
    L['matched'] = L.tag_pct.astype(str).str.split('/').str[1].astype(int)
    L = L[['label', 'term', 'NES', 'ES', 'pval', 'fdr', 'matched', 'lead_genes', 'n_ranked']]
    L = L.merge(meta[['term', 'label', 'category', 'source']].rename(columns={'label': 'set_label'}), on='term', how='left')
    L.to_csv(OUT / 'nes_long.csv', index=False)
    order = [c[0] for c in CONTRASTS]
    W = L.pivot(index='term', columns='label', values='NES')[order]
    F = L.pivot(index='term', columns='label', values='fdr')[order]
    pd.concat({'NES': W, 'FDR': F}, axis=1).to_csv(OUT / 'nes_wide.csv')

    # 要約モジュール図と同じ形式の csv
    name = {c[0]: c[3] for c in CONTRASTS}
    grp = {c[0]: c[4] for c in CONTRASTS}
    ds = {c[0]: c[5] for c in CONTRASTS}
    out = []
    for term, rl, rg in SM.ROWS + SM.DROPPED:
        for key in SUMMARY_COLS:
            nes, fdr = (W.loc[term, key], F.loc[term, key]) if term in W.index else (np.nan, np.nan)
            call = ('n/a' if np.isnan(nes) else 'up' if (fdr < SM.SIG and nes > 0) else 'down' if fdr < SM.SIG
                    else 'trend_up' if (fdr < SM.TREND and nes > 0) else 'trend_down' if fdr < SM.TREND else 'ns')
            out.append(dict(shown=(term, rl, rg) in SM.ROWS, module=rl, term=term, row_group=rg, contrast=key,
                            dataset=ds[key], column_label=name[key], group=grp[key], NES=nes, FDR=fdr, call=call,
                            source='out/prelim/aging2/nes_wide.csv'))
    T = pd.DataFrame(out)
    T.to_csv(OUT / 'summary_modules_values_aging2.csv', index=False)
    print(T.pivot(index='module', columns='contrast', values='call')
          .reindex(index=[r[1] for r in SM.ROWS + SM.DROPPED], columns=SUMMARY_COLS).to_string())

    # GO BP 全体（探索用）
    go = GS.gobp()
    for lb, de, excl, *_ in CONTRASTS:
        if lb not in GOBP:
            continue
        f = OUT / 'gsea' / f'gobp_{lb}.csv'
        if f.exists():
            continue
        r = prerank(ranking(de, excl), go, lb, 15, 500)
        r.to_csv(f, index=False)
        r['FDR q-val'] = r['FDR q-val'].astype(float)
        r['NES'] = r['NES'].astype(float)
        print(lb, 'GO BP:', len(r), 'sets;', int((r['FDR q-val'] < 0.05).sum()), 'FDR<0.05 (up',
              int(((r['FDR q-val'] < 0.05) & (r.NES > 0)).sum()), ')')


if __name__ == '__main__':
    main()

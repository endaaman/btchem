"""4C の prerank GSEA と、要約モジュール図（prelim/summary/summary_modules.py）に列として足せる形の csv。

  uv run python -I prelim/fig4c_v4/gsea.py      # repo root から（de.py の後）

- 案1 Cort_MC3T3（GSE304562, マウス）: prelim/chemo2/genesets.py の collection()（Hallmark 50 + CURATED 17 + 手製 5、
  MSigDB 2025.1.Mm）で prerank。順位 = Wald 統計量（log2FC / lfcSE）、gseapy 1000 perm, seed 0, min 4, max 1000（chemo2 と同じ）。
- 案2 C_DEX / C_DEX_bgadj（GSE276391）: 既存 out/prelim/chemo2/nes_wide.csv の値をそのまま使う（再計算しない。同じ collection・同じ設定）。
- 案2 のヒト版 H_GIOP（GSE276529）と組成補正 2 版（_marrowadj, _bgadj）: 既存 out/prelim/human/de/*.csv の Wald 統計量を、
  要約図のヒト列と同じく prelim/invitro/gsea_iv.py の collection()（同じ term のヒト GMT）・同じ設定で prerank し直す
  （既存 out/prelim/human/enrich/ は Hallmark + SenMayo しか無く、要約図の 8 行が揃わないため）。キャッシュは gsea_human/。

出力（out/prelim/fig4c_v4/）:
  gsea/Cort_MC3T3.csv, gsea_human/<label>.csv, nes_long.csv（Cort_MC3T3 の全 set）
  summary_modules_cols.csv  要約図と同じ列構成（shown,module,term,row_group,contrast,dataset,column_label,group,NES,FDR,call,source）。
                            行は要約図の 8 行（shown=True）＋外した 4 行（shown=False）。
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp

ROOT = Path(__file__).resolve().parents[2]
for p in ('prelim', 'prelim/chemo2', 'prelim/invitro', 'prelim/summary'):
    sys.path.insert(0, str(ROOT / p))
import genesets as GS  # noqa: E402
import gsea_iv  # noqa: E402   collection() だけ借りる
import summary_modules as SM  # noqa: E402   ROWS / DROPPED / SIG / TREND を同じ定義で使う（出力には触らない）

OUT = ROOT / 'out/prelim/fig4c_v4'
CHEMO2 = ROOT / 'out/prelim/chemo2/nes_wide.csv'
HUMAN_DE = ROOT / 'out/prelim/human/de'
PRERANK = dict(permutation_num=1000, min_size=4, max_size=1000, seed=0, threads=8, outdir=None, no_plot=True,
               verbose=False)

# (key, 出どころ, 図のラベル, データセット, 新規計算か)
COLS = [
    ('Cort_MC3T3', 'new_mouse', 'Mouse\nMC3T3-E1 OB\ncorticosterone\n24 h (in vitro)', 'GSE304562'),
    ('C_DEX', 'chemo2', 'Mouse ♂\ntibia\nDEX 5 d', 'GSE276391'),
    ('C_DEX_bgadj', 'chemo2', 'Mouse ♂\ntibia\nDEX 5 d\n(bg-adjusted)', 'GSE276391'),
    ('H_GIOP', 'new_human', 'Human ♀\ncortical bone\nGIOP vs control', 'GSE276529'),
    ('H_GIOP_marrowadj', 'new_human', 'Human ♀\ncortical bone\nGIOP vs control\n(marrow-adjusted)', 'GSE276529'),
    ('H_GIOP_bgadj', 'new_human', 'Human ♀\ncortical bone\nGIOP vs control\n(bg-adjusted)', 'GSE276529'),
]


def wald(df):
    d = df.dropna(subset=['lfcSE', 'gene'])
    d = d[d.lfcSE > 0]
    return (d.log2FoldChange / d.lfcSE).groupby(d.gene.astype(str).str.upper()).mean().sort_values(ascending=False)


def tidy(r):
    L = r.rename(columns={'Term': 'term', 'FDR q-val': 'fdr', 'NOM p-val': 'pval', 'Lead_genes': 'lead_genes',
                          'Tag %': 'tag_pct'})
    for k in ['NES', 'ES', 'fdr', 'pval']:
        L[k] = L[k].astype(float)
    L['matched'] = L.tag_pct.astype(str).str.split('/').str[1].astype(int)
    return L


def mouse_gsea():
    f = OUT / 'gsea/Cort_MC3T3.csv'
    if not f.exists():
        f.parent.mkdir(parents=True, exist_ok=True)
        sets, meta = GS.collection()
        rnk = wald(pd.read_csv(OUT / 'de/Cort_MC3T3.csv'))
        r = gp.prerank(rnk=rnk, gene_sets=sets, **PRERANK).res2d
        r.assign(label='Cort_MC3T3', n_ranked=len(rnk)).to_csv(f, index=False)
        L = tidy(r).merge(meta.rename(columns={'label': 'set_label'}), on='term', how='left')
        L.assign(label='Cort_MC3T3', n_ranked=len(rnk))[
            ['label', 'term', 'set_label', 'category', 'source', 'NES', 'ES', 'pval', 'fdr', 'matched', 'lead_genes',
             'n_ranked']].sort_values('NES', ascending=False).to_csv(OUT / 'nes_long.csv', index=False)
        print('Cort_MC3T3:', len(rnk), 'genes ranked;', (r['FDR q-val'].astype(float) < 0.05).sum(), 'sets FDR<0.05')
    r = pd.read_csv(f)
    return r.set_index('Term')[['NES', 'FDR q-val']].astype(float).rename(columns={'FDR q-val': 'FDR'})


def human_gsea(lb, sets_cache={}):
    d = OUT / 'gsea_human'
    d.mkdir(parents=True, exist_ok=True)
    f = d / f'{lb}.csv'
    if not f.exists():
        if 'sets' not in sets_cache:
            sets_cache['sets'], _ = gsea_iv.collection()
        rnk = wald(pd.read_csv(HUMAN_DE / f'{lb}.csv'))
        r = gp.prerank(rnk=rnk, gene_sets=sets_cache['sets'], **PRERANK).res2d
        r.assign(label=lb, n_ranked=len(rnk)).to_csv(f, index=False)
        print(lb, len(rnk), 'genes ranked;', (r['FDR q-val'].astype(float) < 0.05).sum(), 'sets FDR<0.05')
    r = pd.read_csv(f)
    return r.set_index('Term')[['NES', 'FDR q-val']].astype(float).rename(columns={'FDR q-val': 'FDR'})


def table():
    w = pd.read_csv(CHEMO2, header=[0, 1], index_col=0)
    src = {}
    for key, kind, *_ in COLS:
        if kind == 'new_mouse':
            src[key] = (mouse_gsea(), f'out/prelim/fig4c_v4/gsea/{key}.csv')
        elif kind == 'new_human':
            src[key] = (human_gsea(key), f'out/prelim/fig4c_v4/gsea_human/{key}.csv')
        else:
            src[key] = (pd.DataFrame({'NES': w[('NES', key)], 'FDR': w[('FDR', key)]}).astype(float),
                        str(CHEMO2.relative_to(ROOT)))
    rows = []
    for term, rl, rg in SM.ROWS + SM.DROPPED:
        for key, kind, cl, gse in COLS:
            t, s = src[key]
            nes, fdr = (t.loc[term, 'NES'], t.loc[term, 'FDR']) if term in t.index else (np.nan, np.nan)
            call = ('n/a' if np.isnan(nes) else 'up' if (fdr < SM.SIG and nes > 0) else 'down' if fdr < SM.SIG
                    else 'trend_up' if (fdr < SM.TREND and nes > 0) else 'trend_down' if fdr < SM.TREND else 'ns')
            rows.append(dict(shown=(term, rl, rg) in SM.ROWS, module=rl, term=term, row_group=rg, contrast=key,
                             dataset=gse, column_label=cl.replace('\n', ' '), group='steroid', NES=nes, FDR=fdr,
                             call=call, source=s))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    T = table()
    T.to_csv(OUT / 'summary_modules_cols.csv', index=False)
    pd.set_option('display.width', 250)
    v = T.assign(cell=T.NES.map('{:+.2f}'.format) + ' (' + T.FDR.map('{:.3f}'.format) + ')')
    print(v.pivot(index='module', columns='contrast', values='cell')
          .reindex(index=[r[1] for r in SM.ROWS + SM.DROPPED], columns=[c[0] for c in COLS]).to_string())


if __name__ == '__main__':
    main()

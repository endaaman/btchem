"""GSE318370（GSE296761 のヒト部分: ALL 患者の骨髄吸引液を培養して得た間質, SMART-Seq v4 bulk）。
診断時 vs 寛解導入療法 4 週後のペアで、マウスと同じ軸（脂肪・ECM・骨形成・破骨支持・IFN）を見る。補助解析。

サンプル（SOFT の title と library 名）:
  ペア（診断 / 寛解）: Patient 1, 2, 3, 6, 8 → 主解析（5 ペア）
  ペア（診断 / 無効 No response）: Patient 5, 7 → 感度分析で「寛解導入後」に含める（7 ペア）
  診断 / 再発: Patient 4 → 時期が違うので除外
設計: pydeseq2 で ~ patient + time（対応あり）。低発現フィルタは prelim/de.py と同じ。
GSEA: prelim/invitro/gsea_iv.py の collection()（72 set と同じ term をヒト GMT から）と、panels.py の追加 set をヒト GMT から
（手製パネルはマウスシンボルを大文字化。H2-K1 / H2-Ab1 / Bglap2 はヒトに 1:1 の相手が無いので落とす）。順位 = Wald 統計量。

  uv run python -I prelim/gse296761/human.py      # repo root から
出力: data/seq/GSE318370/samples.csv（SOFT から作った対応表）, out/prelim/gse296761/human/{de_<label>.csv,
      gsea_<label>.csv, gsea_extra_<label>.csv, panel_genes.csv}
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
sys.path.insert(0, str(ROOT / 'prelim/invitro'))
sys.path.insert(0, str(ROOT / 'prelim/gse296761'))
import gsea_iv  # noqa: E402  （collection() / read_gmt() を借りるだけ）
import panels as PN  # noqa: E402

G = ROOT / 'data/seq/GSE318370'
OUT = ROOT / 'out/prelim/gse296761/human'
CONTRASTS = {
    'H318370_Remission': ['Remission'],                       # 5 ペア
    'H318370_PostInduction': ['Remission', 'No response'],    # 7 ペア（感度分析）
}
NOT_IN_HUMAN = {'H2-K1', 'H2-Ab1', 'Bglap2'}


def samples():
    f = G / 'samples.csv'
    if f.exists():
        return pd.read_csv(f, dtype=str)
    rows, cur = [], {}
    with gzip.open(G / 'GSE318370_family.soft.gz', 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE'):
                if cur:
                    rows.append(cur)
                cur = {'gsm': line.split(' = ')[1]}
            elif line.startswith('!Sample_title'):
                cur['title'] = line.split(' = ')[1]
            elif line.startswith('!Sample_description = Library name'):
                cur['sample'] = line.split(': ')[1]
            elif 'time point:' in line:
                cur['time'] = line.split('time point: ')[1]
    rows.append(cur)
    s = pd.DataFrame(rows)
    s['patient'] = s.title.str.extract(r'Patient (\d+)')[0].map(lambda x: f'P{x}')
    s = s[['sample', 'gsm', 'title', 'patient', 'time']]
    s.to_csv(f, index=False)
    return s


def de(label, post):
    c = pd.read_csv(G / 'GSE318370_raw_counts_msc_patient.csv.gz').set_index('gene')
    s = samples()
    pts = s[s.time.isin(post)].patient.unique()
    s = s[s.patient.isin(pts) & s.time.isin(['Diagnosis'] + post)].copy()
    s['cond'] = np.where(s.time == 'Diagnosis', 'ref', 'test')
    s = s.set_index('sample')
    x = c[s.index].T
    lib = x.sum(axis=1)
    cpm = x.div(lib, axis=0) * 1e6
    keep = ((cpm >= 10 / (lib.median() / 1e6)).sum(axis=0) >= s.cond.value_counts().min()) & (x.sum(axis=0) >= 15)
    x = x.loc[:, keep]
    dds = DeseqDataSet(counts=x, metadata=s[['patient', 'cond']], design='~patient + cond', refit_cooks=True,
                       inference=DefaultInference(n_cpus=8), quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=['cond', 'test', 'ref'], inference=DefaultInference(n_cpus=8), quiet=True)
    st.summary()
    r = st.results_df.reset_index()
    r = r.rename(columns={r.columns[0]: 'gene'})[['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj']]
    r.to_csv(OUT / f'de_{label}.csv', index=False)
    print(label, 'pairs', len(pts), 'genes', keep.sum(), 'padj<0.05&|lfc|>=1:',
          ((r.padj < 0.05) & (r.log2FoldChange.abs() >= 1)).sum(),
          'up', ((r.padj < 0.05) & (r.log2FoldChange >= 1)).sum())
    return r


def extra_sets_human():
    allh = {**gsea_iv.read_gmt('c2.all'), **gsea_iv.read_gmt('c5.go.bp')}
    sets = {k: allh[k] for k, _, _ in PN.EXTRA if k in allh}
    for k, v in PN.HAND.items():
        sets[k] = [g.upper() for g in v if g not in NOT_IN_HUMAN]
    return sets


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base, _ = gsea_iv.collection()
    extra = extra_sets_human()
    panel_rows = []
    for label, post in CONTRASTS.items():
        r = de(label, post)
        d = r.dropna(subset=['lfcSE'])
        d = d[d.lfcSE > 0]
        rnk = (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean().sort_values(ascending=False)
        for tag, sets in [('', base), ('extra_', extra)]:
            g = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=4, max_size=1000,
                           seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
            g.assign(label=label, n_ranked=len(rnk)).to_csv(OUT / f'gsea_{tag}{label}.csv', index=False)
        rr = r.set_index('gene')
        for k, gl in PN.PANELS.items():
            for gm in gl:
                if gm in NOT_IN_HUMAN:
                    continue
                gh = gm.upper()
                panel_rows.append(dict(contrast=label, panel=k, gene=gh,
                                       log2FC=rr.log2FoldChange.get(gh, np.nan), padj=rr.padj.get(gh, np.nan),
                                       baseMean=rr.baseMean.get(gh, np.nan)))
    P = pd.DataFrame(panel_rows)
    P.to_csv(OUT / 'panel_genes.csv', index=False)
    pd.set_option('display.width', 200)
    print(P.pivot(index='gene', columns='contrast', values='log2FC').round(2).to_string())


if __name__ == '__main__':
    main()

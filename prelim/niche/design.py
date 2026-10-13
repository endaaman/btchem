"""加齢側 2 データ（GSE130885, GSE133922）と化療側（GSE296758）の設計を SOFT から確定し、DE 用の counts を作る。

  uv run python -I prelim/niche/design.py      # repo root から

- GSE130885: featureCounts（RefSeq, mm9）の raw counts。36 列 = 6 集団 × 若齢/老齢 × 3。
  LepR+ に当たる集団（LEPR+_PDGFRA-_CBM / LEPR+_PDGFRA+_CBM）をマーカーと GSE296758 PBS / GSE133922 CAR との相関で確認する。
- GSE133922: DESeq2 normalized counts のみ。各列を「最小の非ゼロ値 = raw 1」とみなした size factor を掛けて整数に戻せるか確認し、
  戻せた場合は raw counts として使う（全値が整数に戻り、かつ値の最大公約数が 1 であることを条件にする）。
  列名 ADU / AGE / YOU と SOFT の年齢（2 か月 / 2 年 / 2 週）の対応は、PαS 2 か月が SOFT 上 3 匹（replicate 3 欠番）で
  ADU_PAS が 3 列しかないことから ADU = 2 か月、残る YOU = 2 週（早期生後）とした。
- 性別: SOFT に記載が無いので、発現（Xist と Y 染色体遺伝子 Ddx3y / Eif2s3y / Kdm5d / Uty）から推定して併記する。

出力（out/prelim/niche/）: design_samples.csv, design_summary.csv, sex_markers.csv,
  gse130885_population_markers.csv, gse130885_population_similarity.csv, counts/{GSE130885,GSE133922}_counts.tsv.gz
"""
import gzip
import math
import re
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SEQ = ROOT / 'data/seq'
OUT = ROOT / 'out/prelim/niche'
SEX_GENES = ['Xist', 'Ddx3y', 'Eif2s3y', 'Kdm5d', 'Uty']
MARKERS = ['Lepr', 'Cxcl12', 'Kitl', 'Adipoq', 'Angpt1', 'Vcam1', 'Pdgfra', 'Pdgfrb', 'Ly6a', 'Alcam', 'Grem1',
           'Sp7', 'Bglap', 'Col1a1', 'Alpl', 'Pecam1', 'Cdh5', 'Emcn', 'Ptprc', 'Hbb-bs']


def soft_samples(gse):
    """SOFT の Sample ブロックを dict のリストにする（title, characteristics, source, instrument）。"""
    rows, cur = [], None
    with gzip.open(SEQ / gse / f'{gse}_family.soft.gz', 'rt', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE'):
                cur = {'gsm': line.split('=')[1].strip()}
                rows.append(cur)
            elif cur is None:
                continue
            elif line.startswith('!Sample_title'):
                cur['title'] = line.split('=', 1)[1].strip()
            elif line.startswith('!Sample_characteristics_ch1'):
                k, _, v = line.split('=', 1)[1].strip().partition(':')
                cur[k.strip()] = v.strip()
            elif line.startswith('!Sample_instrument_model'):
                cur['instrument'] = line.split('=', 1)[1].strip()
            elif line.startswith('!Sample_platform_id'):
                cur['platform'] = line.split('=', 1)[1].strip()
    return pd.DataFrame(rows)


def cpm(x):
    return x.div(x.sum(axis=0), axis=1) * 1e6


# ---------- GSE130885 ----------
def gse130885():
    s = soft_samples('GSE130885')
    s['sample'] = s.title
    s['population'] = s.title.str.replace(r'_(OLD|YOUNG)_\d$', '', regex=True)
    s['group'] = s.age.str.upper().map({'OLD': 'OLD', 'YOUNG': 'YOUNG'})
    s['age_detail'] = s.group.map({'YOUNG': '8-12 wk (Series overall_design)', 'OLD': '24-25 mo (Series overall_design)'})
    x = pd.read_csv(SEQ / 'GSE130885/GSE130885_Valletta_et_al_RNAseq_stromal_featurecounts.txt.gz', sep='\t', index_col=0)
    assert set(x.columns) == set(s['sample']), 'column/title mismatch'
    x = x.groupby(level=0).sum()  # 重複シンボルは合算
    x = x[s['sample']]
    return s, x


# ---------- GSE133922 ----------
def recover_counts(norm):
    out, info = {}, []
    for c in norm.columns:
        v = norm[c].values
        sf = 1 / v[v > 0].min()
        r = v * sf
        ri = np.round(r).astype(np.int64)
        maxdev = float(np.abs(r - ri).max())
        g = reduce(math.gcd, [int(a) for a in np.unique(ri[ri > 0])[:200]])
        out[c] = ri
        info.append(dict(sample=c, size_factor=1 / sf, max_dev_from_integer=maxdev, gcd=g, total=int(ri.sum())))
    info = pd.DataFrame(info)
    ok = (info.max_dev_from_integer < 1e-3) & (info.gcd == 1)
    print('GSE133922 integer recovery: all ok =', bool(ok.all()), '; max dev =', info.max_dev_from_integer.max())
    return pd.DataFrame(out, index=norm.index), info


def gse133922():
    s = soft_samples('GSE133922')
    norm = pd.read_csv(SEQ / 'GSE133922/GSE133922_ALLsamples.normalizedCounts.txt.gz', sep='\t', index_col=0)
    x, info = recover_counts(norm)
    cols = pd.DataFrame({'sample': x.columns})
    p = cols['sample'].str.extract(r'^(ADU|AGE|YOU|LPS|PIC)_(CAR|PAS|SEC|AEC)_(\d)$')
    cols['group'], cols['population'], cols['rep_in_name'] = p[0], p[1], p[2]
    cols['age_detail'] = cols.group.map({'ADU': '2 mo', 'AGE': '2 y', 'YOU': '2 wk', 'LPS': '2 mo + LPS', 'PIC': '2 mo + poly I:C'})
    # SOFT との数の照合（列名 ↔ GSM の 1 対 1 対応は SOFT に無い）
    s['population'] = s['cell type'].map({'CXCL12-abundant reticular cells': 'CAR', 'PDGFR-alpha+ Sca-1+ mesenchymal cells': 'PAS',
                                          'sinusoidal endothelial cells': 'SEC', 'arterial endothelial cells': 'AEC'})
    s['group'] = np.select([s.agent.str.startswith('LPS'), s.agent.str.startswith('poly'), s.age == '2 years',
                            s.age == '2 months', s.age == '2 weeks'], ['LPS', 'PIC', 'AGE', 'ADU', 'YOU'], 'NA')
    a = s.groupby(['population', 'group']).size().rename('n_soft')
    b = cols.groupby(['population', 'group']).size().rename('n_columns')
    chk = pd.concat([a, b], axis=1)
    print(chk.to_string())
    assert (chk.n_soft == chk.n_columns).all(), 'SOFT と列の群サイズが合わない'
    cols = cols.merge(info, on='sample')
    return cols, x


def sex_table(name, x, s):
    c = cpm(x)
    g = [k for k in SEX_GENES if k in c.index]
    t = c.loc[g].T.round(1)
    t['dataset'] = name
    t['sample'] = t.index
    yg = [k for k in g if k != 'Xist']
    xi = t['Xist'] if 'Xist' in g else pd.Series(np.nan, index=t.index)
    ym = t[yg].max(axis=1) if yg else pd.Series(np.nan, index=t.index)

    def call(xv, yv):
        if np.isnan(yv):  # Y 遺伝子の行が表に無い
            return 'F (Xist+; Y genes absent from table)' if xv > 10 else 'unclear'
        if np.isnan(xv):  # Xist の行が表に無い
            return 'M (Y genes+; Xist absent from table)' if yv > 5 else 'unclear'
        if xv > 10 and yv < 5:
            return 'F'
        if xv < 5 and yv > 5:
            return 'M'
        return 'unclear'
    t['sex_inferred'] = [call(a, b) for a, b in zip(xi, ym)]
    return t.reset_index(drop=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'counts').mkdir(exist_ok=True)
    s1, x1 = gse130885()
    s2, x2 = gse133922()
    x3 = pd.read_csv(SEQ / 'GSE296758/counts.tsv.gz', sep='\t', index_col=0)
    s3 = pd.read_csv(SEQ / 'GSE296758/samples.csv', dtype=str)
    x1.to_csv(OUT / 'counts/GSE130885_counts.tsv.gz', sep='\t')
    x2.to_csv(OUT / 'counts/GSE133922_counts.tsv.gz', sep='\t')

    sex = pd.concat([sex_table('GSE130885', x1, s1), sex_table('GSE133922', x2, s2),
                     sex_table('GSE296758', x3, s3)], ignore_index=True)
    sex.to_csv(OUT / 'sex_markers.csv', index=False)
    print(sex.groupby(['dataset', 'sex_inferred']).size().to_string())

    # GSE130885 の集団マーカー（群平均 CPM）と、GSE296758 PBS・GSE133922 ADU_CAR との log CPM Spearman
    c1 = cpm(x1)
    pm = c1.reindex(MARKERS).T.groupby(s1.set_index('sample').population).mean().T.round(0)
    pm.to_csv(OUT / 'gse130885_population_markers.csv')
    print(pm.to_string())
    ref_mspc = np.log2(cpm(x3)[s3[s3.group == 'PBS']['sample']].mean(axis=1) + 1)
    car = s2[(s2.population == 'CAR') & (s2.group == 'ADU')]['sample']
    ref_car = np.log2(cpm(x2)[car].mean(axis=1) + 1)
    sim = []
    for pop, d in s1.groupby('population'):
        v = np.log2(c1[d['sample']].mean(axis=1) + 1)
        row = dict(population=pop)
        for nm, ref in [('GSE296758_PBS_MSPC', ref_mspc), ('GSE133922_ADU_CAR', ref_car)]:
            j = v.index.intersection(ref.index)
            hv = ref.loc[j].sort_values().index[-3000:]  # 参照側で発現上位 3000 遺伝子で比べる
            row[f'spearman_{nm}'] = v.loc[hv].corr(ref.loc[hv], method='spearman')
        sim.append(row)
    sim = pd.DataFrame(sim)
    sim.to_csv(OUT / 'gse130885_population_similarity.csv', index=False)
    print(sim.round(3).to_string())

    # サンプル表
    s1o = s1.assign(dataset='GSE130885', instrument=s1.instrument)[['dataset', 'gsm', 'sample', 'population', 'group', 'age_detail', 'phenotype', 'strain', 'instrument']]
    s2o = s2.assign(dataset='GSE133922', gsm='(列名と GSM の対応は SOFT に無い)', phenotype='', strain='C57BL/6JRj', instrument='Illumina HiSeq 4000')[
        ['dataset', 'gsm', 'sample', 'population', 'group', 'age_detail', 'phenotype', 'strain', 'instrument', 'size_factor', 'total']]
    s3o = s3.assign(dataset='GSE296758', population='LepR+ MSPC (Lepr-Cre;tdTomato+)', age_detail='8-10 wk at treatment, harvested 8 wk after last dose',
                    phenotype='tdTomato+', strain='Lepr-Cre;tdTomato', instrument='Illumina NovaSeq 6000')[
        ['dataset', 'gsm', 'sample', 'population', 'group', 'age_detail', 'phenotype', 'strain', 'instrument']]
    S = pd.concat([s1o, s2o, s3o], ignore_index=True).merge(sex[['dataset', 'sample', 'sex_inferred']], on=['dataset', 'sample'], how='left')
    S.to_csv(OUT / 'design_samples.csv', index=False)
    summ = S.groupby(['dataset', 'population', 'group', 'age_detail']).agg(n=('sample', 'size'),
                                                                          sex_inferred=('sex_inferred', lambda v: '/'.join(f'{k}:{n}' for k, n in v.value_counts().items()))).reset_index()
    summ.to_csv(OUT / 'design_summary.csv', index=False)
    print(summ.to_string())


if __name__ == '__main__':
    main()

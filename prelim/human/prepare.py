"""ヒト骨 RNA-seq（GEO 補足ファイル）→ 遺伝子シンボル単位の整数 counts と samples.csv。

  uv run python -I prelim/human/prepare.py      # repo root から

出力（prelim/de.py の run() がそのまま読める形）:
  data/seq/GSE141595/counts.tsv.gz, samples.csv        閉経前 15 + 閉経後 placebo 15 + 閉経後 Dmab 15。
                                                       **閉経前ファイルに存在する遺伝子（37,163 ID）だけ**
  data/seq/GSE141595/post30/counts.tsv.gz, samples.csv 閉経後 30 例（placebo / Dmab）の全遺伝子
  data/seq/GSE141610/counts.tsv.gz, samples.csv        遠心のみ（消化前）の閉経後 30 例
  data/seq/GSE72815/counts.tsv.gz, samples.csv         若年 19 / 高齢 19 / 高齢+E 20
  data/seq/GSE276529/counts.tsv.gz, samples.csv        非 GIOP 4 / GIOP 5
  data/seq/GSE141595/human_ensembl2symbol.csv         GSE141595 の GeneID/GeneName 列から作った対応表

処理:
- 同一シンボルは合算。RSEM の非整数 counts は四捨五入。
- GSE72815 の counts は全値が偶数（mate 二重計数）なので /2（GSE222752 と同じ扱い）。
- GSE141595 の閉経前ファイル（2022 年追加分）は chr1〜11, X と chr12 の途中で途切れている（行が欠落）。
  閉経前を含む比較はこの共通遺伝子に限る。
- GSE276529 の Ensembl（バージョン付き）→ シンボルは、バージョンを外して上記対応表で引く。引けないものは ID のまま。
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

SEQ = Path('data/seq')


def soft_samples(gse):
    rows, cur = [], None
    with gzip.open(SEQ / gse / f'{gse}_family.soft.gz', 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE'):
                cur = {'gsm': line.split('= ')[1]}
                rows.append(cur)
            elif cur is None:
                continue
            elif line.startswith('!Sample_title'):
                cur['title'] = line.split('= ', 1)[1]
            elif line.startswith('!Sample_source_name_ch1'):
                cur['source'] = line.split('= ', 1)[1]
            elif line.startswith('!Sample_characteristics_ch1'):
                k, v = line.split('= ', 1)[1].split(': ', 1)
                cur[k.strip().lower()] = v.strip()
    return pd.DataFrame(rows)


def by_symbol(df, sym, cols):
    x = df[cols].copy()
    x.index = df[sym].astype(str).values
    x = x.groupby(level=0).sum().round().astype(int)
    x.index.name = 'gene'
    return x


def write(gse, counts, samples, sub=None):
    d = SEQ / gse / sub if sub else SEQ / gse
    d.mkdir(parents=True, exist_ok=True)
    assert list(counts.columns) == list(samples['sample']), gse
    counts.to_csv(d / 'counts.tsv.gz', sep='\t')
    samples.to_csv(d / 'samples.csv', index=False)
    print(f'{d}: {counts.shape[0]} genes x {counts.shape[1]} samples;',
          samples.group.value_counts().to_dict())


def gse141595():
    post = pd.read_csv(SEQ / 'GSE141595/GSE141595_Osteocyte_enriched_bone_biopsies_gene_count.csv.gz')
    pre = pd.read_csv(SEQ / 'GSE141595/GSE141595_Raw_Data_GSM6246351-GSM6246365_June_13_2022.csv.gz')
    ann = post[['GeneID', 'GeneName', 'Chr']].rename(columns={'GeneID': 'ensembl', 'GeneName': 'symbol'})
    ann.to_csv(SEQ / 'GSE141595/human_ensembl2symbol.csv', index=False)

    s = soft_samples('GSE141595')
    s['sample'] = s.title.str.replace('_Placebo', '', regex=False)
    s['group'] = 'Post'
    s.loc[s.treatment == 'Denosumab', 'group'] = 'Post_Dmab'
    s.loc[s.age == '25-37 years', 'group'] = 'Pre'
    s['batch'] = s.group.map({'Pre': 'b2022', 'Post': 'b2019', 'Post_Dmab': 'b2019'})
    s = s[['sample', 'gsm', 'group', 'batch', 'age', 'treatment']]

    def cols(df):
        c = [x for x in df.columns if x.endswith('_GeneCount')]
        return c, [re.sub(r'^s_(.*)_GeneCount$', r'\1', x) for x in c]

    pc, pn = cols(post)
    rc, rn = cols(pre)
    assert set(pn) | set(rn) == set(s['sample']) and not set(pn) & set(rn)
    # 全遺伝子（閉経後 30 例のみ）
    xp = by_symbol(post.rename(columns=dict(zip(pc, pn))), 'GeneName', pn)
    sp = s[s['sample'].isin(pn)].sort_values(['group', 'sample'])
    write('GSE141595', xp[sp['sample']], sp, sub='post30')
    # 閉経前ファイルにある ID に限定して 45 例
    common = post[post.GeneID.isin(pre.GeneID)]
    m = common.rename(columns=dict(zip(pc, pn)))[['GeneID', 'GeneName'] + pn].merge(
        pre.rename(columns=dict(zip(rc, rn)))[['GeneID'] + rn], on='GeneID', how='inner')
    x = by_symbol(m, 'GeneName', pn + rn)
    sa = s.sort_values(['group', 'sample'])
    write('GSE141595', x[sa['sample']], sa)
    lost = post[~post.GeneID.isin(pre.GeneID)]
    print(f'  pre-file truncation: {len(lost)} / {len(post)} IDs missing in pre file; '
          f'chromosomes present in pre: {sorted(pre.Chr.unique())}')


def gse141610():
    cb = pd.read_csv(SEQ / 'GSE141610/GSE141610_Centrifuged_bone_biopsies_counts.csv.gz')
    c = [x for x in cb.columns if x.endswith('_GeneCount')]
    n = [re.sub(r'^s_(.*)_GeneCount$', r'\1', x) for x in c]
    s = soft_samples('GSE141610')
    s['sample'] = s.title.str.extract(r'\[(.*)\]')[0]
    s['group'] = s.treatment.map({'Placebo': 'Post', 'Denosumab': 'Post_Dmab'})
    s = s[['sample', 'gsm', 'group', 'treatment']].sort_values(['group', 'sample'])
    assert set(n) == set(s['sample'])
    x = by_symbol(cb.rename(columns=dict(zip(c, n))), 'GeneName', n)
    write('GSE141610', x[s['sample']], s)


def gse72815():
    y = pd.read_csv(SEQ / 'GSE72815/GSE72815_YOE_counts.csv.gz', index_col=0)
    c = [x for x in y.columns if x.startswith('s_')]
    assert (y[c] % 2 == 0).all().all()
    y[c] = y[c] // 2
    s = soft_samples('GSE72815')
    s['sample'] = s.title.str.replace('S_', 's_', regex=False)
    s['group'] = s.source.map({'Young Woman': 'Young', 'Old Woman': 'Old',
                               'Old Woman + Estrogen': 'Old_E'})
    s = s[['sample', 'gsm', 'group', 'source']].sort_values(['group', 'sample'])
    assert set(c) == set(s['sample']), set(c) ^ set(s['sample'])
    x = by_symbol(y, 'GeneID', c)
    write('GSE72815', x[s['sample']], s)


def gse276529():
    g = pd.read_csv(SEQ / 'GSE276529/GSE276529_rsem_count.tsv.gz', sep='\t', index_col=0)
    ann = pd.read_csv(SEQ / 'GSE141595/human_ensembl2symbol.csv').drop_duplicates('ensembl')
    m = dict(zip(ann.ensembl, ann.symbol))
    base = g.index.str.replace(r'\.\d+$', '', regex=True)
    base = base.str.replace(r'_PAR_Y$', '', regex=True)
    sym = [m.get(b, i) for b, i in zip(base, g.index)]
    expressed = g.sum(axis=1) >= 10
    unm = sum(1 for b, e in zip(base, expressed) if e and b not in m)
    print(f'  GSE276529: expressed IDs without symbol: {unm}')
    g = g.assign(symbol=sym)
    s = soft_samples('GSE276529')
    # 列名 Ctrl_k / GIOP_k(_2) と title "femur, Non-GIOP, k" / "femur, GIOP, k" を番号で対応させる
    t = s.title.str.extract(r'femur, (Non-GIOP|GIOP), (\d+)')
    s['group'] = t[0].map({'Non-GIOP': 'NonGIOP', 'GIOP': 'GIOP'})
    s['k'] = t[1]
    colmap = {}
    for col in g.columns[:-1]:
        mm = re.match(r'(Ctrl|GIOP)_(\d+)', col)
        colmap[(('NonGIOP' if mm[1] == 'Ctrl' else 'GIOP'), mm[2])] = col
    s['sample'] = [colmap[(a, b)] for a, b in zip(s.group, s.k)]
    s = s[['sample', 'gsm', 'group', 'age']].sort_values(['group', 'sample'])
    x = by_symbol(g, 'symbol', list(s['sample']))
    x = x[s['sample']]
    cpm = x / x.sum() * 1e6
    # 背景リード指標（嗅覚受容体 OR の CPM 合計。マウス GSE276391 と同じ定義・同じ標準化）
    bg = cpm[cpm.index.str.match(r'^OR\d+[A-Z]+\d+')].sum()
    lb = np.log2(bg.values)
    s['bg_or_cpm'] = bg.round(1).values
    s['bg'] = ((lb - lb.mean()) / lb.std()).round(4)
    write('GSE276529', x, s)
    # 感度分析用: bg 列に骨髄・血液混入スコア（HBB, HBA1, MPO, ELANE の log2 CPM 平均の z）を入れた別ディレクトリ
    # （prelim/de.py は 'bg' 列だけを連続共変量として扱うため）
    ms = np.log2(cpm.loc[['HBB', 'HBA1', 'MPO', 'ELANE']] + 1).mean().values
    s2 = s.drop(columns=['bg']).assign(marrow_score=ms.round(4),
                                       bg=((ms - ms.mean()) / ms.std()).round(4))
    write('GSE276529', x, s2, sub='marrowadj')


if __name__ == '__main__':
    gse141595()
    gse141610()
    gse72815()
    gse276529()

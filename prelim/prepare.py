"""GEO から落とした補足ファイルを読み、遺伝子シンボル単位の整数 counts と samples.csv を作る。

出力（各 data/seq/<GSE>/）:
  samples.csv     sample(=counts の列名), gsm, title, group, ...（SOFT から）
  counts.tsv.gz   gene symbol x sample の整数 counts（同一シンボルは合算）
共通:
  data/seq/annot/ensembl2symbol.csv  Ensembl gene ID（版なし）→ MGI symbol

Ensembl→symbol は GTF もネットも使わず、GEO 自体に含まれる対応表で作る:
  1) GSE205237 の count 表（GENCODE M23 = Ensembl 98、id_gene,gene_name）
  2) NIHMS1776347 tS1（Ensembl Id, Gene Id）で補完
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEQ = Path('data/seq')
ANNOT = SEQ / 'annot' / 'ensembl2symbol.csv'


def read_soft(gse):
    """SOFT から sample ごとの属性を dict のリストで返す。"""
    path = SEQ / gse / f'{gse}_family.soft.gz'
    rows, cur = [], None
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE'):
                cur = {'gsm': line.split('=', 1)[1].strip()}
                rows.append(cur)
            elif line.startswith('^') and not line.startswith('^SAMPLE'):
                cur = None
            elif cur is not None and line.startswith('!Sample_'):
                k, _, v = line.partition(' = ')
                k = k[len('!Sample_'):]
                if k == 'characteristics_ch1' and ':' in v:
                    ck, _, cv = v.partition(':')
                    cur[ck.strip().lower().replace(' ', '_')] = cv.strip()
                elif k in ('title', 'source_name_ch1', 'description'):
                    cur[k] = cur[k] + ' | ' + v if k in cur else v
    return pd.DataFrame(rows)


def build_annot():
    a = pd.read_csv(SEQ / 'GSE205237' / 'GSE205237_COUNTS_genes_VENTURAFRA_01.txt.gz',
                    sep='\t', usecols=[0])
    k = a.iloc[:, 0].str.split(',', expand=True)
    m1 = pd.DataFrame({'ensembl': k[0].str.split('.').str[0], 'symbol': k[1],
                       'biotype': k[2], 'source': 'GSE205237 (GENCODE M23)'})
    n = pd.read_excel('data/NIHMS1776347-supplement-tS1.xlsx',
                      sheet_name='All results for 2m-2y', header=None, skiprows=8)
    m2 = pd.DataFrame({'ensembl': n[0].astype(str), 'symbol': n[1].astype(str),
                       'biotype': np.nan, 'source': 'NIHMS1776347 tS1'})
    m2 = m2[m2.ensembl.str.startswith('ENSMUSG') & (m2.symbol != 'nan')]
    m = pd.concat([m1, m2[~m2.ensembl.isin(m1.ensembl)]], ignore_index=True)
    m = m.drop_duplicates('ensembl')
    ANNOT.parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(ANNOT, index=False)
    print(f'annot: {len(m)} ids ({(m.source != m1.source.iloc[0]).sum()} from NIHMS)')
    return dict(zip(m.ensembl, m.symbol))


def to_symbol_counts(df, symbols):
    """df: rows=genes, cols=samples（数値）。symbols: 同じ index のシンボル Series。"""
    df = df.copy()
    df.index = symbols.values
    df = df.groupby(level=0).sum()
    return df.round().astype(np.int64)


def write(gse, counts, samples):
    assert list(counts.columns) == list(samples['sample']), (gse, counts.columns, samples['sample'])
    samples.to_csv(SEQ / gse / 'samples.csv', index=False)
    counts.index.name = 'gene'
    counts.to_csv(SEQ / gse / 'counts.tsv.gz', sep='\t')
    print(f'{gse}: {counts.shape[0]} genes x {counts.shape[1]} samples; '
          f'groups={samples.group.value_counts().to_dict()}')


def gse205237(e2s):
    g = 'GSE205237'
    a = pd.read_csv(SEQ / g / 'GSE205237_COUNTS_genes_VENTURAFRA_01.txt.gz', sep='\t')
    k = a.iloc[:, 0].str.split(',', expand=True)
    c = a.iloc[:, 1:].astype(float)
    counts = to_symbol_counts(c, k[1])
    s = read_soft(g)
    s['sample'] = s['description']
    s['group'] = np.where(s.treatment.str.contains('DoxoR'), 'Dox', 'Vehicle')
    s = s.set_index('sample').loc[list(counts.columns)].rename_axis('sample').reset_index()
    s['note'] = ''
    write(g, counts, s[['sample', 'gsm', 'title', 'group', 'treatment', 'tissue', 'sex', 'note']])


def gse276391(e2s):
    g = 'GSE276391'
    a = pd.read_csv(SEQ / g / 'GSE276391_rsem_count.tsv.gz', sep='\t')
    ids = a.gene_id.str.split('.').str[0]
    sym = ids.map(e2s).fillna(ids)
    c = a.drop(columns='gene_id').astype(float)
    counts = to_symbol_counts(c, sym)
    s = read_soft(g)
    s['key'] = s['description']
    # counts の列名は Vehicle_*_{1,2,4,5,6}、SOFT は Vehicle_*_{1..5}。群内の並び順で対応させる。
    rows = []
    for col in counts.columns:
        prefix, idx = col.rsplit('_', 1)
        rows.append({'sample': col, 'prefix': prefix, 'idx': int(idx)})
    cs = pd.DataFrame(rows)
    cs['rank'] = cs.groupby('prefix').idx.rank().astype(int)
    s['prefix'] = s.key.str.rsplit('_', n=1).str[0]
    s['rank'] = s.key.str.rsplit('_', n=1).str[1].astype(int)
    m = cs.merge(s, on=['prefix', 'rank'], how='left')
    m['note'] = np.where(m['sample'] != m['key'],
                         'counts 列名と SOFT description 不一致: 群内の番号順で対応付け', '')
    m['drug'] = np.where(m.prefix.str.startswith('DEX'), 'DEX', 'Vehicle')
    m['loading'] = np.where(m.prefix.str.endswith('Load'), 'Load', 'Ctrl')
    m['group'] = m.drug + '_' + m.loading
    m = m.set_index('sample').loc[list(counts.columns)].rename_axis('sample').reset_index()
    # gDNA/背景リードの指標: 嗅覚・鋤鼻受容体遺伝子（骨では発現しない）の CPM。z-score した log2 を共変量候補に。
    bg = counts[counts.index.str.match(r'^(Olfr|Vmn1r|Vmn2r)')].sum() / counts.sum() * 1e6
    lb = np.log2(bg.values)
    m['bg_olfr_cpm'] = bg.round(1).values
    m['bg'] = ((lb - lb.mean()) / lb.std()).round(4)
    write(g, counts, m[['sample', 'gsm', 'title', 'group', 'drug', 'loading', 'treatment',
                        'bg_olfr_cpm', 'bg', 'note']])
    unm = (~ids.isin(e2s.keys())) & (c.sum(axis=1) > 0)
    print(f'  {g}: expressed Ensembl IDs without symbol = {unm.sum()} (kept as Ensembl ID)')


def gse222752(e2s):
    g = 'GSE222752'
    d = pd.read_excel(SEQ / g / 'GSE222752_mRNA_counts_and_RPKM.xlsx', header=1)
    cols = [x for x in d.columns if x.endswith('_GeneCount')]
    c = d[cols].astype(float)
    # 全値が偶数: ペアエンドの両 mate を数えている（1 fragment = 2）。fragment 単位に戻す。
    assert (c % 2 == 0).all().all()
    c = c / 2
    c.columns = [x.split('_')[0] for x in cols]
    counts = to_symbol_counts(c, d.GeneID.astype(str))
    s = read_soft(g)
    s['sample'] = s.title.str.replace('_mrna', '', regex=False)
    s['compartment'] = s.tissue_type
    s['group'] = s.treatment.str.replace('+', '_', regex=False)  # Sham / OVX / OVX_ERT
    s = s.set_index('sample').loc[list(counts.columns)].rename_axis('sample').reset_index()
    s['note'] = '各サンプルは 2 匹の RNA を等量プール。counts は原表 /2'
    write(g, counts, s[['sample', 'gsm', 'title', 'source_name_ch1', 'compartment', 'group',
                        'treatment', 'note']])


def gse238244(e2s):
    g = 'GSE238244'
    a = pd.read_csv(SEQ / g / 'GSE238244_RawCounts.txt.gz', sep='\t')
    c = a.drop(columns='Feature').astype(float)
    c.columns = [x.split('_')[0] for x in c.columns]
    counts = to_symbol_counts(c, a.Feature.astype(str))
    s = read_soft(g)
    s['sample'] = s.title.str.split(', ').str[1]
    code = {'s': 'Sham', 'ov': 'OVX_PBS', 'h': 'OVX_HA', 'p': 'OVX_PPA', 'z': 'OVX_ZA'}
    s['group'] = s['sample'].str.extract(r'^([a-z]+)')[0].map(code)
    s = s.set_index('sample').loc[list(counts.columns)].rename_axis('sample').reset_index()
    s['note'] = 'Partek の非整数 counts を四捨五入'
    write(g, counts, s[['sample', 'gsm', 'title', 'group', 'treatment', 'note']])


def gse296757(e2s):
    g = 'GSE296757'
    a = pd.read_csv(SEQ / g / 'GSE296757_raw_counts.csv.gz')
    c = a.drop(columns='gene_name').astype(float)
    counts = to_symbol_counts(c, a.gene_name.astype(str))
    s = read_soft(g)
    s['sample'] = s.description.str.replace('Library name: ', '', regex=False)
    # 群は title（PBS / DOX4X）で決める。GSM8975970 (S1_D1) は title=DOX4X だが
    # treatment 欄が "4X PBS IV" と食い違う。library 名の D も Dox を示すので title を採る。
    s['group'] = np.where(s.title.str.startswith('DOX'), 'DOX4X', 'PBS')
    s['note'] = np.where(s.title.str.startswith('DOX') & s.treatment.str.contains('PBS'),
                         'SOFT の treatment 欄は PBS だが title/library 名は DOX → DOX として扱う', '')
    s = s.set_index('sample').loc[list(counts.columns)].rename_axis('sample').reset_index()
    write(g, counts, s[['sample', 'gsm', 'title', 'group', 'batch', 'treatment', 'note']])


def main():
    e2s = build_annot()
    for fn in (gse205237, gse276391, gse222752, gse238244, gse296757):
        fn(e2s)


if __name__ == '__main__':
    main()

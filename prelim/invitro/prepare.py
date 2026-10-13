"""in vitro ヒト骨系細胞 × DNA 傷害（GEO 補足ファイル）→ シンボル単位の整数 counts と samples.csv。

  uv run python -I prelim/invitro/prepare.py      # repo root から

出力（prelim/de.py の run() がそのまま読める形）:
  data/seq/GSE284666/counts.tsv.gz, samples.csv   NHO / BMMSC × 対照(DMSO 8 日) / エトポシド 8 日 / 放射線 8 日、各 n=4
  data/seq/GSE135323/counts.tsv.gz, samples.csv   TE-Bone（3D 骨芽細胞）の Dox 4.16 µM 2 週 / 溶媒対照（各 n=3）
                                                  ＋ Linsitinib / 対照（使わないが samples に残す）

処理:
- GSE284666: featureCounts（Ensembl v104 ID）。シンボルは data/seq/GSE141595/human_ensembl2symbol.csv で引き、
  引けない ID は ID のまま。同一シンボルは合算。
- GSE135323: TempO-Seq S1500+ のプローブ（<SYMBOL>_<id>）を遺伝子単位に合算。列名は fastq 名で、SOFT の
  Sample_description（fastq 名）と対応させる。
"""
import gzip
from pathlib import Path

import pandas as pd

SEQ = Path('data/seq')


def soft_samples(gse):
    rows, cur = [], None
    with gzip.open(SEQ / gse / f'{gse}_family.soft.gz', 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE'):
                cur = {'gsm': line.split('= ')[1], 'desc': []}
                rows.append(cur)
            elif cur is None:
                continue
            elif line.startswith('!Sample_title'):
                cur['title'] = line.split('= ', 1)[1]
            elif line.startswith('!Sample_characteristics_ch1'):
                k, v = line.split('= ', 1)[1].split(': ', 1)
                cur[k.strip().lower()] = v.strip()
            elif line.startswith('!Sample_description'):
                cur['desc'].append(line.split('= ', 1)[1])
    return pd.DataFrame(rows)


def write(gse, counts, samples):
    d = SEQ / gse
    counts.to_csv(d / 'counts.tsv.gz', sep='\t')
    samples.to_csv(d / 'samples.csv', index=False)
    print(gse, counts.shape, samples.group.value_counts().to_dict())


def gse284666():
    s = soft_samples('GSE284666')
    s = s[s['cell line'].isin(['NHO', 'BMMSC'])].copy()
    files = {f.name.split('_')[0]: f for f in (SEQ / 'GSE284666/raw').glob('GSM*_geneCOUNT.txt.gz')}
    cols = {}
    for gsm in s.gsm:
        x = pd.read_csv(files[gsm], sep='\t', index_col=0)
        cols[gsm] = x.iloc[:, 0]
    c = pd.DataFrame(cols)
    m = pd.read_csv(SEQ / 'GSE141595/human_ensembl2symbol.csv').drop_duplicates('ensembl').set_index('ensembl').symbol
    sym = c.index.map(lambda e: m.get(e, e))
    sym = pd.Index([e if pd.isna(v) else v for e, v in zip(c.index, sym)])
    expressed = c.sum(axis=1) > 0
    print('GSE284666: Ensembl IDs', len(c), 'expressed', expressed.sum(),
          'unmapped(expressed)', (expressed & (sym == c.index)).sum())
    c.index = sym
    c = c.groupby(level=0).sum().astype(int)
    c.index.name = 'gene'
    tr = {'Control proliferating cells': 'Ctrl'}
    s['cell'] = s['cell line']
    s['trt'] = s.treatment.map(lambda t: 'Ctrl' if t in tr else ('Eto' if 'Etoposide' in t else 'IR'))
    s['group'] = s.cell + '_' + s.trt
    s['rep'] = s.title.str.rsplit('_', n=1).str[1]
    s['sample'] = s.gsm
    write('GSE284666', c[s.gsm], s[['sample', 'gsm', 'title', 'cell', 'trt', 'group', 'rep', 'treatment']])


def gse135323():
    s = soft_samples('GSE135323')
    s = s[s['cell type'].str.contains('TE-Bone', na=False)].copy()
    s['fastq'] = s.desc.map(lambda d: next(x for x in d if x.endswith('.fastq.gz')))
    s = s[s.exposure.isin(['Doxorubicin', 'Doxorubicin Control', 'Linsitinib', 'Linsitinib Control'])].copy()
    s['group'] = s.exposure.map({'Doxorubicin': 'Dox', 'Doxorubicin Control': 'DoxCtrl',
                                 'Linsitinib': 'Lins', 'Linsitinib Control': 'LinsCtrl'})
    raw = pd.read_csv(SEQ / 'GSE135323/GSE135323_TCTC_Mapped_Reads.csv.gz', index_col=0)
    print('GSE135323: probes', raw.shape[0], 'samples', raw.shape[1])
    x = raw[s.fastq].copy()
    x.columns = s.gsm.values
    gene = x.index.str.rsplit('_', n=1).str[0]
    print('  genes', gene.nunique(), 'multi-probe genes', (gene.value_counts() > 1).sum())
    x.index = gene
    x = x.groupby(level=0).sum().astype(int)
    x.index.name = 'gene'
    s['sample'] = s.gsm
    write('GSE135323', x, s[['sample', 'gsm', 'title', 'fastq', 'exposure', 'treatment concentration', 'group']])


if __name__ == '__main__':
    gse284666()
    gse135323()

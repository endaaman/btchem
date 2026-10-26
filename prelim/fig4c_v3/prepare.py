"""GSE304562（MC3T3-E1 骨芽細胞 × コルチコステロン 24 h）の counts と samples.csv を作る。

  uv run python -I prelim/fig4c_v3/prepare.py      # repo root から

GEO の補足ファイルは `GSE304562_FPKM_allSamples.txt.gz` 1 本だけで、名前は FPKM だが中身は FPKM ではない:
- 列和が全サンプル ~2,000 万（FPKM なら遺伝子長で割るので ~10^5–10^6 になる）。
- 各列の値はすべて「その列の最小の非ゼロ値」の**整数倍**（全非ゼロ値で |比 − 整数| < 0.01）。
  ＝ DESeq2 型の正規化 counts（count / size factor）。最小値 = 1 / size factor（count 1 の遺伝子が数千個ある）。
- 最小値の逆数（= size factor）の幾何平均は ~1.0（DESeq2 median-of-ratios の性質）、
  復元した列和（1,600–2,560 万）は SRA の read pair 数（2,090–3,490 万）の 73–77 % で、遺伝子への割り当て率として自然。
SOFT の data_processing にも "raw counts for each sample" とある。そこで value / 列最小値 を丸めて整数 counts に戻した。
丸め誤差（|比 − 整数|）の最大値は stdout と data/seq/GSE304562/count_recovery.csv に残す。

Ensembl → symbol は既存の data/seq/annot/ensembl2symbol.csv（GENCODE M23 + NIHMS）。引けない ID は ID のまま残す
（prelim/prepare.py と同じ扱い）。同一シンボルは合算。

出力（data/seq/GSE304562/。既存ファイルは上書きしない）: counts.tsv.gz, samples.csv, count_recovery.csv
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / 'data/seq/GSE304562'
ANNOT = ROOT / 'data/seq/annot/ensembl2symbol.csv'


def read_soft():
    rows, cur = [], None
    with gzip.open(D / 'GSE304562_family.soft.gz', 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE'):
                cur = {'gsm': line.split('=', 1)[1].strip()}
                rows.append(cur)
            elif line.startswith('^'):
                cur = None
            elif cur is not None and line.startswith('!Sample_'):
                k, _, v = line.partition(' = ')
                k = k[len('!Sample_'):]
                if k == 'characteristics_ch1' and ':' in v:
                    ck, _, cv = v.partition(':')
                    cur[ck.strip().lower().replace(' ', '_')] = cv.strip()
                elif k in ('title', 'treatment_protocol_ch1'):
                    cur[k] = v.strip()
    return pd.DataFrame(rows)


def main():
    out = {k: D / k for k in ('counts.tsv.gz', 'samples.csv', 'count_recovery.csv')}
    if all(p.exists() for p in out.values()):
        print('already prepared:', ', '.join(str(p.relative_to(ROOT)) for p in out.values()))
        return
    if any(p.exists() for p in out.values()):
        sys.exit('一部だけ存在する。上書きしないので手で確認すること: ' + str([p for p in out.values() if p.exists()]))

    v = pd.read_csv(D / 'GSE304562_FPKM_allSamples.txt.gz', sep='\t', index_col=0)
    v.columns = [c.strip() for c in v.columns]           # 'Sample 4' など（先頭に空白入りの列名がある）
    rec, cnt = [], {}
    for c in v.columns:
        x = v[c]
        m = x[x > 0].min()
        r = x / m
        err = (r - r.round()).abs()
        cnt[c] = r.round().astype(np.int64)
        rec.append(dict(column=c, min_nonzero=m, size_factor=1 / m, n_nonzero=int((x > 0).sum()),
                        n_count1=int((r.round() == 1).sum()), max_abs_rounding_error=err.max(),
                        recovered_total=int(cnt[c].sum())))
    rec = pd.DataFrame(rec)
    gm = np.exp(np.log(rec.size_factor).mean())
    print(rec.round(4).to_string(index=False))
    print(f'geometric mean of size factors = {gm:.4f}')
    assert rec.max_abs_rounding_error.max() < 0.01, '整数倍になっていない列がある'
    counts = pd.DataFrame(cnt)

    # サンプル: 列 'Sample k' ↔ SOFT title 'Sample k_<grp>_<rep>'
    s = read_soft()
    s['column'] = s.title.str.extract(r'^(Sample \d+)')[0]
    s['group'] = s.treatment.map({'Vehicle': 'Ctrl', 'Corticosterone': 'Cort'})
    assert s.group.notna().all()
    s['rep'] = s.title.str.extract(r'_(\d+)$')[0]
    s['sample'] = s.group + '_' + s.rep
    s = s.set_index('column').loc[[c for c in counts.columns]]   # counts の列順
    s = s.reset_index()
    counts.columns = s['sample'].values

    # symbol
    a = pd.read_csv(ANNOT)
    e2s = dict(zip(a.ensembl, a.symbol))
    sym = pd.Series([e2s.get(i, i) for i in counts.index], index=counts.index)
    expressed = counts.sum(axis=1) > 0
    print(f'Ensembl IDs: {len(counts)}; unmapped {(~counts.index.isin(e2s)).sum()} '
          f'(expressed {(~counts.index.isin(e2s) & expressed).sum()}) -> kept as ID')
    counts.index = sym.values
    counts = counts.groupby(level=0).sum()
    order = sorted(s['sample'], key=lambda z: (z.split('_')[0] != 'Ctrl', z))
    counts = counts[order]
    s = s.set_index('sample').loc[order].reset_index()
    rec['sample'] = rec.column.map(dict(zip(s.column, s['sample'])))

    s[['sample', 'gsm', 'title', 'group', 'column', 'cell_line', 'cell_type', 'treatment', 'batch',
       'treatment_protocol_ch1']].to_csv(out['samples.csv'], index=False)
    counts.index.name = 'gene'
    counts.to_csv(out['counts.tsv.gz'], sep='\t')
    rec.to_csv(out['count_recovery.csv'], index=False)
    print(f'GSE304562: {counts.shape[0]} genes x {counts.shape[1]} samples; libsize',
          counts.sum().to_dict())
    print(s[['sample', 'gsm', 'title', 'group', 'column']].to_string(index=False))


if __name__ == '__main__':
    main()

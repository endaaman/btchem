"""pydeseq2 で差次発現。

  uv run python prelim/de.py            # 登録済み contrast を全部
  uv run python prelim/de.py B_Dox C_DEX

入力: data/seq/<GSE>/counts.tsv.gz, samples.csv（prelim/prepare.py が作る）
出力: out/prelim/de/<label>.csv（gene, baseMean, log2FoldChange, lfcSE, pvalue, padj）
      log2FoldChange は test / ref（例: OVX / Sham）。LFC shrink はしない（MLE）。
低発現フィルタ（edgeR filterByExpr の既定値に準拠、全 contrast 共通）:
  比較に使うサンプルで、CPM >= 10 / (中央値ライブラリサイズ[M]) を満たすサンプルが小さい方の群の n 以上、
  かつ counts 合計 >= 15。
  ※「合計 >= 10」だけだと GSE276391（total RNA pico キット）で Olfr/Gm 等の低カウント背景が
    DEX 側で一斉に上がる artefact が出た（docs/prelim_de.md 参照）。
"""
import sys
from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

sys.path.insert(0, str(Path(__file__).parent))
from contrasts import C, AGING  # noqa: E402

SEQ = Path('data/seq')
OUT = Path('out/prelim/de')
COLS = ['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj']
MIN_COUNT = 10
MIN_TOTAL = 15


def run(c):
    counts = pd.read_csv(SEQ / c['gse'] / 'counts.tsv.gz', sep='\t', index_col=0)
    s = pd.read_csv(SEQ / c['gse'] / 'samples.csv', dtype=str)
    if 'bg' in s:
        s['bg'] = s.bg.astype(float)  # 連続共変量
    for k, v in c['subset'].items():
        s = s[s[k] == v]
    s = s[s.group.isin([c['test'], c['ref']])].copy()
    s['cond'] = (s.group == c['test']).map({True: 'test', False: 'ref'})
    s = s.set_index('sample')
    x = counts[s.index].T  # samples x genes
    lib = x.sum(axis=1)
    cpm_cut = MIN_COUNT / (lib.median() / 1e6)
    cpm = x.div(lib, axis=0) * 1e6
    n_min = s.cond.value_counts().min()
    keep = ((cpm >= cpm_cut).sum(axis=0) >= n_min) & (x.sum(axis=0) >= MIN_TOTAL)
    x = x.loc[:, keep]
    design = '~' + ' + '.join(c['covariates'] + ['cond'])
    dds = DeseqDataSet(counts=x, metadata=s[c['covariates'] + ['cond']], design=design,
                       ref_level=None, refit_cooks=True,
                       inference=DefaultInference(n_cpus=8), quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=['cond', 'test', 'ref'],
                    inference=DefaultInference(n_cpus=8), quiet=True)
    st.summary()
    r = st.results_df.reset_index().rename(columns={'index': 'gene'})
    r = r.rename(columns={r.columns[0]: 'gene'})[COLS]
    n = s.cond.value_counts().to_dict()
    print(f"{c['label']}: {c['test']} (n={n.get('test')}) vs {c['ref']} (n={n.get('ref')}), "
          f"design={design}, genes tested={keep.sum()}, "
          f"padj<0.05&|lfc|>=1: {((r.padj < 0.05) & (r.log2FoldChange.abs() >= 1)).sum()}")
    return r


def aging(label):
    sheet = AGING[label][0]
    n = pd.read_excel('data/NIHMS1776347-supplement-tS1.xlsx', sheet_name=sheet,
                      header=None, skiprows=8)
    n = n.iloc[:, :9]
    n.columns = ['ensembl', 'gene', 'entrez', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat',
                 'pvalue', 'padj']
    n = n[n.ensembl.astype(str).str.startswith('ENSMUSG')]
    n = n[n.gene.notna() & (n.gene.astype(str) != 'nan')]
    n['gene'] = n.gene.astype(str)
    for k in ['baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj']:
        n[k] = pd.to_numeric(n[k], errors='coerce')
    # 同一シンボルの重複は baseMean 最大を残す
    n = n.sort_values('baseMean', ascending=False).drop_duplicates('gene')
    print(f'{label}: imported {len(n)} genes from "{sheet}"')
    return n[COLS]


def main(labels):
    OUT.mkdir(parents=True, exist_ok=True)
    for lb in labels:
        r = aging(lb) if lb in AGING else run(C[lb])
        r.to_csv(OUT / f'{lb}.csv', index=False)


if __name__ == '__main__':
    main(sys.argv[1:] or list(C))

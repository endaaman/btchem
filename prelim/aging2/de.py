"""Figure 4A 案2・案3: pydeseq2 で差次発現（de.py と同じ低発現フィルタ・設定。prelim/aging_f/deutil.py を借りる）。

  uv run python -I prelim/aging2/de.py          # repo root から（先に prepare.py）

contrast（log2FoldChange は test / ref）:
  主（依頼どおり）
    G235_F_age      GSE235942 雌 22 か月 control (n=3) vs 雌 6 か月 control (n=4)、脛骨（骨髄除去）
    G255_MSC_age    GSE255653 老齢 18 か月 control MSC (n=3) vs 若齢 3 か月 control MSC (n=3)、LepR+ 骨髄間質（FPKM から復元した counts）
  同じデータの中の参考（依頼外。比較の解釈のために回した）
    G235_M_age      GSE235942 雄 22 か月 control (n=10) vs 雄 6 か月 control (n=5) — 同じ施設・同じ処理の雄の加齢
    G235_F_age_comp / G235_M_age_comp  上 2 つに組成共変量（筋 = Acta1/Ckm/Myh4、血液 = Hbb-bs/Hba-a1 の log2 CPM 平均を
                    比較に使うサンプル内で z 化）を入れた感度分析（prelim/chemo2/gsea.py の B_Dox_compadj と同じ作り方）
    G235_M6_HLU     雄 6 か月 後肢懸垂 3 週 (n=5) vs 同 control (n=5) — 廃用（加齢ではない）
    G235_M22_HLU    雄 22 か月 後肢懸垂 3 週 (n=3) vs 同 control (n=10) — 同上
    G255_MSC_KO     老齢 Thbs1 KO MSC (n=3) vs 老齢 control MSC (n=3) — 遺伝子欠損（加齢ではない）
出力: out/prelim/aging2/de/<label>.csv（gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj）, de_summary.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim/aging_f'))
import deutil  # noqa: E402

SEQ = ROOT / 'data/seq'
OUT = ROOT / 'out/prelim/aging2/de'

C = {
    'G235_F_age': dict(gse='GSE235942', test='F22_ctrl', ref='F6_ctrl', comp=False),
    'G235_F_age_comp': dict(gse='GSE235942', test='F22_ctrl', ref='F6_ctrl', comp=True),
    'G235_M_age': dict(gse='GSE235942', test='M22_ctrl', ref='M6_ctrl', comp=False),
    'G235_M_age_comp': dict(gse='GSE235942', test='M22_ctrl', ref='M6_ctrl', comp=True),
    'G235_M6_HLU': dict(gse='GSE235942', test='M6_HLU', ref='M6_ctrl', comp=False),
    'G235_M22_HLU': dict(gse='GSE235942', test='M22_HLU', ref='M22_ctrl', comp=False),
    'G255_MSC_age': dict(gse='GSE255653', test='Aged_Control', ref='Young_Control', comp=False),
    'G255_MSC_KO': dict(gse='GSE255653', test='Aged_TSPKO', ref='Aged_Control', comp=False),
}
COMP = {'muscle': ['Acta1', 'Ckm', 'Myh4'], 'blood': ['Hbb-bs', 'Hba-a1']}


def run(lb, c):
    counts = pd.read_csv(SEQ / c['gse'] / 'counts.tsv.gz', sep='\t', index_col=0)
    s = pd.read_csv(SEQ / c['gse'] / 'samples.csv', dtype=str)
    s = s[s.group.isin([c['test'], c['ref']])].copy()
    s['cond'] = np.where(s.group == c['test'], 'test', 'ref')
    s = s.set_index('sample')
    x = counts[s.index].T
    cov = []
    if c['comp']:
        lib = x.sum(axis=1)
        cpm = x.div(lib, axis=0) * 1e6
        for k, g in COMP.items():
            v = np.log2(cpm[g] + 1).mean(axis=1)
            s[k] = ((v - v.mean()) / v.std()).values
        cov = list(COMP)
        print(s[['cond'] + cov].round(2).to_string())
    n = s.cond.value_counts()
    design = '~' + ' + '.join(cov + ['cond'])
    r, ntest = deutil.deseq(x, s[cov + ['cond']], design, ['cond', 'test', 'ref'], int(n.min()))
    up, dn = deutil.ndeg(r)
    up1, dn1 = deutil.ndeg(r, lfc=0.585, p=0.1)
    row = dict(label=lb, gse=c['gse'], test=c['test'], ref=c['ref'], n_test=int(n['test']), n_ref=int(n['ref']),
               design=design, genes_tested=ntest, up=up, down=dn, up_lenient=up1, down_lenient=dn1)
    print(row)
    return r, row


def main(labels):
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for lb in labels:
        r, row = run(lb, C[lb])
        r.to_csv(OUT / f'{lb}.csv', index=False)
        rows.append(row)
    f = OUT / 'de_summary.csv'
    old = pd.read_csv(f) if f.exists() else pd.DataFrame()
    new = pd.DataFrame(rows)
    if len(old):
        new = pd.concat([old[~old.label.isin(new.label)], new], ignore_index=True)
    new = new.set_index('label').reindex([k for k in C if k in set(new.label)]).reset_index()
    new.to_csv(f, index=False)
    print(new.to_string())


if __name__ == '__main__':
    main(sys.argv[1:] or list(C))

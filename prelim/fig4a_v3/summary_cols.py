"""要約モジュール図（prelim/summary/summary_modules.py）に列として足せる形式で、新しい比較の NES / FDR を出す。

  uv run python -I prelim/fig4a_v3/summary_cols.py      # gsea4a.py の後

行は summary_modules.py の ROWS（8 行）+ DROPPED（外した 4 行）をそのまま import して使う（同じ term・同じ判定閾値）。
値は gsea4a.py が in vitro と同じ set 集合・同じ設定（gsea_iv.collection(), 1000 perm, seed 0, min 4, max 1000）で
回した out/prelim/fig4a_v3/gsea/summaryset_<label>.csv から取る。列の書式（shown, module, term, row_group, contrast, dataset,
column_label, group, NES, FDR, call, source）は out/prelim/summary/summary_modules_values.csv と同じ。

出力: out/prelim/fig4a_v3/summary_modules_values_fig4a_v3.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common4a import OUT, ROOT  # noqa: E402
sys.path.insert(0, str(ROOT / 'prelim/summary'))
from summary_modules import ROWS, DROPPED, SIG, TREND  # noqa: E402  （定義を借りるだけ。summary の出力には触らない）

COLS = [  # (contrast, dataset, column_label, group)
    ('H_Dmab', 'GSE141595', 'Human ♀ bone biopsy (digested) denosumab vs placebo, postmeno.', 'denosumab'),
    ('H_Dmab_CB', 'GSE141610', 'Human ♀ bone biopsy (centrifuged) denosumab vs placebo, postmeno.', 'denosumab'),
    ('H_Dmab_marrowadj', 'GSE141595', 'Human ♀ bone biopsy (digested) denosumab vs placebo, marrow-adjusted', 'denosumab'),
    ('H_Dmab_CB_marrowadj', 'GSE141610', 'Human ♀ bone biopsy (centrifuged) denosumab vs placebo, marrow-adjusted', 'denosumab'),
    ('H_PostMeno_procadj', 'GSE141595', 'Human ♀ bone biopsy postmeno. vs premeno., processing-adjusted (heat shock + IEG)', 'aging'),
]


def main():
    rows = []
    for key, gse, lab, grp in COLS:
        f = OUT / 'gsea' / f'summaryset_{key}.csv'
        r = pd.read_csv(f).set_index('Term')
        for term, rl, rg in ROWS + DROPPED:
            nes, fdr = (float(r.loc[term, 'NES']), float(r.loc[term, 'FDR q-val'])) if term in r.index else (np.nan, np.nan)
            call = ('n/a' if np.isnan(nes) else 'up' if (fdr < SIG and nes > 0) else 'down' if fdr < SIG
                    else 'trend_up' if (fdr < TREND and nes > 0) else 'trend_down' if fdr < TREND else 'ns')
            rows.append(dict(shown=(term, rl, rg) in ROWS, module=rl, term=term, row_group=rg, contrast=key, dataset=gse,
                             column_label=lab, group=grp, NES=nes, FDR=fdr, call=call,
                             source=str(f.relative_to(ROOT))))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / 'summary_modules_values_fig4a_v3.csv', index=False)
    print(T.pivot(index='module', columns='contrast', values='NES').round(2).to_string())
    print(T.pivot(index='module', columns='contrast', values='call').to_string())


if __name__ == '__main__':
    main()

"""デノスマブ・閉経後の DE の感度分析（共変量つき）。主の DE は既存の out/prelim/human/de/ をそのまま使う。

  uv run python -I prelim/fig4a_v4/dmab_de.py [label ...]

既存（読むだけ）: H_Dmab（消化後・全遺伝子）, H_Dmab_common, H_Dmab_CB（遠心のみ）, H_PostMeno, H_Old, H_OldE
新規（ここで回す。設定は prelim/de.py と同じ: pydeseq2 0.5.4, MLE, edgeR filterByExpr 相当）:
  H_Dmab_marrowadj     消化後、~ marrow + cond（marrow = HBB/HBA1/MPO/ELANE の log2 CPM 平均の z。H_GIOP_marrowadj と同じ定義）
  H_Dmab_ocadj         消化後、~ oc + cond（oc = 破骨マーカー 7 遺伝子スコアの z）。
                       「破骨細胞が減ったこと」で説明できない変化が残るかを見る（薬の作用経路そのものを補正するので、効果の推定ではなく分解）
  H_Dmab_CB_marrowadj, H_Dmab_CB_ocadj   遠心のみで同じ
  H_PostMeno_procadj   閉経後 vs 閉経前、~ hs + ieg + cond（hs / ieg = 熱ショック・即時初期遺伝子スコアの z。検体処理の差の代理）。
                       ただし閉経前ファイルの欠落で ieg / hs のマーカーは一部しか無い（composition.py の n_genes）
出力: out/prelim/fig4a_v4/de/<label>.csv（gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj）, covariates.csv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common4a import OUT, counts, deseq, scores, zscore, MARROW_GENES, cpm  # noqa: E402

import numpy as np  # noqa: E402

D = OUT / 'de'


def covars(x):
    sc, _ = scores(x)
    c = np.log2(cpm(x) + 1)
    have = [g for g in MARROW_GENES if g in c.index]
    cv = pd.DataFrame({'marrow': zscore(c.loc[have].mean()),
                       'oc': zscore(sc['osteoclast']),
                       'hs': zscore(sc['heat shock (processing)']),
                       'ieg': zscore(sc['immediate early (processing)'])})
    return cv


RUNS = {
    'H_Dmab_marrowadj': ('GSE141595/post30', 'Post_Dmab', 'Post', ['marrow']),
    'H_Dmab_ocadj': ('GSE141595/post30', 'Post_Dmab', 'Post', ['oc']),
    'H_Dmab_CB_marrowadj': ('GSE141610', 'Post_Dmab', 'Post', ['marrow']),
    'H_Dmab_CB_ocadj': ('GSE141610', 'Post_Dmab', 'Post', ['oc']),
    'H_PostMeno_procadj': ('GSE141595', 'Post', 'Pre', ['hs', 'ieg']),
}


def main(labels):
    D.mkdir(parents=True, exist_ok=True)
    cvs = []
    for lb in labels:
        sub, t, r, cov = RUNS[lb]
        x, s = counts(sub)
        cv = covars(x)
        s = s.join(cv)
        sel = s[s.group.isin([t, r])]
        cvs.append(sel[cov + ['group']].assign(label=lb).reset_index())
        res = deseq(x, s, t, r, covariates=cov, label=lb)
        res.to_csv(D / f'{lb}.csv', index=False)
    if cvs:
        f = D / 'covariates.csv'
        old = pd.read_csv(f) if f.exists() else pd.DataFrame()
        new = pd.concat(cvs)
        if len(old):
            old = old[~old.label.isin(new.label.unique())]
        pd.concat([old, new]).to_csv(f, index=False)


if __name__ == '__main__':
    main(sys.argv[1:] or list(RUNS))

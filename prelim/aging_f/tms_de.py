"""GSE132040（Tabula Muris Senis bulk）「Bone」の加齢 DE。主役は雌、雄は性差の参考。

contrast（事前に決めたもの。tms_prepare.py の組成を見て感度分析を足した）:
  TMS_F_old      雌 12–21 か月（n=8） vs 1–6 か月（n=6）             ~cond
  TMS_F_old_comp 同上 + 組成共変量（筋・血液・造血の z スコア）       ~muscle + blood + hem + cond
                 （Bone は骨髄と付着筋を含み、筋の量が個体で 6.7–11.3 log2CPM とばらつくため。B_Dox_compadj と同じ考え方）
  TMS_F_old_ex1  雌 12–21 か月 vs 3–6 か月（1 か月 = 成長期の幼若骨を除く感度分析）
  TMS_F_trend    雌 1–21 か月、月齢を連続変数（10 か月あたりの log2FC）  ~age10
  TMS_M_old      雄 12–21 か月（n=16） vs 1–6 か月（n=11）
  TMS_M_trend    雄 1–21 か月（雌と範囲をそろえる）、~age10
  TMS_M_old_comp 雄の組成補正版
加えて、雌雄それぞれ log CPM と月齢の Spearman（検定遺伝子のみ、BH 補正）を出す。
DE は deutil.deseq（prelim/de.py と同じフィルタ）。連続変数のフィルタは n_min = 全サンプル数の 1/3（四捨五入）とした。

  uv run python -I prelim/aging_f/tms_de.py
出力: out/prelim/aging_f/de/<contrast>.csv, de_summary_tms.csv, out/prelim/aging_f/tms/spearman_{f,m}.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim/aging_f'))
import deutil  # noqa: E402

T = ROOT / 'out/prelim/aging_f/tms'
DE = ROOT / 'out/prelim/aging_f/de'
OLD = [12, 15, 18, 21]
YOUNG = [1, 3, 6]


def z(v):
    return (v - v.mean()) / v.std()


def main():
    DE.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(T / 'samples.csv')
    c = pd.read_csv(T / 'counts_bone_marrow.tsv.gz', sep='\t', index_col=0)
    b = m[m.tissue == 'Bone'].set_index('sample')
    rows = []

    def run(label, s, design, contrast, n_min):
        x = c[s.index].T
        cols = [k for k in ['muscle', 'blood', 'hem', 'age10', 'cond'] if k in s.columns]
        r, ng = deutil.deseq(x, s[cols], design, contrast, n_min)
        r.to_csv(DE / f'{label}.csv', index=False)
        up, dn = deutil.ndeg(r)
        up5, dn5 = deutil.ndeg(r, lfc=0.585)
        d = dict(label=label, design=design, n=len(s), n_test=int((s.get('cond') == 'test').sum()) if 'cond' in s else np.nan,
                 genes_tested=ng, deg_up=up, deg_down=dn, deg_up_lfc0585=up5, deg_down_lfc0585=dn5,
                 padj05_any=int((r.padj < 0.05).sum()))
        rows.append(d)
        print(d)

    for sx, tag in [('f', 'F'), ('m', 'M')]:
        s = b[b.sex == sx].copy()
        s['muscle'] = z(s['comp:skeletal muscle'])
        s['blood'] = z(s['comp:erythroid/blood'])
        s['hem'] = z(s['comp:hematopoietic'])
        # 組成と月齢の関係（記録用）
        print(tag, 'spearman composition vs age:',
              {k: round(spearmanr(s[k], s.age)[0], 2) for k in ['muscle', 'blood', 'hem']})
        oy = s[s.age.isin(OLD + YOUNG)].copy()
        oy['cond'] = np.where(oy.age.isin(OLD), 'test', 'ref')
        nmin = oy.cond.value_counts().min()
        run(f'TMS_{tag}_old', oy, '~cond', ['cond', 'test', 'ref'], nmin)
        run(f'TMS_{tag}_old_comp', oy, '~muscle + blood + hem + cond', ['cond', 'test', 'ref'], nmin)
        if sx == 'f':
            e = oy[oy.age != 1]
            run('TMS_F_old_ex1', e, '~cond', ['cond', 'test', 'ref'], e.cond.value_counts().min())
        t = s[s.age <= 21].copy()
        t['age10'] = t.age / 10.0
        run(f'TMS_{tag}_trend', t, '~age10', 'age10', int(round(len(t) / 3)))
        # Spearman（月齢 vs log CPM）
        x = c[t.index]
        keep = deutil.expr_filter(x.T, int(round(len(t) / 3)))
        l = np.log2(x.loc[keep[keep].index] / x.sum() * 1e6 + 1)
        rho, p = spearmanr(l.T.values, t.age.values[:, None], axis=0)
        rho = rho[-1, :-1]
        p = p[-1, :-1]
        sp = pd.DataFrame({'gene': l.index, 'rho': rho, 'pvalue': p})
        sp['padj'] = multipletests(sp.pvalue.fillna(1), method='fdr_bh')[1]
        sp.to_csv(T / f'spearman_{sx}.csv', index=False)
        print(tag, 'spearman |rho|>=0.6 & padj<0.05:', ((sp.rho.abs() >= 0.6) & (sp.padj < 0.05)).sum(),
              'padj<0.05:', (sp.padj < 0.05).sum())
    pd.DataFrame(rows).to_csv(DE / 'de_summary_tms.csv', index=False)


if __name__ == '__main__':
    main()

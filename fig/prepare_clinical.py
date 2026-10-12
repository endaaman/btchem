"""1 本目のデータから本論文用の臨床テーブル data/clinical/clinical.xlsx を作る（骨のみ）。

入力（data/clinical/op/、1 本目 ~/Sync/Projects/Shimizu/op/data からコピー）:
  enrolled_patients_20260322.xlsx  … 'perioperative chemo'（0/1、43 例が空欄）、閉経、BMFS_Event
  osteoporosis_bone_meta.csv       … BMFS（月）、BMFS_Event、menopause
両方 Exclude==0 の 785 例が 研究ID で 1:1 に対応し、イベント・閉経は完全一致する。
肺・肝（LMFS/HMFS）は 1 本目に無いので列を作らない。清水の拡張データが来たら差し替える。
化療が空欄の例は chemo 列を NaN のまま残し、fig2/fig3 側（km.load_clinical）で除外する。
"""
import os
import pandas as pd

SRC = 'data/clinical/op'
OUT = 'data/clinical/clinical.xlsx'

e = pd.read_excel(os.path.join(SRC, 'enrolled_patients_20260322.xlsx'), sheet_name=0)
e.columns = [str(c).strip() for c in e.columns]
e = e.loc[:, ~e.columns.str.startswith('Unnamed')]
e = e[e['Exclude'] == 0]
c = pd.read_csv(os.path.join(SRC, 'osteoporosis_bone_meta.csv'))
c.columns = [str(c_).strip() for c_ in c.columns]
c = c[c['Exclude'] == 0]

m = e.merge(c[['研究ID', 'menopause', 'BMFS_Event', 'BMFS']], on='研究ID', how='inner', validate='1:1')
assert len(m) == 785, len(m)
ev = [k for k in m.columns if k.startswith('BMFS_Event（')][0]
assert (m[ev] == m['BMFS_Event']).all() and (m['閉経あり=1'] == m['menopause']).all()

out = pd.DataFrame({
    '研究ID': m['研究ID'].astype(int),
    'Age': m['Age'],
    'menopause': m['menopause'].astype(int),
    'chemo': m['perioperative chemo'],           # NaN = 不明（43 例）
    'endocrine': m['endcrine'],
    'Post-Hormone': m['Post-Hormone'],
    'Stage': m['Stage'],
    'BMFS_Event': m['BMFS_Event'].astype(int),
    'BMFS': m['BMFS'],
})
os.makedirs(os.path.dirname(OUT), exist_ok=True)
out.to_excel(OUT, index=False)
print(f'{OUT}: n={len(out)}  chemo missing={int(out.chemo.isna().sum())}')
print(pd.crosstab(out.menopause, out.chemo.fillna(-1)))

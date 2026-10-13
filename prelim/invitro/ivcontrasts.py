"""in vitro ヒト contrast の登録簿（prelim/contrasts.py と同じ形の dict。prelim/de.py の run() にそのまま渡す）。"""
from pathlib import Path

OUT = Path('out/prelim/invitro')
IC = {}


def reg(label, gse, test, ref, title, short, cell, agent):
    IC[label] = dict(label=label, gse=gse, test=test, ref=ref, subset={}, covariates=[],
                     title=title, short=short, cell=cell, agent=agent)


reg('NHO_Eto', 'GSE284666', 'NHO_Eto', 'NHO_Ctrl', 'NHO osteoblasts: etoposide 20 µM 8 d vs DMSO',
    'NHO, etoposide', 'osteoblast (NHO)', 'etoposide')
reg('NHO_IR', 'GSE284666', 'NHO_IR', 'NHO_Ctrl', 'NHO osteoblasts: IR 10 Gy 8 d vs DMSO',
    'NHO, IR 10 Gy', 'osteoblast (NHO)', 'IR')
reg('MSC_Eto', 'GSE284666', 'BMMSC_Eto', 'BMMSC_Ctrl', 'BM-MSC: etoposide 50 µM 8 d vs DMSO',
    'BM-MSC, etoposide', 'BM-MSC', 'etoposide')
reg('MSC_IR', 'GSE284666', 'BMMSC_IR', 'BMMSC_Ctrl', 'BM-MSC: IR 15 Gy 8 d vs DMSO',
    'BM-MSC, IR 15 Gy', 'BM-MSC', 'IR')
reg('TEB_Dox', 'GSE135323', 'Dox', 'DoxCtrl', 'TE-Bone (3D osteoblasts): Dox 4.16 µM 2 wk vs vehicle',
    'TE-Bone, Dox (TempO-Seq)', 'TE-Bone 3D osteoblast', 'doxorubicin')

# 注目遺伝子（依頼の指定。カテゴリ順）
FOCUS = {
    'apoptosis': ['BAX', 'BBC3', 'PMAIP1', 'FAS', 'FASLG', 'CASP3', 'TNFRSF10B'],
    'senescence': ['CDKN1A', 'CDKN2A', 'SERPINE1', 'IGFBP3'],
    'SASP/inflam.': ['IL6', 'IL1B', 'CXCL8', 'CCL2', 'TNF', 'MMP3', 'MMP13'],
    'OC support': ['TNFSF11', 'TNFRSF11B', 'CSF1'],
    'bone formation': ['RUNX2', 'SP7', 'ALPL', 'COL1A1', 'BGLAP', 'SPP1'],
    'ER': ['ESR1', 'PGR', 'GREB1'],
}
# 解釈の補助（依頼外。老化の確認用。volcano のラベルには使わず表だけに出す）
EXTRA = {'proliferation/senescence check': ['MKI67', 'LMNB1', 'GDF15', 'MDM2']}
FOCUS_GENES = [g for v in FOCUS.values() for g in v]
CAT = {g: k for k, v in {**FOCUS, **EXTRA}.items() for g in v}

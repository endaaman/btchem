"""ヒト contrast の登録簿（prelim/contrasts.py と同じ形の dict。prelim/de.py の run() にそのまま渡す）。"""
from pathlib import Path

OUT = Path('out/prelim/human')
HC = {}


def reg(label, gse, test, ref, title, tissue, role='sub', subset=None, covariates=()):
    HC[label] = dict(label=label, gse=gse, test=test, ref=ref, subset=subset or {},
                     covariates=list(covariates), title=title, tissue=tissue, role=role)


reg('H_PostMeno', 'GSE141595', 'Post', 'Pre', 'Postmenopausal vs premenopausal (GSE141595)',
    'iliac biopsy, osteocyte-enriched; Pre は別バッチ（2022 追加）・遺伝子は Pre ファイルの共通分のみ', role='main')
reg('H_Dmab', 'GSE141595/post30', 'Post_Dmab', 'Post', 'Denosumab vs placebo, postmeno (GSE141595)',
    'iliac biopsy, osteocyte-enriched; 全遺伝子', role='main')
reg('H_Dmab_common', 'GSE141595', 'Post_Dmab', 'Post', 'Denosumab vs placebo, common genes only',
    '感度分析: H_PostMeno と同じ遺伝子集合で正規化')
reg('H_Dmab_CB', 'GSE141610', 'Post_Dmab', 'Post', 'Denosumab vs placebo, centrifuged bone (GSE141610)',
    'iliac biopsy, centrifuged (not digested)')
reg('H_Old', 'GSE72815', 'Old', 'Young', 'Old vs young women (GSE72815)',
    'iliac biopsy (whole), 73 vs 30 y', role='main')
reg('H_OldE', 'GSE72815', 'Old_E', 'Old', 'Old+E (3 wk) vs old (GSE72815)', 'iliac biopsy (whole)')
reg('H_GIOP', 'GSE276529', 'GIOP', 'NonGIOP', 'GIOP vs non-GIOP (GSE276529)',
    'femoral cortical bone, women 51-65 y', role='main')
reg('H_GIOP_bgadj', 'GSE276529', 'GIOP', 'NonGIOP', 'GIOP vs non-GIOP, bg-adjusted',
    '感度分析: design ~ bg + cond（bg = OR 遺伝子 CPM の log2 z）', covariates=['bg'])
reg('H_GIOP_marrowadj', 'GSE276529/marrowadj', 'GIOP', 'NonGIOP', 'GIOP vs non-GIOP, marrow-adjusted',
    '感度分析: design ~ bg + cond（この dir の bg 列 = 骨髄/血液スコア z）', covariates=['bg'])

# 横並び表・図で使う主 3 比較
MAIN3 = ['H_PostMeno', 'H_Old', 'H_GIOP']

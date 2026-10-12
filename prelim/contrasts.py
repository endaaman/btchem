"""解析する contrast の登録簿。de.py / volcano.py / enrich.py / overlap.py が共有する。

各 contrast は「比較する 2 群のサンプルだけ」を取り出して pydeseq2 にかける（データセット間で統一）。
"""

C = {}


def reg(label, gse, test, ref, subset=None, covariates=(), axis='', title='', tissue='',
        role='sub'):
    C[label] = dict(label=label, gse=gse, test=test, ref=ref, subset=subset or {},
                    covariates=list(covariates), axis=axis, title=title, tissue=tissue, role=role)


# --- A: 閉経 / エストロゲン欠乏 ---
reg('A_OVX_cort', 'GSE222752', 'OVX', 'Sham', {'compartment': 'Cortical'}, axis='A',
    title='OVX vs Sham (cortical)', tissue='femur cortical bone, ♀, 4 wk post-op, pooled 2 mice/sample',
    role='main')
reg('A_ERT_cort', 'GSE222752', 'OVX_ERT', 'OVX', {'compartment': 'Cortical'}, axis='A',
    title='OVX+ERT vs OVX (cortical)', tissue='femur cortical bone, ♀')
reg('A_OVX_trab', 'GSE222752', 'OVX', 'Sham', {'compartment': 'Trabecular'}, axis='A',
    title='OVX vs Sham (trabecular)', tissue='femur trabecular bone, ♀')
reg('A_ERT_trab', 'GSE222752', 'OVX_ERT', 'OVX', {'compartment': 'Trabecular'}, axis='A',
    title='OVX+ERT vs OVX (trabecular)', tissue='femur trabecular bone, ♀')
reg('A2_OVX_tibia', 'GSE238244', 'OVX_PBS', 'Sham', axis='A',
    title='OVX+PBS vs Sham (whole tibia)', tissue='whole tibia, ♀, 8 wk')

# --- B: 化学療法 ---
reg('B_Dox', 'GSE205237', 'Dox', 'Vehicle', axis='B',
    title='Doxorubicin vs Vehicle', tissue='tibia+femur w/o marrow, ♂, 30 d after 10 mg/kg',
    role='main')
reg('B2_Dox4x_OB', 'GSE296757', 'DOX4X', 'PBS', covariates=['batch'], axis='B',
    title='Dox 4x vs PBS (Col2.3+ osteoblasts)', tissue='sorted osteoblasts, ♀, 8 wk after 4x5 mg/kg')

# --- C: ステロイド ---
reg('C_DEX', 'GSE276391', 'DEX_Ctrl', 'Vehicle_Ctrl', axis='C',
    title='DEX vs Vehicle (no loading)', tissue='tibia, ♂, DEX 1 mg/kg 5 d', role='main')
reg('C_DEX_bgadj', 'GSE276391', 'DEX_Ctrl', 'Vehicle_Ctrl', covariates=['bg'], axis='C',
    title='DEX vs Vehicle (no loading, bg-adjusted)', tissue='tibia, ♂; design ~ bg + cond')
reg('C_DEX_load', 'GSE276391', 'DEX_Load', 'Vehicle_Load', axis='C',
    title='DEX vs Vehicle (loaded)', tissue='tibia, ♂, loaded')

# --- 加齢（DESeq2 済み表をそのまま取り込む。de.py では再計算しない） ---
AGING = {
    'Age_2y': ('All results for 2m-2y', '2 y vs 2 mo (cortical, ♂)'),
    'Age_2.5y': ('All results for 2m-2.5y', '2.5 y vs 2 mo (cortical, ♂)'),
}
for k, (_, t) in AGING.items():
    C[k] = dict(label=k, gse='NIHMS1776347', test=k, ref='2mo', subset={}, covariates=[],
                axis='A-aging', title=t, tissue='cortical bone, ♂', role='sub')

MAIN = {'A': 'A_OVX_cort', 'B': 'B_Dox', 'C': 'C_DEX'}

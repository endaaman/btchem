"""注目遺伝子パネル（事前定義。解析結果を見て足し引きしない）。"""

INFLAM = [
    'Tnfsf11', 'Tnfrsf11b', 'Il6', 'Il1b', 'Il1a', 'Tnf', 'Ccl2', 'Cxcl12',
    'S100a8', 'S100a9', 'Mmp9', 'Mmp13', 'Csf1', 'Sost', 'Dkk1',
]
APOP = [
    'Fas', 'Fasl', 'Casp3', 'Casp8', 'Casp9', 'Bax', 'Bak1', 'Bcl2', 'Bcl2l1',
    'Trp53', 'Cdkn1a', 'Cdkn2a', 'Serpine1', 'Il6', 'Mdm2', 'Pmaip1', 'Bbc3',
]

PANELS = {
    'inflam': INFLAM,   # 炎症・破骨細胞誘導
    'apop': APOP,       # アポトーシス・老化
}

PANEL_LABEL = {
    'inflam': 'Inflammation / osteoclastogenesis',
    'apop': 'Apoptosis / senescence',
}


def panel_of(gene):
    """遺伝子が属するパネル名のリスト（Il6 は両方）。"""
    return [k for k, v in PANELS.items() if gene in v]


ALL_GENES = list(dict.fromkeys(INFLAM + APOP))

"""Figure 4B（GSE296761, Blood 2026 = PMID 42008378）用の事前定義パネルと追加 gene set。解析結果を見て足し引きしない。

清水先生の依頼（10/14）の言い回し「osteoblast における骨基質環境の再編（脂肪形成偏向・ECM 再構築・骨量減少）」を
そのまま 4 つの軸に分け、論文の主張（IFN-γ / CD8+ T、動脈性血管の消失）を 1 軸足した。
- 遺伝子パネル: 遺伝子単位の表と、4 遺伝子以上のものは GSEA の手製 set（HAND_*）にも使う。
- 追加 gene set（MSigDB mouse 2025.1.Mm、prelim/chemo2/genesets.py と同じ GMT）: 72 set の既存コレクション
  （chemo2.genesets.collection()）に無い ECM の下位区分・脂肪/骨芽分化の正負制御・IFN-γ 応答。
  FDR は「追加 set + 手製 set」の小さなコレクションの中で計算する（72 set の FDR とは混ぜない）。
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
import genesets as GS  # noqa: E402

PANELS = {
    'adipo': ['Pparg', 'Cebpa', 'Adipoq', 'Lpl', 'Fabp4', 'Plin1', 'Cd36', 'Cfd', 'Lep'],
    'ecm': ['Col1a1', 'Col1a2', 'Col3a1', 'Col5a1', 'Col5a2', 'Col11a1', 'Fn1', 'Sparc', 'Postn', 'Bgn', 'Dcn',
            'Lum', 'Lox', 'Tnc', 'Eln'],
    'ecm_remodel': ['Mmp2', 'Mmp9', 'Mmp13', 'Mmp14', 'Timp1', 'Timp3'],
    'osteo': ['Runx2', 'Sp7', 'Dlx5', 'Satb2', 'Alpl', 'Ibsp', 'Spp1', 'Bglap', 'Bglap2', 'Ifitm5', 'Dmp1',
              'Phex', 'Sost', 'Wif1'],
    'oc_support': ['Tnfsf11', 'Tnfrsf11b', 'Csf1'],
    'hsc_niche': ['Cxcl12', 'Kitl', 'Lepr', 'Vcam1', 'Angpt1'],
    'ifng': ['Stat1', 'Irf1', 'Gbp2', 'Cxcl9', 'Cxcl10', 'Cd274', 'Psmb8', 'B2m', 'H2-K1', 'H2-Ab1', 'Ciita'],
}
PANEL_LABEL = {
    'adipo': 'Adipogenesis', 'ecm': 'ECM / collagen', 'ecm_remodel': 'ECM remodeling (MMP/TIMP)',
    'osteo': 'Osteogenesis', 'oc_support': 'Osteoclast support', 'hsc_niche': 'HSC niche factors',
    'ifng': 'IFN-γ response',
}
# 手製 set（GSEA / スコア用。oc_support は 3 遺伝子で、向きの違う OPG を含むので set にしない）
HAND = {f'HAND_{k.upper()}': v for k, v in PANELS.items() if k != 'oc_support'}
HAND_LABEL = {f'HAND_{k.upper()}': f'Panel {PANEL_LABEL[k]}' for k in PANELS if k != 'oc_support'}

EXTRA = [
    ('NABA_COLLAGENS', 'NABA collagens', 'ECM'),
    ('NABA_ECM_GLYCOPROTEINS', 'NABA ECM glycoproteins', 'ECM'),
    ('NABA_PROTEOGLYCANS', 'NABA proteoglycans', 'ECM'),
    ('NABA_ECM_REGULATORS', 'NABA ECM regulators', 'ECM'),
    ('GOBP_POSITIVE_REGULATION_OF_FAT_CELL_DIFFERENTIATION', 'GO +reg. fat cell diff.', 'adipo'),
    ('GOBP_NEGATIVE_REGULATION_OF_FAT_CELL_DIFFERENTIATION', 'GO −reg. fat cell diff.', 'adipo'),
    ('GOBP_BROWN_FAT_CELL_DIFFERENTIATION', 'GO brown fat cell diff.', 'adipo'),
    ('GOBP_POSITIVE_REGULATION_OF_OSTEOBLAST_DIFFERENTIATION', 'GO +reg. osteoblast diff.', 'bone formation'),
    ('GOBP_NEGATIVE_REGULATION_OF_OSTEOBLAST_DIFFERENTIATION', 'GO −reg. osteoblast diff.', 'bone formation'),
    ('GOBP_RESPONSE_TO_TYPE_II_INTERFERON', 'GO response to IFN-γ', 'inflammation'),
    ('GOBP_ARTERY_DEVELOPMENT', 'GO artery development', 'vascular'),
]

# 本文の表で軸ごとに並べる set（72 set コレクションの term + 追加 set + 手製）。
AXES = {
    'Adipogenesis': ['HALLMARK_ADIPOGENESIS', 'GOBP_FAT_CELL_DIFFERENTIATION', 'GOBP_POSITIVE_REGULATION_OF_FAT_CELL_DIFFERENTIATION',
                     'GOBP_NEGATIVE_REGULATION_OF_FAT_CELL_DIFFERENTIATION', 'GOBP_BROWN_FAT_CELL_DIFFERENTIATION',
                     'NICHE_ADIPO', 'HAND_ADIPO'],
    'ECM': ['NABA_CORE_MATRISOME', 'NABA_COLLAGENS', 'NABA_ECM_GLYCOPROTEINS', 'NABA_PROTEOGLYCANS', 'NABA_ECM_REGULATORS',
            'GOBP_COLLAGEN_FIBRIL_ORGANIZATION', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'HAND_ECM', 'HAND_ECM_REMODEL'],
    'Bone formation': ['GOBP_OSSIFICATION', 'GOBP_OSTEOBLAST_DIFFERENTIATION', 'GOBP_BONE_MINERALIZATION',
                       'GOBP_POSITIVE_REGULATION_OF_OSTEOBLAST_DIFFERENTIATION', 'GOBP_NEGATIVE_REGULATION_OF_OSTEOBLAST_DIFFERENTIATION',
                       'HAND_OSTEO'],
    'Osteoclast': ['GOBP_POSITIVE_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION', 'GOBP_OSTEOCLAST_DIFFERENTIATION',
                   'GOBP_BONE_RESORPTION', 'NICHE_OC_MARKERS'],
    'IFN / inflammation': ['HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_INTERFERON_ALPHA_RESPONSE',
                           'GOBP_RESPONSE_TO_TYPE_II_INTERFERON', 'HAND_IFNG', 'HALLMARK_INFLAMMATORY_RESPONSE',
                           'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'SAUL_SEN_MAYO', 'HALLMARK_P53_PATHWAY'],
    'Vascular / niche': ['HALLMARK_ANGIOGENESIS', 'GOBP_VASCULATURE_DEVELOPMENT', 'GOBP_SPROUTING_ANGIOGENESIS',
                         'GOBP_ARTERY_DEVELOPMENT', 'HAND_HSC_NICHE', 'NICHE_HSC_RETENTION'],
}

# 要約モジュール図（prelim/summary/summary_modules.py）と同じ 8 行 + 外した 4 行
SUMMARY_ROWS = [
    ('HALLMARK_INFLAMMATORY_RESPONSE', 'Inflammation', 'inflam'),
    ('HALLMARK_INTERFERON_ALPHA_RESPONSE', 'IFN-α response', 'inflam'),
    ('SAUL_SEN_MAYO', 'Senescence (SenMayo)', 'damage'),
    ('HALLMARK_P53_PATHWAY', 'p53 pathway', 'damage'),
    ('GOBP_POSITIVE_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION', 'Osteoclast support', 'niche'),
    ('GOBP_OSSIFICATION', 'Bone formation', 'niche'),
    ('NABA_CORE_MATRISOME', 'Extracellular matrix', 'niche'),
    ('HALLMARK_ANGIOGENESIS', 'Angiogenesis', 'niche'),
]
SUMMARY_DROPPED = [
    ('HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'TNF-α / NF-κB', 'inflam'),
    ('HALLMARK_ESTROGEN_RESPONSE_EARLY', 'Estrogen response', 'niche'),
    ('HALLMARK_APOPTOSIS', 'Apoptosis', 'damage'),
    ('HALLMARK_INTERFERON_GAMMA_RESPONSE', 'IFN-γ response', 'inflam'),
]


def extra_collection(upper=True):
    """追加 set + 手製 set（大文字シンボル）。戻り値: (sets, label dict)"""
    allm = {**GS.read_gmt('m2.all'), **GS.gobp()}
    sets, lab = {}, {}
    for k, l, _ in EXTRA:
        sets[k] = allm[k]
        lab[k] = l
    for k, v in HAND.items():
        sets[k] = [g.upper() for g in v]
        lab[k] = HAND_LABEL[k]
    return sets, lab


def base_collection():
    """chemo2 と同じ 72 set（要約図の FDR と揃えるため、このコレクションだけで回す）。"""
    return GS.collection()


def all_labels():
    _, meta = GS.collection()
    d = dict(zip(meta.term, meta.label))
    d.update({k: l for k, l, _ in EXTRA})
    d.update(HAND_LABEL)
    return d

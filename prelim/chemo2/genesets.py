"""niche 解析用の gene set（MSigDB mouse 2025.1.Mm + 手製 niche パネル）。事前定義。解析結果を見て足し引きしない。

GMT は https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/ から直接取得して
out/prelim/chemo2/msigdb/ に置いた（gseapy.Msigdb は lxml 依存で動かなかったため。依存は増やしていない）。
mouse 版そのものなのでオーソログ変換はしていない。照合は大文字化したシンボルで行う（既存 enrich.py と同じ）。

GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION は mouse 2025.1 に無い（GO 側の改名）ので、
ECM は NABA_CORE_MATRISOME と GOBP_COLLAGEN_FIBRIL_ORGANIZATION で代替した。
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
from genes import INFLAM, APOP  # noqa: E402

GMT = ROOT / 'out/prelim/chemo2/msigdb'
VER = 'v2025.1.Mm'

# MSigDB から取る niche 関連 set: (名前, 短いラベル, カテゴリ)
CURATED = [
    ('SAUL_SEN_MAYO', 'SenMayo', 'senescence/SASP'),
    ('REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP', 'Reactome SASP', 'senescence/SASP'),
    ('REACTOME_CELLULAR_SENESCENCE', 'Reactome cellular senescence', 'senescence/SASP'),
    ('GOBP_POSITIVE_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION', 'GO +reg. osteoclast diff.', 'osteoclast'),
    ('GOBP_OSTEOCLAST_DIFFERENTIATION', 'GO osteoclast diff.', 'osteoclast'),
    ('GOBP_BONE_RESORPTION', 'GO bone resorption', 'osteoclast'),
    ('WP_OSTEOCLAST_SIGNALING', 'WP osteoclast signaling', 'osteoclast'),
    ('NABA_CORE_MATRISOME', 'NABA core matrisome', 'ECM'),
    ('GOBP_COLLAGEN_FIBRIL_ORGANIZATION', 'GO collagen fibril org.', 'ECM'),
    ('GOBP_VASCULATURE_DEVELOPMENT', 'GO vasculature dev.', 'vascular'),
    ('GOBP_SPROUTING_ANGIOGENESIS', 'GO sprouting angiogenesis', 'vascular'),
    ('GOBP_FAT_CELL_DIFFERENTIATION', 'GO fat cell diff.', 'adipo'),
    ('GOBP_OSSIFICATION', 'GO ossification', 'bone formation'),
    ('GOBP_OSTEOBLAST_DIFFERENTIATION', 'GO osteoblast diff.', 'bone formation'),
    ('GOBP_BONE_MINERALIZATION', 'GO bone mineralization', 'bone formation'),
    ('GOBP_MYELOID_LEUKOCYTE_ACTIVATION', 'GO myeloid leukocyte activation', 'myeloid'),
    ('GOBP_HEMATOPOIETIC_STEM_CELL_HOMEOSTASIS', 'GO HSC homeostasis', 'HSC niche'),
]

# Hallmark のうち niche の読みに使うもの（全 50 も別図で出す）
HALLMARK_NICHE = [
    ('HALLMARK_INFLAMMATORY_RESPONSE', 'H inflammatory response', 'inflammation'),
    ('HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'H TNFa/NF-kB', 'inflammation'),
    ('HALLMARK_IL6_JAK_STAT3_SIGNALING', 'H IL6/JAK/STAT3', 'inflammation'),
    ('HALLMARK_INTERFERON_ALPHA_RESPONSE', 'H IFN-a response', 'inflammation'),
    ('HALLMARK_INTERFERON_GAMMA_RESPONSE', 'H IFN-g response', 'inflammation'),
    ('HALLMARK_COMPLEMENT', 'H complement', 'inflammation'),
    ('HALLMARK_ANGIOGENESIS', 'H angiogenesis', 'vascular'),
    ('HALLMARK_HYPOXIA', 'H hypoxia', 'vascular'),
    ('HALLMARK_ADIPOGENESIS', 'H adipogenesis', 'adipo'),
    ('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'H EMT', 'ECM'),
    ('HALLMARK_TGF_BETA_SIGNALING', 'H TGF-b', 'ECM'),
    ('HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'H Wnt/b-catenin', 'bone formation'),
    ('HALLMARK_P53_PATHWAY', 'H p53 pathway', 'senescence/SASP'),
    ('HALLMARK_APOPTOSIS', 'H apoptosis', 'senescence/SASP'),
    ('HALLMARK_E2F_TARGETS', 'H E2F targets', 'proliferation'),
    ('HALLMARK_G2M_CHECKPOINT', 'H G2M checkpoint', 'proliferation'),
    ('HALLMARK_ESTROGEN_RESPONSE_EARLY', 'H estrogen early', 'estrogen'),
    ('HALLMARK_ESTROGEN_RESPONSE_LATE', 'H estrogen late', 'estrogen'),
    ('HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'H OXPHOS', 'composition?'),
    ('HALLMARK_MYOGENESIS', 'H myogenesis', 'composition?'),
]

# 手製 niche パネル。既存 prelim/genes.py（INFLAM / APOP）と重なる遺伝子はそちらから引いて、
# 新しく足すもの（NEW）だけをここに書く。サブパネルは役割で分けた（向きの違う因子を混ぜない）。
_known = set(INFLAM) | set(APOP)
NICHE_PANEL = {
    'NICHE_HSC_RETENTION': ['Cxcl12', 'Kitl', 'Angpt1', 'Lepr', 'Vcam1', 'Spp1'],
    'NICHE_OC_MARKERS': ['Ctsk', 'Acp5', 'Dcstamp', 'Mmp9', 'Mmp13'],
    'NICHE_OC_SUPPORT': ['Tnfsf11', 'Csf1'],          # RANKL, M-CSF（OPG は逆向きなので別に置く）
    'NICHE_OPG': ['Tnfrsf11b'],
    'NICHE_INFLAM_CHEMOKINE': ['Il6', 'Il1b', 'Tnf', 'Ccl2', 'Cxcl1', 'Cxcl2', 'S100a8', 'S100a9'],
    'NICHE_WNT_INHIB': ['Sost', 'Dkk1'],
    'NICHE_ADIPO': ['Pparg', 'Cebpa', 'Adipoq', 'Lpl'],
    'NICHE_SENESCENCE': ['Cdkn1a', 'Cdkn2a', 'Serpine1', 'Igfbp3'],
}
NICHE_LABEL = {
    'NICHE_HSC_RETENTION': 'Panel HSC retention',
    'NICHE_OC_MARKERS': 'Panel osteoclast markers',
    'NICHE_INFLAM_CHEMOKINE': 'Panel inflam./chemokine',
    'NICHE_ADIPO': 'Panel adipocyte',
    'NICHE_SENESCENCE': 'Panel senescence',
}
NEW_GENES = [g for v in NICHE_PANEL.values() for g in v if g not in _known]  # genes.py に無いもの
# GSEA にかけるのは 4 遺伝子以上のサブパネルだけ（2 遺伝子以下は遺伝子単位の表で見る）
NICHE_GSEA = {k: v for k, v in NICHE_PANEL.items() if len(v) >= 4}


def read_gmt(name):
    d = {}
    with open(GMT / f'{name}.{VER}.symbols.gmt') as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            d[x[0]] = [g.upper() for g in x[2:] if g]
    return d


def hallmark():
    return read_gmt('mh.all')


def gobp():
    return read_gmt('m5.go.bp')


def collection():
    """Hallmark 50 + CURATED + 手製パネル（大文字）。戻り値: (sets, meta DataFrame)"""
    import pandas as pd
    h = hallmark()
    allm = {**read_gmt('m2.all'), **gobp()}
    sets, meta = {}, []
    for k, v in h.items():
        sets[k] = v
        lab = dict((a, b) for a, b, _ in HALLMARK_NICHE).get(k, 'H ' + k.replace('HALLMARK_', '').lower().replace('_', ' '))
        cat = dict((a, c) for a, _, c in HALLMARK_NICHE).get(k, 'hallmark (other)')
        meta.append(dict(term=k, label=lab, category=cat, source='Hallmark', niche=k in dict((a, 1) for a, *_ in HALLMARK_NICHE)))
    for k, lab, cat in CURATED:
        if k not in allm:
            raise KeyError(k)
        sets[k] = allm[k]
        meta.append(dict(term=k, label=lab, category=cat, source='MSigDB M2/M5', niche=True))
    for k, v in NICHE_GSEA.items():
        sets[k] = [g.upper() for g in v]
        meta.append(dict(term=k, label=NICHE_LABEL[k], category='hand panel', source='hand', niche=True))
    meta = pd.DataFrame(meta)
    meta['size_total'] = meta.term.map(lambda t: len(sets[t]))
    return sets, meta

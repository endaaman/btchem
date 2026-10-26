"""Fig 4A v3（案1 デノスマブ / 案4 Kaya）共通の部品。

- パス、カウント読み込み、pydeseq2 の DE（prelim/de.py の run() と同じフィルタ・設定を、
  data/seq に samples.csv を書かずに DataFrame で受けられるようにコピーしたもの）
- 速い近似 DE（log2 CPM の平均差と Welch t。split-half と置換に使う）
- gene set（MSigDB 2025.1.Hs の GMT。out/prelim/invitro/msigdb/ にあるものを読むだけ）
- 細胞種マーカー（組成）、即時初期遺伝子・熱ショック（検体処理の痕跡）

既存の prelim/*.py・prelim/human/・prelim/chemo2/・prelim/invitro/ は import して使うだけで変更しない。
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/human'))
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
sys.path.insert(0, str(ROOT / 'prelim/invitro'))

SEQ = ROOT / 'data/seq'
OUT = ROOT / 'out/prelim/fig4a_v3'
HDE = ROOT / 'out/prelim/human/de'          # 既存のヒト DE（読むだけ）
MSIG = ROOT / 'out/prelim/invitro/msigdb'   # MSigDB 2025.1.Hs GMT（読むだけ）
VER = 'v2025.1.Hs'

# de.py と同じ値
MIN_COUNT = 10
MIN_TOTAL = 15
COLS = ['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']


def counts(sub):
    """data/seq/<sub>/counts.tsv.gz と samples.csv。"""
    x = pd.read_csv(SEQ / sub / 'counts.tsv.gz', sep='\t', index_col=0)
    s = pd.read_csv(SEQ / sub / 'samples.csv', dtype=str).set_index('sample')
    return x, s


def keep_genes(x, cond):
    """prelim/de.py と同じ低発現フィルタ（edgeR filterByExpr 既定相当）。x: genes x samples。"""
    lib = x.sum(axis=0)
    cpm_cut = MIN_COUNT / (lib.median() / 1e6)
    cpm = x / lib * 1e6
    n_min = cond.value_counts().min()
    return ((cpm >= cpm_cut).sum(axis=1) >= n_min) & (x.sum(axis=1) >= MIN_TOTAL)


def deseq(x, s, test, ref, covariates=(), label=''):
    """pydeseq2（prelim/de.py の run() と同じ設定）。x: genes x samples, s: index=sample, group 列, 共変量列。"""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    s = s[s.group.isin([test, ref])].copy()
    s['cond'] = np.where(s.group == test, 'test', 'ref')
    x = x[s.index]
    keep = keep_genes(x, s.cond)
    xx = x.loc[keep].T
    cov = list(covariates)
    for c in cov:
        s[c] = s[c].astype(float)
    design = '~' + ' + '.join(cov + ['cond'])
    dds = DeseqDataSet(counts=xx, metadata=s[cov + ['cond']], design=design, refit_cooks=True,
                       inference=DefaultInference(n_cpus=8), quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=['cond', 'test', 'ref'], inference=DefaultInference(n_cpus=8), quiet=True)
    st.summary()
    r = st.results_df.reset_index()
    r = r.rename(columns={r.columns[0]: 'gene'})[COLS]
    n = s.cond.value_counts().to_dict()
    print(f'{label}: {test} (n={n.get("test")}) vs {ref} (n={n.get("ref")}), design={design}, '
          f'tested={int(keep.sum())}, DEG(padj<0.05,|lfc|>=1)={int(((r.padj < 0.05) & (r.log2FoldChange.abs() >= 1)).sum())}')
    return r


def size_factors(x):
    """DESeq2 の median-of-ratios（x: genes x samples の生カウント）。"""
    lx = np.log(x.replace(0, np.nan))
    geo = lx.mean(axis=1)
    ok = geo.notna() & np.isfinite(geo)
    return np.exp((lx.loc[ok].sub(geo[ok], axis=0)).median(axis=0))


def logcpm(x, prior=1.0):
    """median-of-ratios で正規化した log2(count + prior)。"""
    sf = size_factors(x)
    return np.log2(x / sf + prior)


def fast_de(L, a_cols, b_cols):
    """log2 発現行列 L（genes x samples）で a vs b の平均差と Welch t。速い近似（split-half・置換用）。"""
    A, B = L[a_cols].values, L[b_cols].values
    ma, mb = A.mean(1), B.mean(1)
    va, vb = A.var(1, ddof=1), B.var(1, ddof=1)
    se = np.sqrt(va / A.shape[1] + vb / B.shape[1])
    se = np.where(se > 0, se, np.nan)
    return pd.DataFrame({'lfc': ma - mb, 't': (ma - mb) / se}, index=L.index)


def wald_rank(de):
    """DE 表 → GSEA 用の順位（Wald 統計量 = log2FC / lfcSE、大文字シンボル、重複は平均）。既存 gsea_iv.ranking() と同じ。"""
    d = de.dropna(subset=['lfcSE', 'gene'])
    d = d[d.lfcSE > 0]
    return (d.log2FoldChange / d.lfcSE).groupby(d.gene.astype(str).str.upper()).mean().sort_values(ascending=False)


def read_gmt(name):
    d = {}
    with open(MSIG / f'{name}.{VER}.symbols.gmt') as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            d[x[0]] = [g.upper() for g in x[2:] if g]
    return d


def collections():
    """GSEA 用の 3 コレクション（FDR はコレクションごとに計算する）。
    hallmark: Hallmark 50 + SenMayo（既存ヒト解析と同じく Hallmark に足す）
    gobp: GO BP 全体、reactome: C2 の REACTOME_ 全体（サイズ 15–500 は prerank 側で絞る）"""
    h = read_gmt('h.all')
    c2 = read_gmt('c2.all')
    h['SAUL_SEN_MAYO'] = c2['SAUL_SEN_MAYO']
    go = read_gmt('c5.go.bp')
    rea = {k: v for k, v in c2.items() if k.startswith('REACTOME_')}
    return {'hallmark': h, 'gobp': go, 'reactome': rea}


# ---- 細胞種マーカー（組成の目安。マーカーの log2 CPM 平均で見る） ----
MARKERS_H = {
    'osteoclast': ['ACP5', 'CTSK', 'OSCAR', 'DCSTAMP', 'ATP6V0D2', 'SIGLEC15', 'CALCR'],
    'osteoblast': ['COL1A1', 'COL1A2', 'BGLAP', 'SP7', 'IBSP', 'ALPL', 'IFITM5'],
    'osteocyte': ['SOST', 'DMP1', 'PHEX', 'MEPE', 'FGF23'],
    'erythroid': ['HBB', 'HBA1', 'HBA2', 'ALAS2', 'GYPA', 'CA1'],
    'granulocyte': ['MPO', 'ELANE', 'AZU1', 'PRTN3', 'DEFA3', 'LTF', 'CAMP', 'CD177'],
    'mono/macrophage': ['CD14', 'CD68', 'CD163', 'CSF1R', 'MRC1', 'C1QA', 'LYZ'],
    'B cell': ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'VPREB3'],
    'plasma cell': ['JCHAIN', 'MZB1', 'TNFRSF17'],
    'T/NK': ['CD3E', 'CD3D', 'CD2', 'NKG7', 'GZMB'],
    'megakaryocyte': ['PF4', 'PPBP', 'ITGA2B'],
    'endothelial': ['PECAM1', 'CDH5', 'VWF', 'EMCN', 'KDR', 'TEK'],
    'smooth muscle/pericyte': ['ACTA2', 'MYH11', 'CNN1', 'TAGLN', 'RGS5'],
    'skeletal muscle': ['ACTA1', 'CKM', 'MYH1', 'MYH2', 'MYH7', 'MYL1', 'TNNT3'],
    'adipocyte': ['ADIPOQ', 'PLIN1', 'LEP', 'CIDEC', 'GPD1'],
    'LepR+ stroma': ['LEPR', 'CXCL12', 'KITLG', 'VCAM1', 'EBF3'],
    'fibroblast/periosteal': ['PDGFRA', 'DCN', 'LUM', 'COL3A1', 'POSTN', 'PRRX1'],
    'mast cell': ['CPA3', 'TPSAB1', 'MS4A2', 'KIT'],
    'proliferation': ['MKI67', 'TOP2A', 'CENPF'],
    'immediate early (processing)': ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'IER2', 'IER3', 'ATF3', 'DUSP1',
                                     'NR4A1', 'ZFP36', 'KLF2', 'KLF4', 'KLF6', 'PPP1R15A', 'BTG2', 'SOCS3'],
    'heat shock (processing)': ['HSPA1A', 'HSPA1B', 'HSPA6', 'DNAJB1', 'HSPH1', 'HSPB1', 'BAG3'],
}
# 骨髄・血液スコア（H_GIOP_marrowadj と同じ定義: HBB/HBA1/MPO/ELANE の log2 CPM 平均を z 化）
MARROW_GENES = ['HBB', 'HBA1', 'MPO', 'ELANE']
# 検体処理（即時初期・熱ショック）の遺伝子。leading edge に占める割合を見る
IEG = set(MARKERS_H['immediate early (processing)']) | {'EGR2', 'EGR3', 'IER5', 'NR4A2', 'NR4A3', 'FOSL1', 'FOSL2',
                                                       'JUND', 'MAFF', 'CCN1', 'CYR61', 'CCN2', 'CTGF', 'GADD45B',
                                                       'KLF10', 'TRIB1', 'RHOB', 'SGK1', 'CEBPD', 'NFKBIA', 'TNFAIP3'}
HSP = set(MARKERS_H['heat shock (processing)']) | {'HSPA8', 'HSP90AA1', 'HSPD1', 'HSPE1', 'DNAJA1', 'HSPA4L'}

# 図・表で使う主要 set の短い名前
KEY_SETS = {
    'HALLMARK_INFLAMMATORY_RESPONSE': 'Inflammatory response',
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB': 'TNF-α/NF-κB',
    'HALLMARK_IL6_JAK_STAT3_SIGNALING': 'IL-6/JAK/STAT3',
    'HALLMARK_INTERFERON_ALPHA_RESPONSE': 'IFN-α response',
    'HALLMARK_INTERFERON_GAMMA_RESPONSE': 'IFN-γ response',
    'HALLMARK_APOPTOSIS': 'Apoptosis',
    'HALLMARK_P53_PATHWAY': 'p53 pathway',
    'SAUL_SEN_MAYO': 'SenMayo',
}


def cpm(x):
    return x / x.sum() * 1e6


def scores(x, markers=MARKERS_H):
    """カテゴリ別スコア（そのデータにあるマーカーの log2(CPM+1) 平均）。返り値: samples x categories と使った遺伝子。"""
    c = np.log2(cpm(x) + 1)
    out, used = {}, {}
    for k, gs in markers.items():
        have = [g for g in gs if g in c.index]
        if have:
            out[k] = c.loc[have].mean()
            used[k] = have
    return pd.DataFrame(out), used


def zscore(v):
    v = pd.Series(v, dtype=float)
    return (v - v.mean()) / v.std(ddof=1)


def save(fig, name, sub=''):
    import common  # prelim/common.py
    d = OUT / sub if sub else OUT
    d.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(d / f'{name}.{ext}', dpi=common.DPI, bbox_inches='tight', pad_inches=0.02)
    print('saved', d / f'{name}.png|pdf')

"""案4: 現行 Fig 4A（Kaya et al. 2022, NIHMS1776347, PRJNA695408）を「骨全体の環境」の目で読み直す。

  uv run python -I prelim/fig4a_v3/kaya.py

入力は著者の DESeq2 結果表（data/NIHMS1776347-supplement-tS1.xlsx の "All results for 2m-1y / 2m-2y / 2m-2.5y"）だけ。
生データ（SRA）は使っていない（再定量の見積もりは docs/prelim_fig4a_v3.md）。
月齢は論文本文と ENA の sample_alias で確認: 2m = 2 か月（n=3）, 1y = 11 か月, 2y = 23 か月, 2.5y = 30 か月（各 n=4）、雄 C57BL/6、同一コロニー。

やること:
1. 表の限界: 行数、FDR 空欄（DESeq2 の独立フィルタ）と p 値空欄（Cook 外れ値）の数・発現量、注目遺伝子の欠け。
2. 試料の実体: 細胞種マーカーの存在量（baseMean から擬似 CPM）を、他のマウス骨データ（GSE222752 皮質骨、GSE205237 脛骨+大腿骨、
   GSE276391 脛骨丸ごと）の CPM と並べる。マーカーの加齢での変化（組成の変化）。
3. 骨環境の軸（炎症・破骨支持・骨形成/ECM・血管・脂肪化・老化・アポトーシス）: prelim/chemo2/genesets.py の collection()
   （Hallmark 50 + MSigDB 17 + 手製パネル 5、mouse 2025.1）で prerank GSEA を 1y / 2y / 2.5y について回す（chemo2 と同じ設定）。
   加えて、共通の 2m を引き算した近似の「2.5y vs 1y」「2y vs 1y」（統計量 = (lfc_a − lfc_b) / sqrt(se_a² + se_b²)。
   2m 群の共有で実際の分散はこれより小さいので保守的）で、**成熟（2m → 1y）と加齢（1y → 2.5y）を分ける**。
   2 か月齢の雄マウスは骨格成長の途中にあるため、2m を基準にした変化には成長の停止・成熟が混ざる。

出力: out/prelim/fig4a_v3/kaya/{table_limits.csv, composition_abundance.csv, composition_age.csv, gsea_<label>.csv,
      axes_nes.csv, panel_genes.csv, kaya_axes.png|pdf}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common4a import OUT, ROOT, counts, save  # noqa: E402
import common  # noqa: E402
import genesets as GS  # noqa: E402  prelim/chemo2/genesets.py
import genes as PG  # noqa: E402  prelim/genes.py

D = OUT / 'kaya'
XLSX = ROOT / 'data/NIHMS1776347-supplement-tS1.xlsx'
SHEETS = {'Age_1y': 'All results for 2m-1y', 'Age_2y': 'All results for 2m-2y', 'Age_2.5y': 'All results for 2m-2.5y'}
MONTHS = {'Age_1y': '11 mo', 'Age_2y': '23 mo', 'Age_2.5y': '30 mo'}

MARKERS_M = {
    'osteocyte': ['Sost', 'Dmp1', 'Phex', 'Mepe', 'Fgf23'],
    'osteoblast': ['Bglap', 'Bglap2', 'Col1a1', 'Sp7', 'Ibsp', 'Alpl', 'Ifitm5'],
    'osteoclast': ['Acp5', 'Ctsk', 'Oscar', 'Dcstamp', 'Atp6v0d2', 'Calcr', 'Mmp9'],
    'chondrocyte/growth plate': ['Col2a1', 'Acan', 'Col10a1', 'Matn3'],
    'erythroid': ['Hbb-bs', 'Hba-a1', 'Hba-a2', 'Alas2', 'Gypa', 'Car1'],
    'granulocyte': ['Mpo', 'Elane', 'Ltf', 'Camp', 'Ngp', 'Chil3'],
    'mono/macrophage': ['Cd14', 'Cd68', 'Adgre1', 'Csf1r', 'Mrc1', 'C1qa', 'Lyz2'],
    'B cell': ['Cd19', 'Ms4a1', 'Cd79a', 'Cd79b', 'Vpreb3'],
    'plasma cell': ['Jchain', 'Mzb1', 'Tnfrsf17'],
    'T/NK': ['Cd3e', 'Cd3d', 'Cd2', 'Nkg7', 'Gzmb'],
    'endothelial': ['Pecam1', 'Cdh5', 'Emcn', 'Kdr', 'Tek', 'Vwf'],
    'smooth muscle/pericyte': ['Acta2', 'Myh11', 'Cnn1', 'Tagln', 'Rgs5'],
    'skeletal muscle': ['Acta1', 'Ckm', 'Myh1', 'Myh2', 'Myh4', 'Mylpf', 'Tnnt3'],
    'adipocyte': ['Adipoq', 'Plin1', 'Lep', 'Cidec', 'Cfd'],
    'LepR+ stroma': ['Lepr', 'Cxcl12', 'Kitl', 'Vcam1', 'Ebf3'],
    'fibroblast/periosteal': ['Pdgfra', 'Dcn', 'Lum', 'Col3a1', 'Postn', 'Prrx1'],
    'mast cell': ['Cpa3', 'Mcpt4', 'Cma1', 'Kit'],
    'proliferation': ['Mki67', 'Top2a', 'Cenpf'],
    'immediate early': ['Fos', 'Fosb', 'Jun', 'Junb', 'Egr1', 'Ier2', 'Atf3', 'Dusp1', 'Nr4a1', 'Zfp36'],
    'heat shock': ['Hspa1a', 'Hspa1b', 'Dnajb1', 'Hsph1', 'Hspb1'],
}
OTHER = [  # (counts dir, 表示名, サンプルの絞り込み)。対照群だけを使う
    ('GSE222752', 'GSE222752 femur cortical, sham (♀)', lambda s: (s.compartment == 'Cortical') & (s.group == 'Sham')),
    ('GSE205237', 'GSE205237 tibia+femur, marrow removed, vehicle (♂)', lambda s: s.group == 'Vehicle'),
    ('GSE276391', 'GSE276391 whole tibia, vehicle (♂)', lambda s: s.group == 'Vehicle_Ctrl'),
]

# 骨環境の軸 → set（chemo2 collection の term）
AXES = [
    ('Inflammation', 'HALLMARK_INFLAMMATORY_RESPONSE', 'H inflammatory response'),
    ('Inflammation', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'H TNF-α/NF-κB'),
    ('Inflammation', 'HALLMARK_IL6_JAK_STAT3_SIGNALING', 'H IL-6/JAK/STAT3'),
    ('Inflammation', 'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'H IFN-α response'),
    ('Inflammation', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'H IFN-γ response'),
    ('Inflammation', 'GOBP_MYELOID_LEUKOCYTE_ACTIVATION', 'GO myeloid leukocyte activation'),
    ('Inflammation', 'NICHE_INFLAM_CHEMOKINE', 'Panel inflam./chemokine'),
    ('Osteoclast support', 'GOBP_POSITIVE_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION', 'GO +reg. osteoclast diff.'),
    ('Osteoclast support', 'GOBP_BONE_RESORPTION', 'GO bone resorption'),
    ('Osteoclast support', 'NICHE_OC_MARKERS', 'Panel osteoclast markers'),
    ('Bone formation / ECM', 'GOBP_OSSIFICATION', 'GO ossification'),
    ('Bone formation / ECM', 'GOBP_OSTEOBLAST_DIFFERENTIATION', 'GO osteoblast diff.'),
    ('Bone formation / ECM', 'GOBP_BONE_MINERALIZATION', 'GO bone mineralization'),
    ('Bone formation / ECM', 'NABA_CORE_MATRISOME', 'NABA core matrisome'),
    ('Bone formation / ECM', 'GOBP_COLLAGEN_FIBRIL_ORGANIZATION', 'GO collagen fibril org.'),
    ('Vascular', 'HALLMARK_ANGIOGENESIS', 'H angiogenesis'),
    ('Vascular', 'GOBP_VASCULATURE_DEVELOPMENT', 'GO vasculature dev.'),
    ('Adiposity', 'HALLMARK_ADIPOGENESIS', 'H adipogenesis'),
    ('Adiposity', 'GOBP_FAT_CELL_DIFFERENTIATION', 'GO fat cell diff.'),
    ('Adiposity', 'NICHE_ADIPO', 'Panel adipocyte'),
    ('Senescence', 'SAUL_SEN_MAYO', 'SenMayo'),
    ('Senescence', 'REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP', 'Reactome SASP'),
    ('Senescence', 'REACTOME_CELLULAR_SENESCENCE', 'Reactome cellular senescence'),
    ('Apoptosis / DNA damage', 'HALLMARK_APOPTOSIS', 'H apoptosis'),
    ('Apoptosis / DNA damage', 'HALLMARK_P53_PATHWAY', 'H p53 pathway'),
    ('Proliferation', 'HALLMARK_E2F_TARGETS', 'H E2F targets'),
    ('HSC niche', 'NICHE_HSC_RETENTION', 'Panel HSC retention'),
    ('Composition?', 'HALLMARK_MYOGENESIS', 'H myogenesis'),
    ('Composition?', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'H OXPHOS'),
]


def read_sheet(sheet):
    n = pd.read_excel(XLSX, sheet_name=sheet, header=None, skiprows=8).iloc[:, :9]
    n.columns = ['ensembl', 'gene', 'entrez', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    n = n[n.ensembl.astype(str).str.startswith('ENSMUSG')].copy()
    for k in ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']:
        n[k] = pd.to_numeric(n[k], errors='coerce')
    return n


def limits(T):
    rows = []
    panel = PG.ALL_GENES
    for lb, n in T.items():
        sym = n[n.gene.notna() & (n.gene.astype(str) != 'nan')]
        p_na = n.pvalue.isna()
        q_na = n.padj.isna() & ~p_na
        rows.append(dict(label=lb, rows=len(n), with_symbol=len(sym), baseMean_zero=int((n.baseMean == 0).sum()),
                         pvalue_blank=int(p_na.sum()), fdr_blank_only=int(q_na.sum()),
                         fdr_blank_baseMean_median=n.loc[q_na, 'baseMean'].median(),
                         fdr_blank_baseMean_max=n.loc[q_na, 'baseMean'].max(),
                         pvalue_blank_baseMean_median=n.loc[p_na & (n.baseMean > 0), 'baseMean'].median(),
                         deg_fdr01=int((n.padj < 0.1).sum()),
                         deg_strict=int(((n.padj < 0.05) & (n.log2FoldChange.abs() >= 1)).sum()),
                         panel_missing=','.join(g for g in panel if g not in set(sym.gene.astype(str))),
                         panel_fdr_blank=','.join(g for g in panel if g in set(sym.loc[sym.padj.isna(), 'gene'].astype(str)))))
    return pd.DataFrame(rows)


def pseudo_cpm(n):
    s = n.dropna(subset=['gene']).copy()
    s['gene'] = s.gene.astype(str)
    s = s.groupby('gene').baseMean.sum()
    return s / s.sum() * 1e6


def cat_score(cpm_series, markers):
    out = {}
    for k, gs in markers.items():
        have = [g for g in gs if g in cpm_series.index]
        out[k] = (np.log2(cpm_series.reindex(have) + 1).mean(), len(have))
    return out


def main():
    D.mkdir(parents=True, exist_ok=True)
    T = {lb: read_sheet(sh) for lb, sh in SHEETS.items()}
    L = limits(T)
    L.to_csv(D / 'table_limits.csv', index=False)
    print(L.T.to_string())

    # ---- 2. 試料の実体（存在量） ----
    ab = []
    pc = pseudo_cpm(T['Age_2.5y'])  # 2m + 2.5y の平均正規化カウント
    for k, (v, n) in cat_score(pc, MARKERS_M).items():
        ab.append(dict(dataset='Kaya 2022 humerus cortical (♂; baseMean 2m+30mo)', cat=k, log2cpm=v, n_genes=n))
    for sub, name, f in OTHER:
        x, s = counts(sub)
        s = s[f(s)]
        x = x[s.index]
        c = (x / x.sum() * 1e6).median(axis=1)
        for k, (v, n) in cat_score(c, MARKERS_M).items():
            ab.append(dict(dataset=name, cat=k, log2cpm=v, n_genes=n))
    AB = pd.DataFrame(ab)
    AB.to_csv(D / 'composition_abundance.csv', index=False)
    print(AB.pivot(index='cat', columns='dataset', values='log2cpm').round(2).to_string())
    gene_ab = []
    for k, gs in MARKERS_M.items():
        for g in gs:
            gene_ab.append(dict(cat=k, gene=g, kaya_pseudo_cpm=pc.get(g, np.nan)))
    pd.DataFrame(gene_ab).to_csv(D / 'composition_marker_genes_kaya.csv', index=False)

    # ---- 組成の加齢変化（マーカー lfc） ----
    ca = []
    for lb, n in T.items():
        s = n.dropna(subset=['gene']).copy()
        s['gene'] = s.gene.astype(str)
        s = s.drop_duplicates('gene').set_index('gene')
        for k, gs in MARKERS_M.items():
            have = [g for g in gs if g in s.index and s.loc[g, 'baseMean'] > 0]
            q = s.loc[have]
            ca.append(dict(label=lb, cat=k, n_genes=len(have), mean_lfc=q.log2FoldChange.mean(),
                           n_up_fdr01=int(((q.padj < 0.1) & (q.log2FoldChange > 0)).sum()),
                           n_down_fdr01=int(((q.padj < 0.1) & (q.log2FoldChange < 0)).sum()),
                           genes=','.join(f'{g}:{q.loc[g, "log2FoldChange"]:+.2f}' for g in have)))
    CA = pd.DataFrame(ca)
    CA.to_csv(D / 'composition_age.csv', index=False)
    print(CA.pivot(index='cat', columns='label', values='mean_lfc').round(2).to_string())

    # ---- 3. GSEA（成熟 vs 加齢の分解） ----
    sets, meta = GS.collection()
    J = None
    for lb, n in T.items():
        s = n.dropna(subset=['gene', 'lfcSE']).copy()
        s['gene'] = s.gene.astype(str).str.upper()
        s = s[s.lfcSE > 0].groupby('gene')[['log2FoldChange', 'lfcSE']].mean()
        s.columns = [f'lfc_{lb}', f'se_{lb}']
        J = s if J is None else J.join(s, how='inner')
    ranks = {lb: (J[f'lfc_{lb}'] / J[f'se_{lb}']).sort_values(ascending=False) for lb in T}
    for a, b, name in [('Age_2.5y', 'Age_1y', 'Late_2.5y_vs_1y'), ('Age_2y', 'Age_1y', 'Late_2y_vs_1y')]:
        ranks[name] = ((J[f'lfc_{a}'] - J[f'lfc_{b}']) / np.sqrt(J[f'se_{a}'] ** 2 + J[f'se_{b}'] ** 2)).sort_values(ascending=False)
    res = {}
    for lb, rnk in ranks.items():
        f = D / f'gsea_{lb}.csv'
        if not f.exists():
            r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=4, max_size=1000, seed=0, threads=8,
                           outdir=None, no_plot=True, verbose=False).res2d
            r.assign(label=lb, n_ranked=len(rnk)).to_csv(f, index=False)
        r = pd.read_csv(f).set_index('Term')
        res[lb] = r
    rows = []
    for ax_, term, lab in AXES:
        row = dict(axis=ax_, term=term, label=lab)
        for lb, r in res.items():
            row[f'NES_{lb}'] = r.loc[term, 'NES'] if term in r.index else np.nan
            row[f'FDR_{lb}'] = r.loc[term, 'FDR q-val'] if term in r.index else np.nan
        rows.append(row)
    AX = pd.DataFrame(rows)
    AX.to_csv(D / 'axes_nes.csv', index=False)
    print(AX[['axis', 'label'] + [c for c in AX.columns if c.startswith('NES_')]].round(2).to_string(index=False))
    # 既存 chemo2（Age_2y, Age_2.5y; inner join 前の全遺伝子で回した値）との一致を確認
    old = pd.read_csv(ROOT / 'out/prelim/chemo2/nes_wide.csv', header=[0, 1], index_col=0)
    for lb in ['Age_2y', 'Age_2.5y']:
        m = pd.concat([old[('NES', lb)], res[lb].NES], axis=1, join='inner').dropna()
        print(f'check vs chemo2 {lb}: r = {np.corrcoef(m.iloc[:, 0], m.iloc[:, 1])[0, 1]:.3f}, n = {len(m)}')

    # ---- 注目パネル遺伝子 ----
    pg = []
    for g in PG.ALL_GENES + ['Igfbp3', 'Lepr', 'Cxcl12', 'Kitl', 'Spp1', 'Pparg', 'Adipoq', 'Acp5', 'Ctsk']:
        row = dict(gene=g, panel=('inflam' if g in PG.INFLAM else '') + ('apop' if g in PG.APOP else ''))
        for lb in T:
            q = T[lb][T[lb].gene.astype(str) == g]
            if len(q):
                q = q.iloc[0]
                row[f'lfc_{lb}'], row[f'padj_{lb}'], row[f'baseMean_{lb}'] = q.log2FoldChange, q.padj, q.baseMean
        G = g.upper()
        if G in J.index:
            row['late_2.5y_vs_1y_lfc'] = J.loc[G, 'lfc_Age_2.5y'] - J.loc[G, 'lfc_Age_1y']
            row['late_2.5y_vs_1y_z'] = ranks['Late_2.5y_vs_1y'].get(G, np.nan)
        pg.append(row)
    PGT = pd.DataFrame(pg)
    PGT.to_csv(D / 'panel_genes.csv', index=False)
    print(PGT[['gene', 'lfc_Age_1y', 'padj_Age_1y', 'lfc_Age_2y', 'padj_Age_2y', 'lfc_Age_2.5y', 'padj_Age_2.5y',
               'late_2.5y_vs_1y_lfc', 'late_2.5y_vs_1y_z']].round(3).to_string(index=False))
    concordance_aging2(res)
    plot(AX, CA)


def concordance_aging2(res):
    """案2（別担当 prelim/aging2、GSE235942 脛骨 22 vs 6 か月）の NES と、Kaya の各比較の NES の Spearman。
    同じ chemo2 collection（72 set）で回してあるので term で突き合わせる。aging2 の出力は読むだけ（無ければ飛ばす）。"""
    from scipy.stats import spearmanr
    src = ROOT / 'out/prelim/aging2/gsea'
    labs = {'G235_M_age': 'male 22 vs 6 mo', 'G235_M_age_comp': 'male, composition-adjusted', 'G235_F_age': 'female 22 vs 6 mo'}
    if not all((src / f'{k}.csv').exists() for k in labs):
        print('aging2 GSEA not found; skip concordance')
        return
    A = {k: pd.read_csv(src / f'{k}.csv').set_index('Term')['NES'].astype(float) for k in labs}
    rows = []
    for lb, r in res.items():
        row = dict(kaya=lb)
        for k, v in A.items():
            t = r.index.intersection(v.index)
            row[f'rho_{k}'] = spearmanr(r.loc[t, 'NES'].astype(float), v[t])[0]
            row['n_sets'] = len(t)
        rows.append(row)
    C = pd.DataFrame(rows)
    C.to_csv(D / 'concordance_aging2.csv', index=False)
    print(C.round(2).to_string(index=False))


def plot(AX, CA):
    common.setup()
    cols = [('Age_1y', '11 vs 2\nmo'), ('Age_2y', '23 vs 2\nmo'), ('Age_2.5y', '30 vs 2\nmo'),
            ('Late_2y_vs_1y', '23 vs 11\nmo*'), ('Late_2.5y_vs_1y', '30 vs 11\nmo*')]
    fig = plt.figure(figsize=(common.mm(170), common.mm(125)))
    ax = fig.add_axes([0.335, 0.07, 0.30, 0.82])
    W = AX[[f'NES_{c}' for c, _ in cols]].values.astype(float)
    F = AX[[f'FDR_{c}' for c, _ in cols]].values.astype(float)
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    ax.imshow(np.ma.masked_invalid(W), cmap=cmap, vmin=-3, vmax=3, aspect='auto')
    for i in range(W.shape[0]):
        for j in range(W.shape[1]):
            if np.isnan(W[i, j]):
                continue
            if F[i, j] < 0.05:
                ax.text(j, i, f'{W[i, j]:+.1f}'.replace('-', '−'), ha='center', va='center', fontsize=5.5,
                        color='white' if abs(W[i, j]) > 2 else 'black', fontweight='bold')
            elif F[i, j] < 0.25:
                ax.text(j, i, f'{W[i, j]:+.1f}'.replace('-', '−'), ha='center', va='center', fontsize=5, color='#555555')
    ax.set_yticks(range(len(AX)), AX.label, fontsize=6)
    ax.set_xticks(range(len(cols)), [c[1] for c in cols], fontsize=6)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    prev = None
    for i, a in enumerate(AX.axis):
        if prev is not None and a != prev:
            ax.axhline(i - 0.5, c='k', lw=0.5)
        if a != prev:
            fig.text(0.005, ax.transData.transform((0, i - 0.5))[1] / fig.bbox.height, a, fontsize=6, fontweight='bold',
                     va='top', ha='left')
        prev = a
    ax.axvline(2.5, c='k', lw=1.0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title('Prerank NES (bold: FDR < 0.05; grey: < 0.25)', fontsize=7, pad=20)
    fig.text(0.335, 0.015, '* approximate: (lfc_a − lfc_b) / √(se_a² + se_b²) from the shared-2-mo author tables',
             fontsize=5.5, color='#444444')
    fig.text(0.0, 0.985, 'A', fontsize=12, fontweight='bold', va='top')
    # B: 組成マーカーの加齢変化
    ax2 = fig.add_axes([0.83, 0.07, 0.16, 0.82])
    cats = list(MARKERS_M)
    M = CA.pivot(index='cat', columns='label', values='mean_lfc').reindex(index=cats, columns=['Age_1y', 'Age_2y', 'Age_2.5y'])
    ax2.imshow(M.values, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            if not np.isnan(v):
                ax2.text(j, i, f'{v:+.1f}'.replace('-', '−'), ha='center', va='center', fontsize=5.5,
                         color='white' if abs(v) > 1.3 else 'black')
    ax2.set_yticks(range(len(cats)), cats, fontsize=6)
    ax2.set_xticks(range(3), ['11 mo', '23 mo', '30 mo'], fontsize=6)
    ax2.xaxis.tick_top()
    ax2.tick_params(length=0)
    for sp in ax2.spines.values():
        sp.set_visible(False)
    ax2.set_title('Cell-type markers:\nmean log2FC vs 2 mo', fontsize=7, pad=14)
    fig.text(0.665, 0.985, 'B', fontsize=12, fontweight='bold', va='top')
    save(fig, 'kaya_axes', 'kaya')
    plt.close(fig)


if __name__ == '__main__':
    main()

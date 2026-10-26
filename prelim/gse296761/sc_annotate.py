"""GSE296752: sc_cluster.py の Leiden（res 1.0, 39 クラスタ）に細胞種を付け、細胞種の割合を PBS / DOX で比べる。

細胞種の割り当て（ASSIGN）は out/prelim/gse296761/sc/cluster_markers.csv（マーカー平均）と cluster_top_genes.csv
（Wilcoxon 上位）を見て決めた。論文（Blood 2026）の呼び名に合わせ、Lepr+ 間質を
  MSC-Adipo: Lepr/Cxcl12/Kitl 高 + Adipoq/Lpl 高、Wif1/Alpl/Sp7 低（Baccin 2020 の Adipo-CAR に相当）
  MSC-Osteo: Lepr/Cxcl12 あり + Wif1/Alpl/Satb2/Runx2/Spp1 高（Osteo-CAR に相当）
とした。Lepr 低・Spp1/Ibsp/Mmp13 高のクラスタは Osteoprogenitor、Bglap/Col1a1/Ifitm5 高は Osteoblast。
内皮は arteriolar（Gja5/Bmx/Sox17/Sema3g/Ly6a）、capillary（Emcn/Kdr 高・Efnb2 あり・Stab2 弱）、
venular（Vwf/Selp/Ackr1/Nr2f2）、sinusoidal（Stab2/Flt4）、lymphatic（Prox1/Lyve1）に分けた。

低品質クラスタ: クラスタの UMI 中央値 < 1,500（他は 4,000–23,000）。主解析（割合・擬似バルク）からは除き、
親の細胞種に戻した場合を感度分析として並べる（PARENT）。

CD8+ TEM: 標本は系統（CD3/CD4/CD8/B220/CD11b/Gr1/Ter119/CD71）と CD45/cKit を MACS で 2 回除去しているので、
T 細胞はほぼ居ないはず。Cd3e/Cd3d/Cd3g・Cd8a/Cd8b1・Ifng の生 UMI を数えて確かめる（t_cell_check.csv）。

  uv run python -I prelim/gse296761/sc_annotate.py      # repo root から（sc_cluster.py の後）
出力: out/prelim/gse296761/sc/{adata_annot.h5ad, celltype_counts.csv, celltype_counts_sensitivity.csv,
      msc_score_split.csv, t_cell_check.csv}, out/prelim/gse296761/sc_umap.png|pdf, sc_dotplot.png|pdf,
      sc_proportions.png|pdf
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
import common  # noqa: E402

OUT = ROOT / 'out/prelim/gse296761'
SC = OUT / 'sc'
G = ROOT / 'data/seq/GSE296752'

ASSIGN = {
    '6': 'MSC-Adipo', '8': 'MSC-Adipo', '9': 'MSC-Adipo',
    '10': 'MSC-Osteo', '15': 'MSC-Osteo',
    '14': 'Osteoprogenitor',
    '0': 'Osteoblast',
    '2': 'EC arteriolar',
    '11': 'EC capillary',
    '4': 'EC venular',
    '12': 'EC sinusoidal',
    '38': 'EC lymphatic',
    '17': 'Pericyte',
    '16': 'SMC',
    '22': 'Fibroblast', '24': 'Fibroblast', '28': 'Fibroblast', '30': 'Fibroblast',
    '27': 'Chondrocyte',
    '1': 'Schwann/glia', '32': 'Schwann/glia',
    '18': 'Muscle satellite',
    '7': 'Hematopoietic (residual)', '20': 'Hematopoietic (residual)', '23': 'Hematopoietic (residual)',
    '35': 'Hematopoietic (residual)',
    # 低品質（UMI 中央値 < 1,500）
    '3': 'Low quality', '5': 'Low quality', '13': 'Low quality', '21': 'Low quality', '25': 'Low quality',
    '29': 'Low quality', '31': 'Low quality', '34': 'Low quality', '36': 'Low quality', '37': 'Low quality',
    # doublet 様（2 系統のマーカーが同じ細胞に乗る）・小さすぎる
    '19': 'Doublet/unclear', '26': 'Doublet/unclear', '33': 'Doublet/unclear',
}
# 低品質クラスタを親の細胞種に戻すとき（感度分析）。21（巨核/血小板）と 29（分裂中の赤芽球系）は造血、5 は内皮 + 間質の混ざり
PARENT = {'3': 'EC arteriolar', '5': 'Doublet/unclear', '13': 'MSC-Adipo', '21': 'Hematopoietic (residual)',
          '25': 'Fibroblast', '29': 'Hematopoietic (residual)', '31': 'EC sinusoidal', '34': 'MSC-Adipo',
          '36': 'SMC', '37': 'Schwann/glia'}
ORDER = ['MSC-Adipo', 'MSC-Osteo', 'Osteoprogenitor', 'Osteoblast', 'EC arteriolar', 'EC capillary', 'EC venular',
         'EC sinusoidal', 'EC lymphatic', 'Pericyte', 'SMC', 'Fibroblast', 'Chondrocyte', 'Schwann/glia',
         'Muscle satellite', 'Hematopoietic (residual)', 'Doublet/unclear', 'Low quality']
# 依頼で名前の挙がった細胞種（割合の図で強調する）
REQUESTED = ['MSC-Adipo', 'MSC-Osteo', 'Osteoblast', 'EC arteriolar', 'EC capillary', 'EC venular', 'EC sinusoidal',
             'Pericyte']
COLORS = {'PBS': '#7F7F7F', 'DOX': '#C44E52'}

# Lepr+ MSC を 2 群に分ける別法（感度分析）: 細胞ごとの Adipo / Osteo スコアの差で分ける
ADIPO_SIG = ['Adipoq', 'Lpl', 'Kng2', 'Esm1', 'Hp', 'Cxcl12', 'Lepr']
OSTEO_SIG = ['Wif1', 'Alpl', 'Sp7', 'Satb2', 'Tnfrsf19', 'Spp1', 'Runx2']

DOT = ['Lepr', 'Cxcl12', 'Kitl', 'Adipoq', 'Lpl', 'Wif1', 'Alpl', 'Satb2', 'Sp7', 'Spp1', 'Ibsp', 'Bglap', 'Ifitm5',
       'Dmp1', 'Pecam1', 'Cdh5', 'Gja5', 'Bmx', 'Sema3g', 'Emcn', 'Kdr', 'Vwf', 'Ackr1', 'Stab2', 'Prox1', 'Rgs5',
       'Kcnj8', 'Higd1b', 'Acta2', 'Myh11', 'Dcn', 'Clec3b', 'Prg4', 'Col2a1', 'Mpz', 'Pax7', 'Ptprc', 'Hbb-bs', 'Pf4']


def t_cell_check(a):
    rows = []
    for f, s in [('GSM8975917_PBS', 'PBS'), ('GSM8975918_DOX', 'DOX')]:
        r = sc.read_10x_h5(G / f'{f}_filtered_feature_bc_matrix.h5')
        r.var_names_make_unique()
        r.obs_names = [f'{s}_{x}' for x in r.obs_names]
        r = r[r.obs_names.isin(a.obs_names)]
        g = ['Cd3e', 'Cd3d', 'Cd3g', 'Cd8a', 'Cd8b1', 'Cd4', 'Ifng', 'Ptprc', 'Ccl5', 'Cd44', 'Sell', 'Il7r']
        D = pd.DataFrame({x: np.asarray(r[:, x].X.todense()).ravel() for x in g}, index=r.obs_names)
        cd3 = D[['Cd3e', 'Cd3d', 'Cd3g']].sum(axis=1) > 0
        cd8 = D[['Cd8a', 'Cd8b1']].sum(axis=1) > 0
        rows.append(dict(sample=s, cells=len(D), cells_Ptprc_pos=int((D.Ptprc > 0).sum()),
                         cells_Cd3_any=int(cd3.sum()), cells_Cd8_any=int(cd8.sum()),
                         cells_Cd3_and_Cd8=int((cd3 & cd8).sum()),
                         cells_Cd3_Cd8_Cd44_noSell=int((cd3 & cd8 & (D.Cd44 > 0) & (D.Sell == 0)).sum()),
                         cells_Ifng=int((D.Ifng > 0).sum()),
                         **{f'UMI_{x}': int(D[x].sum()) for x in g}))
    t = pd.DataFrame(rows)
    t.to_csv(SC / 't_cell_check.csv', index=False)
    print(t.T.to_string())


def props(lab, samp):
    ct = pd.crosstab(lab, samp).reindex(columns=['PBS', 'DOX']).fillna(0).astype(int)
    for s in ['PBS', 'DOX']:
        ct[f'{s}_%'] = ct[s] / ct[s].sum() * 100
    ct['log2_ratio_DOX_PBS'] = np.log2((ct['DOX_%'] + 1e-9) / (ct['PBS_%'] + 1e-9))
    return ct


def main():
    common.setup()
    a = sc.read_h5ad(SC / 'adata_clustered.h5ad')
    a.obs['leiden'] = a.obs.leiden.astype(str)
    miss = set(a.obs.leiden) - set(ASSIGN)
    assert not miss, miss
    a.obs['celltype'] = pd.Categorical(a.obs.leiden.map(ASSIGN), ORDER)
    a.obs['celltype_sens'] = a.obs.leiden.map(lambda c: PARENT.get(c, ASSIGN[c]))
    a.obs['sample'] = pd.Categorical(a.obs['sample'].astype(str), ['PBS', 'DOX'])

    # 割合: 主解析は「低品質・doublet・残存造血を除いた非血球の細胞」を分母にする
    keep = ~a.obs.celltype.isin(['Low quality', 'Doublet/unclear', 'Hematopoietic (residual)'])
    P = props(a.obs.celltype[keep].astype(str), a.obs['sample'][keep].astype(str)).reindex(
        [c for c in ORDER if c in set(a.obs.celltype[keep])])
    # Lepr+ MSC の中での Adipo / Osteo の比
    msc = P.loc[['MSC-Adipo', 'MSC-Osteo']]
    P.to_csv(SC / 'celltype_counts.csv')
    keep2 = ~a.obs.celltype_sens.isin(['Low quality', 'Doublet/unclear', 'Hematopoietic (residual)'])
    P2 = props(a.obs.celltype_sens[keep2], a.obs['sample'][keep2].astype(str))
    Pall = props(a.obs.celltype.astype(str), a.obs['sample'].astype(str))
    S = pd.concat({'main_excl_lowq': P[['PBS', 'DOX', 'PBS_%', 'DOX_%', 'log2_ratio_DOX_PBS']],
                   'lowq_to_parent': P2[['PBS', 'DOX', 'PBS_%', 'DOX_%', 'log2_ratio_DOX_PBS']],
                   'all_cells_incl_lowq_hem': Pall[['PBS', 'DOX', 'PBS_%', 'DOX_%', 'log2_ratio_DOX_PBS']]}, axis=1)
    S.to_csv(SC / 'celltype_counts_sensitivity.csv')
    pd.set_option('display.width', 250)
    print(S.round(2).to_string())
    print('MSC-Osteo / (Adipo+Osteo):',
          {s: round(msc.loc['MSC-Osteo', s] / msc[s].sum(), 3) for s in ['PBS', 'DOX']})

    # スコアで Lepr+ MSC を分ける（感度分析）
    sc.tl.score_genes(a, [g for g in ADIPO_SIG if g in a.var_names], score_name='adipo_sig', random_state=0)
    sc.tl.score_genes(a, [g for g in OSTEO_SIG if g in a.var_names], score_name='osteo_sig', random_state=0)
    m = a.obs.celltype.isin(['MSC-Adipo', 'MSC-Osteo'])
    d = a.obs.loc[m, ['sample', 'celltype', 'adipo_sig', 'osteo_sig']].copy()
    d['split'] = np.where(d.osteo_sig > d.adipo_sig, 'osteo>adipo', 'adipo>=osteo')
    ms = pd.crosstab(d['sample'], d['split'])
    ms['osteo_frac'] = ms['osteo>adipo'] / ms.sum(axis=1)
    ms['cluster_osteo_frac'] = [msc.loc['MSC-Osteo', s] / msc[s].sum() for s in ms.index]
    ms.to_csv(SC / 'msc_score_split.csv')
    print(ms.round(3).to_string())

    t_cell_check(a)

    # 図: UMAP（細胞種 / 群）
    a.obs['celltype_plot'] = a.obs.celltype.astype(str)
    pal = dict(zip(ORDER, plt.get_cmap('tab20').colors[:len(ORDER)]))
    pal['Low quality'] = '#DDDDDD'
    pal['Doublet/unclear'] = '#BBBBBB'
    fig, axes = plt.subplots(1, 2, figsize=(common.mm(170), common.mm(78)))
    ap = a[np.random.default_rng(0).permutation(a.n_obs)]  # 描画順をランダムに（片方の群が上に重ならないように）
    for ax, (key, title, p) in zip(axes, [('celltype_plot', 'Cell types (Leiden 1.0 → markers)', pal),
                                          ('sample', 'Library (2 mice pooled each)', COLORS)]):
        sc.pl.umap(ap, color=key, ax=ax, show=False, frameon=False, s=2, palette=p, title=title,
                   legend_fontsize=5.5 if key == 'celltype_plot' else 7)
    fig.suptitle('GSE296752 non-hematopoietic BM + bone cells, PBS ×4 vs DOX 5 mg/kg ×4, 8 wk after (♀)', fontsize=8)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'sc_umap.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)

    b = a[a.obs.celltype.isin([c for c in ORDER if c not in ('Low quality', 'Doublet/unclear')])].copy()
    b.obs['celltype'] = b.obs.celltype.cat.remove_unused_categories()
    dp = sc.pl.dotplot(b, [g for g in DOT if g in b.var_names], groupby='celltype', show=False, figsize=(11, 4.6),
                       standard_scale='var')
    fig = dp['mainplot_ax'].figure
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'sc_dotplot.{ext}', dpi=300, bbox_inches='tight')
    plt.close('all')

    # 割合の図（主解析の分母）
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(62)))
    x = np.arange(len(P))
    w = 0.38
    ax.bar(x - w / 2, P['PBS_%'], w, color=COLORS['PBS'], label='PBS ×4')
    ax.bar(x + w / 2, P['DOX_%'], w, color=COLORS['DOX'], label='DOX ×4')
    for i, (c, r) in enumerate(P.iterrows()):
        ax.text(i, max(r['PBS_%'], r['DOX_%']) + 0.4, f"{r.log2_ratio_DOX_PBS:+.1f}".replace('-', '−'),
                ha='center', fontsize=6, color='#C44E52' if r.log2_ratio_DOX_PBS > 0 else '#1F4E79',
                fontweight='bold' if c in REQUESTED else 'normal')
    ax.set_xticks(x, [c + (' *' if c in REQUESTED else '') for c in P.index], rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('% of non-hematopoietic cells')
    ax.legend(frameon=False, fontsize=7)
    ax.set_title('Cell-type composition (n = 1 pooled library per group; numbers = log$_2$ DOX/PBS share). '
                 '* cell types named in the request', fontsize=7, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'sc_proportions.{ext}', dpi=common.DPI, bbox_inches='tight')
    plt.close(fig)

    a.obs = a.obs.drop(columns=['celltype_plot'])
    a.write_h5ad(SC / 'adata_annot.h5ad', compression='gzip')


if __name__ == '__main__':
    main()

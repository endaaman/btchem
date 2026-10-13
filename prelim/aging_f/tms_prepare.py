"""GSE132040（Tabula Muris Senis bulk）から Bone と Marrow を抜き出し、性別・組成を確認する。

  uv run python -I prelim/aging_f/tms_prepare.py
入力: data/seq/GSE132040/GSE132040_190214_*.csv.gz（htseq-count、gencode vM19、全 947 サンプル）、*_MACA_Bulk_metadata.csv.gz
出力（out/prelim/aging_f/tms/）: counts_bone_marrow.tsv.gz, samples.csv, composition.csv, composition.png, qc_pca.png
組成マーカー（事前定義）: 下の COMP。Marrow（洗い流した骨髄）を同じ表に並べて、Bone にどれだけ骨髄が入っているかの目安にする。
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
import common  # noqa: E402

G = ROOT / 'data/seq/GSE132040'
OUT = ROOT / 'out/prelim/aging_f/tms'
COMP = {
    'osteoblast/osteocyte': ['Bglap', 'Bglap2', 'Col1a1', 'Sp7', 'Ifitm5', 'Dmp1', 'Sost', 'Phex'],
    'chondrocyte/growth plate': ['Col2a1', 'Col10a1', 'Acan', 'Ihh'],
    'marrow stroma (CAR)': ['Lepr', 'Cxcl12', 'Kitl', 'Adipoq'],
    'hematopoietic': ['Ptprc', 'Lyz2', 'S100a8', 'Cd74', 'Cd79a', 'Mpo'],
    'erythroid/blood': ['Hbb-bs', 'Hba-a1', 'Alas2'],
    'skeletal muscle': ['Acta1', 'Ckm', 'Myh4', 'Tnnt3'],
    'adipocyte': ['Adipoq', 'Cfd', 'Lep', 'Plin1'],
    'osteoclast': ['Acp5', 'Ctsk', 'Mmp9'],
}
SEXG = ['Xist', 'Ddx3y', 'Eif2s3y', 'Kdm5d', 'Uty']


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(G / 'GSE132040_MACA_Bulk_metadata.csv.gz')
    m['tissue'] = m['source name'].str.rsplit('_', n=1).str[0]
    m = m[m.tissue.isin(['Bone', 'Marrow'])].copy()
    m = m.rename(columns={'Sample name': 'sample', 'characteristics: age': 'age', 'characteristics: sex': 'sex'})
    m['age'] = m.age.astype(int)
    m['plate'] = m['sample'].str.extract(r'(Plate\d)')[0]
    fn = next(G.glob('GSE132040_190214_*.csv.gz'))
    c = pd.read_csv(fn, index_col=0)
    c.columns = c.columns.str.replace('.gencode.vM19', '', regex=False)
    c = c[~c.index.str.startswith('__')]
    miss = set(m['sample']) - set(c.columns)
    print('missing in counts:', miss)
    m = m[m['sample'].isin(c.columns)]
    c = c[m['sample']]
    # 同じシンボルの重複は合計（htseq の gene 行はシンボル）
    c = c.groupby(level=0).sum()
    m['lib'] = c.sum().values
    cpm = c / c.sum() * 1e6
    for g in SEXG:
        m[g + '_CPM'] = cpm.loc[g].values if g in cpm.index else np.nan
    # 雄でも Xist が 5–30 CPM 出る（384 ウェルプレートのサンプル間の持ち込み・index hopping と思われる）ので、
    # Y 4 遺伝子の合計で判定する: 雌 = Xist > 100 かつ Y 合計 < 10、雄 = Y 合計 > 40 かつ Xist < 50
    ysum = m[[g + '_CPM' for g in SEXG[1:]]].sum(axis=1)
    m['Ysum_CPM'] = ysum
    m['sex_expr'] = np.where((m.Xist_CPM > 100) & (ysum < 10), 'f', np.where((ysum > 40) & (m.Xist_CPM < 50), 'm', '?'))
    for k, gl in COMP.items():
        gl = [g for g in gl if g in cpm.index]
        m['comp:' + k] = np.log2(cpm.loc[gl] + 1).mean().values
    m['detected_genes'] = (c > 0).sum().values
    m = m.drop(columns=[x for x in m.columns if x in ['title', 'organism', 'molecule', 'description', 'processed data file ',
                                                       'raw file', 'BioSample', 'Instrument Model',
                                                       'characteristics: developmental stage']])
    m.to_csv(OUT / 'samples.csv', index=False)
    c.to_csv(OUT / 'counts_bone_marrow.tsv.gz', sep='\t')
    pd.set_option('display.width', 250)
    print(pd.crosstab([m.tissue, m.sex], m.sex_expr))
    print(m[m.sex != m.sex_expr][['sample', 'tissue', 'sex', 'age', 'Xist_CPM', 'Ysum_CPM', 'lib']].to_string())
    cc = [k for k in m.columns if k.startswith('comp:')]
    t = m.groupby(['tissue', 'sex', 'age'])[cc + ['lib', 'detected_genes']].mean().round(2)
    t.to_csv(OUT / 'composition.csv')
    print(t.to_string())
    print(pd.crosstab([m.tissue, m.sex, m.age], m.plate))

    common.setup()
    fig, axes = plt.subplots(2, 4, figsize=(common.mm(170), common.mm(85)), sharex=True)
    for ax, k in zip(axes.ravel(), cc):
        for (tis, sx), col, mk in [(('Bone', 'f'), '#C44E52', 'o'), (('Bone', 'm'), '#4A7BA7', 'o'),
                                   (('Marrow', 'f'), '#C44E52', 'x'), (('Marrow', 'm'), '#4A7BA7', 'x')]:
            d = m[(m.tissue == tis) & (m.sex == sx)]
            ax.scatter(d.age + (0.3 if sx == 'm' else -0.3), d[k], s=6, c=col, marker=mk, lw=0.6,
                       label=f'{tis} {"♀" if sx == "f" else "♂"}')
        ax.set_title(k.replace('comp:', ''), fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axes[1]:
        ax.set_xlabel('age (months)', fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel('mean log2(CPM+1)', fontsize=7)
    axes[0, 0].legend(fontsize=5, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / 'composition.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    # PCA（Bone のみ、log CPM、上位分散 2000）
    b = m[m.tissue == 'Bone']
    l = np.log2(cpm[b['sample']] + 1)
    l = l.loc[l.var(axis=1).sort_values().index[-2000:]]
    l = l.sub(l.mean(axis=1), axis=0).T
    u, s, vt = np.linalg.svd(l.values, full_matrices=False)
    pc = u[:, :3] * s[:3]
    ve = s ** 2 / (s ** 2).sum()
    fig, axes = plt.subplots(1, 2, figsize=(common.mm(170), common.mm(70)))
    for ax, (i, j) in zip(axes, [(0, 1), (0, 2)]):
        sca = ax.scatter(pc[:, i], pc[:, j], c=b.age, cmap='viridis', s=14,
                         marker='o', edgecolors=np.where(b.sex == 'f', '#C44E52', '#4A7BA7'), linewidths=0.8)
        ax.set_xlabel(f'PC{i + 1} ({ve[i] * 100:.0f}%)')
        ax.set_ylabel(f'PC{j + 1} ({ve[j] * 100:.0f}%)')
    fig.colorbar(sca, ax=axes, label='age (months)', shrink=0.7)
    axes[0].set_title('GSE132040 Bone (edge: red ♀ / blue ♂)', fontsize=7)
    fig.savefig(OUT / 'qc_pca.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    pcd = pd.DataFrame(pc, columns=['PC1', 'PC2', 'PC3'], index=b['sample'])
    # PC と組成スコアの相関（何が主軸か）
    print(pd.concat([pcd.reset_index(drop=True), b[cc + ['age']].reset_index(drop=True)], axis=1).corr(method='spearman')
          .loc[['PC1', 'PC2', 'PC3'], cc + ['age']].round(2).to_string())


if __name__ == '__main__':
    main()

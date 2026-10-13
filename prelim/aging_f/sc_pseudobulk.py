"""GSE210584: 間質の細胞種を付け（ASSIGN）、個体ごとの擬似バルク → pydeseq2（中年 vs 若齢）。

ASSIGN の根拠は out/prelim/aging_f/sc/GSE210584_cluster_markers_stroma.csv（sc_cluster.py の 2 段階目）:
  CAR (LepR+)   : Lepr 2.8–3.5、Cxcl12 5.1–5.5、Kitl 2.3–3.4、Adipoq 2.3–3.2（クラスタ 2,3,4,5,6,7,8,9,10,15）
  Osteo-lineage : 13 = 骨芽細胞（Bglap 2.1、Col1a1 3.8、Ifitm5 0.64、Sp7 0.43）、
                  14 = osteo-CAR / 骨前駆（Sp7 0.90、Alpl 2.2、Spp1 3.9、Bglap 1.75、Lepr 1.4）
  Endothelial   : 0 = 動脈系（Ly6a 2.5、Cdh5 2.7）、1 = 類洞（Emcn 1.9、Vcam1 2.4）、18 = 動脈系（Ly6a 2.8、Kitl 2.3）
  除外          : 11/12/17（Kit・Mki67 が高い造血前駆の混入）、16/19/20/21/22/25（血球・内皮マーカーとの混在 = ダブレット疑い）、
                  23（周皮細胞/SMC、全 69 細胞で少なすぎ）、24（内皮マーカーが中途半端、233 細胞）
擬似バルク: 個体 × 細胞種で raw UMI を合計。個体ごとの細胞数が MINCELLS 未満の細胞種は DE しない。
DE: deutil.deseq（de.py と同じフィルタ）、~cond（中年 vs 若齢）。

  uv run python -I prelim/aging_f/sc_pseudobulk.py
出力: out/prelim/aging_f/sc/GSE210584_celltype_counts.csv, GSE210584_umap_celltype.png, GSE210584_pb_pca.png,
      out/prelim/aging_f/de/sc210584_<celltype>.csv, de_summary_sc210584.csv, out/prelim/aging_f/sc/GSE210584_pb_<ct>.csv.gz
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
sys.path.insert(0, str(ROOT / 'prelim/aging_f'))
import common  # noqa: E402
import deutil  # noqa: E402

SC = ROOT / 'out/prelim/aging_f/sc'
DE = ROOT / 'out/prelim/aging_f/de'
MINCELLS = 50

ASSIGN = {**{k: 'CAR' for k in ['2', '3', '4', '5', '6', '7', '8', '9', '10', '15']},
          '13': 'Osteo', '14': 'Osteo',
          '0': 'Endothelial', '1': 'Endothelial', '18': 'Endothelial'}


def main():
    DE.mkdir(parents=True, exist_ok=True)
    s = sc.read_h5ad(SC / 'GSE210584_stroma.h5ad')
    s.obs['celltype'] = s.obs['sub'].astype(str).map(ASSIGN).fillna('excluded').astype('category')
    s.obs['subtype'] = s.obs['sub'].astype(str).map({'13': 'Osteoblast', '14': 'Osteo-CAR/progenitor'}).fillna(s.obs.celltype.astype(str))
    n = pd.crosstab(s.obs['sample'].astype(str), s.obs['celltype'].astype(str))
    n2 = pd.crosstab(s.obs['sample'].astype(str), s.obs['subtype'].astype(str))
    n = n.join(n2[['Osteoblast', 'Osteo-CAR/progenitor']])
    n['group'] = s.obs.assign(sample=s.obs['sample'].astype(str)).groupby('sample').group.first().reindex(n.index)
    n.to_csv(SC / 'GSE210584_celltype_counts.csv')
    print(n.to_string())

    common.setup()
    fig, axes = plt.subplots(1, 2, figsize=(common.mm(170), common.mm(75)))
    sc.pl.umap(s, color='celltype', ax=axes[0], show=False, frameon=False, s=3, title='GSE210584 non-hematopoietic')
    sc.pl.umap(s, color='group', ax=axes[1], show=False, frameon=False, s=3, title='Young (n=5) / middle-aged (n=9)',
               palette={'young': '#4A7BA7', 'middle': '#C44E52'})
    fig.savefig(SC / 'GSE210584_umap_celltype.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    mk = ['Lepr', 'Cxcl12', 'Kitl', 'Adipoq', 'Igf1', 'Runx2', 'Sp7', 'Alpl', 'Spp1', 'Bglap', 'Ifitm5', 'Dmp1',
          'Cdh5', 'Pecam1', 'Emcn', 'Ly6a', 'Kit', 'Mki67', 'Ptprc', 'Rgs5', 'Acta2']
    sc.pl.dotplot(s, mk, groupby='sub', show=False)
    plt.savefig(SC / 'GSE210584_dotplot_sub.png', dpi=200, bbox_inches='tight')
    plt.close('all')

    X = s.layers['counts'].tocsr()
    genes = np.array(s.var_names)
    rows = []
    pcs = {}
    for ct in ['CAR', 'Osteo', 'Endothelial']:
        pb = {}
        for smp in n.index:
            idx = np.where((s.obs['sample'] == smp).values & (s.obs.celltype == ct).values)[0]
            pb[smp] = (len(idx), np.asarray(X[idx].sum(axis=0)).ravel())
        x = pd.DataFrame({k: v[1] for k, v in pb.items()}, index=genes).T
        ncell = pd.Series({k: v[0] for k, v in pb.items()})
        x.T.to_csv(SC / f'GSE210584_pb_{ct}.csv.gz')
        meta = pd.DataFrame({'cond': np.where(n.loc[x.index, 'group'] == 'middle', 'test', 'ref')}, index=x.index)
        ok = ncell >= MINCELLS
        if (meta.cond[ok] == 'test').sum() < 3 or (meta.cond[ok] == 'ref').sum() < 3:
            print('skip', ct, ncell.to_dict())
            rows.append(dict(celltype=ct, n_mid=(meta.cond[ok] == 'test').sum(), n_young=(meta.cond[ok] == 'ref').sum(),
                             cells_min=ncell.min(), note='skipped'))
            continue
        x, meta = x[ok], meta[ok]
        r, ng = deutil.deseq(x, meta, '~cond', ['cond', 'test', 'ref'], n_min=meta.cond.value_counts().min())
        r.to_csv(DE / f'sc210584_{ct}.csv', index=False)
        up, dn = deutil.ndeg(r)
        up5, dn5 = deutil.ndeg(r, lfc=0.585)
        rows.append(dict(celltype=ct, n_mid=(meta.cond == 'test').sum(), n_young=(meta.cond == 'ref').sum(),
                         cells_min=ncell[ok].min(), cells_median=ncell[ok].median(), genes_tested=ng,
                         deg_up=up, deg_down=dn, deg_up_lfc0585=up5, deg_down_lfc0585=dn5,
                         padj05_any=int((r.padj < 0.05).sum()), excluded=','.join(ncell.index[~ok])))
        print(rows[-1])
        # 擬似バルクの PCA（log CPM、検定遺伝子）
        keep = deutil.expr_filter(x, meta.cond.value_counts().min())
        l = np.log2(x.loc[:, keep].div(x.sum(axis=1), axis=0) * 1e6 + 1)
        l = l - l.mean()
        u, sv, vt = np.linalg.svd(l.values, full_matrices=False)
        pcs[ct] = (pd.DataFrame(u[:, :2] * sv[:2], index=l.index, columns=['PC1', 'PC2']), sv ** 2 / (sv ** 2).sum())
    pd.DataFrame(rows).to_csv(DE / 'de_summary_sc210584.csv', index=False)
    sensitivity()

    fig, axes = plt.subplots(1, len(pcs), figsize=(common.mm(170), common.mm(55)))
    for ax, (ct, (p, ve)) in zip(np.atleast_1d(axes), pcs.items()):
        g = n.loc[p.index, 'group']
        for grp, col in [('young', '#4A7BA7'), ('middle', '#C44E52')]:
            m = g == grp
            ax.scatter(p.PC1[m], p.PC2[m], c=col, s=14, label=grp)
        for k, (a_, b_) in p.iterrows():
            ax.annotate(k, (a_, b_), fontsize=5, xytext=(2, 2), textcoords='offset points')
        ax.set_xlabel(f'PC1 ({ve[0] * 100:.0f}%)')
        ax.set_ylabel(f'PC2 ({ve[1] * 100:.0f}%)')
        ax.set_title(f'{ct} pseudobulk')
    np.atleast_1d(axes)[0].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(SC / 'GSE210584_pb_pca.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def sensitivity():
    """中年群は擬似バルク PCA で 2 群に割れる（MA2/4/6/8 と MA1/3/5/7/9。CAR・内皮で同じ割れ方）。
    処理日などのバッチの可能性があるので、中年をどちらか片方だけにして若齢と比べ、主解析と向きがそろうかを見る。
    あわせて中年 2 群どうしの DE（何が割れを作っているか）も出す。"""
    A = ['MA2', 'MA4', 'MA6', 'MA8']
    B = ['MA1', 'MA3', 'MA5', 'MA7', 'MA9']
    Y = ['Y1', 'Y2', 'Y3', 'Y4', 'Y5']
    rows = []
    for ct in ['CAR', 'Osteo', 'Endothelial']:
        x = pd.read_csv(SC / f'GSE210584_pb_{ct}.csv.gz', index_col=0).T
        main = pd.read_csv(DE / f'sc210584_{ct}.csv').set_index('gene')
        for name, test, ref in [('MA_2468_vs_Y', A, Y), ('MA_13579_vs_Y', B, Y), ('MA_2468_vs_13579', A, B)]:
            xx = x.loc[test + ref]
            meta = pd.DataFrame({'cond': ['test'] * len(test) + ['ref'] * len(ref)}, index=xx.index)
            r, ng = deutil.deseq(xx, meta, '~cond', ['cond', 'test', 'ref'], n_min=min(len(test), len(ref)))
            r.to_csv(DE / f'sens_sc210584_{ct}_{name}.csv', index=False)
            up, dn = deutil.ndeg(r)
            j = main.join(r.set_index('gene')[['log2FoldChange', 'stat']], rsuffix='_s', how='inner').dropna(subset=['stat', 'stat_s'])
            md = main[(main.padj < 0.05) & (main.log2FoldChange.abs() >= 1)].index
            jj = j.loc[j.index.intersection(md)]
            rows.append(dict(celltype=ct, contrast=name, deg_up=up, deg_down=dn,
                             spearman_stat_vs_main=j[['stat', 'stat_s']].corr(method='spearman').iloc[0, 1],
                             main_deg_same_sign=float((np.sign(jj.log2FoldChange) == np.sign(jj.log2FoldChange_s)).mean()),
                             n_main_deg=len(jj)))
            if name == 'MA_2468_vs_13579':
                s_ = r.dropna(subset=['padj']).sort_values('stat')
                print(ct, 'batch-split top:', list(s_.gene[:15]), list(s_.gene[::-1][:15]))
    d = pd.DataFrame(rows)
    d.to_csv(DE / 'sens_sc210584_split.csv', index=False)
    print(d.round(3).to_string())


if __name__ == '__main__':
    main()

"""in vitro ヒト contrast の prerank GSEA と、マウス生体 Dox（prelim/chemo2）との並置。

  uv run python -I prelim/invitro/gsea_iv.py      # repo root から（de_iv.py の後）

gene set は prelim/chemo2/genesets.py と**同じ term**（Hallmark 50 + CURATED 17 + 手製 niche パネル 5）を、
MSigDB human 2025.1.Hs の GMT から引く（マウス側と同じ集合名で比べるため、足し引きしない）。
GMT は https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/ から直接取得して
out/prelim/invitro/msigdb/ に置いた（h.all / c2.all / c5.go.bp）。
手製パネルはマウスシンボルの大文字化（Il6→IL6 など。パネル内に名前が変わる遺伝子は無い）。
順位: Wald 統計量（log2FC / lfcSE）。gseapy prerank, 1000 perm, seed 0, min_size 4, max_size 1000（chemo2 と同じ）。
TE-Bone（TempO-Seq S1500+）は検定できた遺伝子が ~1,300 しかなく、set との重なりが小さい → 参考値。

出力（out/prelim/invitro/）: gsea/<label>.csv, nes_long.csv, nes_wide.csv, genesets_meta.csv,
  heatmap_vs_mouse.png|pdf|tiff（niche set × [in vitro ヒト | マウス生体 Dox]）, heatmap_hallmark_iv.png|pdf,
  concordance_vs_mouse.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
sys.path.insert(0, str(ROOT / 'prelim/invitro'))
import common  # noqa: E402
import genesets as GS  # noqa: E402
from ivcontrasts import IC  # noqa: E402

OUT = ROOT / 'out/prelim/invitro'
GMT = OUT / 'msigdb'
VER = 'v2025.1.Hs'
MOUSE = ROOT / 'out/prelim/chemo2/nes_wide.csv'
# マウス生体 Dox の列（chemo2 の CONTRASTS の名前）
MOUSE_COLS = [
    ('B_Dox', 'Ms Dox 1x 30d, bone ♂'),
    ('B_Dox_compadj', 'Ms Dox 1x 30d, +comp. adj.'),
    ('B3_Dox4x_MSPC', 'Ms Dox 4x 8wk, LepR+ MSPC ♀'),
    ('B3_Dox1x_MSPC', 'Ms Dox 1x 8wk, LepR+ MSPC ♀'),
    ('B2_Dox4x_OB', 'Ms Dox 4x 8wk, Col2.3+ OB ♀'),
    ('sc_CAR', 'Ms Dox 96h scRNA, CAR (n=1)'),
    ('sc_Osteoprogenitor', 'Ms Dox 96h scRNA, OP (n=1)'),
    ('sc_Osteoblast', 'Ms Dox 96h scRNA, OB (n=1)'),
]
MOUSE_BULK = ['B_Dox', 'B3_Dox4x_MSPC', 'B3_Dox1x_MSPC', 'B2_Dox4x_OB']
IV_NAME = {k: 'Hs ' + v['short'] for k, v in IC.items()}
CATORDER = ['inflammation', 'myeloid', 'osteoclast', 'senescence/SASP', 'HSC niche', 'vascular', 'ECM',
            'bone formation', 'adipo', 'proliferation', 'estrogen', 'composition?', 'hand panel']


def read_gmt(name):
    d = {}
    with open(GMT / f'{name}.{VER}.symbols.gmt') as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            d[x[0]] = [g.upper() for g in x[2:] if g]
    return d


def collection():
    """chemo2/genesets.collection() と同じ term・同じラベル/カテゴリで、ヒト GMT から作る。"""
    h = read_gmt('h.all')
    allh = {**read_gmt('c2.all'), **read_gmt('c5.go.bp')}
    hn = {a: (b, c) for a, b, c in GS.HALLMARK_NICHE}
    sets, meta = {}, []
    for k, v in h.items():
        sets[k] = v
        lab, cat = hn.get(k, ('H ' + k.replace('HALLMARK_', '').lower().replace('_', ' '), 'hallmark (other)'))
        meta.append(dict(term=k, label=lab, category=cat, source='Hallmark', niche=k in hn))
    for k, lab, cat in GS.CURATED:
        sets[k] = allh[k]
        meta.append(dict(term=k, label=lab, category=cat, source='MSigDB C2/C5', niche=True))
    for k, v in GS.NICHE_GSEA.items():
        sets[k] = [g.upper() for g in v]
        meta.append(dict(term=k, label=GS.NICHE_LABEL[k], category='hand panel', source='hand', niche=True))
    meta = pd.DataFrame(meta)
    meta['size_total'] = meta.term.map(lambda t: len(sets[t]))
    return sets, meta


def ranking(lb):
    d = pd.read_csv(OUT / 'de' / f'{lb}.csv').dropna(subset=['lfcSE', 'gene'])
    d = d[d.lfcSE > 0]
    return (d.log2FoldChange / d.lfcSE).groupby(d.gene.str.upper()).mean().sort_values(ascending=False)


def run_gsea(sets):
    (OUT / 'gsea').mkdir(parents=True, exist_ok=True)
    rows = []
    for lb in IC:
        rnk = ranking(lb)
        r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=4, max_size=1000,
                       seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
        r = r.assign(label=lb, n_ranked=len(rnk))
        r.to_csv(OUT / 'gsea' / f'{lb}.csv', index=False)
        rows.append(r)
        print(lb, len(rnk), 'genes;', len(r), 'sets tested;', (r['FDR q-val'].astype(float) < 0.05).sum(), 'FDR<0.05')
    L = pd.concat(rows, ignore_index=True)
    L = L.rename(columns={'Term': 'term', 'FDR q-val': 'fdr', 'NOM p-val': 'pval', 'Lead_genes': 'lead_genes',
                          'Tag %': 'tag_pct'})
    for k in ['NES', 'ES', 'fdr', 'pval']:
        L[k] = L[k].astype(float)
    L['matched'] = L.tag_pct.astype(str).str.split('/').str[1].astype(int)
    return L[['label', 'term', 'NES', 'ES', 'pval', 'fdr', 'matched', 'lead_genes', 'n_ranked']]


def heatmap(W, F, labels, cats, groups, names, fname, title, height_mm):
    common.setup()
    cols = list(W.columns)
    fig, ax = plt.subplots(figsize=(common.mm(170), common.mm(height_mm)))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#EEEEEE')
    im = ax.imshow(np.ma.masked_invalid(W.values.astype(float)), cmap=cmap, vmin=-3, vmax=3, aspect='auto')
    for i in range(W.shape[0]):
        for j in range(W.shape[1]):
            v, f = W.iloc[i, j], F.iloc[i, j]
            if np.isnan(v):
                continue
            if f < 0.05:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=5.5, fontweight='bold',
                        color='white' if abs(v) > 2 else 'black')
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, lw=0.7, ec='black'))
            elif f < 0.25:
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=5.5, color='#333333')
    ax.set_yticks(range(len(labels)), labels, fontsize=6)
    ax.set_xticks(range(len(cols)), [names[c] for c in cols], rotation=60, ha='right', fontsize=6)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    prev = None
    for i, c in enumerate(cats):
        if prev is not None and c != prev:
            ax.axhline(i - 0.5, c='k', lw=0.5)
        prev = c
    for j in range(1, len(cols)):
        if groups[j] != groups[j - 1]:
            ax.axvline(j - 0.5, c='k', lw=1.2)
    cb = fig.colorbar(im, ax=ax, shrink=0.3, pad=0.01)
    cb.set_label('NES', fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=7, loc='left')
    fig.tight_layout()
    for ext in ('png', 'pdf', 'tiff'):
        fig.savefig(OUT / f'{fname}.{ext}', dpi=common.DPI, bbox_inches='tight',
                    **({'pil_kwargs': {'compression': 'tiff_lzw'}} if ext == 'tiff' else {}))
    plt.close(fig)
    print('saved', OUT / f'{fname}.png')


def sgn(w, f):
    return 0 if (np.isnan(w) or f >= 0.05) else int(np.sign(w))


def concordance(W, F, meta):
    """set ごとに「ヒト in vitro（GSE284666 の 4 contrast + TE-Bone）」と「マウス生体 Dox bulk 4 本」の向きを並べる。
    ヒト側の判定: 4 contrast（NHO/MSC × Eto/IR）のうち FDR<0.05 で同符号の数。TE-Bone は別列。
    マウス側: chemo2 の分類と同じ（bulk 4 本のどれかで FDR<0.05、向きは多数決、割れたら mixed）。記述用（検定ではない）。"""
    hs4 = ['NHO_Eto', 'NHO_IR', 'MSC_Eto', 'MSC_IR']
    arrow = {1: '↑', -1: '↓', 0: '·'}
    rows = []
    for t in W.index:
        hs = [sgn(W.loc[t, c], F.loc[t, c]) for c in hs4]
        up, dn = hs.count(1), hs.count(-1)
        h = 1 if up >= 3 and dn == 0 else (-1 if dn >= 3 and up == 0 else (0 if up == dn == 0 else 'mixed'))
        ms = [sgn(W.loc[t, c], F.loc[t, c]) for c in MOUSE_BULK]
        sig = [v for v in ms if v]
        m = 0 if not sig else (sig[0] if all(v == sig[0] for v in sig) else 'mixed')
        if h in (1, -1) and m in (1, -1):
            call = 'same' if h == m else 'OPPOSITE'
        elif h in (1, -1):
            call = 'human in vitro only'
        elif m in (1, -1):
            call = 'mouse in vivo only'
        else:
            call = 'neither/mixed'
        rows.append(dict(term=t, human_consensus=arrow.get(h, '↕'), human_up=up, human_down=dn,
                         **{f'Hs:{c}': arrow[s] for c, s in zip(hs4, hs)},
                         **{'Hs:TEB_Dox': arrow[sgn(W.loc[t, 'TEB_Dox'], F.loc[t, 'TEB_Dox'])]},
                         mouse_bulk=arrow.get(m, '↕'), **{f'Ms:{c}': arrow[s] for c, s in zip(MOUSE_BULK, ms)},
                         call=call))
    d = pd.DataFrame(rows).merge(meta[['term', 'label', 'category', 'niche']], on='term')
    d.to_csv(OUT / 'concordance_vs_mouse.csv', index=False)
    print(d[d.niche].call.value_counts().to_string())
    return d


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sets, meta = collection()
    meta.to_csv(OUT / 'genesets_meta.csv', index=False)
    L = run_gsea(sets)
    L = L.merge(meta[['term', 'label', 'category', 'source', 'niche']].rename(columns={'label': 'set_label'}),
                on='term', how='left')
    L.to_csv(OUT / 'nes_long.csv', index=False)
    iv = list(IC)
    W = L.pivot(index='term', columns='label', values='NES').reindex(columns=iv)
    F = L.pivot(index='term', columns='label', values='fdr').reindex(columns=iv)
    pd.concat({'NES': W, 'FDR': F}, axis=1).to_csv(OUT / 'nes_wide.csv')

    M = pd.read_csv(MOUSE, header=[0, 1], index_col=0)
    mc = [c for c, _ in MOUSE_COLS]
    MW, MF = M['NES'][mc], M['FDR'][mc]
    W2 = W.join(MW, how='outer')
    F2 = F.join(MF, how='outer').fillna(1.0)
    names = {**IV_NAME, **dict(MOUSE_COLS)}
    groups = ['hs-2d'] * 4 + ['hs-3d'] + ['ms-bulk'] * 5 + ['ms-sc'] * 3

    nm = meta[meta.niche].copy()
    nm['co'] = nm.category.map({c: i for i, c in enumerate(CATORDER)})
    nm = nm.sort_values(['co', 'source', 'term'])
    heatmap(W2.loc[nm.term], F2.loc[nm.term], [f'{l}  [{c}]' for l, c in zip(nm.label, nm.category)],
            list(nm.category), groups, names, 'heatmap_vs_mouse',
            'Prerank GSEA NES: human bone cells in vitro (left) vs mouse bone in vivo after Dox (right)\n'
            'bold+box: FDR<0.05; plain: FDR<0.25; grey: not tested. TE-Bone = TempO-Seq panel (~1,300 genes), '
            'scRNA n=1: indicative', 205)
    hm = meta[meta.source == 'Hallmark'].sort_values('term')
    heatmap(W2.loc[hm.term], F2.loc[hm.term], list(hm.label), ['h'] * len(hm), groups, names,
            'heatmap_hallmark_iv', 'Hallmark 50, prerank GSEA NES (bold+box: FDR<0.05; plain: FDR<0.25)', 205)
    concordance(W2, F2, meta)


if __name__ == '__main__':
    main()

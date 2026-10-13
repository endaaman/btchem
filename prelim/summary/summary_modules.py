"""要約モジュール図（下書き）: 8 モジュール × 7 contrast。有意（FDR<0.05）だけ色、0.05–0.25 は枠線、残りは空白。

  uv run python -I prelim/summary/summary_modules.py      # repo root から

値の出どころ（既存の csv は読むだけ）:
- マウス / in vitro: out/prelim/chemo2/nes_wide.csv, out/prelim/aging_f/nes_wide_all.csv, out/prelim/invitro/nes_wide.csv
- ヒト骨生検（H_PostMeno, H_Old）: 既存の out/prelim/human/enrich/gsea_*.csv は Hallmark 50 + SenMayo しか無く、
  GO ossification / GO +reg. osteoclast diff. / NABA core matrisome が無い。そこで out/prelim/human/de/<label>.csv の
  Wald 統計量を、in vitro と**同じ gene set 集合・同じ設定**（prelim/invitro/gsea_iv.py の collection() と prerank 設定:
  1000 perm, seed 0, min 4, max 1000）で回し直し、out/prelim/summary/gsea_human/ に置く。全行をこの再計算版から取る
  （set 集合が違うと FDR が変わるので、行ごとに出どころを混ぜない）。既存版との差は docs/prelim_summary_fig.md に書く。
  図から外した候補（TNF-α/NF-κB, estrogen early, apoptosis, IFN-γ）も csv には shown=False で残す。

出力: out/prelim/summary/summary_modules.{png,pdf}, summary_modules_values.csv（全セルの NES・FDR・判定・出どころ）
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
sys.path.insert(0, str(ROOT / 'prelim/invitro'))
import common  # noqa: E402
import gsea_iv  # noqa: E402  （collection() を借りるだけ。gsea_iv の出力には触らない）

OUT = ROOT / 'out/prelim/summary'
HUMAN_DE = ROOT / 'out/prelim/human/de'
P = ROOT / 'out/prelim'

# 行: (term, 図のラベル, 行グループ)。候補 10 から TNF-α/NF-κB（Inflammation とほぼ同じ形）と
# estrogen early（雄の列で意味が無く、閉経後で ↑ と直感に反する）を外して 8 行にした。値は csv 側には残さない
ROWS = [
    ('HALLMARK_INFLAMMATORY_RESPONSE', 'Inflammation', 'inflam'),
    ('HALLMARK_INTERFERON_ALPHA_RESPONSE', 'IFN-α response', 'inflam'),
    ('SAUL_SEN_MAYO', 'Senescence (SenMayo)', 'damage'),
    ('HALLMARK_P53_PATHWAY', 'p53 pathway', 'damage'),
    ('GOBP_POSITIVE_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION', 'Osteoclast support', 'niche'),
    ('GOBP_OSSIFICATION', 'Bone formation', 'niche'),
    ('NABA_CORE_MATRISOME', 'Extracellular matrix', 'niche'),
    ('HALLMARK_ANGIOGENESIS', 'Angiogenesis', 'niche'),
]

# 列: (key, 出どころ, その csv の列名, 図のラベル, 群)
COLS = [
    ('H_PostMeno', 'human', 'H_PostMeno', 'Human ♀\nbone biopsy\npostmeno.\nvs premeno.', 'aging'),
    ('H_Old', 'human', 'H_Old', 'Human ♀\nbone biopsy\nold vs young', 'aging'),
    ('Age_2.5y', 'chemo2', 'Age_2.5y', 'Mouse ♂\ncortical bone\n30 vs 2 mo', 'aging'),
    ('sc210584_CAR', 'aging_f', 'sc210584_CAR', 'Mouse ♀\nCAR cells\nmid-age\nvs young', 'aging'),
    ('B3_Dox4x_MSPC', 'chemo2', 'B3_Dox4x_MSPC', 'Mouse ♀\nLepR+ MSPC\nDox ×4\n(8 wk)', 'chemo'),
    ('B_Dox_compadj', 'chemo2', 'B_Dox_compadj', 'Mouse ♂\ncortical bone\nDox ×1\n(30 d)*', 'chemo'),
    ('NHO_Eto', 'invitro', 'NHO_Eto', 'Human\nosteoblasts\netoposide\n(in vitro)', 'chemo'),
]
GSE = {'H_PostMeno': 'GSE141595', 'H_Old': 'GSE72815', 'Age_2.5y': 'NIHMS1776347', 'sc210584_CAR': 'GSE210584',
       'B3_Dox4x_MSPC': 'GSE296758', 'B_Dox_compadj': 'GSE205237', 'NHO_Eto': 'GSE284666'}
WIDE = {'chemo2': P / 'chemo2/nes_wide.csv', 'aging_f': P / 'aging_f/nes_wide_all.csv',
        'invitro': P / 'invitro/nes_wide.csv'}
SIG, TREND = 0.05, 0.25


def human_gsea():
    """ヒト骨生検 2 contrast を in vitro と同じ set 集合で prerank。結果はキャッシュ（out/prelim/summary/gsea_human/）。"""
    d = OUT / 'gsea_human'
    d.mkdir(parents=True, exist_ok=True)
    sets = None
    res = {}
    for lb in ['H_PostMeno', 'H_Old']:
        f = d / f'{lb}.csv'
        if not f.exists():
            if sets is None:
                sets, _ = gsea_iv.collection()
            x = pd.read_csv(HUMAN_DE / f'{lb}.csv').dropna(subset=['lfcSE', 'gene'])
            x = x[x.lfcSE > 0]
            rnk = (x.log2FoldChange / x.lfcSE).groupby(x.gene.str.upper()).mean().sort_values(ascending=False)
            r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=4, max_size=1000,
                           seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
            r.assign(label=lb, n_ranked=len(rnk)).to_csv(f, index=False)
            print('ran GSEA', lb, len(rnk), 'genes')
        r = pd.read_csv(f)
        res[lb] = r.set_index('Term')[['NES', 'FDR q-val']].astype(float).rename(columns={'FDR q-val': 'FDR'})
    return res


DROPPED = [
    ('HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'TNF-α / NF-κB', 'inflam'),
    ('HALLMARK_ESTROGEN_RESPONSE_EARLY', 'Estrogen response', 'niche'),
    ('HALLMARK_APOPTOSIS', 'Apoptosis', 'damage'),
    ('HALLMARK_INTERFERON_GAMMA_RESPONSE', 'IFN-γ response', 'inflam'),
]


def table():
    hum = human_gsea()
    wide = {k: pd.read_csv(v, header=[0, 1], index_col=0) for k, v in WIDE.items()}
    rows = []
    for term, rl, rg in ROWS + DROPPED:
        for key, src, col, cl, cg in COLS:
            if src == 'human':
                h = hum[col]
                nes, fdr = (h.loc[term, 'NES'], h.loc[term, 'FDR']) if term in h.index else (np.nan, np.nan)
                srcf = f'out/prelim/summary/gsea_human/{col}.csv'
            else:
                w = wide[src]
                ok = term in w.index
                nes = float(w.loc[term, ('NES', col)]) if ok else np.nan
                fdr = float(w.loc[term, ('FDR', col)]) if ok else np.nan
                srcf = str(WIDE[src].relative_to(ROOT))
            call = ('n/a' if np.isnan(nes) else 'up' if (fdr < SIG and nes > 0) else 'down' if fdr < SIG
                    else 'trend_up' if (fdr < TREND and nes > 0) else 'trend_down' if fdr < TREND else 'ns')
            rows.append(dict(shown=(term, rl, rg) in ROWS, module=rl, term=term, row_group=rg, contrast=key, dataset=GSE[key],
                             column_label=cl.replace('\n', ' '), group=cg, NES=nes, FDR=fdr, call=call, source=srcf))
    return pd.DataFrame(rows)


def plot(T):
    common.setup()
    nr, nc = len(ROWS), len(COLS)
    W, H = common.mm(170), common.mm(80)
    fig = plt.figure(figsize=(W, H))
    # 軸の位置（figure 座標）: 左に行ラベル、上に列ラベルと群の括弧、右に凡例
    L, R, B, Tp = 0.215, 0.79, 0.035, 0.64
    ax = fig.add_axes([L, B, R - L, Tp - B])
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(nr - 0.5, -0.5)
    ax.axis('off')
    cmap = plt.get_cmap('RdBu_r')
    norm = Normalize(-3, 3)
    pad = 0.06
    for _, x in T.iterrows():
        i = [r[1] for r in ROWS].index(x.module)
        j = [c[0] for c in COLS].index(x.contrast)
        box = (j - 0.5 + pad, i - 0.5 + pad)
        s = 1 - 2 * pad
        if x.call in ('up', 'down'):
            ax.add_patch(Rectangle(box, s, s, fc=cmap(norm(x.NES)), ec='none'))
            ax.text(j, i, f'{x.NES:+.1f}'.replace('-', '−'), ha='center', va='center', fontsize=6,
                    color='white' if abs(x.NES) >= 2.0 else '#222222')
        elif x.call.startswith('trend'):
            ec = cmap(norm(1.4 if x.NES > 0 else -1.4))
            ax.add_patch(Rectangle(box, s, s, fc='none', ec=ec, lw=0.45))
        elif x.call == 'n/a':
            ax.text(j, i, '–', ha='center', va='center', fontsize=7, color='#999999')
    # 行グループの区切り
    rg = [r[2] for r in ROWS]
    for i in range(1, nr):
        if rg[i] != rg[i - 1]:
            ax.plot([-0.5, nc - 0.5], [i - 0.5, i - 0.5], c='#888888', lw=0.4, clip_on=False)
    # 群（1st / 2nd insult）の境界
    cg = [c[4] for c in COLS]
    jb = next(j for j in range(1, nc) if cg[j] != cg[j - 1])
    ax.plot([jb - 0.5, jb - 0.5], [-0.5, nr - 0.5], c='black', lw=0.8, clip_on=False)
    # 行ラベル
    for i, (_, rl, _) in enumerate(ROWS):
        ax.text(-0.62, i, rl, ha='right', va='center', fontsize=8)
    # 列ラベル
    for j, c in enumerate(COLS):
        ax.text(j, -0.62, c[3], ha='center', va='bottom', fontsize=6.5, linespacing=1.05)
    # 群の括弧ラベル（列ラベルの上）
    ytop = -3.2
    for (j0, j1, lab) in [(0, jb - 1, 'Menopause / aging  (1st insult)'),
                          (jb, nc - 1, 'Chemotherapy  (2nd insult)')]:
        x0, x1 = j0 - 0.42, j1 + 0.42
        ax.plot([x0, x0, x1, x1], [ytop + 0.18, ytop, ytop, ytop + 0.18], c='black', lw=0.7, clip_on=False)
        ax.text((x0 + x1) / 2, ytop - 0.12, lab, ha='center', va='bottom', fontsize=8.5, fontweight='bold')
    # 凡例（右）
    lx = fig.add_axes([0.815, 0.36, 0.022, 0.26])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=lx)
    cb.set_ticks([-3, 0, 3])
    cb.set_ticklabels(['−3', '0', '+3'])
    cb.ax.tick_params(labelsize=6, length=2, width=0.5)
    cb.outline.set_linewidth(0.4)
    lx.set_title('NES\n(FDR < 0.05)', fontsize=6.5, loc='left', pad=3)
    kx = fig.add_axes([0.80, 0.03, 0.2, 0.28])
    kx.set_xlim(0, 1)
    kx.set_ylim(0, 1)
    kx.axis('off')
    items = [('fill', 'FDR < 0.05'), ('edge', 'FDR 0.05–0.25'), ('blank', 'FDR ≥ 0.25')]
    if (T.call == 'n/a').any():
        items.append(('na', 'not tested'))
    for k, (kind, lab) in enumerate(items):
        y = 0.86 - k * 0.24
        if kind == 'fill':
            kx.add_patch(Rectangle((0.02, y - 0.08), 0.11, 0.16, fc=cmap(norm(2.2)), ec='none'))
        elif kind == 'edge':
            kx.add_patch(Rectangle((0.02, y - 0.08), 0.11, 0.16, fc='none', ec=cmap(norm(1.4)), lw=0.45))
        elif kind == 'blank':
            kx.add_patch(Rectangle((0.02, y - 0.08), 0.11, 0.16, fc='none', ec='#DDDDDD', lw=0.4, ls=':'))
        else:
            kx.text(0.075, y, '–', ha='center', va='center', fontsize=7, color='#999999')
        kx.text(0.18, y, lab, va='center', fontsize=6.5)
    kx.text(0.0, 0.86 - len(items) * 0.24, '* composition-adjusted', va='center', fontsize=6, color='#444444')
    return fig


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    T = table()
    T.to_csv(OUT / 'summary_modules_values.csv', index=False)
    print(T[T.shown].pivot(index='module', columns='contrast', values='call')
          .reindex(index=[r[1] for r in ROWS], columns=[c[0] for c in COLS]).to_string())
    fig = plot(T[T.shown])
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'summary_modules.{ext}', dpi=common.DPI)
        print('saved', OUT / f'summary_modules.{ext}')
    plt.close(fig)


if __name__ == '__main__':
    main()

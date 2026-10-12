"""Figure S2: Hallmark prerank GSEA の NES ヒートマップ（マウス主・補助 + ヒト）。

入力（prelim の GSEA 結果をそのまま読む。再計算しない）:
    マウス  out/prelim/enrich/gsea_{Age_2.5y,A_OVX_cort,B_Dox,C_DEX}.csv
            out/prelim/chemo/gsea/B3_Dox4x_MSPC.csv
    ヒト    out/prelim/human/enrich/gsea_{H_Old,H_PostMeno,H_GIOP}.csv（SenMayo を 51 番目のセットとして含む）
    いずれも Wald 統計量で順位付けした prerank GSEA（gseapy, 1,000 perm）。FDR は各 contrast 内で計算。
行: Hallmark 50 + SenMayo のうち、どれか 1 つの contrast で FDR < 0.05 のもの。NES で階層クラスタリング（平均連結, ユークリッド）。
印: FDR < 0.05 に「*」。SenMayo はマウスで評価していないので空欄（灰）。

実行: uv run python fig/figS2.py -> out/fig/figS2.{png,pdf,tiff}, out/fig/figS2_nes.csv, figS2_fdr.csv
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

COLS = [  # (列見出し, csv, 群)
    ('Aging (mouse)', 'out/prelim/enrich/gsea_Age_2.5y.csv', 'A'),
    ('OVX (mouse)', 'out/prelim/enrich/gsea_A_OVX_cort.csv', 'A'),
    ('Dox (mouse bone)', 'out/prelim/enrich/gsea_B_Dox.csv', 'B'),
    ('Dox ×4 (mouse MSPC)', 'out/prelim/chemo/gsea/B3_Dox4x_MSPC.csv', 'B'),
    ('DEX (mouse)', 'out/prelim/enrich/gsea_C_DEX.csv', 'C'),
    ('Old (human)', 'out/prelim/human/enrich/gsea_H_Old.csv', 'H'),
    ('PostMeno (human)', 'out/prelim/human/enrich/gsea_H_PostMeno.csv', 'H'),
    ('GIOP (human)', 'out/prelim/human/enrich/gsea_H_GIOP.csv', 'H'),
]
FDR = 0.05
VMAX = 3.0
SENMAYO = 'SenMayo (SAUL_SEN_MAYO)'

# 発散カラーマップ: 青（低下）— 中立灰 — 赤（上昇）。両極は fig4 の Down / Up と同じ色。
CMAP = LinearSegmentedColormap.from_list('div', [common.COLOR_LRO, '#F0EFEC', common.COLOR_HRO])
CMAP.set_bad('#FFFFFF')


def load():
    nes, fdr = {}, {}
    for name, f, _ in COLS:
        d = pd.read_csv(f).set_index('Term')
        assert d.index.is_unique, f
        nes[name] = d.NES.astype(float)
        fdr[name] = d['FDR q-val'].astype(float)
    return pd.DataFrame(nes), pd.DataFrame(fdr)


def crosscheck(N):
    """prelim の横並び表（overlap/hallmark_nes.csv）と一致するか。"""
    ref = pd.read_csv('out/prelim/overlap/hallmark_nes.csv', index_col=0)
    m = {'Aging (mouse)': 'Age_2.5y', 'OVX (mouse)': 'A_OVX_cort', 'Dox (mouse bone)': 'B_Dox', 'DEX (mouse)': 'C_DEX'}
    for k, v in m.items():
        dif = (N[k] - ref[v].reindex(N.index)).abs().max()
        print(f'cross-check {k} vs overlap/hallmark_nes.csv[{v}]: max |diff| = {dif:.3g}')


def pretty(t):
    # Enrichr 版 MSigDB_Hallmark_2020 の表記ゆれを直す（表示のみ）
    fix = {SENMAYO: 'SenMayo (senescence)', 'Pperoxisome': 'Peroxisome', 'heme Metabolism': 'Heme Metabolism'}
    t = fix.get(t, t)
    return ' '.join(t.split())


def main():
    common.setup()
    N, F = load()
    crosscheck(N)
    keep = (F < FDR).any(axis=1)
    N, F = N[keep], F[keep]
    # クラスタリング（NaN = SenMayo のマウス列は 0 で代入して距離計算のみに使う）
    Z = linkage(N.fillna(0).values, 'average')
    order = N.index[leaves_list(Z)]
    if SENMAYO in order:  # SenMayo は Hallmark 外なので最下行に分ける
        order = [t for t in order if t != SENMAYO] + [SENMAYO]
    N, F = N.loc[order], F.loc[order]
    N.to_csv(os.path.join(common.OUT_DIR, 'figS2_nes.csv'))
    F.to_csv(os.path.join(common.OUT_DIR, 'figS2_fdr.csv'))
    print(f'rows kept: {len(N)} / {len(keep)}')

    nr, nc = N.shape
    cell_h, cell_w = 3.3, 9.0   # mm
    lab_w, top, bot = 50.0, 24.0, 4.0
    fig_w = common.FULL_WIDTH_MM
    fig_h = top + nr * cell_h + bot
    fig = plt.figure(figsize=(common.mm(fig_w), common.mm(fig_h)))
    ax = fig.add_axes([lab_w / fig_w, bot / fig_h, nc * cell_w / fig_w, nr * cell_h / fig_h])
    im = ax.imshow(np.ma.masked_invalid(N.values), cmap=CMAP, vmin=-VMAX, vmax=VMAX, aspect='auto',
                   interpolation='nearest')
    for i in range(nr):
        for j in range(nc):
            v = N.iat[i, j]
            if np.isnan(v):
                ax.text(j, i, 'n.a.', ha='center', va='center', fontsize=5.5, color='#8C8C8C')
            elif F.iat[i, j] < FDR:
                ax.text(j, i + 0.05, '*', ha='center', va='center', fontsize=8,
                        color='white' if abs(v) > 2.2 else 'black')
    # 白い区切り線（セル間 0.5 pt）
    ax.set_xticks(np.arange(-0.5, nc), minor=True)
    ax.set_yticks(np.arange(-0.5, nr), minor=True)
    ax.grid(which='minor', color='white', lw=0.6)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=0, pad=2)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_yticks(range(nr), [pretty(t) for t in N.index], fontsize=7)
    ax.set_xticks(range(nc), N.columns, rotation=45, ha='left', rotation_mode='anchor', fontsize=7)
    ax.xaxis.tick_top()
    # 群の区切り（マウス加齢/OVX | 化療 | DEX | ヒト）
    groups = [g for _, _, g in COLS]
    for j in range(1, nc):
        if groups[j] != groups[j - 1]:
            ax.axvline(j - 0.5, color='black', lw=0.8)
    if SENMAYO in N.index:
        ax.axhline(nr - 1.5, color='black', lw=0.8)
    # カラーバー
    cax = fig.add_axes([(lab_w + nc * cell_w + 6) / fig_w, (bot + nr * cell_h * 0.55) / fig_h,
                        2.5 / fig_w, nr * cell_h * 0.4 / fig_h])
    cb = fig.colorbar(im, cax=cax, extend='both')
    cb.set_label('NES', fontsize=8)
    cb.ax.tick_params(labelsize=7, width=0.6, length=2)
    cb.outline.set_linewidth(0.5)
    fig.text((lab_w + nc * cell_w + 4) / fig_w, (bot + nr * cell_h * 0.48) / fig_h,
             f'* FDR < {FDR}\nn.a. = not tested\n(SenMayo: human only)', fontsize=7, va='top')
    print(f'figure size: {fig_w} x {fig_h:.1f} mm')
    common.save(fig, 'figS2')
    print(N.round(2).to_string())


if __name__ == '__main__':
    main()

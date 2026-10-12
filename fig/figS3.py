"""Figure S3: 注目遺伝子パネル（prelim/genes.py、事前定義）の log2FC ヒートマップ。

列は Figure S2 と同じ 8 contrast。値は各 DE 表（prelim の出力 csv）の log2FoldChange をそのまま使う（再計算しない）。
ヒトはマウスのシンボルを大文字化して対応付け（Fasl → FASLG, Trp53 → TP53）。
印: padj < 0.05 に「*」。灰 = 検定対象外（低発現フィルタで除外、padj NA、または H_PostMeno の閉経前ファイル欠落）。
Il6 は両パネルに属するが、1 行（炎症・破骨側）にまとめる。

実行: uv run python fig/figS3.py -> out/fig/figS3.{png,pdf,tiff}, out/fig/figS3_lfc.csv, figS3_padj.csv
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import volcano as V  # noqa: E402  (パネル定義・シンボル対応)
from figS2 import CMAP  # noqa: E402

COLS = [  # (列見出し, csv, 種, 群)
    ('Aging (mouse)', 'out/prelim/de/Age_2.5y.csv', 'mouse', 'A'),
    ('OVX (mouse)', 'out/prelim/de/A_OVX_cort.csv', 'mouse', 'A'),
    ('Dox (mouse bone)', 'out/prelim/de/B_Dox.csv', 'mouse', 'B'),
    ('Dox ×4 (mouse MSPC)', 'out/prelim/chemo/de/B3_Dox4x_MSPC.csv', 'mouse', 'B'),
    ('DEX (mouse)', 'out/prelim/de/C_DEX.csv', 'mouse', 'C'),
    ('Old (human)', 'out/prelim/human/de/H_Old.csv', 'human', 'H'),
    ('PostMeno (human)', 'out/prelim/human/de/H_PostMeno.csv', 'human', 'H'),
    ('GIOP (human)', 'out/prelim/human/de/H_GIOP.csv', 'human', 'H'),
]
# prelim の横並び表との突き合わせ用（列名 → panel_table の contrast 名）
REF = {
    'Aging (mouse)': ('out/prelim/panel_table.csv', 'Age_2.5y'),
    'OVX (mouse)': ('out/prelim/panel_table.csv', 'A_OVX_cort'),
    'Dox (mouse bone)': ('out/prelim/panel_table.csv', 'B_Dox'),
    'DEX (mouse)': ('out/prelim/panel_table.csv', 'C_DEX'),
    'Old (human)': ('out/prelim/human/panel_table.csv', 'H_Old'),
    'PostMeno (human)': ('out/prelim/human/panel_table.csv', 'H_PostMeno'),
    'GIOP (human)': ('out/prelim/human/panel_table.csv', 'H_GIOP'),
}
PADJ = 0.05
VMAX = 2.0
GENES = {
    'Inflammation /\nosteoclastogenesis': list(V.G.INFLAM),
    'Apoptosis /\nsenescence': [g for g in V.G.APOP if g not in V.G.INFLAM],
}


def load():
    rows = sum(GENES.values(), [])
    L, P = pd.DataFrame(index=rows), pd.DataFrame(index=rows)
    for name, f, sp, _ in COLS:
        d = pd.read_csv(f).set_index('gene')
        sym = {g: (g if sp == 'mouse' else V.human_symbol(g)) for g in rows}
        L[name] = [d.log2FoldChange.get(sym[g], np.nan) for g in rows]
        P[name] = [d.padj.get(sym[g], np.nan) for g in rows]
    L[P.isna()] = np.nan  # padj が無いもの（独立フィルタ等）は検定対象外として灰
    return L, P


def crosscheck(L, P):
    for name, (f, c) in REF.items():
        t = pd.read_csv(f).set_index('gene')
        mouse = 'mouse' in name
        idx = [g if mouse else V.human_symbol(g) for g in L.index]
        ref_l = t[f'{c}|log2FC'].reindex(idx).values
        ref_p = t[f'{c}|padj'].reindex(idx).values
        ok = ~np.isnan(L[name].values)
        dl = np.nanmax(np.abs(L[name].values[ok] - ref_l[ok])) if ok.any() else 0
        dp = np.nanmax(np.abs(P[name].values[ok] - ref_p[ok])) if ok.any() else 0
        miss = int((ok & np.isnan(ref_l)).sum())
        print(f'cross-check {name} vs {f}[{c}]: max|dLFC|={dl:.3g}, max|dpadj|={dp:.3g}, missing in ref={miss}')


def main():
    common.setup()
    L, P = load()
    crosscheck(L, P)
    L.to_csv(os.path.join(common.OUT_DIR, 'figS3_lfc.csv'))
    P.to_csv(os.path.join(common.OUT_DIR, 'figS3_padj.csv'))
    nr, nc = L.shape
    cell_h, cell_w = 3.6, 9.0
    lab_w, grp_w, top, bot = 18.0, 24.0, 24.0, 4.0
    fig_w = common.FULL_WIDTH_MM
    fig_h = top + nr * cell_h + bot
    fig = plt.figure(figsize=(common.mm(fig_w), common.mm(fig_h)))
    x0 = grp_w + lab_w
    ax = fig.add_axes([x0 / fig_w, bot / fig_h, nc * cell_w / fig_w, nr * cell_h / fig_h])
    im = ax.imshow(np.ma.masked_invalid(L.values), cmap=CMAP, vmin=-VMAX, vmax=VMAX, aspect='auto',
                   interpolation='nearest')
    # 検定対象外は斜線（0 付近の淡色と区別するため、塗りではなくハッチ）
    plt.rcParams['hatch.linewidth'] = 0.4
    for i in range(nr):
        for j in range(nc):
            if np.isnan(L.iat[i, j]):
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor='white', edgecolor='#9A9A9A',
                                       hatch='//////', lw=0, zorder=1))
    for i in range(nr):
        for j in range(nc):
            if P.iat[i, j] < PADJ:
                v = L.iat[i, j]
                ax.text(j, i + 0.05, '*', ha='center', va='center', fontsize=8,
                        color='white' if abs(v) > 1.6 else 'black')
    ax.set_xticks(np.arange(-0.5, nc), minor=True)
    ax.set_yticks(np.arange(-0.5, nr), minor=True)
    ax.grid(which='minor', color='white', lw=0.6)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=0, pad=2)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_yticks(range(nr), L.index, fontsize=7, style='italic')
    ax.set_xticks(range(nc), L.columns, rotation=45, ha='left', rotation_mode='anchor', fontsize=7)
    ax.xaxis.tick_top()
    groups = [g for *_, g in COLS]
    for j in range(1, nc):
        if groups[j] != groups[j - 1]:
            ax.axvline(j - 0.5, color='black', lw=0.8)
    # パネル群の区切りと見出し
    y = -0.5
    for k, v in GENES.items():
        y_end = y + len(v)
        if y_end < nr - 0.5:
            ax.axhline(y_end, color='black', lw=0.8)
        yc = 1 - ((y + y_end) / 2 + 0.5) / nr
        fig.text(2.0 / fig_w, (bot + yc * nr * cell_h) / fig_h, k, fontsize=7.5, va='center', ha='left',
                 fontweight='bold')
        y = y_end
    cax = fig.add_axes([(x0 + nc * cell_w + 6) / fig_w, (bot + nr * cell_h * 0.55) / fig_h,
                        2.5 / fig_w, nr * cell_h * 0.4 / fig_h])
    cb = fig.colorbar(im, cax=cax, extend='both')
    cb.set_label('log$_2$ fold change', fontsize=8)
    cb.ax.tick_params(labelsize=7, width=0.6, length=2)
    cb.outline.set_linewidth(0.5)
    fig.text((x0 + nc * cell_w + 4) / fig_w, (bot + nr * cell_h * 0.48) / fig_h,
             f'* adj. $P$ < {PADJ}\nHatched = not tested\n(low expression /\nnot available)\n\nRows: mouse symbols;\nhuman orthologues\nby upper-case match\n(FASLG, TP53)',
             fontsize=7, va='top')
    print(f'figure size: {fig_w} x {fig_h:.1f} mm')
    common.save(fig, 'figS3')
    print(L.round(2).to_string())


if __name__ == '__main__':
    main()

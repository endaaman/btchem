"""Figure 4: マウス骨のドライ解析 volcano（A 加齢 / B ドキソルビシン / C デキサメタゾン）。

入力（prelim の DE 表をそのまま読む。再計算しない）:
    A out/prelim/de/Age_2.5y.csv   NIHMS1776347 tS1（Kaya 2022 JBMR, PRJNA695408）著者の DESeq2 結果。2.5 y (n=4) vs 2 mo (n=3), ♂ 上腕骨皮質骨
    B out/prelim/de/B_Dox.csv      GSE205237 Dox 10 mg/kg 単回 30 日後 vs vehicle (4 vs 4), ♂ 脛骨+大腿骨（骨髄除去）
    C out/prelim/de/C_DEX.csv      GSE276391 DEX 1 mg/kg 5 日 vs vehicle（非荷重, 5 vs 5）, ♂ 脛骨

実行: uv run python fig/fig4.py  -> out/fig/fig4.{png,pdf,tiff}, out/fig/fig4_labels.csv
"""
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import volcano as V  # noqa: E402

PANELS = [
    ('A', 'out/prelim/de/Age_2.5y.csv', 'Aged vs young cortical bone\n(2.5 y vs 2 mo)'),
    ('B', 'out/prelim/de/B_Dox.csv', 'Doxorubicin vs vehicle\n(bone, day 30)'),
    ('C', 'out/prelim/de/C_DEX.csv', 'Dexamethasone vs vehicle\n(tibia, day 5)'),
]

# レイアウト (mm)
FIG_W = common.FULL_WIDTH_MM
LEFT = 13.0      # y ラベル + パネルラベル
GAP = 12.5       # パネル間（次パネルの y ラベル分）
AX_W = (FIG_W - LEFT - 2 * GAP - 2.0) / 3
AX_H = 44.0
TOP = 11.0       # タイトル 2 行 + パネルラベル
XLAB = 10.0      # x 目盛 + x ラベル
LEG = 16.0       # 下の凡例（2 行 + 注 2 行）
FIG_H = TOP + AX_H + XLAB + LEG


def main():
    common.setup()
    fig = plt.figure(figsize=(common.mm(FIG_W), common.mm(FIG_H)))
    rows = []
    for i, (lab, path, title) in enumerate(PANELS):
        left = LEFT + i * (AX_W + GAP)
        bottom = LEG + XLAB
        ax = fig.add_axes([left / FIG_W, bottom / FIG_H, AX_W / FIG_W, AX_H / FIG_H])
        d = V.load(path)
        up, dn, sp = V.plot(ax, d, title)
        common.panel_label(ax, lab, dx=-11.5 / AX_W, dy=1.0 + 2.5 / AX_H)
        tested = int(d.padj.notna().sum())
        print(f'{lab} {path}: rows={len(d)}, plotted (padj not NaN)={tested}, up={up}, down={dn}')
        for _, r in sp.iterrows():
            rows.append(dict(panel=lab, gene=r.gene, category=r.category, log2FC=round(r.log2FoldChange, 3),
                             padj=f'{r.padj:.2g}', DEG=r.deg, consistent_with_literature=r.consistent))
    h = V.legend_handles()  # 行優先に並べる: 1 行目 Up/Down/NS, 2 行目 パネル記号
    fig.legend(handles=[h[0], h[3], h[1], h[4], h[2]], loc='lower center', ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 6.0 / FIG_H), fontsize=7, handletextpad=0.3, columnspacing=1.0)
    fig.text(0.5, 0.4 / FIG_H, V.LEGEND_NOTE, ha='center', va='bottom', fontsize=6.5, color='#333333')
    lab = pd.DataFrame(rows)
    print(lab.to_string(index=False))
    os.makedirs(common.OUT_DIR, exist_ok=True)
    lab.to_csv(os.path.join(common.OUT_DIR, 'fig4_labels.csv'), index=False)
    print(f'figure size: {FIG_W} x {FIG_H:.1f} mm')
    common.save(fig, 'fig4')


if __name__ == '__main__':
    main()

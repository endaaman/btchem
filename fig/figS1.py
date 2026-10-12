"""Figure S1: 代替・補助データの volcano（Figure 4 と同じ描き方）。

入力（prelim の DE 表をそのまま読む）:
    A out/prelim/de/A_OVX_cort.csv                 GSE222752 OVX vs Sham, ♀ 大腿骨皮質骨, 術後 4 週, 3 vs 3（各 2 匹プール）
    B out/prelim/chemo/de/B3_Dox4x_MSPC.csv        GSE296758 Dox 5 mg/kg ×4 vs PBS, ♀ LepR+ 骨髄 MSPC, 最終投与 8 週後, 4 vs 4
    C out/prelim/human/de/H_Old.csv                GSE72815 高齢 vs 若年女性, 腸骨生検, 19 vs 19
    D out/prelim/human/de/H_PostMeno.csv           GSE141595 閉経後 vs 閉経前, 腸骨生検, 15 vs 15（閉経前は別バッチ・遺伝子は共通分のみ）
    E out/prelim/human/de/H_GIOP.csv               GSE276529 GIOP vs 非 GIOP, 女性大腿骨皮質骨, 5 vs 4

実行: uv run python fig/figS1.py  -> out/fig/figS1.{png,pdf,tiff}, out/fig/figS1_labels.csv
"""
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
import volcano as V  # noqa: E402

PANELS = [
    ('A',  'out/prelim/de/A_OVX_cort.csv', 'mouse', 'Mouse OVX vs sham\n(femoral cortical bone, 4 wk)'),
    ('B', 'out/prelim/chemo/de/B3_Dox4x_MSPC.csv', 'mouse', 'Mouse doxorubicin ×4 vs PBS\n(LepR$^+$ marrow MSPCs, 8 wk)',
     dict(inner={-1: 0.15}, top=0.80)),
    ('C',  'out/prelim/human/de/H_Old.csv', 'human', 'Human old vs young women\n(iliac crest biopsy)'),
    ('D',  'out/prelim/human/de/H_PostMeno.csv', 'human', 'Human postmenopausal vs\npremenopausal (iliac biopsy)'),
    ('E',  'out/prelim/human/de/H_GIOP.csv', 'human', 'Human GIOP vs non-GIOP\n(femoral cortical bone)'),
]

NOTE = ('*Pre-defined panel genes with\nadj. $P$ < 0.05 (any fold change);\nfill = direction of change.\n'
        'Bold label = direction consistent\nwith a pro-inflammatory /\npro-osteoclastogenic (or\n'
        'pro-apoptotic / pro-senescent)\nshift; grey label = opposite\ndirection.')

FIG_W = common.FULL_WIDTH_MM
LEFT = 13.0
GAP = 12.5
AX_W = (FIG_W - LEFT - 2 * GAP - 3.0) / 3
AX_H = 44.0
TOP = 11.0
XLAB = 10.0
ROW = TOP + AX_H + XLAB
FIG_H = 2 * ROW


def main():
    common.setup()
    fig = plt.figure(figsize=(common.mm(FIG_W), common.mm(FIG_H)))
    rows = []
    for i, (lab, path, sp_, title, *kw) in enumerate(PANELS):
        r, c = divmod(i, 3)
        left = LEFT + c * (AX_W + GAP)
        bottom = FIG_H - (r * ROW + TOP + AX_H)
        ax = fig.add_axes([left / FIG_W, bottom / FIG_H, AX_W / FIG_W, AX_H / FIG_H])
        d = V.load(path)
        up, dn, sp = V.plot(ax, d, title, species=sp_, label_kw=kw[0] if kw else None)
        common.panel_label(ax, lab, dx=-11.5 / AX_W, dy=1.0 + 2.5 / AX_H)
        print(f'{lab} {path}: rows={len(d)}, plotted={int(d.padj.notna().sum())}, up={up}, down={dn}')
        for _, x in sp.iterrows():
            rows.append(dict(panel=lab, gene=x.gene, category=x.category, log2FC=round(x.log2FoldChange, 3),
                             padj=f'{x.padj:.2g}', DEG=x.deg, consistent_with_literature=x.consistent))
    # 6 枠目に凡例
    left = LEFT + 2 * (AX_W + GAP)
    lax = fig.add_axes([(left - 6) / FIG_W, XLAB / FIG_H, (AX_W + 8) / FIG_W, AX_H / FIG_H])
    lax.axis('off')
    lax.legend(handles=V.legend_handles(), loc='upper left', frameon=False, fontsize=7,
               handletextpad=0.4, labelspacing=0.6)
    lax.text(0.0, 0.40, NOTE,
             transform=lax.transAxes, fontsize=6.5, va='top', color='#333333')
    lab = pd.DataFrame(rows)
    print(lab.to_string(index=False))
    os.makedirs(common.OUT_DIR, exist_ok=True)
    lab.to_csv(os.path.join(common.OUT_DIR, 'figS1_labels.csv'), index=False)
    print(f'figure size: {FIG_W} x {FIG_H:.1f} mm')
    common.save(fig, 'figS1')


if __name__ == '__main__':
    main()

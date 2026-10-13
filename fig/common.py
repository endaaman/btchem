"""MedComm 共通スタイル（prelim/common.py のコピー。出力先は out/<OUT_VER>/fig）。

投稿規定:
- 図の幅 170 mm (1 段組なら 80 mm)、高さ 210 mm 以内
- フォント 8-10 pt (最低 6 pt)
- パネルラベル A/B/C は大文字・太字
- 線幅 0.25 pt 以上、300 dpi 以上
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

OUT_VER = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'OUT_VER')).read().strip()
OUT_DIR = f'out/{OUT_VER}/fig'  # 出力世代は repo 直下の OUT_VER で指定
FULL_WIDTH_MM = 170
HALF_WIDTH_MM = 80
MAX_HEIGHT_MM = 210
DPI = 600

# rev0 と同じ配色
COLOR_LRO = '#4A7BA7'
COLOR_HRO = '#C44E52'


def mm(v):
    """mm -> inch"""
    return v / 25.4


def setup():
    for name in ['arial.ttf', 'arialbd.ttf', 'ariali.ttf', 'arialbi.ttf']:
        fm.fontManager.addfont(os.path.join('data/fonts', name))
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,
        'axes.titlesize': 9,
        'axes.labelsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'legend.title_fontsize': 8,
        'axes.linewidth': 0.8,
        'lines.linewidth': 1.0,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.6,
        'ytick.minor.width': 0.6,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'axes.spines.right': False,
        'axes.spines.top': False,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'savefig.dpi': DPI,
    })


def panel_label(ax, label, dx=-0.15, dy=1.05):
    """パネルラベル (A, B, C ...) を軸左上に 12 pt 太字で置く（編集部指定: 12 pt, bold, uppercase, top-left）。"""
    ax.text(dx, dy, label, transform=ax.transAxes, fontsize=12,
            fontweight='bold', va='bottom', ha='left')


def save(fig, name, formats=('png', 'pdf', 'tiff')):
    os.makedirs(OUT_DIR, exist_ok=True)
    for ext in formats:
        path = os.path.join(OUT_DIR, f'{name}.{ext}')
        fig.savefig(path, dpi=DPI, bbox_inches='tight', pad_inches=0)
        print('saved', path)
    # 出力サイズの確認（bbox_inches='tight' で規定幅を超えていないか）
    from PIL import Image
    w, h = Image.open(os.path.join(OUT_DIR, f'{name}.png')).size
    w_mm, h_mm = w / DPI * 25.4, h / DPI * 25.4
    flag = '' if w_mm <= FULL_WIDTH_MM + 0.5 and h_mm <= MAX_HEIGHT_MM else '  <-- EXCEEDS LIMIT'
    print(f'output size: {w_mm:.1f} x {h_mm:.1f} mm{flag}')

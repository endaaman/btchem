"""Figure 2 / 3 共通: 臨床データの読み込みと KM パネル描画（周術期化療 (+) vs (−)）。

書式は 1 本目 `,/参考1/rev1/fig3.py` の最終形を踏襲する:
  170 mm 幅・1 行 3 パネル、censor のひげ、At risk / Censored / Events × 10 か月刻み 7 列の表、
  凡例タイトルに log-rank p / Cox 単変量 HR / 95% CI。

入力（臨床 xlsx、清水から未着）: CLINICAL_XLSX。期待する列は COLUMNS（左が xlsx 上の想定列名）。
列名・文字列値は strip してから照合する。到着後に列名が違えば COLUMNS の左側だけ直す。

--mock: 1 本目の csv（MOCK_CSV）から擬似データを作る（make_mock 参照）。数値は実データではない。
"""
import os

import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator
from matplotlib.transforms import blended_transform_factory
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

import common

# ---- 入力 ----
# 本番: 清水の拡張データセット（未着。ファイル名・シート名は到着後に合わせる）
CLINICAL_XLSX = 'data/clinical/clinical.xlsx'
CLINICAL_SHEET = 0
# xlsx 上の列名 -> 内部名。到着後はここ（左側）だけ直す。
COLUMNS = {
    '研究ID': 'id',
    'menopause': 'menopause',      # 0 = pre, 1 = post
    'chemo': 'chemo',              # 周術期化療 0 = なし, 1 = あり
    'BMFS_Event': 'BMFS_Event',    # 骨転移 1 = あり
    'BMFS': 'BMFS',                # 骨転移なし生存期間（月）
    'LMFS_Event': 'LMFS_Event',    # 肺
    'LMFS': 'LMFS',
    'HMFS_Event': 'HMFS_Event',    # 肝
    'HMFS': 'HMFS',
}
# mock: 1 本目のデータ（785 例、同じコホート）
MOCK_CSV = ',/参考1/data/osteoporosis_bone_meta.csv'
MOCK_SEED = 20261012
MOCK_KEEP = {'LMFS': 0.35, 'HMFS': 0.25}  # BMFS のイベントを残す割合（肺・肝の擬似イベント）
MOCK_BANNER = 'MOCK DATA (surrogate grouping / simulated events) — NOT FOR PUBLICATION'

# endpoint key -> (時間列, イベント列, y ラベル)
ENDPOINTS = {
    'bone': ('BMFS', 'BMFS_Event', 'Bone metastasis-free ratio'),
    'lung': ('LMFS', 'LMFS_Event', 'Lung metastasis-free ratio'),
    'liver': ('HMFS', 'HMFS_Event', 'Liver metastasis-free ratio'),
}

# ---- 群 ----
GROUPS = [0, 1]  # chemo の値。描画・表の順は (−) → (+)（1 本目の LRO → HRO に対応）
GROUP_LABEL = {0: 'Chemo (−)', 1: 'Chemo (+)'}
COLORS = {0: common.COLOR_LRO, 1: common.COLOR_HRO}

# ---- 表・軸（1 本目 rev1/fig3.py と同じ値）----
TABLE_TIMES = [0, 10, 20, 30, 40, 50, 60]
TABLE_ROWS = ['At risk', 'Censored', 'Events']
TABLE_FONTSIZE = 8
XMAX = 65.5
TABLE_LABEL_RIGHT = -6.5
TABLE_LABEL_LEFT = -19.5

# ---- レイアウト (mm) ----
FIG_W = common.FULL_WIDTH_MM
PANEL_LEFT = 15.5
PANEL_GAP = 4.0
AX_W = 37.0
AX_H = 42.0
TOP = 7.0
XLABEL_GAP = 13.0
LINE_MM = 4.2
GROUP_GAP_MM = 3.0
BOTTOM = 3.0
BANNER_MM = 5.0  # mock 時のみ上に足す


def _strip_frame(df):
    df.columns = [str(c).strip() for c in df.columns]
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].str.strip()
    return df


def _check(df, endpoints):
    assert df['id'].is_unique, 'duplicate 研究ID'
    need = ['menopause', 'chemo'] + [c for k in endpoints for c in ENDPOINTS[k][:2]]
    na = df[need].isna().sum()
    assert not na.any(), f'missing values:\n{na[na > 0]}'
    for c in ['menopause', 'chemo'] + [ENDPOINTS[k][1] for k in endpoints]:
        df[c] = df[c].astype(int)
        assert df[c].isin([0, 1]).all(), f'{c} not 0/1: {sorted(df[c].unique())}'
    return df


def load_clinical(path=CLINICAL_XLSX, sheet=CLINICAL_SHEET, endpoints=tuple(ENDPOINTS)):
    """本番データ。列の対応は COLUMNS の 1 か所のみ。"""
    df = _strip_frame(pd.read_excel(path, sheet_name=sheet))
    need = {src: dst for src, dst in COLUMNS.items()
            if dst in ('id', 'menopause', 'chemo')
            or any(dst in ENDPOINTS[k][:2] for k in endpoints)}
    missing = [c for c in need if c not in df.columns]
    assert not missing, f'columns not found: {missing}\navailable: {list(df.columns)}'
    df = df[list(need)].rename(columns=need)
    df = df[df['id'].notna()].copy()
    n_na = int(df['chemo'].isna().sum())
    if n_na:
        print(f'*** dropped {n_na} rows with missing chemo (unknown perioperative chemotherapy)')
        df = df[df['chemo'].notna()].copy()
    return _check(df, endpoints)


def make_mock(endpoints=tuple(ENDPOINTS)):
    """1 本目の csv から擬似データ。
    chemo <- 'Cehst xray score'（High=1）。骨 = BMFS そのまま。
    肺・肝 = BMFS のイベントを MOCK_KEEP の割合で残し（seed 固定）、外したイベントは同時点で打ち切り。"""
    df = _strip_frame(pd.read_csv(MOCK_CSV))
    df = df[df['研究ID'].notna()].copy()
    out = pd.DataFrame({
        'id': df['研究ID'].astype(int),
        'menopause': df['menopause'],
        'chemo': df['Cehst xray score'].map({'Low': 0, 'High': 1}),
        'BMFS': df['BMFS'],
        'BMFS_Event': df['BMFS_Event'],
    })
    rng = np.random.default_rng(MOCK_SEED)
    for t, keep in MOCK_KEEP.items():
        out[t] = out['BMFS']
        out[t + '_Event'] = ((out['BMFS_Event'] == 1) & (rng.random(len(out)) < keep)).astype(int)
    return _check(out.reset_index(drop=True), endpoints)


def load(mock, path=CLINICAL_XLSX, endpoints=tuple(ENDPOINTS)):
    if mock:
        print(f'*** MOCK data from {MOCK_CSV} (chemo <- Cehst xray score; lung/liver events simulated)')
        df = make_mock(endpoints)
    else:
        print(f'*** clinical data: {path}')
        df = load_clinical(path, endpoints=endpoints)
    print(f'rows={len(df)}  chemo(+)={int(df.chemo.sum())}  chemo(−)={int((df.chemo == 0).sum())}  '
          f'pre={int((df.menopause == 0).sum())}  post={int((df.menopause == 1).sum())}')
    return df


def figure_height(mock):
    table = 2 * (1 + len(TABLE_ROWS)) * LINE_MM + GROUP_GAP_MM
    return TOP + AX_H + XLABEL_GAP + table + BOTTOM + (BANNER_MM if mock else 0)


def compute_stats(data, T, E):
    neg, pos = data[data.chemo == 0], data[data.chemo == 1]
    lr = logrank_test(neg[T], pos[T], neg[E], pos[E])
    cph = CoxPHFitter().fit(data[[T, E, 'chemo']], T, E)
    s = cph.summary.loc['chemo']
    return dict(
        n_neg=len(neg), ev_neg=int(neg[E].sum()),
        n_pos=len(pos), ev_pos=int(pos[E].sum()),
        logrank_p=lr.p_value,
        hr=s['exp(coef)'], lo=s['exp(coef) lower 95%'], hi=s['exp(coef) upper 95%'],
        wald_p=s['p'],
    )


def at_risk_table(data, T, E, group):
    d = data[data.chemo == group]
    t_, e_ = d[T].values, d[E].values
    rows = {'At risk': [], 'Censored': [], 'Events': []}
    for t in TABLE_TIMES:
        rows['At risk'].append(int((t_ >= t).sum()))
        rows['Censored'].append(int(((t_ < t) & (e_ == 0)).sum()))
        rows['Events'].append(int(((t_ < t) & (e_ == 1)).sum()))
    return {k: rows[k] for k in TABLE_ROWS}


def fmt_p(p):
    return 'p < 0.01' if p < 0.01 else f'p = {p:.2f}'


def plot_km_panel(ax, data, endpoint, title):
    T, E, ylabel = ENDPOINTS[endpoint]
    st = compute_stats(data, T, E)
    print(f'\n=== {title} ({endpoint}) ===')
    print(f"Chemo(−): n={st['n_neg']}, events={st['ev_neg']} | Chemo(+): n={st['n_pos']}, events={st['ev_pos']}")
    print(f"log-rank p={st['logrank_p']:.4f}  Cox HR (+ vs −)={st['hr']:.3f} "
          f"({st['lo']:.3f}-{st['hi']:.3f}), Wald p={st['wald_p']:.4f}")

    for g in GROUPS:
        mask = data.chemo == g
        kmf = KaplanMeierFitter()
        kmf.fit(data[mask][T], data[mask][E], label=f'{GROUP_LABEL[g]} (n={mask.sum()})')
        kmf.plot(ax=ax, ci_show=False, show_censors=True,
                 censor_styles={'ms': 5, 'marker': '|', 'mew': 0.8},
                 color=COLORS[g], linewidth=1.0)

    legend_title = '\n'.join([fmt_p(st['logrank_p']), f"HR {st['hr']:.2f}",
                              f"95% CI {st['lo']:.2f}–{st['hi']:.2f}"])
    ax.set_title(title)
    ax.set_xlabel('Months after surgery')
    ax.set_ylabel(ylabel)
    ax.set_xlim([0, XMAX])
    ax.set_ylim([0, 1])
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.tick_params(axis='both', which='major', direction='out', length=3.5)
    ax.tick_params(axis='both', which='minor', direction='out', length=2)
    leg = ax.legend(loc='lower left', title=legend_title, frameon=True, edgecolor='black',
                    fancybox=False, framealpha=1.0, borderpad=0.6, handlelength=1.8,
                    alignment='left')
    leg.get_title().set_multialignment('left')

    trans = blended_transform_factory(ax.transData, ax.transAxes)
    kw = dict(transform=trans, fontsize=TABLE_FONTSIZE, va='top')
    y_mm = -XLABEL_GAP
    print(f'Number at risk table (t = {TABLE_TIMES}):')
    for g in GROUPS:
        n = int((data.chemo == g).sum())
        ax.text(TABLE_LABEL_LEFT, y_mm / AX_H, f'{GROUP_LABEL[g]} (n={n})', ha='left', **kw)
        y_mm -= LINE_MM
        for label, vals in at_risk_table(data, T, E, g).items():
            ax.text(TABLE_LABEL_RIGHT, y_mm / AX_H, label, ha='right', **kw)
            for t, v in zip(TABLE_TIMES, vals):
                ax.text(t, y_mm / AX_H, str(v), ha='center', **kw)
            print(f'  {GROUP_LABEL[g]} {label:9s}: {vals}')
            y_mm -= LINE_MM
        y_mm -= GROUP_GAP_MM
    return st


def draw(fig, panels, mock):
    """panels: [(label, title, endpoint, data), ...] を 1 行に並べる。統計の行リストを返す。"""
    fig_h = figure_height(mock)
    rows = []
    for i, (label, title, endpoint, data) in enumerate(panels):
        left = i * (PANEL_LEFT + AX_W + PANEL_GAP) + PANEL_LEFT
        bottom = fig_h - TOP - AX_H - (BANNER_MM if mock else 0)
        ax = fig.add_axes([left / FIG_W, bottom / fig_h, AX_W / FIG_W, AX_H / fig_h])
        st = plot_km_panel(ax, data, endpoint, title)
        common.panel_label(ax, label, dx=-14.5 / AX_W, dy=1.0)
        rows.append(dict(panel=label, title=title, endpoint=endpoint, mock=mock, **st))
    if mock:
        fig.text(0.5, 1 - 1.0 / fig_h, MOCK_BANNER, ha='center', va='top',
                 fontsize=9, fontweight='bold', color='#C00000')
    return rows


def write_stats(rows, name):
    os.makedirs(common.OUT_DIR, exist_ok=True)
    path = os.path.join(common.OUT_DIR, f'{name}_stats.csv')
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False, float_format='%.4f')
    print('saved', path)
    print(df.to_string(index=False, float_format=lambda v: f'{v:.3f}'))

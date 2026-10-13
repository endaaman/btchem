"""PXD078161（Izquierdo et al. Geroscience 2026, PMID 42686970）脛骨皮質骨プロテオームで、
年齢・性別ごとの DX（ドキソルビシン）vs 対照の変化を並べる。

  uv run python -I prelim/proteome/proteome.py      # repo root から

入力: data/seq/PXD078161/VG{1-4}_*_Candidates.xlsx の 'All Proteins' シート（Spectronaut の群間比較表）。
  VG1 若齢雌 / VG2 老齢雌 / VG3 若齢雄 / VG4 老齢雄。4 群は別バッチ・別 Spectronaut 解析。
  VG2・VG4 は比較の向きが 'Ctrl / DX' なので符号を反転して DX/Ctrl にそろえる（Numerator/Denominator 列で判定）。
  q 値は Spectronaut の ratio 単位 t 検定（df = # of Ratios − 1。ペプチド × サンプル対を独立とみなす）なので
  個体数 n を反映しない＝過小評価気味。主に log2FC（効果量）で読む。

GSEA: log2FC で prerank（gseapy, 1000 perm, seed 0, min_size 5, max_size 1000）。
  gene set は prelim/chemo2/genesets.py の collection()（Hallmark 50 + niche MSigDB + 手製パネル）に、
  プロテオーム向けに NABA の matrisome 下位分類 6 set を事前に足した（結果を見て足し引きしていない）。
  定量数 125–400 なので、結果は参考値。感度分析として t = log2FC / SE の順位でも回す。

出力: out/prelim/proteome/
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
sys.path.insert(0, str(ROOT / 'prelim/chemo2'))
import common  # noqa: E402
import genesets as GS  # noqa: E402

SEQ = ROOT / 'data/seq/PXD078161'
OUT = ROOT / 'out/prelim/proteome'
common.OUT_DIR = str(OUT)

GROUPS = [  # (key, ファイル, ラベル, 年齢, 性別)
    ('YF', 'VG1_young females_Candidates.xlsx', 'Young ♀ (3 mo)', 'young', 'F'),
    ('AF', 'VG2_Aged_Females_Candidates.xlsx', 'Aged ♀ (20 mo)', 'aged', 'F'),
    ('YM', 'VG3_young_males_Candidates.xlsx', 'Young ♂ (3 mo)', 'young', 'M'),
    ('AM', 'VG4_aged_males_Candidates.xlsx', 'Aged ♂ (20 mo)', 'aged', 'M'),
]
LABEL = {k: lb for k, _, lb, *_ in GROUPS}
LFC_CUT = 0.58   # 著者の閾値 |log2FC| >= 0.58
Q_CUT = 0.05

EXTRA_SETS = [  # プロテオーム（ECM 中心）向けに事前に足す NABA 下位分類
    ('NABA_COLLAGENS', 'NABA collagens', 'ECM'),
    ('NABA_ECM_GLYCOPROTEINS', 'NABA ECM glycoproteins', 'ECM'),
    ('NABA_PROTEOGLYCANS', 'NABA proteoglycans', 'ECM'),
    ('NABA_ECM_REGULATORS', 'NABA ECM regulators', 'ECM'),
    ('NABA_ECM_AFFILIATED', 'NABA ECM-affiliated', 'ECM'),
    ('NABA_SECRETED_FACTORS', 'NABA secreted factors', 'ECM'),
]


# ---------- 読み込み ----------
def read_all_proteins(path):
    raw = pd.read_excel(path, 'All Proteins', header=None)
    title = str(raw.iloc[0, 0])
    note = str(raw.iloc[1, 0]) if raw.shape[0] > 1 else ''
    h = next(i for i in range(15) if (raw.iloc[i] == 'Comparison (group1/group2)').any())
    df = raw.iloc[h + 1:].copy()
    df.columns = raw.iloc[h]
    df = df.loc[:, df.columns.notna()]
    df = df.loc[:, ~df.columns.duplicated()]
    return df, title


def dx_vs_ctrl(df):
    """DX/Ctrl に向きをそろえた 1 タンパク質 1 行の表。"""
    num, den = df['Condition Numerator'].astype(str).str.strip(), df['Condition Denominator'].astype(str).str.strip()
    fwd = (num == 'DX') & (den == 'Ctrl')
    rev = (num == 'Ctrl') & (den == 'DX')
    d = df[fwd | rev].copy()
    sign = np.where(rev[fwd | rev], -1.0, 1.0)
    out = pd.DataFrame({
        'gene': d['Genes'].astype(str).str.split(';').str[0].str.strip(),
        'genes_all': d['Genes'].astype(str),
        'protein': d['ProteinDescriptions'].astype(str),
        'uniprot': d['UniProtIds'].astype(str),
        'log2fc': d['AVG Log2 Ratio'].astype(float) * sign,
        'se': d['Standard Error'].astype(float),
        'pvalue': d['Pvalue'].astype(float),
        'qvalue': d['Qvalue'].astype(float),
        'n_ratios': d['# of Ratios'].astype(float),
        'valid': d['Valid'].astype(str) if 'Valid' in d else 'NA',
        'flipped': sign < 0,
    })
    out = out[out.gene.notna() & (out.gene != 'nan') & (out.gene != '')]
    out = out.sort_values('n_ratios', ascending=False).drop_duplicates('gene')
    out['t'] = out.log2fc / out.se
    return out.set_index('gene')


def sample_key():
    sk = pd.read_excel(SEQ / 'Sample_Key.xlsx', header=1).dropna(how='all')
    sk['Age'] = sk.Age.astype(int)
    return sk


def load():
    D, rows = {}, []
    sk = sample_key()
    for k, f, lb, age, sex in GROUPS:
        df, title = read_all_proteins(SEQ / f)
        d = dx_vs_ctrl(df)
        D[k] = d
        s = sk[(sk.Age == (3 if age == 'young' else 20)) & (sk.Sex == ('Female' if sex == 'F' else 'Male'))]
        nn = s.Condition.value_counts()
        rows.append(dict(
            group=k, label=lb, sheet_title=title,
            n_ctrl_samplekey=int(nn.get('Control', 0)), n_dx_samplekey=int(nn.get('DX', 0)),
            n_dq_samplekey=int(nn.get('DQ', 0)), n_dxdq_samplekey=int(nn.get('DXDQ', 0)),
            proteins_all_comparisons=df['Genes'].nunique(),
            proteins_dx_vs_ctrl=len(d), flipped_sign=bool(d.flipped.any()),
            invalid_rows=int((d.valid == 'False').sum()),
            q_lt_005=int((d.qvalue < Q_CUT).sum()),
            q_lt_005_and_lfc=int(((d.qvalue < Q_CUT) & (d.log2fc.abs() >= LFC_CUT)).sum()),
            lfc_ge_058=int((d.log2fc.abs() >= LFC_CUT).sum()),
            median_abs_lfc=float(d.log2fc.abs().median()),
        ))
    S = pd.DataFrame(rows)
    return D, S


# ---------- タンパク質単位の比較 ----------
def compare(D, a, b):
    """a, b の DX 効果を共通タンパク質で並べ、カテゴリを付ける。"""
    m = D[a][['log2fc', 'qvalue', 'protein']].join(D[b][['log2fc', 'qvalue']], how='inner', lsuffix='_' + a, rsuffix='_' + b)
    la, lb_ = m['log2fc_' + a], m['log2fc_' + b]
    sa = (m['qvalue_' + a] < Q_CUT) & (la.abs() >= LFC_CUT)
    sb = (m['qvalue_' + b] < Q_CUT) & (lb_.abs() >= LFC_CUT)
    ea, eb = la.abs() >= LFC_CUT, lb_.abs() >= LFC_CUT

    def cat(x, y):
        c = np.where(x & y, np.where(np.sign(la) == np.sign(lb_), 'both_same', 'both_opposite'),
                     np.where(x, f'{a}_only', np.where(y, f'{b}_only', 'neither')))
        return c
    m['cat_sig'] = cat(sa, sb)      # q < 0.05 かつ |log2FC| >= 0.58（著者の基準）
    m['cat_effect'] = cat(ea, eb)   # |log2FC| >= 0.58 のみ（探索的）
    r_p = stats.pearsonr(la, lb_)
    r_s = stats.spearmanr(la, lb_)
    st = dict(a=a, b=b, n_common=len(m), pearson=r_p[0], pearson_p=r_p[1], spearman=r_s[0], spearman_p=r_s[1],
              **{f'sig_{k}': int((m.cat_sig == k).sum()) for k in ['both_same', 'both_opposite', f'{a}_only', f'{b}_only']},
              **{f'eff_{k}': int((m.cat_effect == k).sum()) for k in ['both_same', 'both_opposite', f'{a}_only', f'{b}_only']})
    st = {k.replace(f'_{a}_only', '_a_only').replace(f'_{b}_only', '_b_only'): v for k, v in st.items()}
    return m, st


# ---------- GSEA ----------
def collection():
    sets, meta = GS.collection()
    m2 = GS.read_gmt('m2.all')
    extra = []
    for k, lab, cat in EXTRA_SETS:
        sets[k] = m2[k]
        extra.append(dict(term=k, label=lab, category=cat, source='MSigDB M2 (NABA)', niche=True, size_total=len(m2[k])))
    meta = pd.concat([meta, pd.DataFrame(extra)], ignore_index=True)
    return sets, meta


def run_gsea(D, sets, metric):
    rows = []
    for k in D:
        v = D[k][metric].replace([np.inf, -np.inf], np.nan).dropna()
        rnk = pd.Series(v.values, index=v.index.str.upper()).groupby(level=0).first().sort_values(ascending=False)
        r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=5, max_size=1000,
                       seed=0, threads=8, outdir=None, no_plot=True, verbose=False).res2d
        r = r.assign(group=k, metric=metric, n_ranked=len(rnk))
        rows.append(r)
    L = pd.concat(rows, ignore_index=True).rename(columns={
        'Term': 'term', 'FDR q-val': 'fdr', 'NOM p-val': 'pval', 'Lead_genes': 'lead_genes', 'Tag %': 'tag_pct'})
    for c in ['NES', 'ES', 'fdr', 'pval']:
        L[c] = L[c].astype(float)
    L['matched'] = L.tag_pct.astype(str).str.split('/').str[1].astype(int)
    return L[['group', 'metric', 'term', 'NES', 'ES', 'pval', 'fdr', 'matched', 'lead_genes', 'n_ranked']]


# ---------- 図 ----------
def gene_class(genes, sets):
    up = pd.Index(genes).str.upper()
    core = set(sets['NABA_CORE_MATRISOME'])
    gly = set(sets['HALLMARK_GLYCOLYSIS'])
    cls = np.full(len(up), 'other', dtype=object)
    cls[[g in core for g in up]] = 'core matrisome'
    cls[[g in gly and g not in core for g in up]] = 'glycolysis (H)'
    cls[[g.startswith(('IGH', 'IGK', 'IGL', 'JCHAIN', 'HBA', 'HBB', 'KRT')) for g in up]] = 'Ig / keratin'
    return cls


CLS_COLOR = {'core matrisome': '#4A7BA7', 'glycolysis (H)': '#C44E52', 'Ig / keratin': '#E1A33B', 'other': '#C8C8C8'}
LIM = 3.0  # 軸範囲。外れる点は端に寄せて白抜きで描く


def scatter(ax, m, a, b, sets, title):
    x, y = m['log2fc_' + a], m['log2fc_' + b]
    cls = gene_class(m.index, sets)
    xc, yc = x.clip(-LIM + 0.1, LIM - 0.1), y.clip(-LIM + 0.1, LIM - 0.1)
    out = (xc != x) | (yc != y)
    for c in ['other', 'Ig / keratin', 'core matrisome', 'glycolysis (H)']:
        s = (cls == c)
        ax.scatter(xc[s & ~out], yc[s & ~out], s=10, c=CLS_COLOR[c], lw=0.3, edgecolor='white',
                   label=f'{c} ({s.sum()})', zorder=3 if c != 'other' else 2)
        ax.scatter(xc[s & out], yc[s & out], s=10, facecolor='white', edgecolor=CLS_COLOR[c], lw=0.8, zorder=3)
    ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM, LIM)
    ax.axhline(0, c='k', lw=0.5); ax.axvline(0, c='k', lw=0.5)
    ax.plot([-LIM, LIM], [-LIM, LIM], c='#999999', lw=0.5, ls='--')
    for v in (-LFC_CUT, LFC_CUT):
        ax.axhline(v, c='#BBBBBB', lw=0.4, ls=':'); ax.axvline(v, c='#BBBBBB', lw=0.4, ls=':')
    # ラベル: どちらかで q<0.05 かつ |lfc|>=0.58（著者基準）。Ig/keratin は除く
    lab = (m.cat_sig != 'neither') & (cls != 'Ig / keratin')
    for i, g in enumerate(m.index[lab]):
        ax.annotate(g, (xc[g], yc[g]), fontsize=6, xytext=[(4, 3), (4, -8), (-30, 4), (-30, -8)][i % 4], textcoords='offset points', zorder=4)
    rs = stats.spearmanr(x, y)
    ax.text(0.97, 0.03, f'{len(m)} common proteins\nSpearman ρ = {rs[0]:.2f} (P = {rs[1]:.2f})', transform=ax.transAxes,
            va='bottom', ha='right', fontsize=7)
    ax.set_xlabel(f'log2FC DX/Ctrl, {LABEL[a]}')
    ax.set_ylabel(f'log2FC DX/Ctrl, {LABEL[b]}')
    ax.set_title(title, fontsize=9)


HEAT_ROWS = [  # (term, 表示名)  — 年齢別に向きを見たい軸。matched < 5 の group は空欄
    ('HALLMARK_INFLAMMATORY_RESPONSE', 'H inflammatory response'),
    ('HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'H TNFa/NF-kB'),
    ('HALLMARK_IL6_JAK_STAT3_SIGNALING', 'H IL6/JAK/STAT3'),
    ('HALLMARK_INTERFERON_GAMMA_RESPONSE', 'H IFN-g response'),
    ('HALLMARK_COMPLEMENT', 'H complement'),
    ('HALLMARK_COAGULATION', 'H coagulation'),
    ('SAUL_SEN_MAYO', 'SenMayo'),
    ('NABA_CORE_MATRISOME', 'NABA core matrisome'),
    ('NABA_COLLAGENS', 'NABA collagens'),
    ('NABA_ECM_GLYCOPROTEINS', 'NABA ECM glycoproteins'),
    ('NABA_PROTEOGLYCANS', 'NABA proteoglycans'),
    ('NABA_ECM_REGULATORS', 'NABA ECM regulators'),
    ('GOBP_COLLAGEN_FIBRIL_ORGANIZATION', 'GO collagen fibril org.'),
    ('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'H EMT'),
    ('GOBP_OSSIFICATION', 'GO ossification'),
    ('GOBP_OSTEOBLAST_DIFFERENTIATION', 'GO osteoblast diff.'),
    ('GOBP_BONE_MINERALIZATION', 'GO bone mineralization'),
    ('GOBP_VASCULATURE_DEVELOPMENT', 'GO vasculature dev.'),
    ('HALLMARK_ANGIOGENESIS', 'H angiogenesis'),
    ('HALLMARK_HYPOXIA', 'H hypoxia'),
    ('HALLMARK_GLYCOLYSIS', 'H glycolysis'),
    ('HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'H OXPHOS'),
    ('HALLMARK_FATTY_ACID_METABOLISM', 'H fatty acid metab.'),
    ('HALLMARK_ADIPOGENESIS', 'H adipogenesis'),
    ('HALLMARK_MYOGENESIS', 'H myogenesis'),
]
RNA_COLS = [('Age_2.5y', 'RNA aging 2.5y vs 2m, cortical ♂'), ('B3_Dox4x_MSPC', 'RNA Dox 4x, LepR+ MSPC ♀'), ('B_Dox', 'RNA Dox 1x, bone ♂')]


def heat(ax, L, rna, meta_rows):
    keys = [k for k, *_ in GROUPS]
    W = L.pivot(index='term', columns='group', values='NES').reindex(index=[t for t, _ in meta_rows], columns=keys)
    keep = W.notna().any(axis=1)          # 4 群とも 5 タンパク質未満の set は行ごと落とす（キャプションに記載）
    meta_rows = [r for r, k in zip(meta_rows, keep) if k]
    W = W[keep]
    F = L.pivot(index='term', columns='group', values='fdr').reindex_like(W)
    Mt = L.pivot(index='term', columns='group', values='matched').reindex_like(W)
    R = rna.reindex(index=W.index)[[c for c, _ in RNA_COLS]]
    RF = rna_fdr.reindex(index=W.index)[[c for c, _ in RNA_COLS]]
    A = pd.concat([W, pd.DataFrame({'gap': np.nan}, index=W.index), R], axis=1)
    AF = pd.concat([F, pd.DataFrame({'gap': np.nan}, index=W.index), RF], axis=1)
    cmap = plt.get_cmap('RdBu_r').copy(); cmap.set_bad('#FFFFFF')
    im = ax.imshow(np.ma.masked_invalid(A.values.astype(float)), cmap=cmap, vmin=-2.5, vmax=2.5, aspect='auto')
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if A.columns[j] == 'gap':
                continue
            v = A.iat[i, j]
            if pd.isna(v):
                ax.text(j, i, '–' if j < len(keys) else 'n.t.', ha='center', va='center', fontsize=5.5, color='#AAAAAA')
                continue
            f = AF.iat[i, j]
            star = '**' if f < 0.05 else ('*' if f < 0.25 else '')
            ax.text(j, i, f'{v:.1f}{star}', ha='center', va='center', fontsize=6, color='white' if abs(v) > 1.8 else 'black')
    rng = Mt.apply(lambda r: f'{int(r.min())}–{int(r.max())}' if r.notna().any() else '', axis=1)
    ax.set_yticks(range(A.shape[0]), [f'{lb} [{rng.iloc[i]}]' for i, (_, lb) in enumerate(meta_rows)], fontsize=7)
    xl = ['Young ♀', 'Aged ♀', 'Young ♂', 'Aged ♂'] + [''] + [lb for _, lb in RNA_COLS]
    ax.set_xticks(range(A.shape[1]), xl, fontsize=7, rotation=50, ha='left', rotation_mode='anchor')
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(1.5, A.shape[0] - 0.2, 'Proteome DX vs Ctrl', ha='center', va='top', fontsize=7)
    ax.text(len(keys) + 2, A.shape[0] - 0.2, 'RNA (ref.)', ha='center', va='top', fontsize=7)
    return im


rna_fdr = None


def main():
    global rna_fdr
    OUT.mkdir(parents=True, exist_ok=True)
    D, S = load()
    S.to_csv(OUT / 'summary_groups.csv', index=False)
    print(S.to_string())
    long = pd.concat([d.assign(group=k) for k, d in D.items()]).reset_index()
    long.to_csv(OUT / 'dx_vs_ctrl_long.csv', index=False)
    W = pd.concat({k: d.log2fc for k, d in D.items()}, axis=1)
    Q = pd.concat({k: d.qvalue for k, d in D.items()}, axis=1)
    pd.concat({'log2fc': W, 'qvalue': Q}, axis=1).to_csv(OUT / 'dx_vs_ctrl_wide.csv')

    # タンパク質単位: 全ペアの相関、雌雄それぞれの若齢 vs 老齢のカテゴリ
    keys = list(D)
    st = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            m, s = compare(D, keys[i], keys[j])
            st.append(s)
            if (keys[i], keys[j]) in [('YF', 'AF'), ('YM', 'AM')]:
                m.to_csv(OUT / f'compare_{keys[i]}_{keys[j]}.csv')
                for c in ['cat_sig', 'cat_effect']:
                    print(keys[i], keys[j], c, m[c].value_counts().to_dict())
    ST = pd.DataFrame(st)
    ST.to_csv(OUT / 'pairwise_correlation.csv', index=False)
    print(ST[['a', 'b', 'n_common', 'pearson', 'spearman', 'spearman_p']].round(3).to_string())

    # GSEA
    sets, meta = collection()
    meta.to_csv(OUT / 'genesets_meta.csv', index=False)
    L = run_gsea(D, sets, 'log2fc')
    Lt = run_gsea(D, sets, 't')
    LL = pd.concat([L, Lt]).merge(meta[['term', 'label', 'category']], on='term', how='left')
    LL.to_csv(OUT / 'gsea_long.csv', index=False)
    for met, X in [('log2fc', L), ('t', Lt)]:
        NW = X.pivot(index='term', columns='group', values='NES')[keys]
        FW = X.pivot(index='term', columns='group', values='fdr')[keys]
        MW = X.pivot(index='term', columns='group', values='matched')[keys]
        X2 = pd.concat({'NES': NW, 'FDR': FW, 'matched': MW}, axis=1)
        X2.columns = [f'{a}_{b}' for a, b in X2.columns]
        X2.join(meta.set_index('term')[['label', 'category']]).to_csv(OUT / f'gsea_wide_{met}.csv')
    # log2FC と t の NES の一致
    cmp_ = L.merge(Lt, on=['group', 'term'], suffixes=('_lfc', '_t'))
    print('NES log2fc vs t, spearman by group:',
          cmp_.groupby('group').apply(lambda g: stats.spearmanr(g.NES_lfc, g.NES_t)[0]).round(2).to_dict())

    # RNA との突き合わせ（パスウェイ単位、向きの一致）
    nw = pd.read_csv(ROOT / 'out/prelim/chemo2/nes_wide.csv', header=[0, 1], index_col=0)
    rna = nw['NES']; rna_fdr = nw['FDR']
    rna.index.name = rna_fdr.index.name = 'term'
    rows = []
    for t in L.term.unique():
        if t not in rna.index:
            continue
        for k in keys:
            r = L[(L.term == t) & (L.group == k)]
            if r.empty:
                continue
            for rc, _ in RNA_COLS:
                rows.append(dict(term=t, group=k, rna=rc, nes_prot=float(r.NES.iloc[0]), fdr_prot=float(r.fdr.iloc[0]),
                                 matched=int(r.matched.iloc[0]), nes_rna=float(rna.at[t, rc]), fdr_rna=float(rna_fdr.at[t, rc])))
    C = pd.DataFrame(rows)
    C['same_sign'] = np.sign(C.nes_prot) == np.sign(C.nes_rna)
    C = C.merge(meta[['term', 'label', 'category']], on='term', how='left')
    C.to_csv(OUT / 'rna_concordance_long.csv', index=False)
    agg = []
    for (k, rc), g in C.groupby(['group', 'rna']):
        gg = g[(g.fdr_prot < 0.25) & (g.fdr_rna < 0.05)]
        agg.append(dict(group=k, rna=rc, n_sets=len(g), spearman=stats.spearmanr(g.nes_prot, g.nes_rna)[0],
                        same_sign_frac=g.same_sign.mean(), n_both_sig=len(gg), n_both_sig_same=int(gg.same_sign.sum())))
    A = pd.DataFrame(agg)
    A.to_csv(OUT / 'rna_concordance_summary.csv', index=False)
    print(A.round(2).to_string())

    # 図
    common.setup()
    fig = plt.figure(figsize=(common.mm(170), common.mm(120)))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.75)
    ax1 = fig.add_subplot(gs[0, 0]); ax1.set_box_aspect(1)
    m, _ = compare(D, 'YF', 'AF')
    scatter(ax1, m, 'YF', 'AF', sets, 'DX effect, young vs aged females')
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=6.5, frameon=False,
               handletextpad=0.2, columnspacing=0.8)
    common.panel_label(ax1, 'A', dx=-0.25, dy=1.03)
    ax2 = fig.add_subplot(gs[0, 1])
    im = heat(ax2, L, rna, HEAT_ROWS)
    cb = fig.colorbar(im, ax=ax2, fraction=0.05, pad=0.03, shrink=0.35, anchor=(0, 0))
    cb.set_label('NES', fontsize=7); cb.ax.tick_params(labelsize=6)
    common.panel_label(ax2, 'B', dx=-1.1, dy=1.12)
    fig.text(0.0, 0.0, 'A: open circles = clipped at axis limit; labels = q < 0.05 and |log2FC| ≥ 0.58 in either group '
             '(Spectronaut q, not adjusted for animal n).\n'
             'B: prerank GSEA on DX/Ctrl log2FC (min size 5). [a–b] = proteins matched across groups; – = < 5 proteins; '
             'n.t. = set not in RNA run. ** FDR < 0.05, * FDR < 0.25.\n'
             'Inflammatory response, TNFa/NF-kB, IL6/JAK/STAT3, IFN-g and OXPHOS had < 5 proteins in every group. '
             'Each proteome group is a separate batch (125–400 proteins).', fontsize=6, ha='left', va='top')
    common.save(fig, 'fig_proteome_dx_age')
    plt.close(fig)

    # 補助図: 雄の若齢 vs 老齢
    fig, ax = plt.subplots(figsize=(common.mm(80), common.mm(80)))
    m, _ = compare(D, 'YM', 'AM')
    scatter(ax, m, 'YM', 'AM', sets, 'DX effect: young vs aged males')
    ax.legend(loc='lower right', fontsize=5.5, frameon=False, handletextpad=0.2)
    common.save(fig, 'figS_proteome_males')
    plt.close(fig)


if __name__ == '__main__':
    main()

"""Figure 4A 案2・案3 の公開データを counts.tsv.gz / samples.csv にする。

  uv run python -I prelim/aging2/prepare.py          # repo root から

GSE235942（案2。依頼資料本文の「GSE23592」はニワトリの熱ストレスのアレイで別物。URL 側の GSE235942 が該当）:
  GSE235942_raw_counts.txt.gz（BioJupies = kallisto の遺伝子単位の整数 counts、シンボル行）をそのまま使う。
  列名は SOFT の description（A1, Y9, AM1 ...）。

GSE255653（案3）: GEO には FPKM（edgeR の rpkm、小数 2 桁）しか無い。原著の処理は DEE2 ワークフロー
  （STAR --quantMode GeneCounts、GRCm38）→ edgeR で FPKM。そこで FPKM から整数 counts を逆算する:
    count = FPKM × 遺伝子長[bp] × N[サンプルのライブラリサイズ] / 1e9
  遺伝子長は DEE2 の GeneInfo（mmusculus）の "merged"（全アイソフォームのエクソンを合併した長さ）。
    Σ FPKM×merged ≈ 0.98×1e9（mean / median / longest は 0.45–0.83）なので merged が edgeR に渡した長さと判断。
  シンボル → Ensembl は DEE2 のシンボルと MGI の MRK_ENSEMBL.rpt（現行シンボル）の両方から候補を引き、
    候補が複数あるときは全サンプルで整数に最も近くなる長さを選ぶ。長さが決まらない遺伝子（104、うち FPKM>1 は 8）は落とす。
  N はサンプルごとに「短い遺伝子（merged < 900 bp）の count が整数になる」N を格子探索で決める
    （score = mean cos(2π·count)。1 なら全部整数、無相関なら 0）。整数倍 kN も整数になるが、最小のものが真の N。
  検証値（復元の当否）は out/prelim/aging2/qc/GSE255653_reconstruction.csv に出す。

出力:
  data/seq/GSE235942/{counts.tsv.gz, samples.csv}
  data/seq/GSE255653/{counts.tsv.gz, samples.csv}
  out/prelim/aging2/qc/GSE255653_reconstruction.csv, GSE255653_gene_length_source.csv, composition_GSE235942.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'prelim'))
from prepare import read_soft  # noqa: E402  （SOFT の読み取りだけ借りる。prepare.py の main は呼ばない）

SEQ = ROOT / 'data/seq'
QC = ROOT / 'out/prelim/aging2/qc'


def write_new(df, path, **kw):
    """既存ファイルは上書きしない（data/ は Sync 上の実体）。"""
    if path.exists():
        print('exists, skip:', path)
        return
    df.to_csv(path, **kw)
    print('wrote', path)


# ---------------- GSE235942 ----------------
def gse235942():
    d = SEQ / 'GSE235942'
    c = pd.read_csv(d / 'GSE235942_raw_counts.txt.gz', sep='\t', index_col=0)
    assert (c.values % 1 == 0).all()
    s = read_soft('GSE235942')
    s = s.rename(columns={'description': 'sample', 'gender': 'sex_label'})
    s['sex'] = s.sex_label.str[0]
    s['age_m'] = s.age.str.extract(r'(\d+)')[0].astype(int)
    s['trt'] = s.treatment.map({'control': 'ctrl', 'HLU': 'HLU'})
    s['group'] = s.sex + s.age_m.astype(str) + '_' + s.trt
    s = s[['sample', 'gsm', 'title', 'sex', 'age_m', 'trt', 'group', 'tissue']]
    assert set(s['sample']) == set(c.columns)
    # 組成（筋・血液・Y 遺伝子）: 群ごとの比較用に記録しておく
    lib = c.sum()
    cpm = c / lib * 1e6
    comp = pd.DataFrame({'lib_M': lib / 1e6})
    for k, g in [('muscle', ['Acta1', 'Ckm', 'Myh4']), ('blood', ['Hbb-bs', 'Hba-a1']), ('Ychr', ['Ddx3y', 'Eif2s3y', 'Uty']),
                 ('ptprc', ['Ptprc']), ('osteocyte', ['Sost', 'Dmp1']), ('growth_plate', ['Col10a1', 'Acan'])]:
        comp[k] = np.log2(cpm.loc[g] + 1).mean()
    comp = comp.join(s.set_index('sample')[['group']])
    QC.mkdir(parents=True, exist_ok=True)
    comp.round(3).to_csv(QC / 'composition_GSE235942.csv')
    print(comp.groupby('group').mean().round(2).to_string())
    # 性別ラベルと Y 遺伝子の一致
    bad = comp[((comp.Ychr > 2) & s.set_index('sample').sex.eq('F')) | ((comp.Ychr < 2) & s.set_index('sample').sex.eq('M'))]
    print('sex label mismatch:', list(bad.index))
    write_new(c.astype(int), d / 'counts.tsv.gz', sep='\t')
    write_new(s, d / 'samples.csv', index=False)


# ---------------- GSE255653 ----------------
def integerness(x, N):
    """x = FPKM×len/1e9（= count/N）。N 倍した値の整数らしさ（mean cos 2πc）。"""
    c = x * N
    return np.cos(2 * np.pi * c).mean()


def fit_libsize(x):
    """x: 短い遺伝子の count/N。粗い格子 → 細かい格子で N を決める。"""
    N0 = 1 / x.min()   # count = 1 の遺伝子があれば ≈ N（丸め誤差で ±20% 程度ずれる）
    best = None
    for N in np.arange(0.6 * N0, 1.4 * N0, 200.0):
        xx = x[x * N < 2000]
        sc = integerness(xx, N)
        if best is None or sc > best[1]:
            best = (N, sc)
    for N in np.arange(best[0] - 400, best[0] + 400, 1.0):
        xx = x[x * N < 2000]
        sc = integerness(xx, N)
        if sc > best[1]:
            best = (N, sc)
    # 半分の N では奇数 count が半整数になるのでスコアは負になるはず（エイリアスの確認）
    half = integerness(x[x * best[0] < 2000], best[0] / 2)
    return best[0], best[1], half


def gse255653():
    d = SEQ / 'GSE255653'
    f = pd.read_csv(d / 'GSE255653_FPKM_Subseries_4.txt.gz', sep='\t', skiprows=[0], index_col=0)
    f.index.name = 'gene'
    gi = pd.read_csv(d / 'DEE2_mmusculus_GeneInfo.tsv.gz', sep='\t').set_index('GeneID')
    mgi = pd.read_csv(d / 'MRK_ENSEMBL.rpt.gz', sep='\t', header=None, usecols=[1, 5])
    mgi.columns = ['symbol', 'ens']
    cand = {}
    for sym, ens in zip(gi.GeneSymbol, gi.index):
        cand.setdefault(sym, []).append(('DEE2', ens))
    for sym, ens in zip(mgi.symbol, mgi.ens):
        if ens in gi.index and ens not in [e for _, e in cand.get(sym, [])]:
            cand.setdefault(sym, []).append(('MGI', ens))
    # 1) 候補が 1 つの遺伝子だけで N を決める
    uniq = {g: gi.loc[v[0][1], 'merged'] for g, v in cand.items() if len(v) == 1 and g in f.index}
    L = pd.Series(uniq)
    Ns, rows = {}, []
    for s in f.columns:
        x = (f.loc[L.index, s] * L / 1e9)
        sel = (L < 900) & (f.loc[L.index, s] > 0)
        N, sc, half = fit_libsize(x[sel])
        Ns[s] = N
        rows.append(dict(sample=s, libsize_N=round(N), integerness=sc, integerness_at_N_half=half,
                         n_genes_used=int((x[sel] * N < 2000).sum())))
    # 2) 候補が複数の遺伝子: 全サンプルで整数に最も近い長さを選ぶ
    Nv = pd.Series(Ns)[f.columns]
    length, src = {}, {}
    for g in f.index:
        if g not in cand:
            continue
        opts = cand[g]
        if len(opts) == 1:
            length[g], src[g] = gi.loc[opts[0][1], 'merged'], f'{opts[0][0]}:{opts[0][1]}'
            continue
        best = None
        for o, e in opts:
            ln = gi.loc[e, 'merged']
            c = f.loc[g] * ln * Nv / 1e9
            c = c[f.loc[g] > 0]
            sc = np.cos(2 * np.pi * c).mean() if len(c) else 0
            if best is None or sc > best[0] + 1e-9:
                best = (sc, ln, f'{o}:{e}(of {len(opts)})')
        length[g], src[g] = best[1], best[2]
    L = pd.Series(length)
    counts = f.loc[L.index].mul(L, axis=0).mul(Nv, axis=1) / 1e9
    resid = (counts - counts.round()).abs()
    tol = (0.005 * L.values[:, None] * Nv.values[None, :] / 1e9)   # FPKM の小数 2 桁丸めによる最大誤差
    tol = pd.DataFrame(tol, index=L.index, columns=f.columns)
    rc = counts.round().astype(int)
    # 復元の検証: 丸め誤差が小さい（tol < 0.1）遺伝子で、残差 < 0.1 の割合（ランダムなら ~0.2）
    for r in rows:
        s = r['sample']
        m = (tol[s] < 0.1) & (f.loc[L.index, s] > 0)
        r['n_genes_tol<0.1'] = int(m.sum())
        r['frac_resid<0.1_among_tol<0.1'] = float((resid.loc[m, s] < 0.1).mean())
        r['median_tol_counts'] = float(tol[s].median())
        r['sum_FPKMxlen/1e9'] = float((f.loc[L.index, s] * L).sum() / 1e9)
        r['reconstructed_total'] = int(rc[s].sum())
    q = pd.DataFrame(rows)
    QC.mkdir(parents=True, exist_ok=True)
    q.to_csv(QC / 'GSE255653_reconstruction.csv', index=False)
    print(q.round(4).to_string())
    pd.DataFrame({'length_merged': L, 'source': pd.Series(src)}).to_csv(QC / 'GSE255653_gene_length_source.csv')
    print('genes with length:', len(L), '/', len(f), '; dropped (no length):', len(f) - len(L),
          '; of which mean FPKM>1:', int((f.drop(L.index).mean(axis=1) > 1).sum()))
    s = read_soft('GSE255653')
    s = s.rename(columns={'title': 'sample'})
    s['group'] = s['sample'].str.extract(r'^(Young_Control|Aged_Control|Aged_TSPKO)')[0]
    s['age_m'] = s.age.str.extract(r'(\d+)')[0].astype(int)
    s = s[['sample', 'gsm', 'group', 'age_m', 'genotype', 'cell_type', 'description']]
    assert set(s['sample']) == set(rc.columns)
    write_new(rc, d / 'counts.tsv.gz', sep='\t')
    write_new(s, d / 'samples.csv', index=False)


if __name__ == '__main__':
    gse235942()
    gse255653()

---
title: "第 13 章 解析の方法"
date: 2026-10-26
author: 遠田
---

# 解析の方法（素人向けの説明）

## 差次発現解析（DESeq2）

RNA-seq のデータは、遺伝子ごとに「何本の断片（リード）が読まれたか」というカウントの表です。サンプルごとに読んだ総量が違うので、まず総量の差を補正し、そのうえで 2 群の間で遺伝子ごとに「差が偶然で説明できる範囲か」を検定します。
このための標準的な方法が DESeq2 で、本研究では Python 版（PyDESeq2 0.5.4）を全データ共通で使いました。
発現がごく低い遺伝子は検定してもぶれるだけなので、あらかじめ除きました（edgeR の filterByExpr と同じ基準）。比べる 2 群のサンプルだけを取り出して解析しています。
加齢マウスの皮質骨（Kaya 2022）だけは、原著者の DESeq2 の結果表をそのまま使いました。

## 閾値

DEG は全データ共通で「調整済み p < 0.05 かつ |log2FC| ≥ 1（2 倍以上の変化）」としました。2 万個の遺伝子を同時に検定するので、偶然の有意を減らすために Benjamini–Hochberg 法で p 値を補正しています。
結果を見てから閾値を動かすことはしていません。閾値を緩めても（調整済み p < 0.1）、OVX や Dox 単回の DEG が数個〜数十個にとどまることは確かめました（効果が小さいか、匹数が足りない）。

## 注目遺伝子のパネル

仮説に関わる遺伝子を**解析の前に**決めておきました。炎症・破骨 15 個（*Tnfsf11*（RANKL）、*Tnfrsf11b*（OPG）、*Il6*、*Il1b*、*Tnf*、*S100a8/9*、*Mmp9* など）と、アポトーシス・老化 17 個（*Fas*、*Casp3*、*Bax*、*Bcl2*、*Cdkn1a*（p21）、*Cdkn2a*（p16）、*Pmaip1*、*Bbc3* など）です。
結果を見てから都合のよい遺伝子を選ぶのを避けるためです。ヒトのデータでは同じ名前の遺伝子（大文字）に対応させました。

## GSEA（遺伝子セットの濃縮解析）

1 つ 1 つの遺伝子の変化は小さくても、同じ働きの遺伝子がそろって同じ向きに少しずつ動いていれば、それは意味のある変化です。これを拾うのが GSEA です。
全遺伝子を「上がった順」に並べ、ある遺伝子セット（例: 炎症反応の 200 遺伝子）が上の方か下の方に偏っているかを調べます。偏りの強さが **NES**（正なら上昇側、負なら低下側。±1.5 を超えるとかなり強い）、偶然でないかの指標が **FDR** です。FDR < 0.05 を有意としました。
遺伝子セットは公開のデータベース MSigDB（マウス版 2025.1）から、Hallmark 50 セット、niche に関わる 17 セット（SenMayo、破骨分化、骨化、基質、血管など）と、手製の niche パネル 5 つを解析の前に決めて使いました。
**SenMayo** は、骨の加齢から作られた老化・SASP の遺伝子セット（125 遺伝子）で、本研究の「老化」の主な指標です。Reactome の老化・SASP セットは中身の多くがヒストン（増殖の遺伝子）で、増殖が止まると下がってしまうので、老化の指標には使いませんでした。
FDR は解析に含めたセットの数で変わるので、v1（Hallmark 50 セットだけ）と今回（72 セット）では同じ比較でも値が少し違います（例: 加齢マウスの炎症反応 NES は v1 で +2.11、今回 +2.19）。

## 組成の補正

骨の bulk は、筋肉・血液・骨髄がどれだけ混ざったかで結果が大きく変わります。そこで、筋の遺伝子（*Acta1*、*Ckm*、*Myh4*）と血液の遺伝子（ヘモグロビン）の量からサンプルごとの「混入スコア」を作り、それを統計モデルに入れて差し引いた版も計算しました（Dox 単回の骨、ステロイドのヒト骨、TMS の骨）。
これはおおまかな補正で、細胞の割合を正確に推定（デコンボリューション）したものではありません。補正で消える結果は「混入の差」の可能性が高いと判断しました。

## 単一細胞データの擬似バルク

単一細胞のデータでは、細胞を似た発現のグループ（クラスタ）に分け、マーカー遺伝子で CAR 細胞・骨芽細胞などと名前を付けました。
そのうえで、同じ種類の細胞の発現を**個体ごとに合算**して 1 サンプルとし（擬似バルク）、bulk と同じ DESeq2 で比べました。細胞 1 個 1 個を反復として数えると p 値が極端に小さくなってしまうのを避けるためです。各群 1 匹のデータ（GSE289491）は統計にかけず、変化の大きさだけを示しました。

## マウスとヒト、データどうしの比較

異なるデータを比べるときは、遺伝子ごとの log2FC の順位の相関（Spearman ρ）と、遺伝子セットの NES の相関を使いました。遺伝子セットどうしは中身が重なるので、NES の相関は記述用で、検定ではありません。
加齢 × 化療を同じ実験で組んだ RNA データが無いので、**「重なると何が起きるか」（交互作用）は一度も検定していません**。この章の比較はすべて「別々の実験の向きを並べたもの」です。

# Methods 案（英文）

（v1 資料の Methods 案です。）

**Public transcriptomic data and differential expression analysis.** Publicly available bulk RNA-sequencing data of mouse and human bone were obtained from the Gene Expression Omnibus (GEO): GSE205237 (doxorubicin), GSE276391 (dexamethasone), GSE222752 (ovariectomy), GSE296758 (doxorubicin, LepR+ marrow stromal cells), GSE72815 (older vs younger women), GSE141595 (postmenopausal vs premenopausal women) and GSE276529 (glucocorticoid-induced osteoporosis). Gene-level count matrices deposited by the authors were used; counts were summed across identical gene symbols and rounded to integers, and matrices in which both mates of paired-end reads had been counted (GSE222752, GSE72815) were halved. For each comparison, only the two groups being compared were analysed. Genes were filtered with a rule equivalent to the default of edgeR filterByExpr (counts per million ≥ 10/median library size in millions in at least as many samples as the smaller group, and a total count ≥ 15). Differential expression was tested with PyDESeq2 (version 0.5.4; design ~ condition; Wald test; maximum-likelihood log2 fold changes without shrinkage), and P values were adjusted by the Benjamini–Hochberg method. Differentially expressed genes were defined as adjusted P < 0.05 and |log2 fold change| ≥ 1. For the aging comparison (2.5-year-old vs 2-month-old male mice, humeral cortical bone; BioProject PRJNA695408), the DESeq2 results provided by the original authors (Supplementary Table S1 of Kaya et al., J Bone Miner Res 2022) were used without re-analysis. In GSE276391, background reads from non-expressed genes (e.g., olfactory receptors) were more abundant in the dexamethasone group; adding a background index as a covariate did not change the results for the pre-specified gene panel. Pre-specified genes related to inflammation/osteoclastogenesis (Tnfsf11, Tnfrsf11b, Il6, Il1b, Il1a, Tnf, Ccl2, Cxcl12, S100a8, S100a9, Mmp9, Mmp13, Csf1, Sost, Dkk1) and apoptosis/senescence (Fas, Fasl, Casp3, Casp8, Casp9, Bax, Bak1, Bcl2, Bcl2l1, Trp53, Cdkn1a, Cdkn2a, Serpine1, Il6, Mdm2, Pmaip1, Bbc3) were defined before analysis; human orthologues were matched by gene symbol.

**Gene set enrichment analysis.** Genes were ranked by the Wald statistic and analysed by pre-ranked gene set enrichment analysis (GSEApy version 1.3.1; 1,000 permutations; gene set size 15–500) using the MSigDB Hallmark gene sets (MSigDB_Hallmark_2020). For human datasets, the SenMayo senescence gene set (MSigDB C2, SAUL_SEN_MAYO) was added. Gene sets with FDR < 0.05 were considered significant.

（GSE296758 は flowcell と群が完全に重なるため、design に flowcell を入れられませんでした。必要なら Methods に一文足します。）

# v3 の解析の方法と生成物（各解析メモより）

v3 の 4 本の解析メモの前書き・方法・生成物の節を、そのまま移しました。

## 4A 案 1・案 4（`docs/prelim_fig4a_v3.md`）

清水先生の 10/14 依頼資料（`data/20261014_遠田先生解析依頼資料.docx`）の Figure 4A 案1・案4 への回答。
方針は「骨細胞から骨全体の環境へ対象を広げる」。案2・案3（GSE235942 脛骨、GSE255653 LepR+ 間質）は別担当（`prelim/aging2/`、`docs/prelim_aging2.md`）。
再生成は `bash prelim/fig4a_v3/run_all.sh`（~45 分、GSEA は結果があれば飛ばす）。DE の設定・閾値は既存（`docs/prelim_human.md`）と同じ
（pydeseq2 0.5.4・MLE・edgeR filterByExpr 相当・DEG = padj < 0.05 かつ |log2FC| ≥ 1）。

### 生成物

- スクリプト（`prelim/fig4a_v3/`）: `common4a.py`（共通部品: DE・速い近似 DE・gene set・マーカー）, `pairing.py`（同じ人の照合）, `composition.py`, `dmab_de.py`（共変量つき DE）, `gsea4a.py`,
  `dmab_sets.py`, `reversal.py`（`--nes` でセット単位だけ、`--plot` で図だけ）, `inflam_apop.py`, `summary_cols.py`, `kaya.py`, `run_all.sh`
- 既存の `prelim/de.py`・`prelim/human/`・`prelim/chemo2/`・`prelim/invitro/`・`prelim/summary/` は import（または定義を借りる）だけで変更していない。依存の追加なし（statsmodels は既に入っていた）。
- 出力（`out/prelim/fig4a_v3/`）:
  - `pairing/`（GSE141595 × GSE141610 の照合）, `composition/`（組成, 図 `composition.png|pdf`）, `de/`（共変量つき DE と共変量）, `gsea/`（Hallmark / Reactome / GO BP / 要約形式）
  - `dmab/`（DEG, 有意セット, 画分の一致, 破骨スコアとの相関, 図 `dmab.png|pdf`）
  - `reversal/`（`summary.csv`, `genes.csv`, `sets_splithalf.csv`, `sets_nes*.csv`, 帰無分布, 図 `reversal.png|pdf`）
  - `inflam_apop/`（NES 表, サンプル単位の効果量, leading edge, 即時初期遺伝子を除いた NES, 図 `inflam_apop.png|pdf`）
  - `kaya/`（表の限界, 組成, 軸の NES, 注目遺伝子, GSEA, 案2 との一致 `concordance_aging2.csv`, ENA の run 表, 図 `kaya_axes.png|pdf`）
  - `summary_modules_values_fig4a_v3.csv`, `ena_runs_PRJNA594028.tsv`
- データの追加ダウンロードなし（GEO の正規化済みファイル 2 つは中身の確認だけ scratchpad に落とした: 閉経前の正規化ファイルは 13,182 遺伝子で、欠落した染色体の遺伝子も含むが閉経後と別に正規化されており、比較には使わなかった）。

## 4A 案 2・案 3（`docs/prelim_aging2.md`）

前提: 清水先生の 10/14 依頼資料（`data/20261014_遠田先生解析依頼資料.docx`）の Figure 4A 案2・案3。比較の相手は
`docs/prelim_aging_f.md`（雌の加齢）、`docs/prelim_niche.md`（LepR+/CAR の加齢）、`docs/prelim_chemo2.md`（Dox 列・gene set の流儀）、
`docs/prelim_summary_fig.md`（要約モジュール図）、`docs/prelim_human.md`（ヒト骨生検）。
再生成（repo root から）: `uv run python -I prelim/aging2/prepare.py` → `contam_ref.py` → `de.py` → `gsea.py` → `compare.py`。
出力は `out/prelim/aging2/`（図は png / pdf）。

### 2. 方法

- DE: pydeseq2（`prelim/aging_f/deutil.py` = `prelim/de.py` と同じ低発現フィルタ、MLE、shrink なし）。閾値は padj < 0.05 かつ |log2FC| ≥ 1。
- 組成補正版（GSE235942）: `~ muscle + blood + cond`。muscle = *Acta1*/*Ckm*/*Myh4*、blood = *Hbb-bs*/*Hba-a1* の log2 CPM 平均を、比較に使うサンプルの中で z 化（`prelim/chemo2/gsea.py` の B_Dox_compadj と同じ作り方）。
- GSEA: `prelim/chemo2/genesets.collection()`（Hallmark 50 + MSigDB M2/M5 17 + 手製 5 = 72 set、MSigDB 2025.1.Mm）。Wald 統計量で prerank、gseapy 1000 perm、seed 0、min 4、max 1000。
  GO BP 全体（size 15–500）も主な contrast で回した（探索用、`gsea/gobp_*.csv`）。
- 混入の感度分析（`_excl`）: 雌の骨髄の単一細胞 GSE210584（既存の h5ad）で、CAR 細胞より血球または内皮で 4 倍以上高く、かつ 10 CPM 以上の 1,904 遺伝子（血球 1,042 / 内皮 862）を定義した（`prelim/aging2/contam_ref.py`、`ref/contam_genes.csv`）。これを順位から抜いて GSEA をやり直した（DE はやり直していない）。
  10x の 3' データを基準にしているので、肥満細胞・形質細胞に特有の遺伝子は拾いきれない。
- サンプルごとのモジュールスコア（`sample_scores.png` A/B）: 混入しうる遺伝子を除いた set の遺伝子を median-of-ratios で正規化して log2 をとり、データ内で z 化して平均した。群の差が何例で支えられているかを見る記述用。
- 比較: マウス同士はシンボルで、ヒトとは大文字化（*Fasl*→*FASLG*、*Trp53*→*TP53*）で対応づけた（`compare_h.py` と同じ）。
  ヒトの NES は `out/prelim/summary/gsea_human/`（MSigDB Hs の同名 set）。既存の列は既存の csv を読むだけで、再計算していない。

### 8. 出力ファイル

スクリプト（`prelim/aging2/`）: `prepare.py`（counts・samples の作成、FPKM からのカウント復元）、`contam_ref.py`（混入しうる遺伝子の定義）、`de.py`、`gsea.py`、`compare.py`。

データ（`data/seq/`）:
- `GSE235942/`: `GSE235942_raw_counts.txt.gz`、`GSE235942_family.soft.gz`、`counts.tsv.gz`、`samples.csv`
- `GSE255653/`: `GSE255653_FPKM_Subseries_4.txt.gz`、`GSE255653_family.soft.gz`、`DEE2_mmusculus_GeneInfo.tsv.gz`（DEE2 の遺伝子長）、`MRK_ENSEMBL.rpt.gz`（MGI のシンボル → Ensembl）、`counts.tsv.gz`（復元）、`samples.csv`
- GSE23592 は落としていない（ニワトリのデータで無関係）。

結果（`out/prelim/aging2/`）:
- `de/<label>.csv`、`de/de_summary.csv`、`deg_counts.csv`
- `gsea/<label>.csv`（72 set）、`gsea/gobp_<label>.csv`（GO BP 全体）、`nes_long.csv`、`nes_wide.csv`
- `summary_modules_values_aging2.csv`（要約モジュール図と同じ列）、`module_grid.csv|png|pdf`
- `corr_log2fc.csv`、`corr_nes.csv`、`corr_matrix.png|pdf`、`heatmap_niche.png|pdf`
- `volcano_aging2.png|pdf`、`panel_genes.csv`、`heatmap_panel.png|pdf`
- `sample_scores.csv|png|pdf`、`markers_GSE235942.csv`、`markers_GSE255653.csv`
- `qc/composition_GSE235942.csv`、`qc/GSE255653_reconstruction.csv`、`qc/GSE255653_gene_length_source.csv`
- `ref/contam_genes.csv`、`ref/celltype_cpm_GSE210584.csv.gz`

## 4B（`docs/prelim_gse296761.md`）

清水先生の 10/14 依頼「【Figure4B：ドキソルビシン】2026 年 9 月に Blood に掲載されたものを活用しましょう」への回答。
前提: `docs/prelim_chemo.md`（GSE296758 の DE）、`docs/prelim_chemo2.md`（同じ gene set での並置）、v2 資料 `docs/report/shimizu_20261016.md`。
スクリプトは `prelim/gse296761/`、出力は `out/prelim/gse296761/`。再現手順は末尾。

（単一細胞の処理は第 2 章の「5.1 処理」にあります。）

### 再現手順・生成物

```
uv run --with scikit-image python -I prelim/gse296761/sc_cluster.py   # QC・クラスタ（scrublet の自動閾値に scikit-image。pyproject には足していない）
uv run python -I prelim/gse296761/sc_annotate.py                      # 細胞種・割合・T 細胞の確認・UMAP
uv run python -I prelim/gse296761/sc_compare.py                       # 擬似バルク・GSEA・スコア・細胞 DE・bulk との整合
uv run python -I prelim/gse296761/bulk.py                             # GSE296757/58 の軸ごとの整理（既存 DE を読むだけ）
uv run python -I prelim/gse296761/human.py                            # GSE318370
uv run python -I prelim/gse296761/modules.py                          # 横並び図・要約 csv
```

- スクリプト（`prelim/gse296761/`）: `panels.py`（事前定義の遺伝子パネル・追加 set・軸の定義）, `sc_cluster.py`, `sc_annotate.py`, `sc_compare.py`, `bulk.py`, `human.py`, `modules.py`
- データ（新規）: `data/seq/GSE296752/`（SOFT, PBS / DOX の h5）, `data/seq/GSE318370/`（SOFT, counts, `samples.csv`）
- 図（`out/prelim/gse296761/`、png / pdf）: `fig4b_modules`（軸 × データの横並び）, `bulk_panel_heatmap`, `sc_umap`, `sc_dotplot`, `sc_proportions`, `sc_panel_heatmap`, `sc_scores`
- 表: `summary_modules_values_gse296761.csv`, `summary_modules_values_gse296761_existing.csv`,
  `bulk/{panel_genes, axes_nes, gsea_extra_<label>, rankl_opg, rankl_opg_test, purity_cpm}.csv`,
  `sc/{qc, cluster_markers, cluster_top_genes, celltype_counts, celltype_counts_sensitivity, msc_score_split, t_cell_check, panel_pseudobulk, pseudobulk_all_genes(.gz), gsea_<ct>, gsea_extra_<ct>, sc_nes_long, scores_by_celltype, thin_factors, cellde_<ct>, bulk_concordance}.csv`,
  `human/{de_<label>, gsea_<label>, gsea_extra_<label>, panel_genes}.csv`
- 中間: `sc/adata_clustered.h5ad`, `sc/adata_annot.h5ad`（各 ~400 MB）

## 4C（`docs/prelim_fig4c_v3.md`）

前提: 清水先生の 10/14 依頼資料（`data/20261014_遠田先生解析依頼資料.docx`）の 4C（案1 GSE304562「Osteoblast における変化」、案2「現行の Figure 4C の活用」）、
10/13 の打合せ（`data/20261013_meeting.md`: DEX は閉経の前後を問わず動き話を混乱させる → 本編から外す方向）、
これまでの解析（`docs/prelim_de.md`, `docs/prelim_chemo2.md`, `docs/prelim_human.md`, `docs/prelim_aging_f.md`, `docs/prelim_summary_fig.md`）。予備段階の叩き台。

再生成（repo root から）: `bash prelim/fig4c_v3/run_all.sh`（`prepare.py` → `de.py` → `volcano.py` → `gsea.py` → `compare.py`）。出力は `out/prelim/fig4c_v3/`。
閾値・gene set・マーカーは解析前に固定し、結果を見て動かしていない。

### 2. 方法

- **counts の復元**（`prelim/fig4c_v3/prepare.py`）: 上記。Ensembl → シンボルは既存の `data/seq/annot/ensembl2symbol.csv`。引けない 144 ID は ID のまま（`prelim/prepare.py` と同じ扱い）。
- **DE**（`de.py`）: 既存 `prelim/de.py` の `run()` をそのまま呼ぶ（pydeseq2、`~ cond`、MLE、edgeR filterByExpr 相当の低発現フィルタ）。検定 13,868 遺伝子。
- **volcano**（`volcano.py`）: 既存 `prelim/volcano.py` の `plot()`。案2 は既存の `out/prelim/de/C_DEX.csv` を描いただけ。
- **GSEA**（`gsea.py`）: 案1 は `prelim/chemo2/genesets.py` の collection（Hallmark 50 + CURATED 17 + 手製 5、MSigDB 2025.1.Mm）で prerank（Wald 統計量、1000 perm、seed 0、min 4、max 1000 = chemo2 と同じ）。
  案2 は既存 `out/prelim/chemo2/nes_wide.csv` の値をそのまま使った。ヒト GIOP（素・骨髄補正・背景補正）は要約図のヒト列と同じく `prelim/invitro/gsea_iv.py` の collection で回し直した。
- **感度分析**: GR 標的 20 遺伝子（事前定義、`compare.py` の `GR_TARGETS`）を順位から抜いて同じ GSEA を回した（`gsea_noGR.csv`）。
- **骨環境の区画スコア**（`compare.py`）: 区画ごとに事前に決めた目印遺伝子（3–8 個、`COMPART`）の log2(CPM+1) の平均をサンプルごとに出し、DEX − 対照を Welch t 検定、14 区画で BH 補正。
  v1 の炎症パネル遺伝子（*S100a8/9*, *Mmp9*, *Tnf*, *Il1b*）は区画の目印に入れていない（循環を避けるため）。
- **一致の評価**: 遺伝子単位は共通遺伝子の Wald 統計量の Spearman と、共通 DEG（strict: padj < 0.05 & |log2FC| ≥ 1、lenient: padj < 0.1 & |log2FC| ≥ 0.585）の同方向・逆方向（超幾何 P）。
  ヒトとは大文字化したシンボルで照合。モジュール単位は同じ 68–71 set の NES の Spearman と、両方 FDR < 0.05 の set 数。set は遺伝子を共有するので独立ではなく、記述用。

### 11. 生成物

- スクリプト（`prelim/fig4c_v3/`）: `prepare.py`（counts の復元と samples.csv）, `de.py`（DESeq2・QC）, `volcano.py`, `gsea.py`（GSEA と要約図用 csv）, `compare.py`（一致・区画・感度分析・プレビュー）, `run_all.sh`
- データ（`data/seq/GSE304562/`）: `GSE304562_FPKM_allSamples.txt.gz`, `GSE304562_family.soft.gz`（GEO から取得）, `counts.tsv.gz`, `samples.csv`, `count_recovery.csv`（作成）
- 出力（`out/prelim/fig4c_v3/`）:
  - 図: `volcano/Cort_MC3T3.png|pdf`, `volcano/C4_options.png|pdf`（案1 | 案2）, `dex_environment.png|pdf`, `summary_preview.png|pdf`, `qc/pca.png|pdf`
  - 表: `de/Cort_MC3T3.csv`, `deg_counts.csv`, `qc/markers.csv`, `gsea/Cort_MC3T3.csv`, `nes_long.csv`, `gsea_human/H_GIOP*.csv`,
    **`summary_modules_cols.csv`（要約図に足す列）**, `gc_targets.csv`, `panel_genes.csv`, `gene_concordance.csv`, `module_concordance.csv`, `gsea_noGR.csv`,
    `dex_environment.csv`, `dex_environment_scores.csv`, `dex_environment_genes.csv`, `dex_inflam_vs_cells.csv`
- 依存の追加なし（GEO 確認用に一時的に `uv run --with python-docx` を使っただけ。pyproject は触っていない）。

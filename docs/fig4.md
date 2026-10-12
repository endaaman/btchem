# Figure 4 / Fig S1–S3 本番版 作図メモ（2026-10-12）

構成は `docs/prelim_summary.md` の決定どおり（4A 加齢, 4B Dox, 4C DEX、補助は Supplementary）。
入力は prelim の DE / GSEA 表（csv）で、**DE も GSEA も再計算していない**。

- 再生成: `bash fig/all.sh`（個別は `uv run python fig/fig4.py` など）
- スクリプト: `fig/common.py`（`prelim/common.py` のコピー。出力先 `out/fig` と出力サイズの確認だけ追加）、
  `fig/volcano.py`（volcano の共通描画）、`fig/fig4.py`, `fig/figS1.py`, `fig/figS2.py`, `fig/figS3.py`, `fig/all.sh`
- 出力: `out/fig/{fig4,figS1,figS2,figS3}.{png,pdf,tiff}`（600 dpi）、ラベル一覧 `out/fig/fig4_labels.csv`, `figS1_labels.csv`、
  ヒートマップの値 `figS2_nes.csv`, `figS2_fdr.csv`, `figS3_lfc.csv`, `figS3_padj.csv`
- 追加依存なし（adjustText を試したが、最終的にラベルは引き出し線の手動配置にしたので外した）

## 仕上がりサイズ（bbox tight 後の実寸）

| 図 | 幅 × 高さ (mm) | 備考 |
|---|---|---|
| Figure 4 | 167.3 × 76.6 | 1 行 3 パネル + 下に凡例 |
| Fig S1 | 169.5 × 125.4 | 2 行 3 列（6 枠目は凡例） |
| Fig S2 | 140.7 × 188.5 | 51 行 × 8 列 |
| Fig S3 | 141.1 × 131.8 | 30 行 × 8 列 |

フォント Arial（軸ラベル 9 pt、目盛 8 pt、遺伝子ラベル・凡例 7 pt、凡例の注 6.5 pt）、パネルラベル 12 pt 太字、線幅 ≥ 0.4 pt。
点群は PDF 内でラスタ化（600 dpi）、文字・線はベクタ。

## 作図ルール（volcano 共通、`fig/volcano.py`）

- 縦軸 −log10(padj)、横軸 log2FC（x は対称）。padj が NaN の遺伝子（DESeq2 の独立フィルタ / Cook 外れ値。加齢は著者表で FDR 空欄）は**描かない**。
- DEG = padj < 0.05 かつ |log2FC| ≥ 1。上昇 = 赤 `COLOR_HRO`、低下 = 青 `COLOR_LRO`、非有意 = 薄灰。閾値線は破線（|log2FC| = 1, padj = 0.05）。
- 注目遺伝子パネル（`prelim/genes.py`、事前定義）は **padj < 0.05 のものだけ**強調（FC は問わない）。記号: 炎症・破骨 = 丸、アポトーシス・老化 = 菱形（Il6 は両方に属するので丸）。
  塗りは変化の向き。|log2FC| < 1 の遺伝子は閾値線の内側に描かれるので、DEG かどうかは位置で読める。
- 清水指示「既報の炎症・破骨因子の上昇と矛盾しないものを赤く」は、**ラベルを太字（既報の向きと一致）/ 灰色細字（逆向き）**で表現した。
  点の色は有意な変化の向きのまま（有意でないものを有意に見せない）。「既報の向き」の定義（`EXPECTED`）:
  パネルの全遺伝子で上昇が促進方向、ただし **Tnfrsf11b（OPG）、Bcl2、Bcl2l1 は低下が促進方向**。Sost / Dkk1 は上昇を破骨促進側に置いた（Wnt 抑制→RANKL/OPG 比上昇）。
  この割り当ては議論の余地がある（特に Sost/Dkk1, Cdkn1a）ので、清水先生に確認したい。
- 軸は切っていない（指示どおり）。4A は −log10(padj) が 108 の遺伝子（低下側）があるため、パネル遺伝子は下の方に固まる。ラベルは上の空き領域に縦に並べ、引き出し線で結んだ。
  もっと見やすくするなら、拡大インセット（y 0–12）を足す案がある。
- 各パネル上端に DEG 数（Down / Up）。

## Figure 4

| パネル | 入力 | 比較 | n | DEG up / down | 描画した遺伝子数（padj あり） |
|---|---|---|---|---|---|
| A | `out/prelim/de/Age_2.5y.csv` | 2.5 歳 vs 2 か月、♂ C57BL/6、上腕骨皮質骨（骨髄は遠心で除去）| 4 vs 3 | **667 / 968** | 16,467（全 17,011 行） |
| B | `out/prelim/de/B_Dox.csv` | Dox 10 mg/kg i.p. 単回 30 日後 vs vehicle、♂ 15 週、脛骨+大腿骨（骨髄・骨端除去）| 4 vs 4 | **3 / 6** | 17,026 |
| C | `out/prelim/de/C_DEX.csv` | DEX 1 mg/kg s.c. 5 日 vs vehicle（非荷重群）、♂ 3 か月、脛骨 | 5 vs 5 | **518 / 596** | 18,880（全 19,266 行） |

A の出典: Kaya S et al., *J Bone Miner Res* 2022;37(4):804–817（NIHMS1776347, PMC9018503）、Supplementary Table S1 の "All results for 2m-2.5y"。
n（2 mo n=3, 2.5 y n=4）・性・系統・採取法は PMC 本文で確認した。生データは SRA BioProject **PRJNA695408**（GEO accession なし）。

### ラベルした遺伝子（padj < 0.05 のパネル遺伝子）

太字 = 既報の向きと一致、灰 = 逆向き。DEG = |log2FC| ≥ 1 も満たす。

| パネル | 遺伝子 | 区分 | log2FC | padj | DEG | 既報と一致 |
|---|---|---|---|---|---|---|
| A | *Tnf* | 炎症・破骨 | +1.36 | 5.2e-06 | ○ | ○ |
| A | *S100a8* | 炎症・破骨 | +1.06 | 1.1e-04 | ○ | ○ |
| A | *S100a9* | 炎症・破骨 | +0.90 | 0.0026 | | ○ |
| A | *Il1b* | 炎症・破骨 | +1.11 | 0.039 | ○ | ○ |
| A | *Tnfsf11* (RANKL) | 炎症・破骨 | −1.30 | 2.8e-09 | ○ | ×（低下） |
| A | *Fasl* | アポ・老化 | +1.44 | 4.2e-04 | ○ | ○ |
| A | *Pmaip1* (Noxa) | アポ・老化 | +1.14 | 8.8e-04 | ○ | ○ |
| A | *Bbc3* (Puma) | アポ・老化 | +0.72 | 1.5e-04 | | ○ |
| A | *Bak1* | アポ・老化 | +0.50 | 0.0011 | | ○ |
| A | *Bcl2* | アポ・老化 | +0.70 | 1.2e-04 | | ×（抗アポトーシスの上昇） |
| A | *Cdkn1a* (p21) | アポ・老化 | −0.60 | 0.0038 | | ×（低下） |
| B | — | | | | | padj < 0.05 のパネル遺伝子なし（図中に明記） |
| C | *S100a9* | 炎症・破骨 | +1.27 | 1.3e-07 | ○ | ○ |
| C | *Tnf* | 炎症・破骨 | +1.23 | 0.0012 | ○ | ○ |
| C | *Mmp9* | 炎症・破骨 | +0.97 | 5.0e-05 | | ○ |
| C | *Il1b* | 炎症・破骨 | +0.80 | 0.035 | | ○ |
| C | *Tnfrsf11b* (OPG) | 炎症・破骨 | −1.05 | 6.6e-05 | ○ | ○（OPG 低下 = RANKL/OPG 比上昇） |
| C | *Dkk1* | 炎症・破骨 | −0.90 | 9.6e-04 | | × |
| C | *Sost* | 炎症・破骨 | −0.57 | 0.014 | | × |
| C | *Cdkn1a* (p21) | アポ・老化 | −1.00 | 4.0e-04 | | ×（低下） |

読み: 4A は炎症（Tnf, Il1b, S100a8/9）とアポトーシス（Fasl, Noxa, Puma, Bak1）の両方が上がるが RANKL は下がる。4B は何も出ない。
4C は炎症（Tnf, Il1b, S100a9）・Mmp9 上昇と OPG 低下（破骨誘導方向）で、**アポトーシス系は 1 つも padj < 0.05 にならない**（Bak1 0.079 が最接近）。

## Fig S1（2 行 3 列）

| パネル | 入力 | 比較 | n | DEG up / down | ラベル（padj < 0.05 のパネル遺伝子） |
|---|---|---|---|---|---|
| A | `out/prelim/de/A_OVX_cort.csv` | GSE222752 OVX vs sham、♀、大腿骨皮質骨、術後 4 週（各試料 2 匹プール） | 3 vs 3 | 2 / 3 | なし |
| B | `out/prelim/chemo/de/B3_Dox4x_MSPC.csv` | GSE296758 Dox 5 mg/kg i.v. ×4（週 1）vs PBS、♀、LepR+ 骨髄 MSPC、最終投与 8 週後 | 4 vs 4 | 69 / 71 | *Mmp9* −1.91, *Mmp13* −1.27（灰 = 逆向き） |
| C | `out/prelim/human/de/H_Old.csv` | GSE72815 高齢（73 歳）vs 若年（30 歳）女性、腸骨針生検 | 19 vs 19 | 82 / 7 | *CDKN1A* +0.47, *CCL2* +0.80, *TNF* +0.35 |
| D | `out/prelim/human/de/H_PostMeno.csv` | GSE141595 閉経後 vs 閉経前（placebo）、腸骨生検 | 15 vs 15 | 237 / 28 | *CDKN1A* +1.00 (0.999), *IL6* +0.72 |
| E | `out/prelim/human/de/H_GIOP.csv` | GSE276529 GIOP vs 非 GIOP、女性 51–65 歳、大腿骨皮質骨 | 5 vs 4 | 618 / 298 | *TNFRSF11B* −1.00, *CDKN2A* +0.62, *MDM2* −0.39（灰） |

D の CDKN1A は log2FC 0.999 で |log2FC| ≥ 1 をわずかに満たさない（予備表の「+1.00」は丸め）。

## Fig S2（Hallmark NES ヒートマップ）

- 入力: `out/prelim/enrich/gsea_{Age_2.5y,A_OVX_cort,B_Dox,C_DEX}.csv`、`out/prelim/chemo/gsea/B3_Dox4x_MSPC.csv`、
  `out/prelim/human/enrich/gsea_{H_Old,H_PostMeno,H_GIOP}.csv`。マウス 4 列は `out/prelim/overlap/hallmark_nes.csv` と完全一致を確認。
- 列順: Aging, OVX | Dox (bone), Dox ×4 (MSPC) | DEX | Old, PostMeno, GIOP（縦線で群を区切る）。
- 行: **「どれか 1 つで FDR < 0.05」で絞ると 51 セット全部が残った**（どのセットもどこかで有意）。絞り込みが効かないので全行表示。
  行は NES の平均連結クラスタリング順、SenMayo（MSigDB C2 SAUL_SEN_MAYO、ヒトのみ評価）は最下行に分けた（マウス列は n.a.）。
  もっと絞るなら「2 contrast 以上で有意」か、Fig5 の議論に要るセット（Inflammatory Response, TNF-α/NF-κB, IL-6/JAK/STAT3, Apoptosis, p53, E2F, Estrogen Response, EMT, Myogenesis, SenMayo）だけの抜粋版を本文用に作る案がある。
- 列ごとの FDR < 0.05 のセット数: Aging 33, OVX 19, Dox bone 15, Dox×4 MSPC 27, DEX 23, Old 22, PostMeno 33, GIOP 12。
- 色: `COLOR_LRO`（−3）→ 中立灰 #F0EFEC（0）→ `COLOR_HRO`（+3）の発散 3 色、範囲外は端色（カラーバーに矢印）。FDR < 0.05 に「*」。
- Enrichr 版ライブラリの表記ゆれ（"Pperoxisome", "heme Metabolism", 二重スペース）は表示だけ直した。

## Fig S3（パネル遺伝子 log2FC ヒートマップ）

- 入力: Fig S1 / Fig 4 と同じ 8 本の DE 表（列は Fig S2 と同じ）。行はマウス記号、ヒトは大文字化で対応（Fasl → FASLG, Trp53 → TP53）。Il6 は 1 行（炎症・破骨側）。
- 色範囲 ±2（範囲外は端色）、padj < 0.05 に「*」。**検定対象外（低発現フィルタ、padj NA、H_PostMeno の閉経前ファイル欠落）は斜線**（0 付近の淡色と区別するため塗りではなくハッチ）。
- 予備の `out/prelim/panel_table.csv`, `out/prelim/human/panel_table.csv` と 7 列で log2FC・padj とも一致（Dox×4 MSPC は prelim の横並び表にない）。

## 数値の食い違い・注意

1. **`docs/prelim_chemo.md` の Dox×4 MSPC の Mmp9 padj「2e-4」は誤記**。DE 表（`B3_Dox4x_MSPC.csv`）と `bulk_panel.csv` はどちらも padj 5.4e-11（log2FC −1.91 は一致）。図は csv の値。
2. **「tested」の数え方が予備の表と違う**: `docs/prelim_de.md` の tested は Age_2.5y 16,976 / C_DEX 19,254 だが、csv は 17,011 / 19,266 行で、padj が NaN でない（＝volcano に描いた）のは 16,467 / 18,880。
   DEG 数（667/968, 518/596, 3/6 など）は完全に一致するので結論には影響しない。legend には DEG 数だけを書く。
3. 加齢（4A）の padj は著者の DESeq2 表の "FDR" 列。独立フィルタで FDR が空欄の 544 遺伝子は描いていない。著者の閾値・フィルタは本研究と異なる（Methods に明記）。
4. 4A は雄・上腕骨。「閉経」ではなく「加齢」として書く（legend・Limitation）。
5. 4C（GSE276391）は背景リード（嗅覚受容体など）が DEX 群で多い artefact があり、上昇 DEG の約 1/3 は Olfr/Gm 系（`docs/prelim_de.md`）。背景補正版でもパネルの結論は不変。Methods に一文。
6. Fig S1 D（GSE141595）は閉経前がバッチ完全交絡・遺伝子欠落あり。legend に明記。

---

# Figure legends 案（英語）

**Figure 4. Transcriptomic changes in mouse bone associated with aging, chemotherapy and glucocorticoid exposure.** Volcano plots of differential gene expression in (A) cortical bone of aged versus young mice (2.5 years, n = 4 vs 2 months, n = 3; male C57BL/6, humeral cortical bone after marrow removal; differential expression results reported by Kaya et al., BioProject PRJNA695408), (B) bone of mice 30 days after a single intraperitoneal dose of doxorubicin (10 mg/kg) versus vehicle (n = 4 per group; male, tibiae and femora without marrow; GSE205237) and (C) tibiae of mice treated with dexamethasone (1 mg/kg/day subcutaneously for 5 days) versus vehicle (n = 5 per group; male; non-loaded limbs; GSE276391). The x axis shows the log2 fold change (treated or aged versus control) and the y axis the −log10 Benjamini–Hochberg adjusted P value. Dashed lines indicate |log2 fold change| = 1 and adjusted P = 0.05. Red and blue dots denote significantly up- and downregulated genes (adjusted P < 0.05 and |log2 fold change| ≥ 1), with their numbers shown at the top of each panel; grey dots denote non-significant genes. Pre-specified genes related to inflammation/osteoclastogenesis (circles) and apoptosis/senescence (diamonds) are labelled when adjusted P < 0.05, irrespective of fold change; bold labels indicate a change in the direction reported to promote inflammation, osteoclastogenesis, apoptosis or senescence (for Tnfrsf11b [OPG], Bcl2 and Bcl2l1, a decrease), and grey labels the opposite direction. No pre-specified gene reached adjusted P < 0.05 in (B). Genes without an adjusted P value (filtered by independent filtering) are not shown.

**Figure S1. Transcriptomic changes in additional mouse and human bone datasets.** Volcano plots, drawn as in Figure 4, of (A) femoral cortical bone of ovariectomized versus sham-operated mice 4 weeks after surgery (n = 3 per group, each sample pooled from two mice; female; GSE222752), (B) LepR+ bone marrow mesenchymal stromal/progenitor cells 8 weeks after four weekly intravenous doses of doxorubicin (5 mg/kg) versus PBS (n = 4 per group; female; GSE296758), (C) iliac crest bone biopsies from older versus younger women (73 vs 30 years on average; n = 19 per group; GSE72815), (D) iliac bone biopsies from postmenopausal versus premenopausal women (n = 15 per group; GSE141595; premenopausal samples were sequenced in a separate batch and only genes present in both data files were analysed) and (E) femoral cortical bone from women with glucocorticoid-induced osteoporosis (GIOP) versus without (n = 5 vs 4; GSE276529). Human genes were matched to the mouse gene panel by orthologous symbols. GIOP, glucocorticoid-induced osteoporosis.

**Figure S2. Hallmark pathway enrichment across mouse and human bone datasets.** Normalized enrichment scores (NES) of the MSigDB Hallmark gene sets from pre-ranked gene set enrichment analysis (genes ranked by the Wald statistic) for each comparison shown in Figure 4 and Figure S1. Red indicates enrichment among genes upregulated in the aged, treated or disease group, and blue enrichment among downregulated genes. Asterisks indicate false discovery rate (FDR) < 0.05. The SenMayo senescence gene set was evaluated in human datasets only (n.a., not analysed). Rows are ordered by hierarchical clustering of NES. Aging, 2.5 years vs 2 months; OVX, ovariectomy vs sham; Dox, doxorubicin vs vehicle; MSPC, mesenchymal stromal/progenitor cells; DEX, dexamethasone vs vehicle; Old, older vs younger women; PostMeno, postmenopausal vs premenopausal women; GIOP, glucocorticoid-induced osteoporosis vs control.

**Figure S3. Expression changes of pre-specified inflammation/osteoclastogenesis and apoptosis/senescence genes.** Heatmap of log2 fold changes of the pre-specified gene panel in each comparison (columns as in Figure S2). Asterisks indicate adjusted P < 0.05. Hatched cells indicate genes not tested because of low expression, unavailable adjusted P values or absence from the data file. Mouse gene symbols are shown; human orthologues were matched by symbol (FASLG for Fasl, TP53 for Trp53).

# Methods 案（RNA-seq 解析）

**Public transcriptomic data and differential expression analysis.** Publicly available bulk RNA-sequencing data of mouse and human bone were obtained from the Gene Expression Omnibus (GEO): GSE205237 (doxorubicin), GSE276391 (dexamethasone), GSE222752 (ovariectomy), GSE296758 (doxorubicin, LepR+ marrow stromal cells), GSE72815 (older vs younger women), GSE141595 (postmenopausal vs premenopausal women) and GSE276529 (glucocorticoid-induced osteoporosis). Gene-level count matrices deposited by the authors were used; counts were summed across identical gene symbols and rounded to integers, and matrices in which both mates of paired-end reads had been counted (GSE222752, GSE72815) were halved. For each comparison, only the two groups being compared were analysed. Genes were filtered with a rule equivalent to the default of edgeR filterByExpr (counts per million ≥ 10/median library size in millions in at least as many samples as the smaller group, and a total count ≥ 15). Differential expression was tested with PyDESeq2 (version 0.5.4; design ~ condition; Wald test; maximum-likelihood log2 fold changes without shrinkage), and P values were adjusted by the Benjamini–Hochberg method. Differentially expressed genes were defined as adjusted P < 0.05 and |log2 fold change| ≥ 1. For the aging comparison (2.5-year-old vs 2-month-old male mice, humeral cortical bone; BioProject PRJNA695408), the DESeq2 results provided by the original authors (Supplementary Table S1 of Kaya et al., J Bone Miner Res 2022) were used without re-analysis. In GSE276391, background reads from non-expressed genes (e.g., olfactory receptors) were more abundant in the dexamethasone group; adding a background index as a covariate did not change the results for the pre-specified gene panel. Pre-specified genes related to inflammation/osteoclastogenesis (Tnfsf11, Tnfrsf11b, Il6, Il1b, Il1a, Tnf, Ccl2, Cxcl12, S100a8, S100a9, Mmp9, Mmp13, Csf1, Sost, Dkk1) and apoptosis/senescence (Fas, Fasl, Casp3, Casp8, Casp9, Bax, Bak1, Bcl2, Bcl2l1, Trp53, Cdkn1a, Cdkn2a, Serpine1, Il6, Mdm2, Pmaip1, Bbc3) were defined before analysis; human orthologues were matched by gene symbol.

**Gene set enrichment analysis.** Genes were ranked by the Wald statistic and analysed by pre-ranked gene set enrichment analysis (GSEApy version 1.3.1; 1,000 permutations; gene set size 15–500) using the MSigDB Hallmark gene sets (MSigDB_Hallmark_2020). For human datasets, the SenMayo senescence gene set (MSigDB C2, SAUL_SEN_MAYO) was added. Gene sets with FDR < 0.05 were considered significant.

（GSE296758 の design は `~ cond`（flowcell が群と完全交絡のため共変量に入れられない）。必要なら Methods に一文足す。）

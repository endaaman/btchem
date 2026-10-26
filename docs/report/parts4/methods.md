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

---
title: "Figure 4 / Supplementary Figure S1–S3 — 公開 RNA-seq データの解析（dry 解析）"
date: 2026-10-12
author: 遠田
---

# 先に結論

ご指示どおり、公開されているマウス・ヒトの骨の RNA-seq データで、仮説の各要素（閉経・化療・ステロイド）を偏りなく調べました。
閾値はすべてのデータで共通（調整済み p < 0.05 かつ |log2 fold change| ≥ 1）、注目遺伝子のパネル（炎症・破骨系 15 遺伝子、アポトーシス・老化系 17 遺伝子）は解析前に決めておきました。

| 仮説の要素 | データで見えたこと |
|---|---|
| 閉経後の骨（1st hit）＝ 炎症・破骨 | **加齢マウスの骨では炎症（Tnf, Il1b, S100a8/9）が上昇**。ヒト高齢・閉経後の骨では **老化（p21 = CDKN1A、SenMayo）** が上昇。ただし RANKL（*Tnfsf11*）はマウス加齢でむしろ低下し、ヒトでは動かない。OVX マウスは転写変化がほとんど無い |
| 化療（2nd hit）＝ アポトーシス | **どのデータでも、アポトーシス・老化の遺伝子は有意に上昇しなかった**。Dox 単回 30 日後の骨は変化自体がほぼ無く、雌の反復投与では**骨形成系と ER 標的の低下**が見えた |
| ステロイド（2nd hit）＝ アポトーシス | **アポトーシス系は有意に上昇せず**、代わりに **炎症因子の上昇と OPG の低下（RANKL/OPG 比の上昇）** がはっきり出た。OPG 低下はヒト GIOP でも一致 |

まとめると、**「化療・ステロイドが骨細胞のアポトーシスを起こす」という部分は、今回調べた公開データの転写レベルでは確認できませんでした**。
一方で、「閉経後（加齢）の骨は炎症・老化の状態にあり、ステロイドが RANKL/OPG のバランスを破骨側に傾け、化療が骨形成を抑える」という形であれば、データと矛盾なく説明できます。
Figure 5（概念図）はこの形に沿った修正案を作りました（別資料）。

なお、アポトーシスは主にタンパク質レベルで進む過程（カスパーゼの切断など）で、転写データでは捉えにくいという事情もあります。
化療・ステロイドによる骨細胞アポトーシスは組織学的には既報があります（Weinstein 1998、Wasi 2024 など。文献資料を参照）。
そのため、**論文では「既報の知見」として Discussion で触れ、本研究のデータの主張からは外す**のが安全だと考えます。

# Figure 4

![Figure 4（`out/fig/fig4.png`）](fig/fig4.png){width=100%}

ご指示の構成どおり、in vivo マウス骨の volcano plot を A / B / C の 3 枚にしました。

| パネル | データ | 比較 | n | DEG（上昇 / 低下） |
|---|---|---|---|---|
| A | Kaya et al. 2022（NIHMS1776347、BioProject PRJNA695408） | 2.5 歳 vs 2 か月、**雄** C57BL/6、上腕骨皮質骨（骨髄除去） | 4 vs 3 | 667 / 968 |
| B | GSE205237 | Dox 10 mg/kg 腹腔内 単回、30 日後 vs vehicle、雄 15 週、脛骨 + 大腿骨（骨髄除去） | 4 vs 4 | 3 / 6 |
| C | GSE276391 | DEX 1 mg/kg 皮下 5 日 vs vehicle、雄 3 か月、脛骨 | 5 vs 5 | 518 / 596 |

## 図の読み方

- 横軸 = 変化の大きさ（log2 fold change、右が上昇）、縦軸 = 有意さ（−log10 調整済み p、上ほど有意）。
- 赤 = 有意に上昇、青 = 有意に低下、灰 = 有意でない。破線は閾値（|log2FC| = 1、調整済み p = 0.05）。
- 注目遺伝子のうち調整済み p < 0.05 のものに名前を付けました（丸 = 炎症・破骨、菱形 = アポトーシス・老化）。
- ご指示の「既報の炎症・破骨因子の上昇と矛盾しないものを赤く」は、**ラベルの太字（既報の向きと一致）/ 灰色（逆向き）**で表しました。
  点の色は実際の変化の向きのままにしています（有意でないものを有意に見せないため）。
- 縦軸は切っていません（ご指示どおり）。4A は −log10 p が 108 に達する遺伝子があるため、注目遺伝子は下の方に固まります。見づらければ拡大図（インセット）を足せます。

## パネルごとの結果（`out/fig/fig4_labels.csv`）

**4A 加齢（2.5 歳 vs 2 か月）** — 炎症とアポトーシスの両方が上がる、唯一のマウスデータでした。

| 区分 | 遺伝子 | log2FC | 調整済み p | 既報の向き |
|---|---|---:|---:|---|
| 炎症・破骨 | *Tnf* | +1.36 | 5.2e-06 | 一致 |
| | *S100a8* | +1.06 | 1.1e-04 | 一致 |
| | *S100a9* | +0.90 | 0.0026 | 一致 |
| | *Il1b* | +1.11 | 0.039 | 一致 |
| | *Tnfsf11*（RANKL） | -1.30 | 2.8e-09 | 逆（低下） |
| アポトーシス・老化 | *Fasl* | +1.44 | 4.2e-04 | 一致 |
| | *Pmaip1*（Noxa） | +1.14 | 8.8e-04 | 一致 |
| | *Bbc3*（Puma） | +0.72 | 1.5e-04 | 一致 |
| | *Bak1* | +0.50 | 0.0011 | 一致 |
| | *Bcl2* | +0.70 | 1.2e-04 | 逆（抗アポトーシスの上昇） |
| | *Cdkn1a*（p21） | -0.60 | 0.0038 | 逆（低下） |

ただし、Hallmark の Apoptosis 遺伝子セット全体では有意な上昇はありません（NES -1.13、FDR 0.30、`out/fig/figS2_nes.csv`）。上がったのはパスウェイとしては炎症系（Inflammatory Response NES +2.11、IL-6/JAK/STAT3 +2.40）です。

**4B Dox 単回（30 日後）** — 変化がほとんどありません（DEG 9 個）。注目遺伝子で調整済み p < 0.05 のものは 1 つもありませんでした（図中に明記）。
元論文（Pimenta-Lopes et al. *Mol Ther* 2023、PMID 37143324）が報告する老化マーカー（p21、p16）の上昇は qPCR の結果で、RNA-seq の結果ではありません。
同じ閾値（FC の閾値なし）で数え直すと DEG は 29 個で、元論文の 31 個とほぼ一致するので、解析自体は元論文と矛盾していません。

**4C DEX 5 日** — 大きく動きます（DEG 1,114 個）。

| 区分 | 遺伝子 | log2FC | 調整済み p | 既報の向き |
|---|---|---:|---:|---|
| 炎症・破骨 | *S100a9* | +1.27 | 1.3e-07 | 一致 |
| | *Tnf* | +1.23 | 0.0012 | 一致 |
| | *Mmp9* | +0.97 | 5.0e-05 | 一致 |
| | *Il1b* | +0.80 | 0.035 | 一致 |
| | *Tnfrsf11b*（OPG） | -1.05 | 6.6e-05 | 一致（OPG 低下 = RANKL/OPG 比の上昇） |
| | *Dkk1* | -0.90 | 9.6e-04 | 逆 |
| | *Sost* | -0.57 | 0.014 | 逆 |
| アポトーシス・老化 | *Cdkn1a*（p21） | -1.00 | 4.0e-04 | 逆（低下） |

アポトーシス系で調整済み p < 0.05 になった遺伝子はありません（最も近いのは *Bak1* +0.30、調整済み p 0.079）。

# Supplementary Figure S1 — 追加のマウス・ヒトデータ

![Figure S1（`out/fig/figS1.png`）](fig/figS1.png){width=100%}

| パネル | データ | 比較 | n | DEG（上昇 / 低下） | 注目遺伝子（調整済み p < 0.05） |
|---|---|---|---|---|---|
| A | GSE222752 | OVX vs sham、雌、大腿骨皮質骨、術後 4 週 | 3 vs 3 | 2 / 3 | なし |
| B | GSE296758 | Dox 5 mg/kg 静注 × 4（週 1）、8 週後 vs PBS、雌、LepR+ 骨髄 MSPC | 4 vs 4 | 69 / 71 | *Mmp9* -1.91、*Mmp13* -1.27（逆向き） |
| C | GSE72815 | 高齢（平均 73 歳）vs 若年（30 歳）女性、腸骨生検 | 19 vs 19 | 82 / 7 | *CDKN1A* +0.47、*CCL2* +0.80、*TNF* +0.35 |
| D | GSE141595 | 閉経後 vs 閉経前（placebo）、腸骨生検 | 15 vs 15 | 237 / 28 | *CDKN1A* +1.00、*IL6* +0.72 |
| E | GSE276529 | GIOP vs 非 GIOP、女性 51–65 歳、大腿骨皮質骨 | 5 vs 4 | 618 / 298 | *TNFRSF11B*（OPG）-1.00、*CDKN2A* +0.62、*MDM2* -0.39（逆向き） |

（数値は `out/fig/figS1_labels.csv`。D の *CDKN1A* は正確には +0.999 で、|log2FC| ≥ 1 をわずかに満たしません。）

- **A（OVX）**: OVX が効いていること自体は確認できます（ER 標的の *Greb1* が低下）が、皮質骨の転写変化は小さく、炎症・破骨の遺伝子は 1 つも有意になりませんでした。そのため、閉経軸の主パネル（4A）には加齢データを使いました。
- **B（雌・反復投与の Dox）**: 骨形成系（*Spp1* -1.58、*Dmp1* -2.00）とエストロゲン標的（*Pgr* -2.73、*Greb1* -1.43）が低下しました。アポトーシス・老化の遺伝子は動きません（*Cdkn1a* -0.02、*Bax* +0.12、*Fas* +0.17）。
- **C・D（ヒト高齢・閉経後）**: どちらも老化マーカー p21（*CDKN1A*）と SenMayo（骨で作られた老化遺伝子セット）が上昇。2 つのデータの変化は互いに似ています（Spearman ρ = 0.56）。D は閉経前群が別バッチで、遺伝子の一部が欠けているため、単独では根拠にしません。
- **E（ヒト GIOP）**: マウス DEX と同じ論文のヒト版です。**OPG の低下がマウスとヒトで一致**しました。炎症系はヒトでは有意になりませんでした。

# Supplementary Figure S2 — パスウェイ（Hallmark）の比較

![Figure S2（`out/fig/figS2.png`）](fig/figS2.png){width=100%}

全 8 比較の Hallmark 遺伝子セットの濃縮スコア（NES）を並べたヒートマップです（赤 = 上昇側、青 = 低下側、* = FDR < 0.05）。主な行の値（`out/fig/figS2_nes.csv`, `figS2_fdr.csv`）:

| 遺伝子セット | 加齢 | OVX | Dox 骨 | Dox×4 MSPC | DEX | ヒト高齢 | ヒト閉経後 | ヒト GIOP |
|------------------------|------:|------:|------:|------:|------:|------:|------:|------:|
| Inflammatory Response | **+2.11** | **-1.76** | **-1.84** | **+1.42** | **+1.75** | -1.23 | **+1.81** | +0.83 |
| TNF-α / NF-κB | **+1.58** | -1.01 | **-2.17** | +1.08 | +1.09 | -1.10 | **+2.35** | -0.76 |
| IL-6 / JAK / STAT3 | **+2.40** | **-1.77** | -1.17 | -1.33 | **+2.15** | -1.08 | **+1.62** | +0.93 |
| Apoptosis | -1.13 | -1.14 | -1.30 | -1.22 | -1.30 | +1.19 | **+2.31** | -1.17 |
| p53 Pathway | **+1.37** | -1.16 | +1.09 | +1.07 | -1.34 | +1.41 | **+2.14** | +0.74 |
| SenMayo（ヒトのみ） | — | — | — | — | — | **+1.65** | **+1.72** | +0.85 |

太字は FDR < 0.05。

- Apoptosis が有意に上がったのは**ヒト閉経後 vs 閉経前（GSE141595）だけ**で、これは化療でもステロイドでもなく、しかも群とバッチが完全に重なっているデータです。化療・ステロイドのデータでは、Apoptosis はむしろ負の方向（有意でない）でした。
- 炎症系は、加齢・DEX で上がり、OVX・Dox 単回で下がり、データによって向きがばらばらです。

# Supplementary Figure S3 — 注目遺伝子の一覧

![Figure S3（`out/fig/figS3.png`）](fig/figS3.png){width=100%}

全 8 比較で、注目遺伝子パネルの log2FC を並べたヒートマップです（* = 調整済み p < 0.05、斜線 = 低発現などで検定対象外）。
1 枚で「どの比較でアポトーシス系が動かなかったか」を確認できます。

# アポトーシスについて — 調べたデータの一覧

化療・ステロイドによる骨のアポトーシスを探して、今回見たデータは次のとおりです。

| データ | 曝露・組織 | アポトーシス・老化の遺伝子 | 代わりに出たもの |
|---|---|---|---|
| GSE205237（4B） | Dox 単回 30 日後、雄、骨 | 有意なものなし（*Cdkn2a* +0.63、調整済み p 0.96） | 変化自体がほぼ無い。炎症系パスウェイは低下 |
| GSE296758（S1B） | Dox × 4、8 週後、雌、骨髄 MSPC | 有意なものなし | 骨形成系（*Spp1*, *Dmp1*）と ER 標的（*Pgr*, *Greb1*）の低下 |
| GSE296757 | Dox × 4、8 週後、雌、骨芽細胞 | 有意なものなし（Hallmark Apoptosis NES +1.44、FDR 0.064 で境界） | DEG 11 個のみ |
| GSE289491 | Dox 96 時間、雌、骨の細胞の scRNA（各群 1 匹） | p53 標的の一部（*Phlda3*, *Bax*）が上昇方向。p21 はむしろ低下、SASP も低い（n = 1 のため統計なし） | インターフェロン応答、脂肪化、増殖低下 |
| GSE276391（4C） | DEX 5 日、雄、脛骨 | 有意なものなし（*Bak1* +0.30、調整済み p 0.079） | 炎症因子の上昇、OPG 低下 |
| GSE276529（S1E） | ヒト GIOP、女性、皮質骨 | *CDKN2A* +0.62（p 0.035）、*BAX* +0.47（p 0.051）で境界的。補正の仕方で変わる | OPG 低下（マウスと一致） |

（数値は `docs/prelim_de.md`, `docs/prelim_chemo.md`, `docs/prelim_human.md` と `out/fig/` の csv。）

参考として、ヒト胚性幹細胞から分化させた骨芽細胞の培養に 5-FU・メトトレキサートを 7 日間直接当てたデータ（GSE324491）では、p53 → p21 / p16・*BAX* / *PUMA* の上昇がはっきり出ました。
**解析の手順そのものはアポトーシス・DNA 損傷応答を検出できる**ことの確認になります（ただし培養系・高濃度なので、患者の骨には当てはめられません）。

また、**化療を受けた乳がん患者の骨**、**乳がん骨転移モデルに化療をかけた骨**、**タキサン・シクロホスファミドを骨に当てたデータ**は、公開データには見つかりませんでした（`docs/datasets.md`）。

# 図に関する注意点

1. **4A の加齢データは雄**です（ご提案の加齢データそのもの）。閉経の「代理」としての加齢であり、閉経そのものではないので、本文・legend では「aging」と書き、Limitation に明記します（遠田の判断で、雄であることは「仕方ない」として Limitation に書く方針）。
2. 4A の調整済み p は原著者の解析結果（DESeq2）をそのまま使っており、閾値・フィルタが他のパネルと少し違います（Methods に明記）。
3. 4C（GSE276391）は DEX 群で背景のリード（嗅覚受容体など、骨で発現しない遺伝子）が多いという技術的な偏りがあり、上昇遺伝子の約 1/3 がこの系統です。背景を補正した解析でも注目遺伝子の結論（*Tnf*, *Il1b*, *S100a9*, *Mmp9* 上昇、OPG 低下）は変わりませんでした。
4. マウスとヒトの加齢データは遺伝子全体では**逆相関**しました（Spearman ρ -0.30 前後）。組織（マウス皮質骨 vs ヒト腸骨生検 = 骨髄が主体）と細胞組成が大きく違うためと考えられ、「マウスとヒトで同じシグネチャ」とは書けません。
5. 各データで組織・性・系統・投与法がばらばらで、腫瘍を持たない動物のデータです。Fig 4 は「仮説を立てるためのデータ（hypothesis-generating）」と位置づけます。

# 清水先生に確認したいこと

1. **「既報の向き」の割り当て**: *Sost* / *Dkk1* は上昇を破骨促進側（Wnt 抑制 → RANKL/OPG 比上昇）、*Cdkn1a* は上昇を老化側としました。違和感があれば教えてください（太字・灰色が入れ替わります）。
2. **アポトーシスの扱い**: 上記のとおりデータでは出なかったため、Figure 5 の骨格からは外し、Discussion で既報として触れる形をご提案します。
3. **4B の構成**: 現在は GSE205237（雄・単回）です。雌・反復投与の GSE296758 を 4B に上げ、GSE205237 を Supplementary に回す案もあります。
4. ヒト閉経前後データ（GSE141595）を本番で使う場合、生データからの再解析（欠けている遺伝子の回復）をするかどうか。

# Limitation 案（英文）

> Our transcriptomic analyses have several limitations. First, the datasets were generated independently and differ in species,
> sex, strain, age, tissue, drug, dose and timing; none was derived from tumor-bearing animals or from patients receiving
> perioperative chemotherapy, and to our knowledge no public transcriptome of bone exposed to chemotherapy in the context of
> breast cancer, or of bone from patients treated with chemotherapy, exists. These analyses are therefore hypothesis-generating.
> Second, the aging dataset used as a model of the postmenopausal bone (Figure 4A) was obtained from male mice, so it reflects
> aging rather than estrogen deficiency per se; ovariectomy itself produced only minimal transcriptional changes in cortical bone
> (Figure S1A). Third, the doxorubicin and dexamethasone datasets in Figure 4 were also derived from male mice, and the
> doxorubicin regimen (single dose, 30 days) differs from clinical anthracycline–taxane regimens. Fourth, apoptosis is largely
> regulated post-translationally, and bulk transcriptomes of bone are dominated by marrow, muscle and vascular components;
> the absence of an apoptotic signature therefore does not exclude osteocyte apoptosis, which has been reported histologically
> after doxorubicin and glucocorticoid exposure. Fifth, the human premenopausal and postmenopausal samples (GSE141595) were
> processed in different batches, and the human aging dataset reflects combined effects of aging and menopause.

# Methods 案（英文）

**Public transcriptomic data and differential expression analysis.** Publicly available bulk RNA-sequencing data of mouse and human bone were obtained from the Gene Expression Omnibus (GEO): GSE205237 (doxorubicin), GSE276391 (dexamethasone), GSE222752 (ovariectomy), GSE296758 (doxorubicin, LepR+ marrow stromal cells), GSE72815 (older vs younger women), GSE141595 (postmenopausal vs premenopausal women) and GSE276529 (glucocorticoid-induced osteoporosis). Gene-level count matrices deposited by the authors were used; counts were summed across identical gene symbols and rounded to integers, and matrices in which both mates of paired-end reads had been counted (GSE222752, GSE72815) were halved. For each comparison, only the two groups being compared were analysed. Genes were filtered with a rule equivalent to the default of edgeR filterByExpr (counts per million ≥ 10/median library size in millions in at least as many samples as the smaller group, and a total count ≥ 15). Differential expression was tested with PyDESeq2 (version 0.5.4; design ~ condition; Wald test; maximum-likelihood log2 fold changes without shrinkage), and P values were adjusted by the Benjamini–Hochberg method. Differentially expressed genes were defined as adjusted P < 0.05 and |log2 fold change| ≥ 1. For the aging comparison (2.5-year-old vs 2-month-old male mice, humeral cortical bone; BioProject PRJNA695408), the DESeq2 results provided by the original authors (Supplementary Table S1 of Kaya et al., J Bone Miner Res 2022) were used without re-analysis. In GSE276391, background reads from non-expressed genes (e.g., olfactory receptors) were more abundant in the dexamethasone group; adding a background index as a covariate did not change the results for the pre-specified gene panel. Pre-specified genes related to inflammation/osteoclastogenesis (Tnfsf11, Tnfrsf11b, Il6, Il1b, Il1a, Tnf, Ccl2, Cxcl12, S100a8, S100a9, Mmp9, Mmp13, Csf1, Sost, Dkk1) and apoptosis/senescence (Fas, Fasl, Casp3, Casp8, Casp9, Bax, Bak1, Bcl2, Bcl2l1, Trp53, Cdkn1a, Cdkn2a, Serpine1, Il6, Mdm2, Pmaip1, Bbc3) were defined before analysis; human orthologues were matched by gene symbol.

**Gene set enrichment analysis.** Genes were ranked by the Wald statistic and analysed by pre-ranked gene set enrichment analysis (GSEApy version 1.3.1; 1,000 permutations; gene set size 15–500) using the MSigDB Hallmark gene sets (MSigDB_Hallmark_2020). For human datasets, the SenMayo senescence gene set (MSigDB C2, SAUL_SEN_MAYO) was added. Gene sets with FDR < 0.05 were considered significant.

（GSE296758 は flowcell と群が完全に重なるため、design に flowcell を入れられませんでした。必要なら Methods に一文足します。）

# Figure legend 案（英文）

**Figure 4. Transcriptomic changes in mouse bone associated with aging, chemotherapy and glucocorticoid exposure.** Volcano plots of differential gene expression in (A) cortical bone of aged versus young mice (2.5 years, n = 4 vs 2 months, n = 3; male C57BL/6, humeral cortical bone after marrow removal; differential expression results reported by Kaya et al., BioProject PRJNA695408), (B) bone of mice 30 days after a single intraperitoneal dose of doxorubicin (10 mg/kg) versus vehicle (n = 4 per group; male, tibiae and femora without marrow; GSE205237) and (C) tibiae of mice treated with dexamethasone (1 mg/kg/day subcutaneously for 5 days) versus vehicle (n = 5 per group; male; non-loaded limbs; GSE276391). The x axis shows the log2 fold change (treated or aged versus control) and the y axis the −log10 Benjamini–Hochberg adjusted P value. Dashed lines indicate |log2 fold change| = 1 and adjusted P = 0.05. Red and blue dots denote significantly up- and downregulated genes (adjusted P < 0.05 and |log2 fold change| ≥ 1), with their numbers shown at the top of each panel; grey dots denote non-significant genes. Pre-specified genes related to inflammation/osteoclastogenesis (circles) and apoptosis/senescence (diamonds) are labelled when adjusted P < 0.05, irrespective of fold change; bold labels indicate a change in the direction reported to promote inflammation, osteoclastogenesis, apoptosis or senescence (for Tnfrsf11b [OPG], Bcl2 and Bcl2l1, a decrease), and grey labels the opposite direction. No pre-specified gene reached adjusted P < 0.05 in (B). Genes without an adjusted P value (filtered by independent filtering) are not shown.

**Figure S1. Transcriptomic changes in additional mouse and human bone datasets.** Volcano plots, drawn as in Figure 4, of (A) femoral cortical bone of ovariectomized versus sham-operated mice 4 weeks after surgery (n = 3 per group, each sample pooled from two mice; female; GSE222752), (B) LepR+ bone marrow mesenchymal stromal/progenitor cells 8 weeks after four weekly intravenous doses of doxorubicin (5 mg/kg) versus PBS (n = 4 per group; female; GSE296758), (C) iliac crest bone biopsies from older versus younger women (73 vs 30 years on average; n = 19 per group; GSE72815), (D) iliac bone biopsies from postmenopausal versus premenopausal women (n = 15 per group; GSE141595; premenopausal samples were sequenced in a separate batch and only genes present in both data files were analysed) and (E) femoral cortical bone from women with glucocorticoid-induced osteoporosis (GIOP) versus without (n = 5 vs 4; GSE276529). Human genes were matched to the mouse gene panel by orthologous symbols. GIOP, glucocorticoid-induced osteoporosis.

**Figure S2. Hallmark pathway enrichment across mouse and human bone datasets.** Normalized enrichment scores (NES) of the MSigDB Hallmark gene sets from pre-ranked gene set enrichment analysis (genes ranked by the Wald statistic) for each comparison shown in Figure 4 and Figure S1. Red indicates enrichment among genes upregulated in the aged, treated or disease group, and blue enrichment among downregulated genes. Asterisks indicate false discovery rate (FDR) < 0.05. The SenMayo senescence gene set was evaluated in human datasets only (n.a., not analysed). Rows are ordered by hierarchical clustering of NES. Aging, 2.5 years vs 2 months; OVX, ovariectomy vs sham; Dox, doxorubicin vs vehicle; MSPC, mesenchymal stromal/progenitor cells; DEX, dexamethasone vs vehicle; Old, older vs younger women; PostMeno, postmenopausal vs premenopausal women; GIOP, glucocorticoid-induced osteoporosis vs control.

**Figure S3. Expression changes of pre-specified inflammation/osteoclastogenesis and apoptosis/senescence genes.** Heatmap of log2 fold changes of the pre-specified gene panel in each comparison (columns as in Figure S2). Asterisks indicate adjusted P < 0.05. Hatched cells indicate genes not tested because of low expression, unavailable adjusted P values or absence from the data file. Mouse gene symbols are shown; human orthologues were matched by symbol (FASLG for Fasl, TP53 for Trp53).

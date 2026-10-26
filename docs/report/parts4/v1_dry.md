---
title: "第 12 章 v1 の Figure 4・Supplementary Figure S1–S3（再掲）"
date: 2026-10-26
author: 遠田
---

v1（10/12）資料の dry 解析の章（「Figure 4 / Supplementary Figure S1–S3 — 公開 RNA-seq データの解析（dry 解析）」）のうち、各案に振り分けなかった部分（横断の結論・図・表・注意点・Limitation・legend・確認事項）をそのまま再掲します。
Figure 4 の図と 4A は第 1 章 案 4、4B は第 2 章の関連データ（GSE205237）、4C は第 3 章 案 2 に、S1 の各パネルの説明は各案に移しました。Methods 案は第 13 章です。

# 先に結論（v1）

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

# Supplementary Figure S1 — 追加のマウス・ヒトデータ

![Figure S1（`out/v1/fig/figS1.png`）](v1/fig/figS1.png){width=100%}

| パネル | データ | 比較 | n | DEG（上昇 / 低下） | 注目遺伝子（調整済み p < 0.05） |
|---|---|---|---|---|---|
| A | GSE222752 | OVX vs sham、雌、大腿骨皮質骨、術後 4 週 | 3 vs 3 | 2 / 3 | なし |
| B | GSE296758 | Dox 5 mg/kg 静注 × 4（週 1）、8 週後 vs PBS、雌、LepR+ 骨髄 MSPC | 4 vs 4 | 69 / 71 | *Mmp9* -1.91、*Mmp13* -1.27（逆向き） |
| C | GSE72815 | 高齢（平均 73 歳）vs 若年（30 歳）女性、腸骨生検 | 19 vs 19 | 82 / 7 | *CDKN1A* +0.47、*CCL2* +0.80、*TNF* +0.35 |
| D | GSE141595 | 閉経後 vs 閉経前（placebo）、腸骨生検 | 15 vs 15 | 237 / 28 | *CDKN1A* +1.00、*IL6* +0.72 |
| E | GSE276529 | GIOP vs 非 GIOP、女性 51–65 歳、大腿骨皮質骨 | 5 vs 4 | 618 / 298 | *TNFRSF11B*（OPG）-1.00、*CDKN2A* +0.62、*MDM2* -0.39（逆向き） |

（数値は `out/v1/fig/figS1_labels.csv`。D の *CDKN1A* は正確には +0.999 で、|log2FC| ≥ 1 をわずかに満たしません。）

（各パネルの説明は、A → 第 1 章の関連データ（GSE222752）、B → 第 2 章 bulk、C・D → 第 1 章 案 1、E → 第 3 章 案 2 に移しました。）

# Supplementary Figure S2 — パスウェイ（Hallmark）の比較

![Figure S2（`out/v1/fig/figS2.png`）](v1/fig/figS2.png){width=100%}

全 8 比較の Hallmark 遺伝子セットの濃縮スコア（NES）を並べたヒートマップです（赤 = 上昇側、青 = 低下側、* = FDR < 0.05）。主な行の値（`out/v1/fig/figS2_nes.csv`, `figS2_fdr.csv`）:

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

![Figure S3（`out/v1/fig/figS3.png`）](v1/fig/figS3.png){width=100%}

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

（数値は `docs/prelim_de.md`, `docs/prelim_chemo.md`, `docs/prelim_human.md` と `out/v1/fig/` の csv。）

参考として、ヒト胚性幹細胞から分化させた骨芽細胞の培養に 5-FU・メトトレキサートを 7 日間直接当てたデータ（GSE324491）では、p53 → p21 / p16・*BAX* / *PUMA* の上昇がはっきり出ました。
**解析の手順そのものはアポトーシス・DNA 損傷応答を検出できる**ことの確認になります（ただし培養系・高濃度なので、患者の骨には当てはめられません）。

また、**化療を受けた乳がん患者の骨**、**乳がん骨転移モデルに化療をかけた骨**、**タキサン・シクロホスファミドを骨に当てたデータ**は、公開データには見つかりませんでした（`docs/datasets.md`）。

# 図に関する注意点（v1、5）

5. 各データで組織・性・系統・投与法がばらばらで、腫瘍を持たない動物のデータです。Fig 4 は「仮説を立てるためのデータ（hypothesis-generating）」と位置づけます。

（1・2・4 は第 1 章 案 4、3 は第 3 章 案 2 に移しました。）

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

# Figure legend 案（英文、Supplementary）

**Figure S1. Transcriptomic changes in additional mouse and human bone datasets.** Volcano plots, drawn as in Figure 4, of (A) femoral cortical bone of ovariectomized versus sham-operated mice 4 weeks after surgery (n = 3 per group, each sample pooled from two mice; female; GSE222752), (B) LepR+ bone marrow mesenchymal stromal/progenitor cells 8 weeks after four weekly intravenous doses of doxorubicin (5 mg/kg) versus PBS (n = 4 per group; female; GSE296758), (C) iliac crest bone biopsies from older versus younger women (73 vs 30 years on average; n = 19 per group; GSE72815), (D) iliac bone biopsies from postmenopausal versus premenopausal women (n = 15 per group; GSE141595; premenopausal samples were sequenced in a separate batch and only genes present in both data files were analysed) and (E) femoral cortical bone from women with glucocorticoid-induced osteoporosis (GIOP) versus without (n = 5 vs 4; GSE276529). Human genes were matched to the mouse gene panel by orthologous symbols. GIOP, glucocorticoid-induced osteoporosis.

**Figure S2. Hallmark pathway enrichment across mouse and human bone datasets.** Normalized enrichment scores (NES) of the MSigDB Hallmark gene sets from pre-ranked gene set enrichment analysis (genes ranked by the Wald statistic) for each comparison shown in Figure 4 and Figure S1. Red indicates enrichment among genes upregulated in the aged, treated or disease group, and blue enrichment among downregulated genes. Asterisks indicate false discovery rate (FDR) < 0.05. The SenMayo senescence gene set was evaluated in human datasets only (n.a., not analysed). Rows are ordered by hierarchical clustering of NES. Aging, 2.5 years vs 2 months; OVX, ovariectomy vs sham; Dox, doxorubicin vs vehicle; MSPC, mesenchymal stromal/progenitor cells; DEX, dexamethasone vs vehicle; Old, older vs younger women; PostMeno, postmenopausal vs premenopausal women; GIOP, glucocorticoid-induced osteoporosis vs control.

**Figure S3. Expression changes of pre-specified inflammation/osteoclastogenesis and apoptosis/senescence genes.** Heatmap of log2 fold changes of the pre-specified gene panel in each comparison (columns as in Figure S2). Asterisks indicate adjusted P < 0.05. Hatched cells indicate genes not tested because of low expression, unavailable adjusted P values or absence from the data file. Mouse gene symbols are shown; human orthologues were matched by symbol (FASLG for Fasl, TP53 for Trp53).

（Figure 4 の legend 案は第 1 章 案 4 の「図の候補」にあります。）

# 清水先生に確認したいこと

1. **アポトーシスの扱い**: 上記のとおりデータでは出なかったため、Figure 5 の骨格からは外し、Discussion で既報として触れる形をご提案します。（v1・Figure 4 の章の確認 2）

（確認 1 は第 1 章 案 4、3 は第 2 章 bulk、4 は第 1 章 案 1 に移しました。）

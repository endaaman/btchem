---
title: "仮説の経路オプション — 閉経/加齢と化療を骨転移に繋ぐ"
date: 2026-10-13
author: 遠田
---

# この資料について

dry 解析の到達点は「**閉経/加齢（A）と化療（B）は、骨で別のことを起こす**」です。

- **A（閉経/加齢）**: ヒト閉経後・高齢の骨生検で老化（SenMayo、p21）と p53 が上がる。マウス雄の加齢骨では炎症が上がる。
- **B（化療）**: 生体マウスでは炎症・老化は上がらず、雌 LepR+ 間質で**骨形成・ECM・血管の低下と ER 応答の低下**が出る。ヒト骨芽細胞 in vitro では**アポトーシス・老化・SASP** が強く上がる。
- A と B の掛け合わせ（交互作用）を見た公開データは無い。

これを骨転移（休眠 DTC の再活性化、破骨細胞の活性化）に繋ぐ経路は 1 本ではなく、**どの経路を主にするかは、橋渡しする先行研究の質と適合性で決まります**。
ここでは経路を 5 つに分け、それぞれについて (1) 流れ、(2) 本データのどの結果が乗るか、(3) 橋渡し文献の表、(4) 未検証の飛躍、(5) 反証を並べました。
**決め打ちはしていません**。最後に判断の観点と、経路別の Discussion 段落案（英文）を付けました。

表の記号:

- **適合性**: ◎ = 本仮説の一段をほぼそのまま支える / ○ = 支えるが条件（臓器・種・薬剤）がずれる / △ = 間接的・総説・方向の一部のみ。
- **系**: 種・細胞・モデルと in vivo / in vitro。文献の要点はすべて抄録の記述から（本文は読んでいない）。PMID は `docs/literature.md`・`literature2.md`・`datasets3_lit.md` で E-utilities により確認済みのものだけを使った。

# 経路の一覧

| 経路 | 流れ（要約） | 本データの乗り方 | 文献の厚さ | 臨床（BMA）との整合 |
|---|---|---|---|---|
| 1. 老化・SASP → 破骨 | 老化骨芽細胞/骨細胞の SASP（IL-6）→ 破骨 → DTC 再活性化 | A の老化 ◎、B は in vitro のみ | 厚い（骨・乳がんの原著あり） | ○（破骨を止める BMA で下流を切れる） |
| 2. 骨形成/血管 niche の喪失 | 休眠維持シグナル（GDF10/TGF-β2、E-selectin/SDF-1、TSP-1）が消える → 休眠解除 | B の niche 低下 ◎（マウス）、ヒトは逆向き | 各段はあるが繋いだ実験は無い | △（BMA で niche は戻らない） |
| 3. 間質 ER 応答の消失 | ストローマ傷害 + エストロゲン欠乏 → 分泌性老化 → ER+ 休眠細胞の再増殖 | B の ER 応答 ↓（マウス雌のみ） | 薄い（in vitro 1 本が中核） | △ |
| 4. RANKL/OPG 軸（清水原案） | 閉経・化療・ステロイド → RANKL/OPG ↑ → 破骨 → DTC 再活性化 | ほぼ乗らない（DEX の OPG ↓ のみ） | 臨床・前臨床とも最も厚い | ◎ |
| 5. その他 | IFN、NET、宿主反応一般、ステロイドの腫瘍直接作用、血管 pericyte | IFN ↑ は乗るが向きが逆 | 臓器が肺・肝のものが多い | — |

# 経路 1: 老化した骨芽細胞/骨細胞の SASP → 破骨活性化 → DTC 再活性化

## 流れ

- 1st: 閉経後・加齢の骨 → 骨芽細胞系・骨細胞の老化（p21、SenMayo）
- 2nd: 化療（DNA 傷害）→ 骨芽細胞系の老化をさらに誘導（p53 → 増殖停止 → SASP: IL-6, IL-1β, CXCL8, CCL2）
- → SASP が局所の破骨細胞形成を促す → 骨吸収 → 基質由来の増殖因子の放出 → 休眠 DTC の再活性化 → 骨転移
- 介入点: 破骨（デノスマブ/BP）は下流を切れる。上流（老化）は senolytics が候補だが臨床には無い。

## 本データのどの結果が乗るか

| 結果 | 値 | 乗り方 |
|---|---|---|
| ヒト閉経後 vs 閉経前の骨生検（GSE141595） | SenMayo NES +1.76（FDR 0.002）、p53 +2.17、CDKN1A +1.00（padj 4e-07）、IL6 +0.72（padj 0.002） | 1st の老化 ◎（ただしバッチ交絡） |
| ヒト高齢 vs 若年（GSE72815） | SenMayo +1.65（FDR 0.011）、CDKN1A +0.47（padj 0.0006） | 1st の老化の再現 ◎。炎症パスウェイは上がらない |
| ヒト骨芽細胞・MSC × エトポシド in vitro（GSE284666） | SenMayo 1.8–2.7、p53 1.9–2.7（FDR < 0.01）、IL6 +2.7 / +1.3、IL1B +5.9 / +5.4、CXCL8 +3.8 / +6.0（NHO / MSC） | 2nd の老化・SASP ◎（in vitro） |
| マウス生体の Dox（3 データ） | SenMayo は上がらない（-1.38 / -1.41） | **乗らない**（時間・希釈で説明する必要） |
| 雌マウス AC 2 週後の皮質骨 qPCR（PMID 40389522 の個体値を再計算） | Cdkn1a × 2.1、Cdkn2a × 1.9、Mmp13 × 2.2（p < 0.05）、Il6 × 0.73 | 2nd の老化は急性期に生体で出る ○。IL-6 は上がらない |
| 要約図の Osteoclast support 行 | ほぼ空白（有意はマウス雄加齢の +1.60 のみ） | 出口（破骨）は**データでは支えられない** |
| ヒト デノスマブ vs placebo（GSE141595） | 破骨・骨形成マーカーは ↓、CDKN1A・IL6・SenMayo は戻らない（有意でない） | RANKL を止めても老化シグネチャは残る、という読みは可能（効果量小・n = 15 で「戻らない」の証明ではない） |

## 橋渡し文献

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 適合性 | 弱点 |
|---|---|---|---|---|---|
| 26725121 | Cell Rep 2016（Luo） | マウス骨、乳がん、in vivo | 老化骨芽細胞を誘導するだけで乳がんの骨コロニー形成が増える。老化骨芽細胞は石灰化せず、SASP の IL-6 で局所の破骨形成を増やす。IL-6 中和で両方が止まる | ◎ | 老化は実験的に誘導したもので、閉経・化療による老化ではない |
| 28825716 | Nat Med 2017（Farr） | 加齢マウス、in vivo + 馴化培地 | 老化細胞の除去で骨量・骨質が改善、機序は骨吸収の低下。老化細胞の馴化培地は石灰化を阻害し破骨形成を促す | ◎（1st 側） | 転移は扱っていない |
| 39312185 | Cancer Res 2024（Kaur） | マウス・患者の骨転移、scRNA-seq | 骨転移巣の骨細胞が早期老化し破骨誘導性 SASP を獲得。senolytics で骨吸収抑制・骨量保持 | ○ | 老化は腫瘍が来た**後**の変化で、転移の前段（niche）ではない |
| 31932453 | Cancer Res 2020（Yao） | 化療マウス、in vivo | 化療後の骨減少は老化細胞・SASP による（エストロゲン非依存）。p38MAPK-MK2 阻害・senolytics で保護 | ◎（2nd 側） | 転移は扱っていない |
| 41469373 | Nat Commun 2025 | 雌マウス 12 週、Doxo / PTX、PTX 下の骨転移モデル（PyMT-BO1 心腔内接種） | 化療による骨量減少は、老化した CAR 細胞・骨髄脂肪細胞の RANKL を介する。骨転移モデルの値もある | ◎（2nd → 破骨 → 転移を 1 本で持つ） | 若齢のみ。PTX の骨転移モデルの結論は抄録からは確認していない（Source Data に値がある） |
| 27979832 | Cancer Discov 2017（Demaria） | マウス、in vivo | 化療誘発の老化細胞と SASP が有害事象を増やし、再発・転移を促進。除去で改善 | ○ | 骨に特化していない |
| 39147353 | Bone 2024（Wasi） | 成熟マウス（50 週）、Dox × 6 回、in vivo | Dox が骨細胞アポトーシスを増やしリモデリングを乱す。TUNEL 個体値: 非投与 9.9 % → 低用量 47 % → 高用量 81 % | ○ | 老化ではなくアポトーシス。若齢群が無い |
| 42686970 | GeroScience 2026（Izquierdo） | マウス、Dox、性・年齢別 | Dox で老化/炎症遺伝子 ↑、効果は雌・高齢で強い。D+Q で部分回復 | ○（A × B の唯一の近傍） | 転移は扱っていない |
| 22137794 | Cancer Cell 2011 | 乳がん骨転移休眠モデル | VCAM-1 が α4β1+ 破骨前駆細胞を呼び、休眠微小転移を骨溶解病変に転換 | ○（出口） | 起点は腫瘍側の VCAM-1 |
| 26632274 | Nat Commun 2015（Lawson） | マウス骨髄腫 | 破骨によるリモデリングが休眠細胞を再活性化 | ○（出口） | 骨髄腫 |
| 24687923 / 25603921 | Clin Cancer Res 2014 / Int J Cancer 2015（Ottewell） | OVX マウス | OVX だけで骨内の休眠 DTC が増殖、ZOL / OPG-Fc で阻止。骨外は無関係 | ◎（出口と閉経後特異性） | 老化は扱っていない |

## 経路として成立するための未検証の飛躍

1. **転写の老化シグネチャ ≠ SASP の分泌**。ヒト生検の SenMayo ↑ が局所の IL-6 などのタンパク量に結びつくかは見ていない。IL-6 は v1 のほとんどのデータで発現が低く測れない。
2. **化療の老化は in vitro と急性期だけ**。生体マウスの Dox（30 日〜8 週後）では老化は出ない。in vitro のエトポシドは臨床レジメンではない。
3. 「老化 → 破骨 → DTC 再活性化」を**閉経後 × 化療の宿主**で示した実験は無い（Luo 2016 は誘導した老化、Kaur 2024 は転移後）。
4. **老化が加算される（A の老化に B の老化が上乗せされる）**ことは、どのデータでも検定できない。

## 反証・不利な文献

| PMID | 雑誌・年 | 内容 | 扱い |
|---|---|---|---|
| 40389522 | Sci Rep 2025（Ali） | AC-T 相当で骨量減少は持続するが、senolytic で防げず、後期には老化マーカーも検出されない | 化療の老化は一過性。経路 1 を唯一の機序にしない |
| 30913313 | J Bone Miner Res 2019 | 若年成体の OVX/ORCH では老化マーカーも SASP も増えない。エストロゲン欠乏と老化は独立の機序 | 「閉経が老化を速める」とは書けない。1st は「加齢」の老化 |
| 38457059 | Stem Cell Rev Rep 2024 | 老化促進マウスの骨髄 MSC の SASP は乳がん細胞の viability を下げ、dormancy を誘導する | 老化間質 = 転移促進、とは限らない |
| 35650435 | Nature 2022 | 加齢の効果は臓器依存（肺は休眠脱出を促し、皮膚は抑える） | 「加齢 = 一律に転移促進」とは書けない |

# 経路 2: 骨形成/血管 niche の喪失 → 休眠維持シグナルの喪失 → 休眠解除

## 流れ

- 休眠 DTC は骨の niche（分化した骨芽細胞、血管周囲）から**静止シグナル**を受けて眠っている（GDF10 / TGF-β2、E-selectin / SDF-1、TSP-1）
- 1st: 加齢で骨の血管 niche（type H 血管、pericyte）と骨形成プログラムが減る
- 2nd: 化療が骨形成・ECM・血管 niche をさらに削る（本データ）
- → 静止シグナルが失われる → 休眠解除 → 骨転移
- 介入点: 破骨を止める BMA では niche は戻らない。この経路を主にすると「なぜ D-CARE が陰性だったか」の説明には使えるが、BMA 提案とは直接繋がらない。

## 本データのどの結果が乗るか

| 結果 | 値 | 乗り方 |
|---|---|---|
| 雌 LepR+ 間質 Dox × 4（GSE296758） | ossification -2.13、osteoblast diff. -1.93、core matrisome -2.50、collagen fibril -2.80、angiogenesis -2.16、sprouting -1.58 | 2nd の niche 低下 ◎ |
| 加齢 LepR+ / CAR と Dox の並置（`docs/prelim_niche.md`） | core matrisome・collagen fibril・EMT は 4 contrast すべてで ↓（NES -1.6〜-2.9）。*Dmp1* ↓ が 3 contrast で共通 | A と B が**同じ向き**で niche を削る ◎（加算の形） |
| HSC 保持の niche 因子 | *Cxcl12*、*Kitl*、*Angpt1*、*Vcam1*、*Lepr* は加齢でも化療でも動かない | **SDF-1（Cxcl12）経由の繋ぎは乗らない** |
| マウス雄加齢皮質骨 | ossification -2.36、matrisome -3.37、angiogenesis -2.35 | 1st の niche 低下 ○（雄） |
| ヒト生検（閉経後・高齢） | ossification +2.10 / +1.47、matrisome +2.70 / +2.58、angiogenesis +1.85 / +1.69 | **逆向き**（組成差の疑い） |
| Dox 単回の骨（組成補正版） | ossification +1.63、matrisome +1.88 | 逆向き |
| 3D 骨（TE-Bone）× Dox | COL1A1 -4.0、core matrisome -2.5（FDR < 0.001） | 2nd の基質低下 ○（細胞毒性の可能性） |

## 橋渡し文献

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 適合性 | 弱点 |
|---|---|---|---|---|---|
| 29514796 | Cancer Res 2018 | 前立腺がん + マウス大腿骨、in vitro + in vivo | **分化した**骨芽細胞の分泌因子（GDF10 / TGF-β2）が TGFβRIII-p38MAPK-RB を介して DTC を静止させる。未分化骨芽細胞では起きない | ◎（骨形成 ↓ → 休眠維持の喪失の分子根拠） | 前立腺がん。乳がんでの検証は無い |
| 27225183 | Sci Transl Med 2016（Price） | マウス異種移植 + 患者骨髄 IHC | 休眠乳がん DTC は E-selectin / SDF-1 に富む perisinusoidal 血管周囲に局在。CXCR4 阻害で休眠細胞が動員される | ○（休眠の場所 = 血管 niche） | 本データで *Cxcl12* は動かない |
| 23728425 | Nat Cell Biol 2013（Ghajar） | マウス、乳がん DTC | 血管周囲 niche の TSP-1 が休眠を維持し、内皮の活性化・sprouting で破れる | ○ | 「sprouting で破れる」は、本データの sprouting ↓ と**向きが逆**になりうる |
| 27074508 | Nature 2016 | マウス骨、加齢 | 加齢で type H 血管・動脈・PDGFRβ+ 血管周囲細胞が強く減る | ◎（1st 側の血管 niche 低下） | 転移は扱っていない |
| 28956197 | Angiogenesis 2017 | 総説（ラット/マウス） | 化療・放射線は骨髄 sinusoidal endothelium を傷害する | △（2nd 側の血管障害） | 総説 |
| 39139448 | Front Cell Dev Biol 2024 | 雄マウス、若齢 vs 成体 × busulfan | 骨髄血管の変化は若齢では回復し、成体では 6 週後も残る | ○（年齢 × 化療の唯一の新しい実験） | 雄、閉経後相当の年齢ではない、移植前処置。データ寄託なし |
| 30664790 | Nat Cell Biol 2019（Carlson） | マウス、乳がん DTC | DTC は血管周囲 niche にいることで化療から守られる。integrin / vWF / VCAM-1 を切ると化療感受性になり骨転移が予防される | △ | 主張は「niche が DTC を守る」で、「化療が niche を壊す」ではない |
| 41465319 | Int J Mol Sci 2025 | 総説 | 骨の dormancy escape: 破骨による TGF-β/IGF-1 の遊離、**骨芽細胞由来の静止シグナルの喪失**、骨髄脂肪細胞の支援 | △（枠組み） | 総説 |
| 27761369 | J Bone Oncol 2016 | 総説 | 骨の血管微小環境は休眠を支え、再活性化も指示する | △ | 総説 |

## 経路として成立するための未検証の飛躍

1. **「骨芽細胞/血管 niche が休眠を維持する」と「化療が niche を傷める」は、それぞれ別の論文にあるが、両者を繋いだ実験は無い**（`docs/literature2.md` (d) 2）。
2. 本データの niche 低下は**ソートに混ざる骨芽細胞の量（組成）**でも同じ形になる。細胞内の転写変化か、集団の割合の変化かは区別できない。
3. 休眠維持因子のうち本データで測れる *Cxcl12* は動かない。GDF10 / TGF-β2 / TSP-1 の個別の値はこの資料では確認していない（要確認）。
4. ヒト生検では向きが逆（↑）。ヒトで支えるデータが無い。
5. 休眠維持の喪失が**閉経後だけ**で効く理由を、この経路だけでは説明しにくい（niche 低下は化療単独でも起きる）。

## 反証・不利な文献

| PMID | 雑誌・年 | 内容 | 扱い |
|---|---|---|---|
| 31292293 | JCI Insight 2019（Kusumbe） | 加齢骨では pericyte が減って増殖性の secretome になるが、放射線・化療は type H 血管を作り替え pericyte を増やし、**静止促進** secretome を作る | **化療が niche を「守る」側**。本経路の 2nd と向きが逆。ただし「A と B を同じ臓器・同じ軸で扱った唯一の論文」でもある |
| 25600338 | Cancer Cell 2015（Wang） | osteogenic niche が播種直後の乳がん細胞の初期増殖を**促す** | 骨形成 niche の喪失は、むしろ定着を減らす方向にも読める |

# 経路 3: ストローマ傷害 + エストロゲン欠乏 → 分泌性老化 → ER+ 休眠細胞の再増殖

## 流れ

- 1st: 閉経（エストロゲン欠乏）
- 2nd: 化療が骨髄間質の ER 応答を消す（本データ: *Pgr*、*Greb1* ↓）＋ 間質を傷害する
- → 間質の分泌性老化（IL-6 / IL-8 / TGF-β1）→ ER+ の休眠細胞が再増殖（ERα 低下・間葉様）
- 本研究の対象（ER+）と直接関係する唯一の経路。

## 本データのどの結果が乗るか

| 結果 | 値 | 乗り方 |
|---|---|---|
| 雌 LepR+ 間質 Dox × 4（GSE296758） | Hallmark estrogen early -1.60（FDR 0.014）、*Pgr* -2.73、*Greb1* -1.43 | 2nd の間質 ER 応答 ↓ ○ |
| 加齢側の比較 | 加齢 CAR は雄、加齢 LepR+ では *Pgr* が発現フィルタで落ちる。雌の中年 CAR では *Esr1* -0.01、*Greb1* +0.58（有意でない） | A 側では ER 応答の低下は**見えない** |
| ヒト閉経後 vs 閉経前 | estrogen early +2.16（直感に反して ↑） | 乗らない（バッチ交絡を疑う） |
| in vitro | ESR1/PGR/GREB1 は培養でほとんど発現しない | 検証できない |

## 橋渡し文献

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 適合性 | 弱点 |
|---|---|---|---|---|---|
| 30119678 | Cell Commun Signal 2018（Tivari） | ヒト MCF-7 休眠モデル + マウス骨髄間質、in vitro | 酸化・低酸素・**エストロゲン除去（フルベストラント）**で間質に分泌性老化を誘導すると、休眠 ER+ 細胞が再増殖。rhIL-6 / IL-8 / TGF-β1 単独でも覚醒。休眠細胞は ERα 低・N-cadherin/SLUG 高 | ◎（本仮説に最も近い） | in vitro のみ。1 本だけ |
| 33878299 | Dev Cell 2021（Bado） | ER+ 乳がん骨微小転移、マウス | 骨の osteogenic niche が**腫瘍細胞側**の ER の発現・活性を一過性・可逆的に下げる（EZH2）。内分泌抵抗性と basal/stem 様形質 | ○ | 間質ではなく腫瘍側の ER |
| 28368409 | Oncogene 2017 | 総説 | ER は乳がんの driver であり骨の恒常性でも中心。ER を標的にする治療は骨微小環境にも作用する | △ | 総説 |
| 24451543 | Nature 2014 | マウス HSC | エストロゲンが ERα を介して HSC の自己複製を増やす | △（骨髄の幹細胞は ER 応答性） | DTC ではない |
| 24687923 / 25603921 | Clin Cancer Res 2014 / Int J Cancer 2015（Ottewell） | OVX マウス | OVX のみで骨内の休眠 DTC が増殖 | ○（1st 側） | 機序は破骨依存で、間質 ER ではない |
| 34733240 | Front Endocrinol 2021（Canuas-Landero） | マウス、エストラジオール濃度の作り分け | ZOL の骨転移抑制は低エストラジオール下でのみ | ○ | ER 応答そのものは見ていない |

## 経路として成立するための未検証の飛躍

1. **「骨髄間質の ERα が DTC の休眠を制御する」を示した実験（conditional KO など）は見つからなかった**（`docs/literature2.md` (d) 4）。本データの「化療後の間質で ER 標的が消える」は**先行文献の無い所見**として書くことになる。
2. 本データの ER 応答 ↓ は、卵巣が保たれた若い雌マウスでの結果。**Dox による卵巣機能低下で血中エストロゲンが下がった結果なのか、間質の細胞内在の変化なのかは区別していない**（未検証）。
3. **臨床との整合に問題がある**: 閉経後の女性はすでにエストロゲンが低いので、「化療がエストロゲン応答を消す」効果は閉経前の方が大きいはず。ところが臨床では閉経前で化療の効果が無い（HR 1.24, n.s.）。化療は閉経前女性で E2 を下げる（PMID 25923354: 閉経前 E2 149 → 7 pg/mL）。この経路を主にするなら、査読の R9（化療性無月経）への答えが一層必要になる。
4. Tivari 2018 は in vitro の 1 本だけで、追試を確認していない。

## 反証・不利な文献

- 直接の反証文献は見つかっていない（先行研究自体がほぼ無い）。
- PMID 25923354（PLoS One 2015）: 術前化療の前後で、骨吸収マーカー ↑・骨形成マーカー ↓ は**閉経状態によらない**（閉経前 β-CTX +36 %、閉経後 +32 %）。化療の骨への作用に閉経後特異性が無いことを示す臨床データで、経路 3 と経路 4 の両方に不利。

# 経路 4: RANKL/OPG 軸（清水原案）

## 流れ

- 1st: 閉経 → 骨髄の RANKL ↑、TNF 産生 T 細胞 ↑、骨細胞アポトーシス
- 2nd: 化療・ステロイド → 骨細胞アポトーシス → 周囲の生存骨細胞が RANKL を出す、ステロイドで OPG ↓
- → RANKL/OPG 比 ↑ → 破骨活性化 → 骨吸収 → DTC 再活性化
- 介入点: デノスマブ / BP が**直接**止める。臨床の提案（化療前の BMA）と最も素直に繋がる。

## 本データのどの結果が乗るか

| 結果 | 値 | 乗り方 |
|---|---|---|
| マウス雄加齢皮質骨 | *Tnfsf11*（RANKL）-1.30（padj 2.8e-09） | **逆向き** |
| ヒト閉経後・高齢 | RANKL は動かない | 乗らない |
| マウス DEX（v1 の 4C） | *Tnfrsf11b*（OPG）-1.05（padj 6.6e-05） | 2nd（ステロイド）◎。ただし DEX は Supplementary に下げる方向 |
| ヒト GIOP | OPG -1.00（padj 0.003）、骨髄混入で補正すると -0.58（padj 0.22） | ○（不安定） |
| LepR+ / CAR の加齢と化療 | *Tnfsf11*・*Tnfrsf11b*・*Adipoq* は動かない | 乗らない |
| in vitro（ヒト骨芽細胞・MSC） | OPG は 4 contrast すべてで ↑（+0.5〜+1.9）、RANKL はほぼ発現なし | **逆向き**（OPG 優位） |
| 雌マウス AC / TAC の骨 qPCR（PMID 40389522 の個体値を再計算） | AC 2 週後 Tnfsf11 × 0.98、TAC 7 週後の椎体 × 0.57 | 乗らない |
| 要約図の Osteoclast support 行 | ほぼ空白 | 乗らない |

## 橋渡し文献

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 適合性 | 弱点 |
|---|---|---|---|---|---|
| 26211824 | Lancet 2015（EBCTCG） | 臨床、18,766 例の IPD メタ解析 | 術後 BP は**閉経後でのみ**骨再発 RR 0.72、乳がん死 0.82 | ◎（閉経後特異性の臨床の鏡像） | BP の効果であって化療の作用ではない |
| 25035292 | Lancet Oncol 2014（AZURE 最終） | 臨床 RCT | 閉経後 5 年超の群でのみ IDFS HR 0.77 | ◎ | 同上 |
| 30795951 | Lancet Oncol 2019（ABCSG-18） | 臨床 RCT、閉経後 HR+、AI 併用 | デノスマブで DFS HR 0.82（多重性未補正の記述的解析） | ○ | AI 併用集団 |
| 24687923 / 25603921 | Clin Cancer Res 2014 / Int J Cancer 2015（Ottewell） | OVX マウス | OVX で骨内の休眠 DTC が破骨依存に増殖、ZOL / OPG-Fc で阻止 | ◎ | 化療は扱っていない |
| 12697741 | J Clin Invest 2003 | ヒト閉経直後の女性 | 骨髄細胞表面の RANKL が上昇、骨吸収と相関 | ◎（1st 側、ヒト） | 本データの生検では RANKL は動かない |
| 21909105 / 21909103 | Nat Med 2011（Nakashima / Xiong） | マウス条件付き KO | 骨細胞が骨の主要な RANKL 供給源 | ○（収束点） | 化療・閉経は扱っていない |
| 22342796 | Bone 2012（Kennedy） | 骨の微小損傷モデル | アポトーシス骨細胞の近傍の生存骨細胞が RANKL を出す | ○（アポトーシス → RANKL） | 微小損傷モデル |
| 21771887 | Endocrinology 2011（Weinstein） | マウス、ステロイド | OPG がステロイド誘発の骨細胞アポトーシスを防ぐ | ○ | ステロイド側のみ |
| 41469373 | Nat Commun 2025 | 雌マウス 12 週、Doxo / PTX | 化療による骨量減少は老化 CAR 細胞・骨髄脂肪細胞の RANKL を介する | ◎（化療 → RANKL の最新の直接証拠） | 若齢。**細胞分画**で見た RANKL で、全骨では出ない可能性 |
| 40102423 | Bone Res 2025 | 雌マウス 3 か月、OVX、Adipoq-CreER × Rankl fl/fl | OVX による骨量減少には骨髄脂肪系前駆細胞（MALP）の RANKL が必要 | ◎（閉経 → RANKL を同じ細胞で） | 41469373 と並べると「閉経と化療が脂肪系 niche 細胞の RANKL で重なる」と言えるが、両者を同時に扱った実験は無い |
| 25575458 | Mol Cell Endocrinol 2015（Oostra） | 臨床、補助化療後 | 骨減少と血中 OPG の関連 | ○ | 閉経状態別は未確認 |

## 経路として成立するための未検証の飛躍

1. **本データの転写レベルでは RANKL/OPG 比は上がらない**（DEX を除く）。図とデータの矛盾を避けるなら、この経路は「既報と臨床の根拠に立つ経路」として書き、本データの役割は「否定はしないが支持もしない」になる。
2. 41469373 / 40102423 の「脂肪系 niche 細胞の RANKL」は細胞分画で出る話。手元の LepR+ / CAR の bulk では *Tnfsf11* は動かなかった。
3. 化療の骨細胞アポトーシス（Wasi 2024 の TUNEL）→ RANKL の連結は、化療では示されていない（Kennedy 2012 は微小損傷）。

## 反証・不利な文献

| PMID | 雑誌・年 | 内容 | 扱い |
|---|---|---|---|
| 31806543 | Lancet Oncol 2020（D-CARE） | 高リスク早期乳がん 4,509 例、標準治療 + デノスマブで骨転移無発生生存は改善せず（HR 0.97） | 最大の反証。閉経前後が混在。骨の健康は改善（PMID 34185259） |
| 42138086 | J Clin Invest 2026（Lin） | ER+ 骨転移の約 40 % はデノスマブ下でも進行（RANKL 非依存の破骨活性化） | RANKL 阻害が唯一の手ではない |
| 25923354 | PLoS One 2015 | 化療前後の骨代謝マーカーの変化は閉経状態によらない | 経路 3 と同じ |

# 経路 5: その他

## 5a. IFN 応答（本データで化療と加齢に共通の唯一の軸、ただし向きが逆）

- 本データ: IFN-α / IFN-γ 応答 ↑ は bulk の Dox 4 データすべてと急性期 CAR で同じ向き（NES +1.3〜+2.5）。加齢でも ↑（+2.6 / +2.7）。OVX では ↓、DEX では動かない。
- **文献上、I 型 IFN はほぼすべて「転移を抑える」側**。この軸を 2-hit の支持には使えず、「宿主の防御反応で、他の軸を相殺しうる」と書くのが正確。

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 適合性 | 弱点 |
|---|---|---|---|---|---|
| 30546090 | Oncogene 2019 | ER- 乳がん、Dox / MTX 後 | 化療を生き延びた細胞で IRF7 / IFN-β / IFNAR 軸が持続し、**免疫学的 dormancy を作る**。IFNAR 阻害で解除 | △（向きが逆） | 腫瘍細胞側 |
| 26198985 | Cancer Immunol Res 2015 | Ifnar1-/- マウス | 宿主の I 型 IFN の**喪失**が骨転移を加速 | △（向きが逆） | — |
| 22820642 | Nat Med 2012 | 乳がん骨転移 | 腫瘍側の Irf7 経路の沈黙が骨転移を促す | △（向きが逆） | 腫瘍細胞側 |
| 11117749 / 11961557 | Nature 2000 / 2002 | マウス骨 | IFN-γ・IFN-β はいずれも破骨細胞分化を抑える | △（向きが逆） | — |
| 37512126 | Medicina 2023 | マウス局所照射 + MLO-Y4 骨細胞 | 照射が STING-TBK1 を介して骨微小環境を壊し老化を起こし、転移に好都合な環境を作る | ○ | 放射線で、化療ではない |

## 5b. 化療の宿主反応・NET による休眠解除（臓器は主に肺・肝）

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 適合性 | 弱点 |
|---|---|---|---|---|---|
| 38161426 | iScience 2024 | 腫瘍の無い雌 BALB/c を AC で前処置 → NRG マウス + E2 ペレットに ER+ ZR-75-1 を静注 | 化療で前処置した宿主では ER+ 細胞の肝・肺・**骨**転移が増える（骨は 2 切片 / 匹、n = 3 / 群） | ○（「化療 → 骨転移 ↑」を示した唯一の近傍） | 若齢・免疫不全、骨だけでなく肺・肝も増える（臨床の骨特異性と合わない）、n が小さい |
| 40614736 | Cancer Cell 2025（He） | マウス、DormTracer | Dox・シスプラチンが休眠乳がん DTC を目覚めさせ肺転移を増やす（NET 経由） | ○ | 肺 |
| 42696464 | Cancer Res 2026（Shi） | マウス、膵がん | 化療が老化肝細胞を介した NET で肝の休眠 DTC を再活性化 | △ | 肝・膵がん |
| 30262472 | Science 2018（Albrengues） | マウス、肺 | 持続炎症 → NET → laminin 切断 → 休眠細胞の覚醒 | △（A の炎症の出口として） | 肺 |
| 37699010 | PLoS Biol 2023 | ドセタキセル、腫瘍間質オルガノイド・同系乳がんモデル | タキサンが間質を傷害 → IL-6 / G-CSF → 休眠がんが目覚める | ○ | 乳腺で、骨ではない |
| 21975929 | Cancer Res 2011（Daenen） | マウス | 化療で前処置した宿主（薬は消失済み）で転移が増える | ○（時間差の説明） | 肺 |

## 5c. ステロイドの腫瘍への直接作用

- PMID 30867597（Nature 2019、Obradovic）: 遠隔転移巣でグルココルチコイド受容体が活性化するとコロニー形成が増え生存が短縮（ROR1 経由）。臨床でのステロイド使用に注意を促す。
- DEX を本編から外す場合でも、Discussion で「支持療法ステロイドは骨の niche と腫瘍の両方に作用しうる co-factor」として残せる。臨床データに累積ステロイド量は無い（と思われる）ので、トーンは抑える。

## 5d. 血管 pericyte（Kusumbe 2019）

- 経路 2 の反証欄に記載。化療が骨の pericyte を増やし**静止促進** secretome を作る、という報告で、本経路群とは**逆向き**。ただし「加齢では静止維持能が落ちる」の半分は経路 2 の 1st を支える。

# 経路を組み合わせる場合

単独の経路ではどれも 1st と 2nd の両方を支えきれないので、2 本を組み合わせる形も考えられます。`docs/literature2.md` の Discussion 案は **経路 1（1st の炎症・老化が破骨方向の信号を出す）＋ 経路 2（2nd が休眠維持の信号を取り去る）** の組み合わせです。

| 組み合わせ | 1st の役割 | 2nd の役割 | 長所 | 短所 |
|---|---|---|---|---|
| 1 + 2（literature2 の案） | 炎症・老化 → 破骨方向の信号を**足す** | 骨形成・血管 niche → 休眠維持の信号を**取り去る** | 「A と B は別のことを起こす」という dry の到達点にそのまま合う | 経路 2 の文献の繋ぎが弱い。Kusumbe 2019 が逆向き |
| 1 + 4 | 老化・SASP | RANKL/OPG（文献） | 出口が破骨に一本化され、BMA の提案に直結 | 2nd の RANKL は本データで出ない |
| 1 のみ（A も B も老化） | 加齢の老化 | 化療誘発老化 | 単純。in vitro が陽性対照として効く | 生体マウスで化療の老化が出ない。Ali 2025 の反証 |
| 4 のみ（清水原案） | RANKL ↑ | アポトーシス → RANKL | 臨床と最も整合 | 図とデータが矛盾する |

# 判断の観点

**どの経路を主にするかは、清水先生と ken で判断します。** 判断の観点として次を挙げます。

1. **橋渡し文献の質**: 骨・乳がん・in vivo の原著があるか（経路 1 は Luo 2016、経路 4 は Ottewell・EBCTCG がある。経路 2 は前立腺がんの原著と総説、経路 3 は in vitro 1 本）。
2. **本データとの整合**: 図（Fig4）とデータが矛盾しないか。経路 4 は本データで出ない。経路 2 はマウスでは出るがヒトでは逆。経路 1 はヒトの 1st と in vitro の 2nd で出る。
3. **臨床との整合**: 閉経後でだけ、骨でだけ、という臨床所見を説明できるか。BMA（デノスマブ / BP）の提案に繋がるか。経路 4 > 経路 1 > 経路 2・3。
4. **査読者の予想反応**: D-CARE（経路 4・1）、化療の老化は一過性（Ali 2025、経路 1）、niche の向きが逆（Kusumbe 2019、経路 2）、先行研究なし（経路 3）、肺・肝の話を骨に持ち込む（経路 5）。
5. **Fig5 の描きやすさ**: 経路 1・4 は「1st + 2nd → 破骨 → DTC」と一本で描ける。経路 2 は「信号が消える」を図にする必要がある。組み合わせ（1 + 2）は 2 本の矢印を並べる形になる。
6. **どの Figure 4 構成案と組むか**: 案 A（ヒト中心）は経路 1、案 B（マウス生体中心）は経路 2、案 C はどれでも（`fig_options` の章を参照）。

# Discussion 段落案（英文、経路別）

`docs/literature2.md` (c) の段落を経路別に書き分けました。どの版も最後に「交互作用は示していない」ことを明記しています。

## 版 1: 老化・SASP → 破骨（経路 1、＋経路 3 を 1 文）

> Our transcriptomic survey indicates that postmenopausal, ageing bone is characterised by senescence-associated programmes (SenMayo and p53, together with increased CDKN1A), and that a topoisomerase II inhibitor induces the same programmes, with a pronounced senescence-associated secretory phenotype (SASP), in human osteoblasts in vitro. Senescence of bone-resident cells has repeatedly been linked to osteoclastic bone resorption and to metastatic outgrowth. Senescent osteoblasts alone are sufficient to increase breast cancer colonisation of bone through IL-6-driven osteoclastogenesis (PMID 26725121); clearance of senescent cells in aged mice reduces bone resorption (PMID 28825716); chemotherapy-induced bone loss is driven by senescent cells independently of oestrogen (PMID 31932453); and osteocytes in bone metastases acquire a pro-osteoclastogenic SASP (PMID 39312185). Stromal secretory senescence, including that induced by oestrogen withdrawal, also reawakens dormant ER+ breast cancer cells in vitro (PMID 30119678). A plausible synthesis is that menopause and chemotherapy each add to a senescent, pro-resorptive bone environment in which osteoclast activity releases dormant tumour cells (PMID 24687923, 26632274). We note, however, that senescence signatures were not detected in doxorubicin-exposed mouse bone weeks after treatment, consistent with reports that chemotherapy-induced skeletal senescence is transient (PMID 40389522). No study has examined both insults together in bone, and we present this as a hypothesis rather than an interaction we have demonstrated.

## 版 2: 骨形成/血管 niche の喪失（経路 2）

> Our transcriptomic survey argues that menopause/ageing and chemotherapy do not converge on a single bone programme: ageing bone is dominated by inflammation and senescence, whereas doxorubicin-exposed murine bone-marrow stroma loses its osteogenic, matrix and angiogenic programmes without inflammatory induction. Dormancy in bone is actively imposed by the very compartments that chemotherapy erodes. Differentiated osteoblasts secrete GDF10 and TGF-beta2 that enforce tumour-cell quiescence (PMID 29514796), dormant breast cancer cells reside in E-selectin- and SDF-1-rich perisinusoidal niches (PMID 27225183), and endothelial thrombospondin-1 sustains quiescence (PMID 23728425). Bone vascular niches decline with age (PMID 27074508), and chemotherapy injures marrow sinusoidal endothelium (PMID 28956197). A plausible synthesis is that chemotherapy removes osteoblastic and vascular signals that keep disseminated tumour cells quiescent, in a bone that has already lost part of these signals with age. Two caveats apply. First, therapy has also been reported to expand bone pericytes and to create a quiescence-promoting niche (PMID 31292293), the opposite of our interpretation. Second, the osteogenic and vascular programmes moved in the opposite direction in human bone biopsies, possibly reflecting tissue composition. No study has linked chemotherapy-induced loss of these niches to dormancy exit, and we therefore present this as a hypothesis.

## 版 3: RANKL/OPG 軸を臨床の根拠に置く（経路 4、本データは中立と書く）

> The strongest clinical and experimental evidence connecting the bone microenvironment to dormancy exit concerns osteoclast activity. Ovariectomy alone induces the outgrowth of dormant tumour cells in bone, but not elsewhere, and this is prevented by zoledronic acid or RANKL blockade (PMID 24687923, 25603921); adjuvant bisphosphonates reduce bone recurrence only in postmenopausal women (PMID 26211824, 25035292). Recent work shows that bone loss after oestrogen deficiency requires RANKL from marrow adipogenic progenitors (PMID 40102423), and that chemotherapy-induced bone loss is mediated by RANKL from senescent CXCL12-abundant reticular cells and marrow adipocytes (PMID 41469373), suggesting a cellular point at which the two insults could converge. Our bulk transcriptomes, however, did not show an increase in RANKL or in the RANKL/OPG ratio after ageing or doxorubicin, and osteoprotegerin increased in osteoblasts exposed to DNA damage in vitro; only glucocorticoid exposure reduced osteoprotegerin. Convergence on RANKL may therefore be restricted to specific stromal subsets that bulk profiling cannot resolve. Denosumab given with standard therapy did not improve bone metastasis-free survival in a mixed pre- and postmenopausal population (PMID 31806543), and a substantial fraction of ER+ bone metastases progress through RANKL-independent osteoclast activation (PMID 42138086). We therefore propose a targeted trial in postmenopausal women starting before chemotherapy as a hypothesis, not as a recommendation.

# 清水先生に確認したいこと

1. 主にする経路（または組み合わせ）をどれにするか。候補として、dry の到達点に最も合うのは **1 + 2**、臨床の提案に最も合うのは **1 + 4** です。
2. 経路 4（RANKL/OPG）を主にしない場合でも、BMA の提案（化療前のデノスマブ / BP）を Discussion の結論に残すか。残すなら「下流の破骨を止める」という位置づけでよいか。
3. 経路 3（間質の ER 応答の消失）は先行研究がほぼ無い所見です。新規性として前面に出すか、Discussion の 1 文に留めるか。
4. 経路 5b の iScience 2024（PMID 38161426、化療で前処置した宿主で ER+ 細胞の骨転移が増える）を引用するか。骨特異性が無い点は臨床所見と合いません。
5. 各経路の橋渡し文献のうち、本文を取り寄せて確認すべきもの（例: Luo 2016、Tivari 2018、PMID 29514796、Nat Commun 2025）の優先順位。現状はすべて抄録だけに基づいています。

---
title: "第 6 章 Figure 4 / Figure 5 の構成オプション（v2）"
date: 2026-10-26
author: 遠田
---

10/14 のご依頼との関係: ご依頼の構成（4A = 加齢の 4 案、4B = *Blood* 2026 の GSE296761、4C = 骨芽細胞の GSE304562 または現行の DEX）は、下の 3 案では**案 B（マウス生体中心、4B を GSE296758 に差し替え）に最も近く**、4A 案 1（GSE141595）は案 A の 4A と同じデータです。
一方、3 案はどれも 4C を要約モジュール図にし DEX を Supplementary に下げていますが、ご依頼では 4C を骨芽細胞（案 1）か現行の DEX（案 2）にするので、4C についてはご依頼の方を背骨にしました（第 3 章）。以下は 10/16 時点の検討として、中身を変えずに残します。v3 の解析結果を受けた構成の案は、第 1–3 章の「まとめと推奨（v3）」にあります。

（この章について: v2 では下の「この章について」の節は PDF に載っていませんでした。v2 の原稿（`parts2/fig_options.md`）の前書きで、内容は変えていません。）


# この章について

10/13 の打ち合わせで、Figure 4 は「4C（デキサメタゾン）を本編から外し、閉経/加齢（4A）と化療（4B）の 2 本を主軸にする」方向になりました。
その後、ヒト骨芽細胞の in vitro データ（第三波）と、全データを同じ遺伝子セットで並べた要約図（モジュール図）を作りました。
これを踏まえ、Figure 4 / Figure 5 の組み方を 3 案に整理しました。**どれかに決め打ちはせず、長所・短所・査読リスク・追加作業を並べています。**
Figure 5（概念図）の中身は、別章「仮説の経路オプション」でどの経路を主にするかによって変わります。この章では、各案で Figure 5 と GA の**どこが**変わるかまでを書きます。

前提となる dry の到達点（詳細は v1 資料と `docs/prelim_*.md`）:

- **閉経/加齢（1st insult 側）**: ヒト閉経後・高齢の骨生検で老化（SenMayo、p21 = CDKN1A）と p53 が上がる。マウス雄の加齢皮質骨では炎症（Tnf, Il1b, S100a8/9）が上がる。
- **化療（2nd insult 側）**: 生体マウスでは炎症・老化・アポトーシスは上がらない。雌 LepR+ 間質（Dox × 4）で骨形成・ECM・血管の低下と ER 応答の低下。ヒト骨芽細胞に DNA 傷害（エトポシド）を与えると p53・アポトーシス・老化・SASP がそろって上がる。
- **つまり、閉経/加齢と化療は「同じこと」ではなく「別のこと」を起こしている**。2 つの掛け合わせ（交互作用）を直接見たデータは公開されていない（第一〜三波の調査で見つからず）。

# 3 案の概要

| | 案 A: hypothesis-generating（ヒト中心） | 案 B: 清水案①を維持（マウス生体中心） | 案 C: dry を Supplementary へ |
|---|---|---|---|
| 4A | ヒト閉経後 vs 閉経前の骨生検（GSE141595）。高齢 vs 若年（GSE72815）を再現として添える | マウス雄の加齢皮質骨 2.5 歳 vs 2 か月（Kaya 2022、v1 の 4A） | なし（S に移す） |
| 4B | ヒト骨芽細胞 × エトポシド in vitro（GSE284666） | マウス Dox（v1 の GSE205237 単回、または差し替え候補 GSE296758 の LepR+ 間質 Dox × 4） | なし（S に移す） |
| 4C | 要約モジュール図（マウス生体は対照の列として入る） | 要約モジュール図 | — |
| 本編の Fig4 | 上の 3 パネル | 上の 3 パネル | 概念図（現 Fig5）を Fig4 に繰り上げ |
| volcano | Supplementary | 4A・4B は本編、他は Supplementary | すべて Supplementary |
| DEX（v1 の 4C） | Supplementary | Supplementary | Supplementary |
| 論文中の位置づけ | "hypothesis-generating" と明記 | 清水仮説の検証として提示（出なかった結果も書く） | 「公開データの探索的解析」として補足扱い |

要約モジュール図（下書き、`out/prelim/summary/summary_modules.png`）:

![要約モジュール図（下書き）。FDR < 0.05 のみ色付き、FDR 0.05-0.25 は枠線、FDR ≥ 0.25 は空白](prelim/summary/summary_modules.png){width=100%}

この図は**どの案でも 4C の候補**です。ただし下書きの自己評価として、現在の列の並びでは「2-hit」は図から読めません（後述）。

# 案 A: hypothesis-generating（ヒト中心）

## 構成

- **4A** ヒト女性の骨生検: 閉経後 vs 閉経前（GSE141595、15 vs 15）。同じ方向を高齢 vs 若年（GSE72815、19 vs 19、同一バッチ）で再現（log2FC の Spearman ρ = 0.56）。
- **4B** ヒト骨芽細胞（NHO）× エトポシド（Dox と同じ TOP2 阻害）in vitro（GSE284666）。p53・アポトーシス・老化・SASP の陽性対照。
- **4C** 要約モジュール図。マウス生体（加齢・Dox）は「ヒトの所見がマウス生体では同じ形に出ない」ことを示す対照の列として入る。
- volcano（v1 の 4A/4B/4C、ヒト GIOP）は Supplementary。

## 長所

- **ヒトのデータが主役**になる。臨床コホート（Fig1–3）との距離が最も近い。
- 4A は老化（SenMayo NES +1.72〜1.76, FDR ≤ 0.003; 高齢側 +1.65, FDR 0.011）と p53（+2.14〜2.17）が出ており、**「閉経後の骨は老化側にある」が陽性の結果として言える**。
- 4B は「化療で骨系細胞が p53 → 老化 → SASP に入る」が教科書どおり強く出る（4 contrast すべてで p53 NES 1.9–2.7、SenMayo 1.8–2.7、FDR < 0.01）。v1 で問題になった「化療でアポトーシスが出ない」を、**生体と細胞で何が違うか（時間・希釈）**という形で説明できる。
- 「hypothesis-generating」と最初から宣言するので、交互作用を示していないことを査読で突かれても立場が一貫する。

## 短所

- **4A の GSE141595 は閉経前 15 例が別バッチ**（2022 年追加・別処理・ライブラリ深度 4,200 万 vs 7,300 万）で、群とバッチが完全に交絡している。さらに閉経前のファイルが途中で切れており、遺伝子の約 4 割（chr12 後半〜22, Y, M）が比較できない。単独では根拠にならない。
  （v3 注: この「ライブラリ深度 4,200 万 vs 7,300 万」は、閉経前のファイルが途中で切れていたことによる見かけでした。共通の遺伝子で数えると閉経前 4,140 万、閉経後 4,230 万で同じです（SRA のリード数も 6,050 万 vs 5,760 万で同じ）。別バッチであること自体は、FASTQ の形式の違いから v3 でも確かめています（第 1 章 案 1）。）
- 4A の比較は 25–37 歳 vs 60–79 歳で、**閉経単独ではなく「閉経後の加齢骨」**。高齢 vs 若年（30 vs 73 歳）も同じ。
- ヒト生検の炎症は一貫しない（閉経後では Inflammatory Response +1.81 だが、高齢 vs 若年では -1.23 で有意でない）。強く言えるのは老化で、炎症は TNF・CCL2 の遺伝子レベルにとどまる。
- ヒト生検とマウス加齢は遺伝子全体で**逆相関**（ρ = -0.30 / -0.33）。骨形成・ECM・血管の行は、ヒトで ↑、マウスで ↓ に割れる（組織と組成の違いが大きい）。
- 4B の薬剤は**エトポシド（乳がんの周術期レジメンではない）**、20–50 µM・8 日連続は老化誘導の条件で臨床より強い。1 ドナー由来と思われる細胞の培養反復。ER 系（ESR1/PGR/GREB1）は培養で発現しない。RANKL も発現しない。
- 要約図の自己評価: 列を分けているのは insult の種類より「種・試料」で、1st / 2nd の括弧と色の並びが対応しない（下記「4C の作り直し案」）。

## 査読リスク

| 想定される指摘 | 重さ | 応答の方向 |
|---|---|---|
| 4A の閉経前は別バッチで交絡している | 高 | 高齢 vs 若年（同一バッチ）で同じ方向が出る（ρ = 0.56）ことを主に置く。SRA から再アライメントして欠落遺伝子を回復（バッチ交絡は残る） |
| 4B は乳がん化療ではなく、濃度も臨床と違う | 中 | TOP2 阻害で Dox と機序が同じことを根拠に「陽性対照」と位置づける。3D 骨（TE-Bone）× Dox の結果（p53 標的 ↑、基質・骨形成 ↓）を Supplementary に添える |
| in vitro で出たものがなぜ生体マウスで出ないのか | 中 | 採取時期（30 日〜8 週後）、組織での希釈、死細胞の除去と再増殖（マウス MSPC の E2F ↑）で説明。in vivo の組織レベルの既報（Wasi 2024 の骨細胞 TUNEL 個体値）を Supplementary で示す案もある（`docs/datasets3_lit.md` 案 E） |
| 交互作用を示していない | 高（どの案でも） | 最初から hypothesis-generating と宣言し、Limitation に書く |

## 必要な追加作業

1. **GSE141595 の再アライメント**: SRA（SRP234994）から 45 例をそろえて処理し直す。閉経前ファイルの欠落遺伝子を回復する（バッチ交絡は残る）。
2. 4B のパネル化: `heatmap_vs_mouse` を絞った版（p53 / apoptosis / SenMayo / TNF-α / inflammatory / E2F / osteoblast diff. の 7 行 × ヒト列 + マウス bulk 列）が候補。
3. 4C の作り直し（下記）。
4. volcano を Supplementary 用に組み直し、図番号を振り直す（v1 の S1–S3 の対応も更新）。
5. （任意）Wasi 2024・Nat Commun 2025（PMID 41469373）の公開個体値の再プロット。図に再利用するならライセンス確認が要る。

## Limitation の書き方（英文案）

> Our transcriptomic analyses are hypothesis-generating. The comparison of postmenopausal with premenopausal bone biopsies is confounded by batch, because the premenopausal samples were processed separately; we therefore present the age-matched-batch comparison of older with younger women as the primary replication. Both comparisons contrast postmenopausal, older women with younger women and cannot separate menopause from ageing. The chemotherapy arm relies on human osteoblasts exposed to etoposide in vitro, a topoisomerase II inhibitor that shares its target with doxorubicin but is not part of perioperative breast cancer regimens, at concentrations that exceed clinical exposure. No public dataset combines menopause or ageing with chemotherapy in bone, and the interaction between the two insults was not tested.

# 案 B: 清水案①を維持（マウス生体中心）

## 構成

- **4A** マウス雄の加齢皮質骨（2.5 歳 vs 2 か月、Kaya 2022、4 vs 3）。v1 の 4A のまま。
- **4B** マウス Dox。v1 は GSE205237（雄 15 週、Dox 10 mg/kg 単回、30 日後、4 vs 4、DEG 9）。10/13 の打ち合わせで、より適切に動くデータを再探索することになった。
  現時点の差し替え候補は **GSE296758（雌 LepR+ 間質、Dox × 4、8 週後）**: 骨形成（ossification -2.13）、ECM（core matrisome -2.50）、血管（angiogenesis -2.16）、ER 応答（estrogen early -1.60, FDR 0.014; *Pgr* -2.73）が下がる。
  第三波（`docs/datasets3_lit.md`）でも「加齢/OVX × 化療」の骨のトランスクリプトームは GEO 以外を含めて見つからなかった。
- **4C** 要約モジュール図。
- DEX は Supplementary。

## 長所

- **清水先生の当初の構成と 10/13 の合意に最も近い**（4A・4B を主軸、DEX を外す）。
- 生体（in vivo）の骨のデータで揃う。
- 4A は炎症（Inflammatory Response +2.11、IL-6/JAK/STAT3 +2.40）と破骨支持（GO positive regulation of osteoclast differentiation +1.60, FDR 0.014）が出る唯一のデータで、Fig5 の「出口（破骨）」に直接つながる。
- 4B を GSE296758 にすれば、雌で「化療は炎症を足さずに骨形成・血管 niche を削る」が volcano / モジュールとして見える。加齢 LepR+ / CAR と同じ向き（ECM・骨形成 ↓）が 4 contrast すべてで再現する（`docs/prelim_niche.md`）。

## 短所

- **4A が雄**。「閉経」の代理にならず、「加齢」としか書けない。雌の加齢データ（GSE210584 の中年 CAR、TMS の雌の骨）で置き換えを試したが、4A のストーリー（炎症・老化 ↑、骨形成 ↓）とは骨形成の向きが合わず、TMS 雌は DEG 7 個で volcano にならなかった（`docs/prelim_aging_f.md`）。
- 4B を v1 の GSE205237 のままにすると**ほぼ空の volcano**（DEG 9 個、注目遺伝子で有意なものなし）。DEG の大半は組織組成（筋・血液）の差。
- 4B を GSE296758 にすると、**Dox × 4 vs PBS が flowcell と完全交絡**している。性・組織・時期が 4A と揃わない（4A は雄の皮質骨、4B は雌のソート細胞）。
- 加齢と化療が共通に下げる骨形成・ECM は、**ソートに混ざる骨芽細胞の量（組成）でも同じ形になる**。
- マウス加齢とヒト生検は遺伝子全体で逆相関（ρ = -0.30）。「マウスとヒトで同じシグネチャ」とは書けない。

## 査読リスク

| 想定される指摘 | 重さ | 応答の方向 |
|---|---|---|
| 閉経の話なのに 4A が雄 | 高 | 「加齢」と書き、閉経の代理としての限界を Limitation に。OVX（雌）では転写変化がほぼ無かったことも書く |
| 4B の Dox で何も起きていない（GSE205237） | 高 | GSE296758 に差し替えるか、「化療単独では加齢型の炎症 niche は再現されない」を負の結果として明示する |
| 4A と 4B で性・組織・時期がばらばら | 中 | 同じ遺伝子セットで並べた 4C で比較の単位をそろえる。交互作用は検定していないと明記 |
| マウスの正常骨から転移を語るのは飛躍 | 中（どの案でも） | hypothesis-generating と明記 |

## 必要な追加作業

1. **4B の差し替え先の決定**: GSE296758（Dox × 4 LepR+）で volcano を作る（解析は済、図は未作成）。他に「閉経後（加齢）+ 化療」で動いている公開データは、第三波までの調査では見つかっていない。
2. 4A の legend・Limitation に雄であることを明記。
3. 4C の作り直し（下記）。
4. Supplementary の再編（DEX、ヒト GIOP、雌の加齢 CAR、in vitro）。

## Limitation の書き方（英文案）

> The ageing dataset was obtained from male mice and therefore models skeletal ageing rather than menopause; ovariectomy alone produced few transcriptional changes in bone. The chemotherapy datasets differ from the ageing dataset in sex, tissue compartment and time after exposure, and the doxorubicin-exposed stromal dataset is confounded with sequencing run. Shared decreases in bone-formation and matrix programmes may partly reflect the proportion of osteoblast-lineage cells in the samples. No dataset combines ageing or oestrogen deficiency with chemotherapy, and the interaction was not tested.

# 案 C: dry を Supplementary に下げる

## 構成

- 本編は**臨床の Fig1–3（study design、KM）＋ 概念図**。概念図は現 Fig5 を Fig4 に繰り上げる。
- 公開データの解析（加齢・Dox・DEX・ヒト生検・in vitro・要約図）はすべて Supplementary。
- 本文の Results では dry を 1 段落にまとめ、概念図の各要素に「Supplementary Fig. Sx で整合的な変化を確認した」と添える。

## 長所

- **査読で dry を突かれても論文の主張（臨床の関連）が揺らがない**。臨床結果が主役であることが明確になる。
- 追加作業が最も少ない。
- 「出なかった結果」（化療でアポトーシスが出ない、RANKL が動かない）も Supplementary に正直に置ける。

## 短所

- 機序の裏付けが弱く見え、**「なぜ閉経後でだけ、骨でだけか」に本編で答えられない**。概念図が文献だけで支えられることになる。
- 機序の図が本編に無いと、査読で「機序の裏付けが無い」と評価される可能性がある（MedComm の方針・査読傾向は未確認）。
- ご依頼資料の構成（Fig4 = dry）から最も遠い。

## 査読リスク

| 想定される指摘 | 重さ | 応答の方向 |
|---|---|---|
| 概念図の根拠が文献だけ | 中 | Supplementary の dry を本文で引く。必要ならリバイスで本編に上げる |
| 機序が無いので臨床の関連が交絡では | 中 | 臓器特異性（骨のみ、肺・肝は無関連）と閉経後特異性を主な反論にする |

## 必要な追加作業

1. 図番号の振り直し（Fig5 → Fig4、v1 の Fig4/S1–S3 → Supplementary）。
2. 概念図の文言を、Supplementary に置いた結果だけで支えられる範囲に弱める。

## Limitation の書き方（英文案）

> The proposed model is supported by published experimental work and by an exploratory analysis of public bone transcriptomes (Supplementary Figures), none of which combined menopause or ageing with chemotherapy or included metastatic tumour cells. The model should therefore be regarded as a hypothesis for prospective and experimental testing.

# 3 案の比較（まとめ）

| 観点 | 案 A | 案 B | 案 C |
|---|---|---|---|
| 清水先生の当初案・10/13 の合意との近さ | △（構成は変わるが「4A・4B 主軸、DEX を外す」は守る） | ◎ | △ |
| 陽性の結果が本編に出るか | ◎（ヒトの老化、in vitro の p53/老化） | ○（4A は出る、4B は差し替え次第） | — |
| 臨床（閉経後女性）との距離 | ◎（ヒト女性） | △（雄マウス・若齢マウス） | — |
| 交絡・データの弱点 | バッチ交絡（4A）、エトポシド・in vitro（4B） | 雄（4A）、空の volcano または flowcell 交絡（4B） | 本編には無い |
| 査読で dry を突かれる度合い | 中 | 中〜高 | 低 |
| 追加作業 | 多（再アライメント、4B パネル、4C 作り直し） | 中（4B 差し替え、4C 作り直し） | 少 |
| Fig5 の描きやすさ | 経路 1（老化・SASP）と相性が良い | 経路 2（niche 喪失）と相性が良い | どの経路でも描けるが根拠は文献のみ |

# 4C（要約モジュール図）の作り直し案

どの案でも 4C に使う場合、下書きのままでは主張が伝わりません（`docs/prelim_summary_fig.md` の自己評価）。

- **列を分けているのは insult の種類ではなく種・試料**。骨形成・ECM・血管の 3 行で、ヒト 3 列が赤、マウスの加齢骨と MSPC Dox が青に割れる。読者が最初に受け取るのは「ヒトとマウスで niche 軸の向きが逆」になる。
- 老化（SenMayo, p53）は「ヒト ↑、マウス in vivo は出ない」で、これが図から一目で読める唯一の対比。
- Osteoclast support の行はほぼ空白（有意はマウス雄加齢の +1.60 のみ）。2-hit の出口（破骨）を図が支えない。

提案（未実施）:

1. 列を**「ヒト / マウス」で先に分け**、その中で加齢と化療を並べる。マウス内では「加齢と MSPC Dox が niche 軸 ↓ でそろう、炎症は加齢だけ」が読める。in vitro は別枠（参考）にする。
2. IFN-α を **IFN-γ に替える**と、化療 3 列（MSPC +1.75 / Dox 骨 +1.91 / NHO +2.25）とマウス加齢・閉経後がそろって ↑ になる（ヒト高齢 vs 若年は ↑ にならない）。ただし IFN ↑ は文献上「転移を抑える」側に働く報告がほとんどで、仮説の支持には使いにくい（別章の経路 5）。
3. 行を 6 本（Inflammation / IFN / SenMayo / Bone formation / ECM / Angiogenesis）に絞り、Osteoclast support は「出ない」ことを本文で書いて図から外す。

# Figure 5 と GA は各案でどう変わるか

## 現状（v1）

![Figure 5 修正版（v1、`design/fig5_revised.png`）](v1/design/fig5_revised.png){width=85%}

v1 の修正版は 2nd hit を **RANKL/OPG 比 ↑（OPG ↓、ステロイド）、炎症因子 ↑（ステロイド）、骨形成 ↓（化療）** としていました。
**このうち RANKL/OPG と炎症の 2 項目は DEX（v1 の 4C）が根拠**です。DEX を Supplementary に下げると、2nd hit の本編の根拠は「骨形成 ↓」だけになります。
したがって、**どの案を取っても Fig5 の 2nd hit の中身は書き直しが必要**です。清水原案版（2nd hit = 骨細胞アポトーシス）は、生体の転写データでは引き続き支持されません（in vitro では出る）。

## 案ごとの Fig5 の中身（たたき台）

| | 案 A | 案 B | 案 C |
|---|---|---|---|
| 1st insult の見出し | Postmenopausal (aged) bone | Aged bone | Postmenopausal bone |
| 1st insult の中身 | ↑ Cellular senescence (p21, SenMayo), ↑ p53 | ↑ Inflammatory cytokines (TNF, IL-1β, S100A8/A9), ↓ bone formation | 文献ベース（RANKL ↑、炎症、骨細胞アポトーシス: PMID 12697741, 11086024, 9284757） |
| 2nd insult の見出し | Perioperative chemotherapy | Perioperative chemotherapy | Perioperative chemotherapy (+ steroids) |
| 2nd insult の中身 | ↑ p53 / senescence / SASP in osteoblast lineage (in vitro) | ↓ Osteogenic, matrix and vascular niche, ↓ stromal oestrogen response | 文献ベース（骨細胞アポトーシス・老化: PMID 39147353, 31932453） |
| Data 行 | Fig. 4A / Fig. 4B, 4C | Fig. 4A / Fig. 4B, 4C | Supplementary のみ |
| 右カラム（niche → 骨転移） | 別章の経路の選択で決まる（経路 1 なら SASP → 破骨 → DTC 再活性化） | 別章の経路の選択で決まる（経路 2 なら休眠維持シグナルの喪失 → DTC 再活性化） | 同左 |
| ステロイド | 図から外すか、Supplementary を根拠に小さく残す | 同左 | 文献ベースで残せる |

- どの案でも、**「1st と 2nd が同じ経路に収束する」とは描かず、「別々の変化が、それぞれ独立に DTC の休眠解除に結びつく（と既報が示す）」と描く**のがデータに忠実です（`docs/literature2.md` の Discussion 案と同じ立場）。
- "two-hit" の言い換え（"sequential two-insult" / "priming-and-trigger"）は v1 の確認事項のままです。ただし今回の dry では 1st が「下地を作り」2nd が「引き金を引く」という順序の証拠は無いので、**"priming-and-trigger" より "sequential two-insult" の方が事実に近い**と考えます。

## GA の文言

| | 案 A | 案 B | 案 C |
|---|---|---|---|
| 中央の Hit 1 | Aged bone: senescence | Aged bone: inflammation | Postmenopausal bone |
| 中央の Hit 2 | Chemo: osteoblast senescence | Chemo: loss of osteogenic niche | Chemo (+ steroid) |
| 中央の出口 | 経路の選択に合わせる（例: → DTC reactivation → bone metastasis） | 同左 | 同左 |
| 左（KM・HR） | 変わらず（HR は Table 3 の確定値待ち） | 同左 | 同左 |
| 右（Proposal） | "BMA before chemo" は D-CARE を踏まえ "hypothesis" 寄りに | 同左 | 同左 |

v1 の GA の Hit 2 "Chemo + steroid: RANKL/OPG ↑" は DEX が根拠なので、DEX を本編から外すなら変更が必要です。

# 清水先生に確認したいこと

1. Figure 4 の構成は **案 A / B / C** のどれで進めるか（または組み合わせ。例: 案 B の構成で 4B を GSE296758 に差し替え、案 A の in vitro を Supplementary に置く）。
2. 案 A の場合、GSE141595 の再アライメント（SRA から 45 例）まで行うか。
3. 案 B の場合、4B を GSE296758（雌 LepR+ 間質、Dox × 4、flowcell 交絡あり）に差し替えてよいか。GSE205237（ほぼ空の volcano）を残すか。
4. DEX を Supplementary に下げた結果、Fig5 の 2nd hit から RANKL/OPG とステロイドを外すか、Supplementary を根拠に小さく残すか。
5. 論文全体を **"hypothesis-generating"** と明記する方針でよいか。
6. 用語は "sequential two-insult" でよいか（"priming-and-trigger" は順序の証拠が無いため控える案）。

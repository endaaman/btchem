---
title: "第 7 章 Figure 5（概念図）と Graphical abstract の下書き"
date: 2026-10-26
author: 遠田
---

現行の Figure 5（概念図）と Graphical abstract（GA）の下書きです。各案で Figure 5 と GA がどう変わるかは、第 6 章の「Figure 5 と GA は各案でどう変わるか」にあります。


# 下書きについて

Figure 5（2 段階の hit による骨転移の概念図）と Graphical abstract（GA）の下書きです。drawio で作成し、後で BioRender のイラストに差し替える前提です（現段階では**位置と文言**を固めることが目的）。

Figure 5 は 2 つの版を作りました。

- **清水原案版**: ご依頼資料の文言（2nd hit = 骨細胞アポトーシス）をそのまま図にしたもの。比較のための版です。
- **修正版**: Figure 4 / S1–S3 の解析結果に合わせて中身を差し替えたもの。**こちらを推します。**

構図（3 カラム: 1st hit ＋ 2nd hit → pre-metastatic niche / 骨転移、右下に「臨床で介入できる / できない」の 2 箱）はどちらも同じで、ご依頼資料の案に沿っています。

# 修正版（推奨）

![Figure 5 修正版](design/fig5_revised.png){width=100%}

# 清水原案版（比較用）

![Figure 5 清水原案版](design/fig5_shimizu.png){width=100%}

# 2 つの版の違い

| | 清水原案版 | 修正版（データに基づく） | 修正の根拠（Figure 4 資料） |
|---|---|---|---|
| 1st hit の名前 | Menopause (estrogen depletion) | Postmenopausal (aged) bone | OVX（エストロゲン欠乏）単独では骨の転写変化がほぼ無く、炎症が出たのは加齢データ |
| 1st hit の中身 | 炎症性サイトカイン（IL-6 など）↑、RANKL ↑ | 炎症性サイトカイン（IL-1β, TNF, S100A8/A9）↑、細胞老化（p21, SenMayo）↑ | マウス加齢で *Tnf*, *Il1b*, *S100a8/9* ↑。ヒト高齢・閉経後で p21・SenMayo ↑。**RANKL はマウス加齢で低下、ヒトでは変化なし**。IL-6 はほとんどのデータで発現が低く測れない |
| 2nd hit の名前 | Chemotherapy + steroids | Perioperative chemotherapy + supportive-care steroids | 同じ（表現のみ） |
| 2nd hit の中身 | 骨細胞アポトーシス ↑（FAS/FASL, caspase-3, BCL-2） | RANKL/OPG 比 ↑（OPG 低下、ステロイド）、炎症因子 ↑（TNF, IL-1β, S100A9, MMP9、ステロイド）、骨形成 ↓（*Spp1*, *Dmp1*、化療） | **化療・ステロイドのどのデータでもアポトーシス遺伝子は有意に上がらなかった**。DEX で OPG 低下と炎症因子上昇（OPG 低下はヒト GIOP でも一致）。雌の反復 Dox で骨形成系の低下 |
| niche 〜 骨転移 | 破骨細胞活性化 → 骨吸収 → DTC 再活性化 → 骨転移 | 同じ（骨吸収で基質から増殖因子が放出される、を追記） | 文献（Ottewell 2014/2015、Lawson 2015） |
| 臨床で介入できる | デノスマブ、BP | RANKL 阻害（デノスマブ）、BP | 同じ |
| 介入できない | アポトーシス因子、炎症性サイトカイン（IL-6） | 炎症、老化、骨形成の抑制 | 2nd hit の中身に合わせた |

（DEX のマウス脛骨で上がる炎症因子の遺伝子（*S100a9*・*Tnf*・*Mmp9* など）は、骨髄の細胞構成が好中球・単球の側に傾いた変化と読むのが正確です（第 3 章 案 2）。2nd hit の中身は、第 6 章の各案と第 3 章の結果に合わせて組み直します。）

# 修正版を推す理由

1. **図とデータが矛盾しない**: 原案版の "Data: Fig. 4B–C" の行は、Figure 4 で示すデータ（アポトーシスは出ない、RANKL は上がらない）と食い違います。査読者は図と Figure 4 を見比べるので、そのまま使うと最初に指摘される箇所になります。
2. **論理の筋は原案と変わらない**: 「閉経後の骨が下地を作り、化療（＋ステロイド）が重なって破骨細胞が活性化し、休眠中の DTC が目覚める」という骨格はそのままです。変えたのは 2nd hit の「中身」を、アポトーシスから **RANKL/OPG 比・炎症・骨形成抑制**に置き換えた点だけです。出口（破骨細胞 → 骨吸収）は同じなので、デノスマブ / BP を提案する結論にもつながります。
3. **むしろ介入の論理が強くなる**: 2nd hit が「RANKL/OPG 比の上昇」であれば、RANKL 阻害（デノスマブ）がそこを直接止める、という説明が素直に成り立ちます。
4. **アポトーシスは捨てるのではなく、Discussion へ**: 化療・ステロイドによる骨細胞アポトーシスは組織学的な既報（Weinstein 1998、Wasi 2024 など）があるので、「既報ではこういう機序も報告されている」と Discussion で触れるのが安全です。

# 用語: "two-hit" について

"two-hit" は Knudson のがん抑制遺伝子仮説（2 つのアレルの不活化）を強く連想させるため、査読で誤解・指摘を受ける可能性があります（文献資料 7 節）。
本文では **"sequential two-insult model"** または **"priming-and-trigger model"**（1st hit = priming、2nd hit = trigger）と定義してから使うことをお勧めします。
図の "1st hit / 2nd hit" のラベルも、決まれば "Priming / Trigger" または "First insult / Second insult" に揃えます。

# Figure 5 legend 案（修正版、英文）

> **Figure 5. Proposed sequential two-insult (priming-and-trigger) model of chemotherapy-associated bone metastasis in postmenopausal women.**
> The first insult is the aged, postmenopausal bone microenvironment, characterized by increased inflammatory cytokines
> (IL-1β, TNF, S100A8/A9) and cellular senescence (p21, SenMayo signature) (Figure 4A; Supplementary Figure S1C, D). The second insult is perioperative
> chemotherapy with supportive-care glucocorticoids, which decreases osteoprotegerin (OPG) and thereby increases the
> RANKL/OPG ratio, induces inflammatory factors (Figure 4C; Supplementary Figure S1E), and suppresses the bone-formation program (*Spp1*, *Dmp1*)
> (Supplementary Figure S1B). Together, these changes promote osteoclast activation and bone resorption, releasing
> matrix-derived growth factors that may reactivate disseminated tumor cells (DTC), forming a pre-metastatic niche in bone.
> Of these steps, only RANKL-driven osteoclast activity is currently targetable in the clinic (denosumab or
> bisphosphonates), whereas inflammation, senescence and suppressed bone formation are not.

（図中の "Data: Suppl. Fig. [ ]" は空欄のままです。上の legend の対応（S1B–E）で確定してよければ図も埋めます。
2nd hit の Data 行は "Fig. 4B–C" としていますが、4B は「持続的な変化なし」の図なので、"Fig. 4C; Suppl. Fig. S1B, E" に絞る案もあります。）

（legend の "induces inflammatory factors (Figure 4C; Supplementary Figure S1E)" も、2nd hit を組み直すときに直します。）

# Graphical abstract（60 × 50 mm）

![Graphical abstract 下書き](design/ga.png){width=60%}

- 左「Postmenopausal women」: KM の**模式図**（実データの曲線ではありません）と **HR 1.72（multivariable）**。
- 中「Two-hit model」: Hit 1 = Aged bone: inflammation、Hit 2 = Chemo + steroid: RANKL/OPG ↑ → 破骨細胞活性化 → 骨転移（修正版の要約）。
- 右「Proposal」: 化療前の BMA（デノスマブまたは BP）、対象 = 周術期化療を受ける閉経後 ER+ 早期乳がん。

注意:

- **HR 1.72 は抄録の閉経後・多変量の値**です。Table 3 の確定値が出たら更新します（ご依頼資料のメモでは閉経後の p 値が変動中とのことでした）。
  参考: Figure 3 の閉経後の**単変量** HR は 2.46（95% CI 1.37–4.40）です。
- 中央の見出し "Two-hit model" は、上の用語の結論に合わせて "Priming and trigger" などに変えます。
- KM は Figure 3 確定後に実際の曲線の縮小版に差し替えることもできます。
- 右の提案は、D-CARE 試験（化療と同時のデノスマブで骨転移無発生生存は改善せず）があるため、本文では「推奨」ではなく「検証すべき仮説」として書くのが安全です（文献資料を参照）。GA の文言もそれに合わせて弱めるか、ご相談させてください。

# 清水先生に確認したいこと

1. Figure 5 は修正版の方向（2nd hit を骨細胞アポトーシスから、データで出た変化に置き換える）で進めてよいか。2nd hit の中身は第 6 章の確認とあわせて決めます。
2. GA の HR を何にするか（Table 3 の閉経後多変量 HR の確定値）。
3. GA の "Proposal" の書き方（D-CARE を踏まえて、"hypothesis" 寄りの表現にするか）。

# 手描き系の図（Fig1 / Fig5 / Graphical abstract）の下書き（2026-10-12）

drawio 下書き。後で BioRender のイラストに差し替える前提で、**位置と文言**を固めるのが目的。
配色・箱の書式は 1本目 Fig1（`,/参考1/data/fig1/fig1_rev1.drawio`）を踏襲:
灰 `#F2F2F2`（コホート）/ 緑 `#EAFFEB` + 見出し `#005200` / 桃 `#F8E1E1`（悪い側 = 化療あり・2nd hit）/
青灰 `#DCE3E8`（化療なし・1st hit）/ 主アウトカム赤 `#C44E52`・副アウトカム青 `#4A7BA7`。Arial、角丸箱、黒枠。
Fig1 の女性アイコンは 1本目の drawio に埋め込まれている画像をそのまま使った。

## ファイル

| ソース | 出力 | 最終サイズ |
|---|---|---|
| `design/fig1.drawio` | `out/design/fig1.{png,pdf}` | 170 × 143.2 mm |
| `design/fig5_revised.drawio` | `out/design/fig5_revised.{png,pdf}` | 170 × 111.5 mm |
| `design/fig5_shimizu.drawio` | `out/design/fig5_shimizu.{png,pdf}` | 170 × 111.5 mm |
| `design/ga.drawio` | `out/design/ga.{png,pdf}` | 60 × 49.7 mm |

`design/render.sh [name ...]` で一括書き出し（png は `-s 4 -b 20`、pdf は `--crop -b 20`）。

## サイズとフォント換算

**drawio の 2 px = 最終 1 pt** で作図した。各ページに不可視の枠（id `frame`、線・塗りなし）を置き、
その寸法が最終サイズに一致する:

- 170 mm = 482 pt = **964 px**（Fig1, Fig5）
- 60 mm × 50 mm = 170 × 141.7 pt = **340 × 282 px**（GA、高さ 49.7 mm）

よって fontSize の px ÷ 2 が最終 pt:

| px | pt | 使い所 |
|---|---|---|
| 20 | 10 | Fig1/Fig5 の見出し |
| 18 | 9 | Fig1 の N・Outcomes の行頭、Fig5 の箇条見出し |
| 17 | 8.5 | 本文 |
| 16 | 8 | 補足（Excluded の箇条、Data: 行、介入点の中身、骨・細胞ラベル） |
| 12–14 | 6–7 | GA のみ（最小 12 px = 6 pt） |

線幅: 枠 1.5 px = 0.75 pt、GA の枠 1 px = 0.5 pt、矢印 2 px = 1 pt（すべて 0.25 pt 以上）。
単一パネルなのでパネルラベル（A, B…）は付けていない。

注意（書き出し）: `-b 20` の余白は png に含まれる（Fig1/Fig5 は 3958 px 幅 = 964×4 + 余白）。
**入稿で幅 170 mm に貼るときは余白を除いた版を使う**こと（余白込みで 170 mm に縮めると 2.6% 小さくなり、
8 pt が 7.8 pt になる）。`drawio ... -b 0` で書き出すと 964 px × 4 = 3856(+2) px で、170 mm で 576 dpi。
PDF は drawio が 1 px = 0.72 pt で出す（Fig1 695 × 585 pt）ので、170 mm 幅に縮小して使う
（1本目と同じく pypdf でページ縮小すれば文字はベクタのまま）。

## Figure 1（study design）

縦の本流（Initial screening → Multicenter retrospective cohort）+ 右に Excluded、下で化療 (+)/(−) に分岐し、
それぞれ閉経前/後に分ける。最後に Outcomes（1本目と同じ「Primary/Secondary」の色付き行）。
1本目は縦長 80 mm 幅だったが、今回は 2×2 の層別が入るので 170 mm 幅の横長にした。

数値の出典: 862 / 77 / 785 / 97.8% / 325 / 460 / 除外理由 / 統計の共変量 = 依頼資料の Fig1 案（AI 作図）。
785・2007–2020・ER+・BMA 非使用 = 抄録。

### プレースホルダ

| 場所 | 中身 |
|---|---|
| 化療 (+) の n | `n = [ ]` |
| 化療 (−) の n | `n = [ ]` |
| 化療 (+) × 閉経前 / 閉経後 | `n = [ ]` ×2 |
| 化療 (−) × 閉経前 / 閉経後 | `n = [ ]` ×2 |

### 清水に確認すること

- 除外 77 例の内訳の n（図には理由だけ載せた）。
- 抄録の方法では「リンパ節、肺、肝臓、骨」の 5 年無転移生存を比較しているが、依頼のとおり Secondary は
  「骨・肺・肝」にした。リンパ節を入れるか。
- 依頼資料の Fig1 案は「Perioperative use of BMA for baseline osteoporosis」、抄録は「術後補助療法として BMA を使用していない」。
  図は "Perioperative BMA use" とした。
- Statistics の共変量（年齢、BMI、stage、RO score、RT、内分泌療法）は依頼資料の Fig1 案のまま。Table 2/3 と一致させる。

### Legend 案

> **Figure 1. Study design.** Women who underwent surgery for estrogen receptor-positive, HER2-negative,
> clinical stage 2–3 early breast cancer (eBC) between 2007 and 2020 at participating institutions were screened
> (N = 862). After exclusion of male patients, bilateral breast cancer, patients without a perioperative chest
> radiograph and those receiving perioperative bone-modifying agents (BMA) (n = 77), 785 women were enrolled
> (5-year follow-up rate, 97.8%). Patients were stratified by perioperative chemotherapy and menopausal status.
> The primary endpoint was 5-year bone metastasis-free survival (BMFS); secondary endpoints were organ-specific
> (bone, lung and liver) metastasis-free survival. Survival was compared using the Kaplan–Meier method with the
> log-rank test, and multivariable Cox proportional hazards models were adjusted for age, body mass index, stage,
> risk of osteoporosis (RO) score, radiotherapy (RT) and endocrine therapy.

## Figure 5（2-hit 概念図）

3 カラム（依頼資料の AI 案の構図）: **1st hit**（青灰）＋ **2nd hit**（桃）→ **Pre-metastatic niche / Bone metastasis**（灰）。
各カラムの上部に drawio 図形の簡略イラスト（骨 = 黄土色の矩形、骨細胞 = 楕円 + 樹状突起、破骨細胞 = 多核の紫楕円、
DTC = 橙の小円クラスタ、薬 = カプセル）、下部に文言。右カラムの下に「Clinically targetable」（緑）と
「Not targetable at present」（灰・赤字見出し）の 2 箱。イラストは BioRender に差し替える前提。

### fig5_revised（解析結果に基づく版。`docs/prelim_summary.md`「Fig5」節の論理）

| カラム | 文言 |
|---|---|
| 1st hit: Postmenopausal (aged) bone | ↑ Inflammatory cytokines: IL-1β, TNF, S100A8/A9 / ↑ Cellular senescence: p21 (CDKN1A), SenMayo signature / Data: Fig. 4A; Suppl. Fig. [ ] |
| 2nd hit: Perioperative chemotherapy + supportive-care steroids | ↑ RANKL/OPG ratio: ↓ OPG (*Tnfrsf11b*) — steroid / ↑ Inflammatory factors: TNF, IL-1β, S100A9, MMP9 — steroid / ↓ Bone formation: *Spp1*, *Dmp1* — chemotherapy / Data: Fig. 4B–C; Suppl. Fig. [ ] |
| Pre-metastatic niche → Bone metastasis | Osteoclast activation → Bone resorption: release of matrix-derived growth factors → DTC reactivation → Bone metastasis |
| 介入点 | Targetable: RANKL inhibition (denosumab), bisphosphonates (osteoclasts) / Not targetable: inflammation, senescence, suppressed bone formation |

- 骨細胞の周りの赤い小点 = 老化（SASP）を示す記号。
- 2nd hit のイラストのラベルは "Osteocyte / osteoblast lineage"（*Spp1*/*Dmp1* の低下は Dox×4 の MSPC = 骨芽細胞系、
  Fig. 4B の皮質骨では持続的変化なし）。
- *Spp1*/*Dmp1* の出典は Suppl.（GSE296758）で、4B 自体は「持続的変化なし」の負の結果。Data 行を "Fig. 4B–C" と
  まとめているが、Fig4 確定後に 4C + Suppl. へ絞るか判断。
- 遺伝子名の表記ゆれ: 概念図なのでタンパク名（IL-1β, TNF…）が基本、マウス遺伝子でしか言えないもの（*Tnfrsf11b*, *Spp1*, *Dmp1*）は
  イタリックの遺伝子名にした。

### fig5_shimizu（清水原案の文言を残した版）

| カラム | 文言 |
|---|---|
| 1st hit: Menopause (estrogen depletion) | ↑ Inflammatory cytokines: IL-6 and others / ↑ RANKL: osteoclast induction / Data: Fig. 4A |
| 2nd hit: chemotherapy + steroids | ↑ Osteocyte apoptosis: FAS/FASL, caspase-3, BCL-2 / Data: Fig. 4B–C |
| Niche | Osteoclast activation → Bone resorption: pre-metastatic niche formation → DTC reactivation → Bone metastasis |
| 介入点 | Targetable: denosumab, BP / Not targetable: apoptotic factors, inflammatory cytokines (IL-6) |

**注意**: この版の 2nd hit（アポトーシス）と 1st hit の RANKL↑ は今回の解析では支持されていない
（`docs/prelim_summary.md`: 化療・ステロイドでアポトーシスはどのデータでも出ない、加齢骨では *Tnfsf11* はむしろ↓、
OVX では炎症も出ない）。図の "Data:" 行は清水案の想定を書いただけで、そのまま使うと図とデータが矛盾する。
比較・説明用の版として作った。

### Legend 案（revised 版）

> **Figure 5. Proposed two-hit model of chemotherapy-associated bone metastasis in postmenopausal women.**
> The first hit is the aged, postmenopausal bone microenvironment, characterized by increased inflammatory cytokines
> (IL-1β, TNF, S100A8/A9) and cellular senescence (p21, SenMayo signature) (Figure 4A). The second hit is perioperative
> chemotherapy with supportive-care glucocorticoids, which decreases osteoprotegerin (OPG) and thereby increases the
> RANKL/OPG ratio, induces inflammatory factors (Figure 4C), and suppresses the bone-formation program (*Spp1*, *Dmp1*)
> (Supplementary Figure [ ]). Together, these changes promote osteoclast activation and bone resorption, releasing
> matrix-derived growth factors that may reactivate disseminated tumor cells (DTC), forming a pre-metastatic niche in bone.
> Of these steps, only RANKL-driven osteoclast activity is currently targetable in the clinic (denosumab or
> bisphosphonates), whereas inflammation, senescence and suppressed bone formation are not.

## Graphical abstract（60 × 50 mm）

3 パネル（各 108 px = 19 mm 幅）:
- 左 "Postmenopausal women": KM の模式（階段曲線。Chemo (−) 青が上、Chemo (+) 赤が下、軸ラベル BMFS / Years、目盛り・数値なし）、
  下に **HR 1.72 (multivariable)**（抄録の閉経後多変量 HR）。
- 中 "Two-hit model": Hit 1 Aged bone: inflammation（青灰）+ Hit 2 Chemo + steroid: RANKL/OPG ↑（桃）→ Osteoclast activation → bone metastasis（revised 版の要約）。
- 右 "Proposal": **BMA before chemo (denosumab or BP)**（緑箱）+ "for postmenopausal women with ER+ eBC receiving perioperative chemotherapy"。

KM は模式で実データの曲線ではない（Fig2/3 確定後に実曲線の縮小版へ差し替えも可）。最小フォント 6 pt（12 px）。
HR 1.72 は Fig3/Table3 の確定値に合わせて更新すること（依頼資料のメモでは閉経後 p 値が変動中）。

## 作り方メモ

drawio XML は Python の生成スクリプト（セッションの scratchpad、repo 外）で書き出した。以後の修正は drawio を直接編集する前提
（生成器は保存していない）。

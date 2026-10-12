---
title: "Figure 1 — Study design（研究デザイン図）"
date: 2026-10-12
author: 遠田
---

# この資料について

Figure 1（研究デザインのフロー図）の下書きと、図に載せる症例数の確定値をまとめたものです。
図は drawio で作成した下書きで、最終的には BioRender のイラストに差し替える前提です（今回は**位置と文言**を固めることが目的）。
配色・箱の書式は 1 本目の Fig1 を踏襲しています。

![Figure 1 下書き（`out/design/fig1.png`）](design/fig1.png){width=100%}

# 図の構成

- 縦の本流: Initial screening（N = 862）→ Multicenter retrospective cohort（N = 785）。右に Excluded（n = 77）。
- 下で周術期化療 (+) / (−) に分岐し、それぞれを閉経前 / 閉経後に分ける（2 × 2 の層別）。
- 最後に Outcomes: Primary = 5 年骨転移無発生生存（BMFS、赤）、Secondary = 臓器別無転移生存（骨・肺・肝、青）、Statistics（KM + log-rank、多変量 Cox）。
- 1 本目は縦長 80 mm 幅でしたが、2 × 2 の層別が入るため今回は 170 mm 幅の横長にしました（仕上がり 170 × 143.2 mm）。

# 症例数（確定値）

化療の有無と閉経状態の n は、1 本目の登録症例データ（`perioperative chemo` 列）から集計した値です。

| 区分 | 閉経前 | 閉経後 | 計 |
|---|---:|---:|---:|
| 登録（Enrolled） | 325 | 460 | 785 |
| 化療情報が空欄のため除外 | 23 | 20 | 43 |
| **解析対象** | **302** | **440** | **742** |
| うち 化療 (−) | 81 | 217 | 298 |
| うち 化療 (+) | 221 | 223 | 444 |

- 862 / 77 / 785 / 97.8% / 325 / 460 / 除外理由 / 共変量は、ご依頼資料の Fig1 案と抄録から取りました。
- 解析対象 742 例と 2 × 2 の n は、Figure 3 の KM（`out/fig/fig3_stats.csv`）の n と一致しています。

## 図への反映状況

- 依頼時点では化療別の n は `[ ]`（空欄）でしたが、**現在の下書き（上図）には上表の解析対象の値（444 / 298 / 221 / 223 / 81 / 217）が入っています。**
- ただし、**化療情報が空欄で除外した 43 例（閉経前 23 / 閉経後 20）は図にまだ示されていません**。図では 785 例からそのまま分岐しているため、
  785 → 742 の差が読めません。「Chemotherapy status unknown (n = 43)」の箱を足すなど、**図は後で更新**します。
- 43 例の化療情報が先生の手元で補えるなら、解析対象が 785 例に戻り、図と KM の n を差し替えます。

# Figure legend 案（英文）

> **Figure 1. Study design.** Women who underwent surgery for estrogen receptor-positive, HER2-negative,
> clinical stage 2–3 early breast cancer (eBC) between 2007 and 2020 at participating institutions were screened
> (N = 862). After exclusion of male patients, bilateral breast cancer, patients without a perioperative chest
> radiograph and those receiving perioperative bone-modifying agents (BMA) (n = 77), 785 women were enrolled
> (5-year follow-up rate, 97.8%). After exclusion of 43 women with unknown perioperative chemotherapy status,
> 742 women were stratified by perioperative chemotherapy (with, n = 444; without, n = 298) and menopausal status.
> The primary endpoint was 5-year bone metastasis-free survival (BMFS); secondary endpoints were organ-specific
> (bone, lung and liver) metastasis-free survival. Survival was compared using the Kaplan–Meier method with the
> log-rank test, and multivariable Cox proportional hazards models were adjusted for age, body mass index, stage,
> risk of osteoporosis (RO) score, radiotherapy (RT) and endocrine therapy.

（43 例の扱いが決まったら、"After exclusion of 43 women ..." の文を調整します。）

# 清水先生に確認したいこと

1. **化療情報が空欄の 43 例**（閉経前 23 / 閉経後 20）: 化療の有無を補えるか。補えない場合、図に「除外」として明示してよいか。
2. **除外 77 例の内訳の n**: 図には理由（男性、両側乳がん、周術期の胸部 X 線なし、周術期 BMA 使用）だけを載せています。理由ごとの n をいただければ追記します。
3. **Secondary endpoint にリンパ節を入れるか**: 抄録の方法では「リンパ節、肺、肝臓、骨」の 5 年無転移生存を比較していますが、ご依頼どおり図は「骨・肺・肝」にしました。
4. **BMA の表記**: 依頼資料の Fig1 案は "Perioperative use of BMA for baseline osteoporosis"、抄録は「術後補助療法として BMA を使用していない」です。図は "Perioperative BMA use" としました。
5. **化療の定義**: 群分けは「周術期化療」でよいか（抄録には「術後補助化学療法 (AC)」という表現も出てきます）。今回の n は 1 本目データの `perioperative chemo` 列によるものです。
6. **Statistics の共変量**（年齢、BMI、stage、RO score、RT、内分泌療法）: 依頼資料の Fig1 案のままです。Table 2 / 3 の共変量と一致させてください。

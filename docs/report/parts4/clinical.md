---
title: "第 4 章 臨床パート（Figure 1–3）— v1 からの変更点と現状"
date: 2026-10-26
author: 遠田
---

# v1（10/12 資料）からの変更

- Fig. 3（骨転移の KM）は 1 本目のコホートデータで実データ描画済み（変更なし）。Fig. 2（肺・肝）は拡張データ未着のまま。
- Fig. 1 の n は化療 (+) 444 / (-) 298（計 742）で埋めた。化療不明 43 例の扱いは未決。
- 以下は v1 資料の該当章をそのまま再掲する。


# Figure 1 — Study design（v1 再掲）

## 下書きの図（v1 の前書きより）

Figure 1（研究デザインのフロー図）の下書きと、図に載せる症例数の確定値をまとめたものです。
図は drawio で作成した下書きで、最終的には BioRender のイラストに差し替える前提です（今回は**位置と文言**を固めることが目的）。
配色・箱の書式は 1 本目の Fig1 を踏襲しています。

![Figure 1 下書き（`out/v2/design/fig1.png`）](v2/design/fig1.png){width=100%}

## 図の構成

- 縦の本流: Initial screening（N = 862）→ Multicenter retrospective cohort（N = 785）。右に Excluded（n = 77）。
- 下で周術期化療 (+) / (−) に分岐し、それぞれを閉経前 / 閉経後に分ける（2 × 2 の層別）。
- 最後に Outcomes: Primary = 5 年骨転移無発生生存（BMFS、赤）、Secondary = 臓器別無転移生存（骨・肺・肝、青）、Statistics（KM + log-rank、多変量 Cox）。
- 1 本目は縦長 80 mm 幅でしたが、2 × 2 の層別が入るため今回は 170 mm 幅の横長にしました（仕上がり 170 × 143.2 mm）。

## 症例数（確定値）

化療の有無と閉経状態の n は、1 本目の登録症例データ（`perioperative chemo` 列）から集計した値です。

| 区分 | 閉経前 | 閉経後 | 計 |
|---|---:|---:|---:|
| 登録（Enrolled） | 325 | 460 | 785 |
| 化療情報が空欄のため除外 | 23 | 20 | 43 |
| **解析対象** | **302** | **440** | **742** |
| うち 化療 (−) | 81 | 217 | 298 |
| うち 化療 (+) | 221 | 223 | 444 |

- 862 / 77 / 785 / 97.8% / 325 / 460 / 除外理由 / 共変量は、ご依頼資料の Fig1 案と抄録から取りました。
- 解析対象 742 例と 2 × 2 の n は、Figure 3 の KM（`out/v2/fig/fig3_stats.csv`）の n と一致しています。

### 図への反映状況

- 依頼時点では化療別の n は `[ ]`（空欄）でしたが、**現在の下書き（上図）には上表の解析対象の値（444 / 298 / 221 / 223 / 81 / 217）が入っています。**
- ただし、**化療情報が空欄で除外した 43 例（閉経前 23 / 閉経後 20）は図にまだ示されていません**。図では 785 例からそのまま分岐しているため、
  785 → 742 の差が読めません。「Chemotherapy status unknown (n = 43)」の箱を足すなど、**図は後で更新**します。
- 43 例の化療情報が先生の手元で補えるなら、解析対象が 785 例に戻り、図と KM の n を差し替えます。

## Figure legend 案（英文）

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


# Figure 2 / 3（v1 再掲）

## 現状のまとめ

| 図 | 内容 | 状態 |
|---|---|---|
| Figure 2 | 骨・肺・肝の無転移生存（KM、全体） | **未作成**。肺・肝のイベントを含む拡張データをお待ちしています |
| Figure 3 | 骨転移無発生生存（BMFS）を閉経前後で層別 | **実データで描画済み**（1 本目のデータを使用） |

- 図の書式は 1 本目の Fig3 と同じにしました（170 mm 幅・1 行 3 パネル、打ち切りのひげ、At risk / Censored / Events の表、凡例に log-rank p・HR・95% CI）。
- 色は 化療 (+) = 赤、化療 (−) = 青です。
- Figure 2 用のプログラムは用意済みで、データが届けば列名を合わせるだけで同じ書式の図が出ます。

## Figure 3（実データ）

![Figure 3（`out/v2/fig/fig3.png`）](v2/fig/fig3.png){width=100%}

### データ

- 1 本目の登録症例ファイルに `perioperative chemo`（0/1）の列があり、BMFS の月数・イベント・閉経状態のデータと研究 ID で 785 例が 1 対 1 で対応しました（イベントと閉経状態は完全に一致）。
- 化療欄が空欄の 43 例（閉経前 23 / 閉経後 20）を除外し、**解析対象は 742 例**です（Figure 1 の資料を参照）。

### 結果（`out/v2/fig/fig3_stats.csv`）

HR は化療 (+) vs (−) の**単変量** Cox 比例ハザードモデル。p は log-rank 検定（括弧内は Cox の Wald 検定）。

| パネル | 化療 (−) n（イベント） | 化療 (+) n（イベント） | HR (95% CI) | log-rank p（Wald p） |
|---|---|---|---|---|
| A 全体 | 298 (25) | 444 (69) | 1.89 (1.20–2.99) | 0.006 (0.006) |
| B 閉経前 | 81 (9) | 221 (30) | 1.24 (0.59–2.61) | 0.57 (0.57) |
| C 閉経後 | 217 (16) | 223 (39) | 2.46 (1.37–4.40) | 0.002 (0.002) |

- **閉経後でのみ有意**という抄録の所見と、向きが一致しました。閉経前では HR 1.24 で有意差はありません。
- 全体の HR 1.89 は抄録の HR 2.05 と一致しません。考えられる理由は次の 2 点です（下の確認事項 1）。
  - 抄録の値が多変量の可能性（本図は単変量）。
  - 解析対象が 742 例（化療不明 43 例を除外）である。
- Table 3 の閉経後多変量 HR 1.72 は多変量なので、本図の単変量 HR 2.46 と一致しなくて正常です。

## Figure 2（未作成）

- 依頼内容: 単変量 KM、骨・肺・肝の 3 パネル（脳は不要）。
- 1 本目のデータには肺・肝のイベント列が無いため、**拡張データをお待ちしています**。
- 必要な列（名前は違っていても構いません）: 研究 ID / 閉経状態（0/1）/ 周術期化療（0/1）/ 骨・肺・肝それぞれの「イベント（0/1）」と「無転移生存期間（月）」。
- 届いたら、骨のパネル（2A）が Figure 3A と一致することを確認したうえで作図します。
- 作図の動作確認用に擬似データ版も作ってありますが、**実データではないため、本資料には図も数値も載せていません**。

## Figure legend 案（英文）

**Figure 2.** Kaplan–Meier curves of metastasis-free survival according to perioperative chemotherapy in the analysed women (n = [ ]).
(A) Bone metastasis-free survival. (B) Lung metastasis-free survival. (C) Liver metastasis-free survival.
Red, women who received perioperative chemotherapy (Chemo (+)); blue, women who did not (Chemo (−)).
Vertical ticks indicate censored observations. P values were calculated with the log-rank test; hazard ratios (HR) and 95% confidence
intervals (CI) for Chemo (+) versus Chemo (−) were estimated with univariate Cox proportional hazards models.
The numbers at risk, censored and events at each time point are shown below each panel.

**Figure 3.** Kaplan–Meier curves of bone metastasis-free survival according to perioperative chemotherapy, stratified by menopausal status.
(A) All analysed women (n = 742). (B) Premenopausal women (n = 302). (C) Postmenopausal women (n = 440).
Forty-three women with unknown perioperative chemotherapy status were excluded.
Red, Chemo (+); blue, Chemo (−). Vertical ticks indicate censored observations. P values were calculated with the log-rank test;
HR and 95% CI for Chemo (+) versus Chemo (−) were estimated with univariate Cox proportional hazards models.
The numbers at risk, censored and events at each time point are shown below each panel.

（n は化療不明 43 例の扱いが決まったら確定します。Figure 2 の n は拡張データ到着後に記入します。）


# 清水先生に確認したいこと

1. **化療情報が空欄の 43 例**（閉経前 23 / 閉経後 20）: 化療の有無を補えるか。補えない場合、図に「除外」として明示してよいか。
2. **除外 77 例の内訳の n**: 図には理由（男性、両側乳がん、周術期の胸部 X 線なし、周術期 BMA 使用）だけを載せています。理由ごとの n をいただければ追記します。
3. **Secondary endpoint にリンパ節を入れるか**: 抄録の方法では「リンパ節、肺、肝臓、骨」の 5 年無転移生存を比較していますが、ご依頼どおり図は「骨・肺・肝」にしました。
4. **BMA の表記**: 依頼資料の Fig1 案は "Perioperative use of BMA for baseline osteoporosis"、抄録は「術後補助療法として BMA を使用していない」です。図は "Perioperative BMA use" としました。
5. **化療の定義**: 群分けは「周術期化療」でよいか（抄録には「術後補助化学療法 (AC)」という表現も出てきます）。今回の n は 1 本目データの `perioperative chemo` 列によるものです。
6. **Statistics の共変量**（年齢、BMI、stage、RO score、RT、内分泌療法）: 依頼資料の Fig1 案のままです。Table 2 / 3 の共変量と一致させてください。


1. **抄録の骨転移 HR 2.05 (95% CI 1.27–3.32)** は単変量ですか、多変量ですか。また、何例で計算した値ですか。
   単変量・785 例であれば、43 例の化療情報が補えれば本図の HR と一致するはずです。
2. **化療情報が空欄の 43 例**（閉経前 23 / 閉経後 20）の化療の有無。
3. **Table 3 のメモとの食い違い**: 依頼資料の Table 3 のメモには「現時点 閉経後 p = 0.98、閉経前 p = 0.03」とあり、抄録（閉経後で有意）と本図（閉経後で有意、閉経前 p = 0.57）とは**逆**です。データの版の違いか、記載ミスかを教えてください。
4. **群の定義**: 「周術期化療」と「術後補助化療 (AC)」のどちらで分けるか（値が変わります）。
5. **抄録の p 値**は log-rank と Cox の Wald のどちらですか（stats の表には両方を出してあります）。
6. Figure 2 用の拡張データ（肺・肝のイベントと期間）の送付時期。

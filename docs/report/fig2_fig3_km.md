---
title: "Figure 2 / 3 — 周術期化療と転移無発生生存（Kaplan–Meier）"
date: 2026-10-12
author: 遠田
---

# 現状のまとめ

| 図 | 内容 | 状態 |
|---|---|---|
| Figure 2 | 骨・肺・肝の無転移生存（KM、全体） | **未作成**。肺・肝のイベントを含む拡張データをお待ちしています |
| Figure 3 | 骨転移無発生生存（BMFS）を閉経前後で層別 | **実データで描画済み**（1 本目のデータを使用） |

- 図の書式は 1 本目の Fig3 と同じにしました（170 mm 幅・1 行 3 パネル、打ち切りのひげ、At risk / Censored / Events の表、凡例に log-rank p・HR・95% CI）。
- 色は 化療 (+) = 赤、化療 (−) = 青です。
- Figure 2 用のプログラムは用意済みで、データが届けば列名を合わせるだけで同じ書式の図が出ます。

# Figure 3（実データ）

![Figure 3（`out/fig/fig3.png`）](fig/fig3.png){width=100%}

## データ

- 1 本目の登録症例ファイルに `perioperative chemo`（0/1）の列があり、BMFS の月数・イベント・閉経状態のデータと研究 ID で 785 例が 1 対 1 で対応しました（イベントと閉経状態は完全に一致）。
- 化療欄が空欄の 43 例（閉経前 23 / 閉経後 20）を除外し、**解析対象は 742 例**です（Figure 1 の資料を参照）。

## 結果（`out/fig/fig3_stats.csv`）

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

# Figure 2（未作成）

- 依頼内容: 単変量 KM、骨・肺・肝の 3 パネル（脳は不要）。
- 1 本目のデータには肺・肝のイベント列が無いため、**拡張データをお待ちしています**。
- 必要な列（名前は違っていても構いません）: 研究 ID / 閉経状態（0/1）/ 周術期化療（0/1）/ 骨・肺・肝それぞれの「イベント（0/1）」と「無転移生存期間（月）」。
- 届いたら、骨のパネル（2A）が Figure 3A と一致することを確認したうえで作図します。
- 作図の動作確認用に擬似データ版も作ってありますが、**実データではないため、本資料には図も数値も載せていません**。

# Figure legend 案（英文）

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

1. **抄録の骨転移 HR 2.05 (95% CI 1.27–3.32)** は単変量ですか、多変量ですか。また、何例で計算した値ですか。
   単変量・785 例であれば、43 例の化療情報が補えれば本図の HR と一致するはずです。
2. **化療情報が空欄の 43 例**（閉経前 23 / 閉経後 20）の化療の有無。
3. **Table 3 のメモとの食い違い**: 依頼資料の Table 3 のメモには「現時点 閉経後 p = 0.98、閉経前 p = 0.03」とあり、抄録（閉経後で有意）と本図（閉経後で有意、閉経前 p = 0.57）とは**逆**です。データの版の違いか、記載ミスかを教えてください。
4. **群の定義**: 「周術期化療」と「術後補助化療 (AC)」のどちらで分けるか（値が変わります）。
5. **抄録の p 値**は log-rank と Cox の Wald のどちらですか（stats の表には両方を出してあります）。
6. Figure 2 用の拡張データ（肺・肝のイベントと期間）の送付時期。

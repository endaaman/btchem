# Figure 2 / 3（臨床 KM、周術期化療 (+) vs (−)）作図メモ（2026-10-12）

> **現状の出力はすべて mock（擬似データ）。** 清水からの臨床データ（拡張データセット）が未着のため、
> 1 本目のデータで群分けを代用して描いている。`out/fig/*_mock*` の数値は**実データの結果ではない**ので、
> 本文・抄録・メールに引用しないこと。

- 依頼（清水 2026-09-30）: Fig2 = 単変量 KM、骨・肺・肝（脳は不要）。Fig3 = 骨の KM を閉経前後で層別。いずれも「一本目みたいに作る」。
- スクリプト: `fig/km.py`（読み込み・統計・描画の共通部）、`fig/fig2.py`, `fig/fig3.py`。`fig/all.sh` にも追記済み（xlsx が無い間は自動で `--mock`）
- 書式は 1 本目 `,/参考1/rev1/fig3.py` の最終形と同じ: 170 mm 幅・1 行 3 パネル、軸 37 × 42 mm、censor のひげ、
  At risk / Censored / Events × 10 か月刻み 7 列の表、凡例タイトルに log-rank p / HR / 95% CI（**全パネル**に載せる）、
  y 0.2 刻み・x 20 刻み・xlim 0–65.5、Arial 8/9 pt、パネルラベル 12 pt 太字。
- 色: 化療 (+) = `COLOR_HRO` `#C44E52`（赤）、化療 (−) = `COLOR_LRO` `#4A7BA7`（青）。描画・表の順は (−) → (+)。
- 統計（lifelines 0.30.3）: log-rank は `logrank_test`、HR は `CoxPHFitter` 単変量（共変量 chemo のみ、化療 (+) vs (−)、Efron）。
  p 表示は 1 本目と同じ（< 0.01 は `p < 0.01`、それ以外は小数 2 桁）。凡例の p は log-rank。

## 実行

```
uv run python fig/fig2.py --mock   # -> out/fig/fig2_mock.{png,pdf,tiff}, out/fig/fig2_mock_stats.csv
uv run python fig/fig3.py --mock   # -> out/fig/fig3_mock.{png,pdf,tiff}, out/fig/fig3_mock_stats.csv
# 臨床 xlsx 到着後（--mock を外すだけ）
uv run python fig/fig2.py          # -> out/fig/fig2.{png,pdf,tiff}, out/fig/fig2_stats.csv
uv run python fig/fig3.py          # -> out/fig/fig3.{png,pdf,tiff}, out/fig/fig3_stats.csv
uv run python fig/fig2.py --xlsx data/clinical/<別名>.xlsx   # パスだけ差し替える場合
```

mock の図には上端に赤字で `MOCK DATA (surrogate grouping / simulated events) — NOT FOR PUBLICATION` が入る
（本番では出ない）。出力名も `*_mock` に分けてあり、本番名 `fig2` / `fig3` は本データでしか作られない。
stats csv の列: `panel, title, endpoint, mock, n_neg, ev_neg, n_pos, ev_pos, logrank_p, hr, lo, hi, wald_p`
（neg = 化療 (−)、pos = 化療 (+)）。Number at risk 表の値は標準出力に出る。

仕上がり: 165.9 × 101.0 mm（mock バナー込み。本番はバナー分 5 mm 低くなる）。1 本目 fig3（170 × 102 mm 相当）と目視で並べて、
軸・表・凡例の位置とフォントがそろっていることを確認した。凡例は「Chemo (−) (n=360)」の分だけ 1 本目より横に広いが、軸内に収まる。

## 期待する入力（`fig/km.py` 冒頭の定数）

- `CLINICAL_XLSX = 'data/clinical/clinical.xlsx'`（仮名）、`CLINICAL_SHEET = 0`（先頭シート）
- 列の対応は `COLUMNS`（xlsx 上の列名 → 内部名）の 1 か所だけ。到着後に名前が違えば左側を直す:

| 想定列名 | 内容 | 値 |
|---|---|---|
| `研究ID` | 患者 ID | 一意 |
| `menopause` | 閉経状態 | 0 = 閉経前, 1 = 閉経後 |
| `chemo` | 周術期化療 | 0 = なし, 1 = あり |
| `BMFS_Event` / `BMFS` | 骨転移 / 骨転移なし生存期間（月） | 0/1 / 数値 |
| `LMFS_Event` / `LMFS` | 肺 | 同上 |
| `HMFS_Event` / `HMFS` | 肝 | 同上 |

- 列名・文字列値は前後の空白を除いてから照合する（1 本目の xlsx には `'High '` のような末尾空白があった）。
- `研究ID` が空の行（1 本目の csv にある末尾の空行 216 行など）は捨てる。必要な列に欠損があるか、0/1 列に他の値があれば assert で止まる。
- fig3 は骨の列しか読まないので、肺・肝の列が無くても動く。
- chemo が 'Yes'/'No' などの文字列、または術前/術後で列が分かれている場合は `load_clinical` に変換を 1 行足す
  （抄録には「周術期化学療法」と「術後補助化学療法 (AC)」の両方が出てくる。**どちらを群にするかは清水に確認**）。

## mock の中身（`make_mock`）

入力 `,/参考1/data/osteoporosis_bone_meta.csv`（1 本目、785 例）。

- chemo の代用 = `Cehst xray score`（High = 1, Low = 0。1 本目 rev0 の LRO/HRO）。**化療とは無関係**
- 骨 = 1 本目の BMFS / BMFS_Event そのまま
- 肺・肝 = BMFS のイベントを seed 固定（`MOCK_SEED = 20261012`）で間引いた擬似イベント（残す割合 肺 0.35 / 肝 0.25）。
  外したイベントは同じ時点で打ち切り扱い → 打ち切りのひげが実データより多い

## 到着後の確認手順

1. xlsx を `data/clinical/` に置き、列名を見て `COLUMNS` / `CLINICAL_XLSX` を合わせる。n = 785、閉経前 325 / 閉経後 460 か
   （1 本目と同じコホートなら一致するはず）。`研究ID` で 1 本目の csv と join し、menopause / BMFS / BMFS_Event が一致するかも見る。
2. `--mock` を外して fig2 / fig3 を出し、`out/fig/fig2_stats.csv`, `fig3_stats.csv` を抄録と照合する:

| 抄録（2026-09-30） | 対応する図の値 | 備考 |
|---|---|---|
| 骨 HR 2.05 (95% CI 1.27–3.32), p < 0.01 | Fig2 A（= Fig3 A）の HR / CI / log-rank p | 単変量なので一致するはず |
| 肺・肝は関連なし | Fig2 B, C の p | |
| 閉経後で AC と骨転移の関連 p < 0.001（層別） | Fig3 C の log-rank p | 層別単変量 |
| 閉経後 HR 1.72 (1.01–2.92), p = 0.04 | **図には出ない** | 多変量（Table 3、清水担当）。Fig3 C の単変量 HR とは一致しないのが正常 |

3. 食い違いがあれば、図は計算値のまま清水に確認する（1 本目 rev1 と同じ方針）。照合の際の注意:
   - 抄録は「閉経後で有意」だが、依頼資料の Table 3 メモは「現時点 閉経後 p = 0.98、閉経前 p = 0.03」と**逆**になっている。
     データ版の違いか記載ミスか、Fig3 の結果を出したところで清水に確認する。
   - 群が「周術期化療」か「術後補助化療 (AC)」かで値が変わる（上記）。
   - 抄録の p が log-rank か Cox の Wald かは不明。stats csv には両方出してある。

## Figure legend 案（英語、1 本目の書式）

**Figure 2.** Kaplan–Meier curves of metastasis-free survival according to perioperative chemotherapy in the enrolled women (n = 785).
(A) Bone metastasis-free survival. (B) Lung metastasis-free survival. (C) Liver metastasis-free survival.
Red, women who received perioperative chemotherapy (Chemo (+)); blue, women who did not (Chemo (−)).
Vertical ticks indicate censored observations. P values were calculated with the log-rank test; hazard ratios (HR) and 95% confidence
intervals (CI) for Chemo (+) versus Chemo (−) were estimated with univariate Cox proportional hazards models.
The numbers at risk, censored and events at each time point are shown below each panel.

**Figure 3.** Kaplan–Meier curves of bone metastasis-free survival according to perioperative chemotherapy, stratified by menopausal status.
(A) Enrolled women (n = 785). (B) Premenopausal women (n = 325). (C) Postmenopausal women (n = 460).
Red, Chemo (+); blue, Chemo (−). Vertical ticks indicate censored observations. P values were calculated with the log-rank test;
HR and 95% CI for Chemo (+) versus Chemo (−) were estimated with univariate Cox proportional hazards models.
The numbers at risk, censored and events at each time point are shown below each panel.

（n は 1 本目のコホート値。到着データで確かめてから確定する。）

## 2026-10-12 追記: 1 本目データから骨の KM を実データで描いた

- 1 本目 `~/Sync/Projects/Shimizu/op/data/fig4/enrolled_patients_20260322.xlsx` に `perioperative chemo`（0/1）があった。
  `osteoporosis_bone_meta.csv`（BMFS 月数）と 研究ID で 785 例が 1:1 対応、イベント・閉経は完全一致。
- `data/clinical/op/` にコピーし、`fig/prepare_clinical.py` で `data/clinical/clinical.xlsx` を作成（骨のみ、肺・肝は無し）。
- 化療欄が空欄の 43 例（閉経前 23 / 閉経後 20）は `km.load_clinical` で除外 → 解析 742 例。
- `uv run python fig/fig3.py` → `out/fig/fig3.*`, `fig3_stats.csv`:
  All HR 1.89 (1.20–2.99) p<0.01 / Pre HR 1.24 (0.59–2.61) p=0.57 / Post HR 2.46 (1.37–4.40) p<0.01（Cox 単変量）。
  抄録の「閉経後のみ有意」と向きが一致。抄録の HR 2.05 / Table3 の 1.72 は多変量なので一致しなくてよい。
- Fig2（肺・肝）は清水の拡張データ待ちのまま。43 例の化療情報も清水に確認。

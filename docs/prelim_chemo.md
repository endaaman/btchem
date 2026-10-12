# Fig4B 再検討: 化療 × 骨の急性期・雌のデータ（2026-10-12）

前提は `docs/prelim_de.md`（GSE205237 は DEG 9 個、GSE296757 は 11 個で、アポトーシス・老化・炎症のどれも出なかった）。
ここでは「急性期・雌・乳がんレジメン」のデータで 4B を救えるかを見た。予備段階の叩き台。

## 結論（先に）

- **急性期（Dox 96 h, 雌, GSE289491 scRNA）でも、アポトーシス・老化パネルは「はっきり上がる」とは言えない。** 各群 1 匹なので統計は出せない。
  擬似バルクで見ると、CAR 細胞で p53 標的の一部（Phlda3 +1.05, Bax +0.66）が上がるが、**Cdkn1a（p21）はむしろ −0.51、Il6 −2.90、Serpine1 −0.78、Ccl2 −3.99、Tnfsf11（RANKL）−0.63 で、老化・SASP・炎症はどれも Dox で低い側**。
  Cdkn2a（p16）はどの細胞種でもほぼ検出されない（CPM ≈ 1）。
  著者が GEO に置いた Partek の全細胞比較表でも、Cdkn1a は +0.04、Cdkn2a −0.55、Il6 −0.85、Mmp3 −0.45（Dox / Veh, log2）で、**GEO の summary にある「老化マーカー・SASP の上昇」を、公開されている処理済みデータからは再現できなかった**（論文の図は別の解析・qPCR・組織染色に基づく可能性がある。本文は未確認）。
- **GSE296758（雌, LepR+ 骨髄 MSPC, Dox 8 週後）は DEG が出る**（5 mg/kg×4: 140 個、10 mg/kg×1: 24 個）。ただし中身は**骨形成系の低下**
  （Spp1 −1.58, Dmp1 −2.00, Fam20c −1.40, Mmp13 −1.27, Mmp9 −1.91, Ngfr −3.65）と**エストロゲン標的の低下**（Pgr −2.73, Greb1 −1.43; Hallmark Estrogen Response Early NES −1.58, FDR 0.012）、
  増殖系の上昇（E2F NES +3.24, G2-M +2.69）で、**アポトーシス・老化パネルは全滅**（Cdkn1a −0.02, Bax +0.12, Fas +0.17; Hallmark Apoptosis NES −1.22, FDR 0.19）。
- したがって「化療 = 骨のアポトーシス」は、急性期 scRNA でも慢性期 bulk でも転写レベルでは支持されない。出た結果を言い換えると
  「**Dox は急性期に骨髄間質で p53 の一部標的（Phlda3/Bax）を動かすが、老化/SASP の転写は出ない。数週後には骨形成系（骨芽・骨細胞系の遺伝子）と ER 標的の低下が残る**」。
- 4B の推奨は末尾。要旨: **(c) 複数を並べる**（パネル遺伝子ヒートマップ）を主、volcano を 1 枚出すなら GSE296758 Dox×4（雌・臨床寄りのレジメン・DEG あり）。
  GSE205237 単独の volcano（a）は「出なかった」図にしかならない。急性期 scRNA への差し替え（b）は n=1 なので主図には向かない。
- 急性期 bulk（投与後 ≤ 7 日、骨系細胞、乳がんレジメン）の公開データは見つからなかった。近いのはヒト hESC 骨芽分化培養に 7 日曝露した GSE324491 だけで、5-FU・MTX では p53 / p21 / p16 / Bax / Puma が明確に上がる（in vitro の陽性対照）。CPA（プロドラッグ）では上がらない。
- GSE205237 で元論文と DEG 数が合わない件は**閾値の違い**で説明できる。論文は FDR < 0.05 のみ（FC 閾値なし）+ SVA で 31 個、予備解析を同条件で数えると 29 個で、Cfd↑ も一致する。
  論文が主張する老化マーカー（p21 / p16 / Csf2）の上昇は **qRT-PCR の結果で、RNA-seq の結果ではない**。

## 1. GSE289491（雌 12 週, Dox 10 mg/kg 96 h, 骨常在細胞 scRNA, 各 n=1）

### データの中身

- SOFT: 大腿骨+脛骨を砕いてコラゲナーゼ 2 回消化 → MACS/FACS で **CD45+ / CD71+ / Ter119+ / CD31+ を除去**した細胞を 10x 3' v3。各群 1 ライブラリ（Veh = GSM8792668, Doxo = GSM8792669）。性別は SOFT 上 Female。
  PMID 41469373（"Chemotherapy-induced CAR cell senescence disrupts bone homeostasis to induce bone loss"）。
- **CD45 除去のため破骨細胞は含まれない**（Acp5/Ctsk/Oscar/Dcstamp/Mmp9 を持つクラスタなし）。**骨細胞クラスタもない**（Sost ≈ 0、Dmp1/Phex は 43 細胞の小クラスタにわずか）。
- `GSE289491_Processed_data_GSA.xlsx` は Partek の「全細胞 Veh vs Doxo」の 1 シート（21,905 遺伝子; P 値, FDR, 比, LSMean）。細胞を反復として扱った検定（擬似反復）なので、P 値は個体差を含まない。
- **RAW.tar の `*filtered_feature_bc_matrix.txt` は生 counts ではない。** 値の段階と頻度から逆算すると、各細胞で `u = log(CPM+1)` を作った後にそれを再度 CPM 化して `log2(x+1)` した二重正規化になっていた。
  最頻段（UMI=1）と次の段（UMI=2）から細胞ごとに CPM 係数を解いて UMI 数に戻せる（`prelim/chemo/sc_counts.py`）。99.3% の細胞で逆変換値が整数（同一シンボル重複の遺伝子は 0.5 刻み）かつ総和が一致。
  残り 33 細胞も誤差 < 0.6% で、そのまま使った（`data/seq/GSE289491/qc_inverse.csv`）。
- 細胞数が **両群とも正確に 2,488** で、Partek 側で同数にダウンサンプルされている可能性がある（細胞種の割合比較は参考程度）。

### QC

| | Veh | Dox |
|---|---|---|
| 細胞数 | 2,488 | 2,488 |
| UMI 中央値 | 8,694 | 7,348 |
| 検出遺伝子中央値 | 2,728 | 2,341 |
| mt% 中央値 | 2.5 | 3.3 |

- 統合 → 正規化（1e4, log1p）→ HVG 2,000 → PCA 30 → Leiden(0.8) → 14 クラスタ → マーカー（`out/prelim/chemo/sc/cluster_marker_means.csv`, `dotplot_celltype.png`）で 9 細胞種に統合。バッチ補正はしない（1 ライブラリずつなので処置効果ごと消える）。
- 細胞種（Veh / Dox 細胞数）: CAR (Lepr+Cxcl12hiKitl+Adipoq+) 378/429、骨芽細胞 (Bglap/Ifitm5) 193/116、骨前駆 (Col1a1/Alpl/Ibsp) 267/228、肥大軟骨/骨系 (Col10a1/Dmp1/Mmp13) 24/19、軟骨細胞 411/466、Prg4+ 線維芽 638/857、Ly6a+ 線維芽 345/219、腱/靭帯線維芽 (Tnmd) 130/89、内皮 102/65。
  関節・骨膜由来の線維芽・軟骨細胞が過半を占める（骨端ごと消化したため）。
- **UMAP 上で CAR は Veh/Dox で分かれる**（他の細胞種は混ざる）。CAR が Dox に最も反応した細胞種、という論文の主張とは整合する。
- **技術差の痕跡**: Veh 側で全細胞種にヘモグロビン（Hbb-bs: Dox で −6〜−8）、好中球（Ngp, S100a8/9, Camp）、血小板（Pf4, Ppbp）の ambient RNA が多い。Dox の骨髄抑制で血球が減った結果として自然だが、
  「Veh で高い」遺伝子の一部は混入由来の可能性がある。また CAR の Veh 側で Cxcl1/Cxcl2/Ccl2/Il6/Atf3/Nfkbia（解離ストレス・NF-κB 系）が高く、Fkbp5/Dbp/Tsc22d3 が多くの細胞種で Dox 側に高い（ストレス性グルココルチコイドか採取時刻の差）。n=1 ではこれらを処置効果と区別できない。

### 結果（擬似バルク, Dox vs Veh）

方法: 細胞種ごとに UMI を合算 → CPM → log2((CPM_Dox+1)/(CPM_Veh+1))。括弧内は**同じ 1 匹から細胞を復元抽出した 500 回のブートストラップ 2.5–97.5%**で、細胞サンプリングの揺らぎだけを表す（個体差は含まないので、これが 0 をまたがなくても「有意」ではない）。
pctile は同じ細胞種の全遺伝子（平均 CPM ≥ 10）の log2FC 分布での順位（0 = 最も Dox で低い、100 = 最も高い）。全細胞種・全パネルは `out/prelim/chemo/sc/panel_pseudobulk.csv`。

**CAR 細胞（Veh 378 / Dox 429 細胞）**

| gene | CPM Veh / Dox | log2FC (boot) | 発現細胞 % Veh / Dox | pctile |
|---|---|---|---|---|
| Fas | 19 / 25 | +0.34 (−0.1, +0.8) | 12 / 16 | 77 |
| Casp3 | 34 / 27 | −0.31 (−0.8, +0.2) | 17 / 20 | 18 |
| Bax | 100 / 158 | +0.66 (+0.4, +0.9) | 46 / 59 | 93 |
| Bbc3 | 20 / 26 | +0.38 (−0.2, +0.9) | 12 / 16 | 80 |
| Pmaip1 | 1 / 1 | — | 1 / 0 | 低発現 |
| Cdkn1a | 341 / 238 | −0.51 (−0.9, −0.1) | 54 / 46 | 9 |
| Cdkn2a | 1 / 1 | — | 1 / 1 | 低発現 |
| Serpine1 | 388 / 226 | −0.78 (−1.3, −0.2) | 36 / 24 | 4 |
| Il6 | 79 / 10 | −2.90 (−4.3, −1.8) | 8 / 3 | 0 |
| Tnfsf11 | 143 / 92 | −0.63 (−1.0, −0.2) | 36 / 40 | 7 |
| Il1b | 3 / 6 | — | 2 / 1 | 低発現 |
| Ccl2 | 2840 / 178 | −3.99 (−5.0, −3.2) | 36 / 12 | 0 |
| Cxcl12 | 20860 / 24506 | +0.23 (+0.1, +0.3) | 100 / 100 | 69 |
| (p53) Phlda3 | 93 / 193 | +1.05 (+0.7, +1.4) | 34 / 60 | 98 |
| (p53) Mdm2 | 39 / 50 | +0.34 (−0.1, +0.8) | 23 / 26 | 78 |
| Mki67 | 17 / 14 | −0.27 (−1.8, +1.3) | 5 / 4 | 20 |

**他の細胞種（log2FC のみ。— は両群とも CPM < 10 程度で評価不能）**

| gene | 骨芽細胞 (193/116) | 骨前駆 (267/228) | 軟骨細胞 (411/466) |
|---|---|---|---|
| Fas | — | −0.17 | −0.32 |
| Casp3 | +0.23 | −0.93 | −0.33 |
| Bax | +0.65 | +0.06 | +0.55 |
| Bbc3 | +0.31 | −0.81 | −0.72 |
| Pmaip1 | — | — | — |
| Cdkn1a | +0.26 | −0.03 | −0.32 |
| Cdkn2a | — | — | — |
| Serpine1 | −1.23 | +1.26 | +0.43 |
| Il6 | — | — | −2.51 |
| Tnfsf11 | — | +0.77（CPM 6→11） | — |
| Il1b | — | −6.10（CPM 68→0、数細胞） | — |
| Ccl2 | — | +0.37 | −1.86 |
| Cxcl12 | −0.24 | +0.12 | −0.88 |
| Phlda3 | +0.92 | −0.07 | +0.54 |
| Mki67 | −0.89 | −1.92 | −1.18 |

- 事前定義の p53 標的スコア（19 遺伝子）・SASP スコア（18 遺伝子）の細胞平均（`scores.csv`）: p53 は CAR −0.074→−0.051、骨芽 −0.068→−0.044 とわずかに上がるだけ。**SASP は CAR で 0.42→0.28 と Dox で低い**。
- Dox 側で目立って上がるのは、CAR の**インターフェロン応答**（Ifit1 +2.25, Ifit3 +1.98, Usp18 +1.77, Oasl2 +1.51, Irf7 +1.15）と**脂肪細胞化**（Fabp4 +2.76, Plin1 +2.73, Cd36 +2.51, Car3 +1.71）、全体で**増殖マーカー（Mki67/Top2a）の低下**。
  Dox の骨髄抑制・脂肪化として既知の方向で、「急性期に CAR が反応している」こと自体は見える。
- 著者の Partek 表（全細胞, 擬似反復）: Cdkn1a +0.04、Cdkn2a −0.55、Il6 −0.85、Mmp3 −0.45、Serpine1 +0.24、Spp1 +0.28、Cstb +0.03、Mki67 −1.29、Bax +0.26、Fas −0.35、Casp3 −0.15（`out/prelim/chemo/sc/partek_panel.csv`）。

### 限界

- **各群 1 匹**。個体差が推定できないので、上の数字はすべて記述統計。どの変化も「1 匹と 1 匹の差」以上のことは言えない。
- 生 counts ではなく正規化値からの逆変換（整合性は上記のとおり確認済み）。
- 破骨細胞・骨細胞はデータに含まれない（CD45 除去・消化条件）。Dox が直接骨細胞死を起こすかはこのデータでは見られない。
- ambient RNA と解離ストレスの量が群で違う（上記）。

## 2. GSE296758（雌, Lepr-Cre;tdTomato の骨髄 MSPC, Dox 投与 8 週後, bulk SMART-Seq v4）

### データ

- SOFT: タイトルはすべて "female"。Xist 発現あり（counts 表に Y 染色体遺伝子の行自体がない）。若年成体（8–10 週）に PBS×4 / Dox 10 mg/kg i.v.×1 / Dox 5 mg/kg i.v.×4（週 1）、最終投与 8 週後に採取。GRCm39, STAR + featureCounts の raw counts（gene_name 単位）。
- 群: PBS 4（C1–C4）/ DOX1X 3 / DOX4X 4。**GSM8975978 は title が "DOX bio rep 1 female" だが、library 名 C4・treatment "4X PBS IV"** で、他の Dox 試料の title は必ず DOX1X/DOX4X を含む → PBS とした。
  PCA でも C4 は PBS 側の端に来る（`out/prelim/chemo/qc/GSE296758_pca.png`）。念のため C4 を除いた感度分析も実施。
- **flowcell**: C1–C4 と DOX1X_1/2 が HGLGLDSX7_L2、DOX1X_3 と DOX4X_1–4 が HGLNTDSX7_L1 で、**DOX4X vs PBS は flowcell と完全交絡**するため design に入れられない（`~ cond`）。
  ただし同じ flowcell の DOX1X_3 は PCA で PBS 側、別 flowcell の DOX1X_1/2 は Dox 側に来るので、PC1 を flowcell が支配しているわけではない。library ID は連番で、調製バッチは同じとみられる。
- 組成（CPM）: Lepr 4,000–7,000、Cxcl12 ~50,000、Kitl ~1,100 で LepR+ MSPC として妥当。Ptprc 12–35、Hbb-bs 15–86 と血球混入は少ないが Dox 群でやや多い。

### 結果（pydeseq2, prelim/de.py の run をそのまま使用）

| contrast | tested | padj<0.05, abs(lfc)≥1 (up/down) | padj<0.1, abs(lfc)≥0.585 |
|---|---|---|---|
| B3_Dox4x_MSPC（4x vs PBS, 4/4） | 16,643 | **69 / 71** | 211 / 205 |
| B3_Dox4x_MSPC_noC4（C4 除外, 4/3） | 16,930 | 51 / 97 | 200 / 270 |
| B3_Dox1x_MSPC（1x vs PBS, 3/4） | 17,052 | 12 / 12 | 27 / 36 |

パネル（log2FC (padj)、Dox×4 / Dox×1）。全体は `out/prelim/chemo/bulk_panel.csv`。

| gene | Dox×4 | Dox×1 |
|---|---|---|
| Fas | +0.17 (0.61) | +0.01 (1.0) |
| Casp3 | +0.34 (0.20) | +0.30 (0.85) |
| Bax | +0.12 (0.62) | +0.12 (0.96) |
| Bbc3 | −0.23 (0.59) | +0.10 (0.97) |
| Pmaip1 | +0.27 (0.78) | −0.07 (1.0) |
| Cdkn1a | −0.02 (0.98) | +0.13 (0.98) |
| Cdkn2a | 低発現 | 低発現 |
| Serpine1 | +0.04 (0.99) | +0.99 (0.86) |
| Il6 | +0.10 (0.93) | +0.11 (0.98) |
| Tnfsf11 | −0.52 (0.16) | −0.18 (0.96) |
| Il1b | −0.71 (0.43) | −0.26 (0.97) |
| Ccl2 | +0.33 (0.60) | +0.44 (0.89) |
| Cxcl12 | −0.10 (0.68) | −0.01 (1.0) |
| Mmp9 | **−1.91 (2e-4)** | **−1.26 (0.020)** |
| Mmp13 | **−1.27 (<1e-4)** | **−0.87 (0.002)** |
| Spp1 | **−1.58 (<1e-4)** | **−0.99 (1e-4)** |

- 下がった遺伝子: Ngfr −3.65, Pgr −2.73, Wnt7b −2.66, Dmp1 −2.00, Mmp9, Spp1, Fam20c −1.40, Greb1 −1.43, Ibsp −0.63, Runx2 −0.26（padj 0.017）。**骨芽・骨細胞系への分化の低下**（論文 summary の "loss of osteogenesis from MSPCs" と整合）と **ER 標的（Pgr, Greb1）の低下**。
- 上がった遺伝子: Fabp4 +1.42, Lpl +0.47（脂肪化）、B 細胞・赤血球系（Cd79a, Ms4a1, Hbb-bt）= 血球混入の増加。
- Hallmark GSEA（FDR < 0.05 抜粋）: Dox×4 で up E2F (+3.24), G2-M (+2.69), Myc, DNA Repair, IFN-γ/α, Inflammatory Response (+1.42, FDR 0.033) / down EMT (−2.93), Protein Secretion, Coagulation, Angiogenesis, TGF-β, **Estrogen Response Early (−1.58, FDR 0.012)**。
  Apoptosis −1.22 (FDR 0.19)、p53 Pathway +1.07 (0.35) は動かない。Dox×1 もほぼ同じ向き（ER Early −1.58, FDR 0.013）。
- C4 除外でも結論は同じ（Pgr −2.41, Greb1 −1.35, Dmp1 −2.32, Spp1 −1.71; Cdkn1a −0.09, Bax +0.18, Fas +0.22 はいずれも n.s.）。
- 他の化療データとの log2FC 相関（Spearman）は弱い: GSE296757 骨芽細胞（同じ研究室, 雌）0.17、GSE205237 −0.12、GSE289491 の CAR 擬似バルク 0.13。
- **ER 標的の低下**は、Dox による卵巣機能低下（化療誘発性閉経）を示唆し、本論文の「化療と閉経」の両軸をつなぐ所見になりうる。ただし同じ研究室の GSE296757（骨芽細胞）では Greb1 は +1.39 (n.s.)、ER Early は NES −1.37 (FDR 0.093) と弱く、再現は不十分。血中 E2 や卵巣のデータはこの GEO にはない。

## 3. 急性期 bulk の追加探索 と GSE205237 元論文の照合

### 検索

E-utilities（db=gds の統合検索式＋薬剤 8 語 × 骨関連語 17 語の総当たり、db=sra、db=bioproject）と ArrayExpress で検索した。
薬剤は doxorubicin / adriamycin / cyclophosphamide / 4-HC / mafosfamide / paclitaxel / docetaxel、骨関連語は bone / osteoblast / osteocyte / MLO-Y4 / IDG-SW3 / MC3T3 / calvaria / BMSC / MSC / osteoclast / RAW264.7 など。
`docs/datasets.md` で既に評価済みの 11 本は除いた。ヒットした series は SOFT を読んで確認した。

**条件（乳がんレジメンの薬剤・骨系細胞・投与後 ≤ 7 日・counts あり・n ≥ 2）をすべて満たすデータは見つからなかった。** 部分的に合うもの:

| accession | 中身 | 時間 | n | 処理済み | 評価 |
|---|---|---|---|---|---|
| **GSE324491** | ヒト H9 hESC を骨芽方向へ分化させた培養に、9 種の毒性物質を IC50 で連続曝露。うち化療薬は **シクロホスファミド (CPA), 5-FU, MTX** | day 0–7（day 7 回収） | 各 3 | raw counts（featureCounts） | ○ 唯一 DE を回せた。in vitro・胚性細胞・IC50 の発生毒性モデルなので、生体の骨とは距離がある |
| GSE290573 | マウス RAW264.7（破骨前駆）に Dox | 24 h | 5 / 5 | FPKM のみ（counts なし） | △ 骨系ではなく破骨前駆。濃度表記が SOFT 内で 10 nM / 100 nM と食い違う |
| GSE17839 / GSE17842 | ヒト骨芽細胞株に melphalan 50 µg/mL（または melphalan 処理した間質の上清） | 6 h | 4 / 4 | Agilent 2 色の log 比 | △ 急性の骨芽細胞だが、薬剤が乳がんレジメン外 |
| GSE284666 / GSE285826 | ヒト BMMSC・骨芽細胞を etoposide で老化誘導（Dox は血管細胞のみ） | 8 日 | 4 | RAW | △ 参考 |

見つからなかった組み合わせ: Dox / CPA / タキサンを骨芽細胞・骨細胞（MC3T3, MLO-Y4, IDG-SW3, 頭蓋冠）に 7 日以内で当てたデータ、4-HC・マフォスファミド（CPA の活性体）× 骨系細胞、
タキサン × 正常な骨系細胞、in vivo で 7 日以内に骨を採ったデータ（GSE289491 を除く）。

### GSE324491 の DE（同じパイプライン）

- 対照列が 6 本（UNT_I–III と UNT_IA–IIIA）あるが GSM は 3 本。PCA の座標がほぼ一致し、ライブラリサイズが約半分なので、**A 系列は同じライブラリの再シーケンス（技術的反復）**とみられる。
  主解析は UNT_I–III、感度分析で UNT_IA–IIIA を使い、DEG 数は同程度（CPA 1,731 → 1,479、5-FU 3,589 → 3,104、MTX 7,380 → 6,832）。
- シンボルは GEO の xlsx の hgnc_symbol 列を使った（26,578 遺伝子で txt の counts と完全一致）。パネルと照合するため頭文字だけ大文字に変換（CDKN1A → Cdkn1a）。TP53 は Tp53 の表記になる。

| gene | CPA | 5-FU | MTX |
|---|---|---|---|
| Cdkn1a | +0.58 (0.083) | **+2.09** | **+3.31** |
| Cdkn2a | +0.36 (0.76) | **+2.34** | **+6.00** |
| Bax | +0.33 (0.50) | **+0.83** | **+1.74** |
| Bbc3 | −0.28 (0.58) | **+0.79** | **+1.83** |
| Fas | −0.03 (0.98) | +0.58 (0.12) | **+1.96** |
| Casp3 | −0.18 (0.75) | −0.16 (0.55) | **+1.32** |
| Phlda3 / Gdf15 / Eda2r | n.s. | **+1.76 / +1.73 / +1.35** | **+2.53 / +3.07 / +0.84** |
| Mki67 | −0.27 (0.77) | **−0.97** | **−2.60** |
| Sp7 | **−1.75** | **−3.31** | **−4.92** |
| Tnfsf11 / Tnfrsf11b | **−1.82 / −2.60** | **−1.90 / −2.85** | −0.81 / **−1.22** |
| Serpine1 | **−2.22** | **−1.19** | **+3.05** |

太字は padj < 0.05。Hallmark GSEA: p53 Pathway は 5-FU で NES +2.09、MTX で +2.77、Apoptosis は +1.52 / +2.20 で、いずれも FDR < 0.02。
CPA は p53 −1.18 (0.26)、Apoptosis −1.54 (0.017)。CPA は肝で活性化されるプロドラッグで、培養系ではほとんど活性体にならない。
CPA でも DEG が 1,731 個出て Sp7 や RANKL/OPG が下がるのは、p53 経路を通らない分化阻害（IC50 で曝露しているため）とみられる。

→ **活性のある代謝拮抗薬を骨芽系の細胞に直接当てると、7 日で p53 → p21/p16・Bax/Puma の典型的な DNA 損傷応答と、増殖・骨芽分化（Sp7）の低下がはっきり出る。** これは in vitro の陽性対照として扱える。
ただし Dox ではなく、胚性細胞の分化培養で、濃度も IC50 なので、乳がん患者の骨に外挿できるデータではない。

### GSE205237 元論文との不一致の説明

論文: Pimenta-Lopes C, …, Ventura F. *Mol Ther* 2023;31(8):2507–2523（PMID 37143324, PMC10422003; "Inhibition of C5AR1 impairs osteoclast mobilization and prevents bone loss"）。PMC の本文を読んで確認した。
Table S1（DEG 一覧）は PMC のボット対策で取得できず、未確認。

- 解析条件（Methods）: "DESeq2/1.18 using the default parameters. Unwanted sources of variation were removed with surrogate variable analysis using the sva R package. Genes with FDR < 5% were considered significant."
  → **FC の閾値なし**、SVA で補正、低発現フィルタの記載なし。
- 結果（Results）: "RNA-seq analysis showed **31 protein-coding genes** differentially expressed (FDR < 0.05)"。本文で名前を挙げた RNA-seq の DEG は **Cfd（上昇）だけ**。
  Discussion には "we did not find alterations in the expression levels of Csf1, RANKL (Tnfsf11), or OPG (Tnfrsf11b)" とある。C3 / C4b / Cfh は本文に出てこない。
- **老化マーカー（p21, p16, SASP の Csf2）の上昇と Cfd / Cfb の上昇は qRT-PCR の結果**（Fig 1B, 3D, n = 4–8）で、RNA-seq の結果ではない。
- 予備解析を同じ閾値（padj < 0.05、FC 閾値なし）で数え直すと **29 個（上昇 7 / 下降 22）**で、論文の 31 個とほぼ同じ。**Cfd も +1.32（padj 0.045）で論文と一致**する。
  予備解析の 9 個は |log2FC| ≥ 1 を課した数。**つまり不一致は閾値の違いで説明でき、データ解析そのものは論文と矛盾しない。**
  一方、Cdkn1a は −0.24 (0.97)、Cdkn2a は +0.63 (0.96) で、論文の qPCR が示す老化マーカーの上昇は RNA-seq には現れていない。これは論文自身の RNA-seq とも食い違っていない（論文も RNA-seq で老化を主張していない）。
- 残る差の候補は SVA（予備解析では使っていない）と低発現フィルタ。DEG 数がほぼ一致しているので、ここは追わなかった。

## 4. 4B の推奨

並べた図は `out/prelim/chemo/panel_heatmap.png`（左が GSE289491 の擬似バルク、右が bulk DESeq2）。

| 案 | 中身 | 長所 | 短所 |
|---|---|---|---|
| (a) GSE205237 のまま | 雄、Dox 単回、30 日後の骨 | 骨細胞を含む骨組織。元論文と DEG 数が一致し、解析の妥当性は示せる | DEG 9（FC 閾値なしでも 29）。アポトーシス・老化パネルは全滅。雄 |
| (b) 急性期 scRNA（GSE289491）に差し替え | 雌、Dox 96 h、CAR などの細胞種別 | 急性期・雌・細胞種別に見られる | **各群 n=1 で volcano（P 値）が描けない**。老化・SASP はむしろ低い側。破骨細胞と骨細胞が含まれない |
| (b') GSE296758 Dox×4 に差し替え | 雌、臨床寄りの反復投与、LepR+ MSPC、8 週後 | n=4/4 で DEG 140。骨形成系と ER 標的の低下という読める結果 | 骨細胞ではなく骨髄の間質前駆。アポトーシスは出ない。flowcell と群が交絡 |
| (c) 複数を並べる | パネル遺伝子のヒートマップ（+ volcano 1 枚） | 「急性期でも慢性期でもアポトーシスは出ない」を正直に一枚で示せる。in vitro の陽性対照（hESC-OB の 5-FU）と対比できる | 図が説明的になる。データごとに組織・性・時期が違う |

**推奨: (c)。** 主パネルは化療データを並べたパネル遺伝子ヒートマップ。volcano を 1 枚残すなら GSE296758 Dox×4（雌・反復投与・DEG あり）を 4B に、GSE205237 は Supplementary に回す。
本文は「Dox は骨髄間質に急性の p53 応答の一部（Phlda3 / Bax）と IFN 応答、増殖低下を起こすが、老化・SASP・アポトーシスの転写は生体の骨では確認できなかった。
数週後には骨形成系遺伝子と ER 標的の低下が残る」とし、清水仮説の「化療 = アポトーシス」は転写データからは支持されないと書く。

- (b) を主図にするのは勧めない。n=1 であることに加え、著者自身の公開表でも老化マーカーが上がっていない。
- 判断を ken / 清水に仰ぎたい点:
  1. ER 標的の低下（化療誘発性の卵巣機能低下を示唆）を、閉経軸（4A）とつなぐ所見として使うか。同じ研究室の骨芽細胞（GSE296757）では再現が弱い。
  2. in vitro の hESC-OB（5-FU で p53 / p21 / p16 が明確に上がる）を「活性のある化療薬に直接曝露すれば出る」陽性対照として Supplementary に入れるか。
     入れる場合も、Dox ではなく、胚性細胞・IC50 である点は明記が要る。

## 生成物

- スクリプト（`prelim/chemo/`）: `sc_counts.py`（GSE289491 の逆変換）, `sc_cluster.py`（クラスタリング）, `sc_pseudobulk.py`（細胞種付与・擬似バルク・スコア・Partek 表抽出）,
  `bulk.py`（GSE296758 / GSE324491 の counts 作成・QC・DE・volcano・GSEA・パネル表）, `panel_heatmap.py`（化療データ横並びのパネルヒートマップ）。`sc_peek*.py` は値の形式を調べた使い捨て。
- データ: `data/seq/GSE289491/`（SOFT, xlsx, RAW.tar, `raw/`, `counts_{Veh,Doxo}.npz`, `cells_*.txt`, `genes.txt`, `qc_inverse.csv`）, `data/seq/GSE296758/`, `data/seq/GSE324491/`（SOFT, 補足ファイル, `counts.tsv.gz`, `samples.csv`）
- 出力（`out/prelim/chemo/`）: `sc/`（umap.png, dotplot_celltype.png, ma_*.png, panel_pseudobulk.csv, pseudobulk_all_genes.csv.gz, celltype_counts.csv, scores.csv, partek_panel.csv, cluster_*.csv, adata_unlabeled.h5ad）,
  `de/`, `volcano/`, `gsea/`, `qc/`, `bulk_panel.csv`, `bulk_deg_counts.csv`, `panel_heatmap.png|pdf|csv`
- 依存追加: `uv add scanpy igraph leidenalg`（pyproject.toml / uv.lock が変わった）

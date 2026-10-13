# 化療 × niche モジュール: 加齢・OVX・DEX と同じ gene set で並べる（2026-10-13）

前提: `docs/prelim_de.md`（GSE205237 Dox 単回は DEG 9 個）、`docs/prelim_chemo.md`（GSE296758 / GSE289491）、
`data/20261013_meeting.md`（4A/4B を主軸に、4C の DEX は本編から外す方向）。予備段階の叩き台。
再生成: `uv run python -I prelim/chemo2/gsea.py` → `uv run python -I prelim/chemo2/dox_deg.py`（repo root から）。

## 結論（先に）

1. **Dox 単回（GSE205237）の DEG 9 個は「ノイズ」と切ってよいものが大半だが、全部ではない。**
   - B 細胞・形質細胞（Igkv4-70, Ighg3）とマクロファージ（Timd4, Cd163）の 4 個、および Ccdc148 / Gm44861 は、
     所属するどの set も FDR < 0.05 で動いておらず（孤立ヒット）、サンプル間で血液・マクロファージ量と相関する（r = 0.5–0.8）。
     Ig 2 個は組成（筋・血液）を共変量に入れると有意でなくなる（padj 0.14 / 0.15）。**→ 組成の揺れで説明できる孤立ヒット**。
     Timd4 / Cd163 / Ccdc148 / Gm44861 は補正後も padj < 0.05 で残るが、モジュールの裏付けはない。
   - **Cyp2e1↑ は孤立ヒットではない**: GO lipid modification（NES +2.00, FDR 0.010）、fatty acid metabolic process（+1.96, 0.015）、
     Hallmark xenobiotic metabolism（+1.55, 0.026）、GO xenobiotic metabolic process（+1.83, 0.037）の leading edge にあり、組成補正後も残る（padj 9e-5）。
     骨髄脂肪細胞系の上昇（Hallmark adipogenesis +1.89）と同じ方向。
   - Cilp↑ は NABA core matrisome（+1.97, FDR < 0.001）の leading edge。Cfd↑ は Hallmark myogenesis（+2.53）の leading edge で、
     組成補正で有意でなくなる（padj 0.25）＝筋混入量の影響を否定できない。
   - **「DEG はノイズ」より正確には「DEG はほぼ組成差。モジュール単位では一貫した動きがあり、それも半分は組成差」**。
     組成補正前は Dox の NES ベクトルが加齢と逆相関（Spearman −0.39）するが、補正後は −0.08 でほぼ無相関になる。
2. **化療で動く niche モジュール**（FDR < 0.05、全化療データに共通なのは 1 つだけ）:
   - **IFN-α / IFN-γ 応答 ↑ は、bulk の Dox 4 データすべて（単回骨・組成補正後も・MSPC×4・MSPC×1・OB×4）と急性期 CAR で同じ向き**（NES +1.3〜+2.5）。
     加齢でも ↑（+2.6 / +2.7）、OVX では逆に ↓、DEX では動かない（−0.8; bg 補正版で IFN-γ のみ +1.40）。**化療と加齢に共通で、DEX には無い唯一の niche 系モジュール**。
   - LepR+ MSPC（GSE296758, ♀, 8 週後）: 骨形成（ossification −2.13, osteoblast diff. −1.93, mineralization −1.95）、ECM（matrisome −2.50, collagen fibril −2.80）、
     血管（angiogenesis −2.16, vasculature dev. −1.88, sprouting −1.58）、EMT・TGF-β ↓、E2F / G2M ↑、ER 応答 ↓。Dox×1 でもほぼ同じ。
   - 単回骨（GSE205237）: 炎症（inflammatory response −1.81, TNF-α/NF-κB −2.20）↓、手製 HSC retention パネル ↓（−1.68）、E2F/G2M ↓、
     骨形成・ECM はむしろ ↑（+1.6〜+2.0）、破骨マーカーパネル ↑（+1.57、組成補正後 +1.88）。組成補正後も炎症↓・TNF↓・ECM↑・ossification↑ は残る。
   - **老化/SASP 系（SenMayo, Reactome SASP, 手製 senescence パネル）はどの化療 bulk でも ↑ にならない**（SenMayo は −1.3〜−1.4 で弱く ↓）。
     Reactome cellular senescence だけ MSPC で ↑（+1.58 / +2.02; 中身はヒストン・増殖関連が多い）。
3. **加齢と化療の向き**: データで割れる。
   - **MSPC の Dox（×4 / ×1）は加齢とほぼ同じ向き**（niche + Hallmark 72 set の NES の Spearman 0.69 / 0.76）。両方で有意かつ同方向 25 / 24 set、逆向き 2 / 0。
     共通: 骨形成・ECM・血管・EMT ↓、E2F/G2M・IFN ↑。**違い: 炎症（IL6/JAK/STAT3, 手製炎症パネル）・SenMayo は MSPC では ↑ にならない**（加齢で ↑ の炎症は化療に無い）。
   - **単回骨の Dox は加齢と逆向き**（−0.39）。炎症・TNF・補体・myeloid activation ↓、骨形成・ECM・myogenesis ↑ と、加齢の主要軸がすべて反転。
     ただし組成補正で −0.08 まで弱まるので、逆向きの多くは筋・血液の量の差。
   - 交互作用（加齢 × 化療）はどのデータにも無いので検定していない。表はモジュールの向きの並置にとどまる。
4. **Dox 単独が「対照」として働くと言える根拠**（= 化療単独では加齢型の niche は作れない、という主張の根拠）:
   - 言える: 若年骨への Dox（単回 30 日、×4 / ×1 の 8 週後）で、**加齢で上がる炎症モジュール（Inflammatory response, IL6/JAK/STAT3, TNF-α, 手製炎症パネル）・myeloid activation・破骨分化の GO set は、どの Dox bulk でも有意に ↑ しない**（単回骨ではむしろ ↓、組成補正後も ↓）。
     老化/SASP も上がらない。「化療単独では加齢骨の炎症性 niche は再現されない」は 4 データで一貫。
   - 言えない: 「化療単独では何も起きない」。MSPC では骨形成・ECM・血管が加齢と同じ向きで落ち、IFN 応答は全データで上がる。
     **化療単独は "炎症なしで骨形成 niche を削る + IFN を上げる" 方向に働く**、が正確。臨床の「加齢背景 + 化療で骨転移」に乗せるなら、
     加齢（炎症・破骨方向）と化療（骨形成・血管 niche 低下）が**別の軸を削る**という 2-hit の言い方になり、「単独では動かない」とは書けない。
5. **DEX を加齢の差分として使う価値（追加依頼）**: **低い。Supplementary 行きを推奨。** 詳細は §4。
   DEX の NES ベクトルは加齢と Spearman 0.87（bg 補正版でも 0.78）で、加齢で有意な 43 set のうち 29 set で DEX も同方向に有意。
   加齢から DEX を差し引くと加齢シグナルの大半が消え、残るのは IFN-α/γ・Myc・DNA repair・heme metabolism などに限られる。
   しかも重なりの中心（ECM・骨形成・myogenesis ↓、myeloid ↑）は組織組成でも同じ形になる軸なので、「どんな侵襲でも出る汎用炎症」と「組成差」を分離できない。

## 方法

- gene set（`prelim/chemo2/genesets.py`, 事前定義）: MSigDB mouse **2025.1.Mm** の GMT を Broad から直接取得（`out/prelim/chemo2/msigdb/`）。
  mouse 版なのでオーソログ変換なし（照合は大文字シンボル）。gseapy の `Msigdb` クラスは lxml が無く動かなかったので使っていない（依存は増やしていない）。
  - Hallmark 50 全部 + niche 関連 17 set（SenMayo = `SAUL_SEN_MAYO`, Reactome SASP / cellular senescence, GO +reg. osteoclast diff. / osteoclast diff. / bone resorption,
    `WP_OSTEOCLAST_SIGNALING`, `NABA_CORE_MATRISOME`, GO collagen fibril org., GO vasculature dev. / sprouting angiogenesis, GO fat cell diff., GO ossification / osteoblast diff. / bone mineralization,
    GO myeloid leukocyte activation, GO HSC homeostasis）。
  - `GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION` は mouse 2025.1 に存在しない（GO 側の改名）ため、ECM は NABA core matrisome と collagen fibril organization で代替。
  - 手製 niche パネル（役割別サブパネル。OPG のように向きが逆の因子は混ぜない）: HSC retention（Cxcl12, Kitl, Angpt1, Lepr, Vcam1, Spp1）、
    osteoclast markers（Ctsk, Acp5, Dcstamp, Mmp9, Mmp13）、inflam./chemokine（Il6, Il1b, Tnf, Ccl2, Cxcl1, Cxcl2, S100a8, S100a9）、
    adipocyte（Pparg, Cebpa, Adipoq, Lpl）、senescence（Cdkn1a, Cdkn2a, Serpine1, Igfbp3）。RANKL/M-CSF、OPG、Sost/Dkk1 は 1–2 遺伝子なので遺伝子単位の図だけ。
    `prelim/genes.py` と重なる遺伝子はそちらと同じ名前で使い、新規は Kitl, Angpt1, Lepr, Vcam1, Spp1, Ctsk, Acp5, Dcstamp, Cxcl1, Cxcl2, Pparg, Cebpa, Adipoq, Lpl, Igfbp3。
- 順位: bulk は既存 DE 表の Wald 統計量（log2FC / lfcSE、既存 `enrich.py` と同じ）。scRNA（GSE289491, n=1）は細胞種別擬似バルクの log2FC（平均 CPM ≥ 10）で、統計量ではない。
- gseapy prerank（1000 perm, seed 0, min_size 4, max_size 1000）。FDR は 72 set のコレクション全体で計算（既存の Hallmark 単独 GSEA とは FDR が少し変わる）。
- 組成補正版 **B_Dox_compadj**: GSE205237 を `~ muscle + blood + cond` で再計算（muscle = Acta1/Ckm/Myh4、blood = Hbb-bs/Hba-a1 の log2 CPM 平均を z 化。de.py と同じフィルタ）。n = 4/4 に共変量 2 つなので残差自由度は 4。
- 問い 1: B_Dox の順位で GO BP 全体（size 15–500, 3,953 set）の prerank も回し（FDR < 0.05 は 172 set）、DEG ごとに所属 set・leading edge・向きを照合。

| 列 | データ | 性 | 中身 |
|---|---|---|---|
| Age_2.5y / Age_2y | NIHMS1776347 | ♂ | 皮質骨, 著者 DESeq2 |
| A_OVX_cort / trab | GSE222752 | ♀ | OVX vs Sham |
| B_Dox (+compadj) | GSE205237 | ♂ | 骨（骨髄除去）, Dox 単回 30 日 |
| B3_Dox4x / Dox1x_MSPC | GSE296758 | ♀ | LepR+ MSPC, 最終投与 8 週後 |
| B2_Dox4x_OB | GSE296757 | ♀ | Col2.3+ 骨芽細胞, 8 週後 |
| sc_CAR / Osteoprogenitor / Osteoblast | GSE289491 | ♀ | Dox 96 h scRNA 擬似バルク, n=1 |
| C_DEX / C_DEX_bgadj | GSE276391 | ♂ | 脛骨, DEX 5 日（背景リード補正版も） |

## 1. Dox 単回の DEG 9 個（`out/prelim/chemo2/dox_deg/deg_summary.csv`）

| gene | log2FC | 所属 set 数 | うち FDR<0.05 | 同方向かつ leading edge | 組成補正後 log2FC (padj) | r 血液 / r マクロファージ | 判定 |
|---|---|---|---|---|---|---|---|
| Cyp2e1 | +1.86 | 18 | 4 | 4（lipid modification, fatty acid metab., xenobiotic ×2） | +1.75 (9e-5) | −0.40 / −0.15 | **モジュールあり・組成に依らない** |
| Cfd | +1.32 | 12 | 2 | 1（Hallmark myogenesis）; coagulation とは逆 | +1.39 (0.25) | −0.12 / +0.25 | 筋混入と区別できない |
| Cilp | +1.01 | 14 | 1 | 1（NABA core matrisome） | +0.83 (0.021) | −0.49 / −0.54 | ECM モジュールと一致 |
| Gm44861 | −1.05 | 0 | 0 | 0 | −1.10 (0.021) | +0.47 / +0.48 | 孤立 |
| Cd163 | −1.19 | 2 | 0 | 0（acute inflam. response FDR 0.24 の LE） | −0.96 (0.015) | +0.71 / +0.74 | 孤立（マクロファージ量） |
| Ccdc148 | −1.41 | 0 | 0 | 0 | −1.38 (0.017) | +0.51 / +0.48 | 孤立 |
| Ighg3 | −1.56 | 12 | 0 | 0（humoral immunity 系 FDR 0.05–0.2） | −1.04 (0.15) | +0.76 / +0.80 | 孤立（形質細胞・血液量） |
| Timd4 | −2.07 | 4 | 0 | 0（phagocytosis engulfment FDR 0.085） | −1.82 (0.034) | +0.61 / +0.59 | 孤立（常在マクロファージ） |
| Igkv4-70 | −4.15 | 0 | 0 | 0 | −3.90 (0.14) | +0.25 / +0.08 | 孤立（形質細胞, 1 クローン） |

- 「同じ set の他遺伝子」の統計量平均（`deg_membership.csv` の `mean_stat_others`）も、FDR<0.05 の set を持たない 6 個では 0 付近〜−0.4 で、DEG だけが突出している。
- r は 8 サンプルの log2 CPM と組成スコアの Pearson（n=8 なので目安）。
- 結論は冒頭 1 のとおり。下がった 6 個は免疫細胞量（マクロファージ・形質細胞）のマーカーで、骨細胞の転写変化としては読めない。

## 2. 化療で動く niche モジュール（`heatmap_niche.png`, `nes_wide.csv`）

NES（** FDR<0.05）。抜粋。

| set | Age 2.5y | DEX | Dox 単回骨 | 同 組成補正 | MSPC ×4 | MSPC ×1 | OB ×4 | CAR 96h (n=1) |
|---|---|---|---|---|---|---|---|---|
| H inflammatory response | +2.19** | +1.76** | −1.81** | −1.71** | +1.38 | +1.28 | −0.99 | −2.02 |
| H IL6/JAK/STAT3 | +2.39** | +2.19** | −1.17 | +1.06 | −1.46** | −1.21 | +0.96 | −2.05 |
| H TNFα/NF-κB | +1.63** | +1.14 | −2.20** | −2.75** | +1.05 | +1.39 | +1.75** | −2.86 |
| H IFN-α response | +2.59** | −0.83 | +1.96** | +2.45** | +1.42 | +1.38 | +2.16** | +2.10 |
| H IFN-γ response | +2.75** | −0.79 | +1.34 | +1.91** | +1.75** | +1.60** | +2.01** | +1.40 |
| GO myeloid leukocyte activation | +2.66** | +2.21** | −1.32 | +1.23 | +1.16 | +1.19 | −1.09 | −1.60 |
| GO +reg. osteoclast diff. | +1.60** | +1.50** | −1.00 | +1.18 | +0.88 | +0.94 | −0.93 | −1.19 |
| Panel inflam./chemokine | +1.72** | +1.73** | −1.17 | −1.30 | +0.55 | +1.57** | — | −2.10 |
| Panel osteoclast markers | +0.62 | +1.25 | +1.57** | +1.88** | −1.57** | −1.41 | +0.78 | — |
| SenMayo | −1.28 | −1.10 | −1.34 | −1.41 | −1.38 | −1.37 | +0.85 | −2.60 |
| Panel HSC retention | +0.74 | +0.52 | −1.68** | −1.58 | −1.38 | −1.38 | −0.87 | +0.71 |
| GO ossification | −2.36** | −2.50** | +1.64** | +1.63** | −2.13** | −1.58** | −1.34 | −1.64 |
| NABA core matrisome | −3.37** | −2.74** | +1.97** | +1.88** | −2.50** | −2.14** | −1.07 | −2.01 |
| H angiogenesis | −2.35** | −1.79** | +0.77 | −1.01 | −2.16** | −2.02** | −1.49 | −1.99 |
| GO vasculature dev. | −2.21** | −1.89** | +1.31 | +1.16 | −1.88** | −1.97** | −1.58** | −1.80 |
| H adipogenesis | −1.38 | −1.83** | +1.89** | +0.88 | +1.24 | −0.90 | +0.95 | +1.83 |
| H E2F targets | +2.91** | +1.39 | −2.53** | −0.99 | +3.27** | +3.22** | +0.87 | +0.85 |
| H myogenesis | −2.40** | −2.88** | +2.53** | +1.72** | −1.85** | −1.82** | −0.97 | −1.30 |

- 単回骨の adipogenesis ↑・E2F ↓・OXPHOS ↑ は組成補正で FDR ≥ 0.05 に落ちる（組成依存）。炎症↓・TNF↓・IFN↑・ECM↑・ossification↑・破骨マーカー↑ は残る。
- scRNA 列は n=1 の log2FC 順位で、FDR は参考。Osteoblast（116/193 細胞）は低カウントで順位の同値が多く（65.8%）、ほぼ何も出ない。
- 遺伝子単位（`heatmap_panel_genes.png`）: MSPC で Spp1 / Mmp9 / Mmp13 ↓、単回骨で Acp5 ↑（組成補正後）、加齢で Il1b / Tnf / S100a8/9 ↑・Tnfsf11 ↓・Igfbp3 ↓。
  Cxcl12 / Kitl / Lepr / Vcam1 / Angpt1 はどの化療データでも動かない（HSC retention パネルの ↓ は個々の遺伝子では有意でない）。

## 3. 加齢と化療の向き（`concordance.csv`, `concordance_counts.csv`）

72 set の NES ベクトルの Spearman（vs Age_2.5y）と、両方 FDR<0.05 の set 数。set 同士は遺伝子を共有するので独立ではなく、記述用。

| 化療 contrast | Spearman | 同方向 | 逆方向 | 加齢のみ | 化療のみ |
|---|---|---|---|---|---|
| MSPC ×1 | +0.76 | 24 | 0 | 19 | 4 |
| MSPC ×4 | +0.69 | 25 | 2 | 16 | 8 |
| OB ×4 | +0.42 | 9 | 2 | 31 | 3 |
| CAR 96h (n=1) | +0.09 | 14 | 10 | 19 | 12 |
| Dox 単回骨 | −0.39 | 2 | 13 | 28 | 5 |
| Dox 単回骨 組成補正 | −0.08 | 3 | 6 | 34 | 3 |

- **加齢と同方向**（複数の化療データで）: 骨形成・ECM・血管・EMT・TGF-β ↓、IFN-α/γ ↑、E2F/G2M・Myc・DNA repair ↑（MSPC）。
- **加齢と逆方向**: 炎症（Inflammatory response, TNF-α, 補体）は単回骨・CAR で ↓、IL6/JAK/STAT3 は MSPC×4 で ↓。破骨 GO set は OB で ↓。
- 加齢の炎症 ↑ を共有する化療データは無い（MSPC×1 の手製炎症パネル +1.57 が唯一の例外で、Hallmark では FDR ≥ 0.05）。
- MSPC と加齢の一致は、片方がソートした LepR+ 細胞、もう片方が皮質骨 bulk なので、同じ細胞で同じことが起きているとは限らない。
  骨形成・ECM の ↓ は「骨形成系の細胞の割合が減る」でも同じ形になる。

## 4. DEX を加齢の差分として使えるか（`classify_aging_dex_chemo.csv`）

各 set を加齢（Age_2.5y）・DEX（C_DEX）・化療（bulk Dox 4 本のどれか）の FDR<0.05 と向きで分類（化療の向きは有意なものの多数決、割れたら mixed）。
DEX_bgadj で同じ判定になるかも列 `dex_robust` に残した。

| 分類 | set 数 | 主な中身 |
|---|---|---|
| 加齢∩DEX（同方向）∩化療（同方向） | 12 | EMT・TGF-β ↓, 血管（angiogenesis, vasculature dev., sprouting, hypoxia）↓, bone mineralization ↓, **手製炎症パネル ↑**（化療は MSPC×1 のみ） |
| 加齢∩DEX（同方向）∩化療（逆方向） | 4 | **Inflammatory response ↑, IL6/JAK/STAT3 ↑**（化療は ↓）, OXPHOS ↓, fatty acid metab. ↓（化療は ↑） |
| 加齢∩DEX（同方向）∩化療（mixed） | 7 | collagen fibril, matrisome, ossification, osteoblast diff., myogenesis ↓, allograft rejection・G2M ↑（単回骨と MSPC で向きが逆） |
| **加齢∩DEX のみ**（汎用炎症/侵襲） | 6 | **myeloid leukocyte activation ↑, +reg. osteoclast diff. ↑**, ROS ↑, HSC homeostasis ↓, fat cell diff. ↓, Wnt/β-catenin ↓ |
| **加齢∩化療**（DEX なし） | 12 | **IFN-α / IFN-γ ↑**, Myc v1/v2 ↑, DNA repair ↑, coagulation ↓, KRAS up ↓, apical surface ↓（同方向 8）; complement（逆）; E2F・TNF-α・heme（mixed） |
| 加齢のみ | 2 | PI3K/AKT/mTOR ↑, GO osteoclast diff. ↑（いずれも DEX_bgadj では ↑ になる → 実質 加齢∩DEX） |
| DEX のみ | 3 | Hallmark apoptosis ↓, UPR ↓, pancreas β ↓ |
| DEX∩化療 | 3 | adipogenesis（DEX ↓, 化療 ↑）, androgen response ↓, hedgehog ↓ |
| 化療のみ | 12 | ER early ↓, cholesterol homeostasis ↓, xenobiotic・bile acid ↑（単回骨）, 手製 HSC retention ↓・adipocyte ↑, Reactome cellular senescence ↑（MSPC）ほか |
| なし | 11 | SenMayo, Reactome SASP, p53 pathway, WP osteoclast signaling, 手製 senescence ほか |

判断材料:
- **DEX は加齢とほぼ同じ形**（NES Spearman 0.87）。加齢で有意な 43 set のうち、DEX と同方向に有意なのが 29、DEX と逆向きは 0。
  「加齢 − DEX」で残るのは **加齢∩化療の行（IFN-α/γ, Myc, DNA repair）と、加齢のみの 2 set（bg 補正版では DEX にも出る）だけ**。
- 加齢と DEX が共有するのは、炎症・myeloid・破骨誘導 ↑ と、ECM・骨形成・myogenesis ↓。後者は組織組成（筋・結合組織の混入量、骨髄比率）の差でも同じ形になり、
  `docs/prelim_de.md` の問題点 1 で指摘した組成の軸と区別できない。DEX データは背景リード artefact も抱える。
- したがって「DEX は加齢から差し引く汎用炎症の差分」として使うと、**加齢の炎症・破骨シグナルをほぼ全部消してしまい**、残りは IFN/増殖系だけになる。
  逆に「加齢の炎症・破骨方向は、DEX という別の侵襲でも同じ形で出る非特異的な応答」と書く材料にはなるが、それは 4A の主張を弱める方向。
- **推奨: DEX は Supplementary**（会議の方針と一致）。差分として使うなら、加齢∩化療で DEX に無い **IFN-α/γ 応答** を「加齢と化療が共有し、ステロイドには無い軸」として示す用途に限る。

## 根拠が弱い点

- **組成**: bulk の骨はどれも筋・血液・骨髄の混入量で NES が大きく動く。組成補正は GSE205237 だけ、マーカー 2 軸・残差自由度 4 の粗い補正。加齢（著者 DESeq2 表）と DEX は補正していない。
  加齢∩DEX∩化療の ECM・骨形成・myogenesis の重なりは組成で説明できる可能性が高い。
- **IFN 応答**: bulk 骨では免疫細胞量の差でも上がる。ただしソートした MSPC・OB・CAR でも上がるので、細胞内在の応答の可能性はある（Dox の DNA 損傷 → I 型 IFN は既知の経路として知られるが、文献は未確認・要引用）。OVX では逆向き。
- 性・組織・時期がデータごとに違う（加齢 ♂ 皮質骨、MSPC ♀ ソート細胞、単回骨 ♂）。「加齢と同方向」は異なる細胞集団どうしの比較。
- scRNA は n=1、順位は log2FC で同値が多い（gseapy の重複警告: CAR 19.5%, Osteoprogenitor 34.2%, Osteoblast 65.8%）。FDR は参考値。
- 手製パネルは 4–8 遺伝子の GSEA で、FDR の安定性は低い。低発現で 4 遺伝子未満になった contrast は空欄（senescence パネルはほとんど空欄）。
- FDR は 72 set のコレクションで計算しており、既存の Hallmark 単独 GSEA（`out/prelim/enrich/`）とは値が少し違う（例: B_Dox の TNF-α/NF-κB は −2.17 → −2.20）。
- 交互作用（加齢 × 化療）は検定していない。加齢マウスに化療を当てた公開データは無いので、向きの並置以上のことは言えない。

## 生成物

- スクリプト（`prelim/chemo2/`）: `genesets.py`（gene set 定義・GMT 読み込み）, `gsea.py`（組成補正 DE・GSEA・ヒートマップ・一致表・分類表・パネル遺伝子図）, `dox_deg.py`（問い 1）
- 出力（`out/prelim/chemo2/`）:
  - 図: `heatmap_niche.png|pdf`（niche 42 set × 14 contrast）, `heatmap_hallmark.png|pdf`（Hallmark 50）, `heatmap_panel_genes.png|pdf`（手製パネル遺伝子の log2FC）
  - 表: `nes_long.csv`, `nes_wide.csv`, `gsea/<label>.csv`, `concordance.csv`, `concordance_counts.csv`, `concordance_nes_spearman.csv`,
    `classify_aging_dex_chemo.csv`, `panel_genes.csv`, `genesets_meta.csv`, `de/B_Dox_compadj.csv`,
    `dox_deg/deg_summary.csv`, `dox_deg/deg_membership.csv`, `dox_deg/gsea_gobp_B_Dox.csv`
  - gene set: `msigdb/{mh.all,m2.all,m5.go.bp}.v2025.1.Mm.symbols.gmt`

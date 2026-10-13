# LepR+ / CAR 間質で加齢と化療を並べる（案 B、2026-10-13）

前提: `docs/datasets2.md`（案 B）、`docs/prelim_chemo.md`（GSE296758）、`docs/prelim_chemo2.md`（gene set・GSEA の流儀）。予備段階の叩き台。
再生成: `uv run python -I prelim/niche/design.py` → `prelim/niche/de.py` → `prelim/niche/compare.py`（repo root から）。

## 結論（先に）

1. **加齢と化療が同じ向きに動かすのは「ECM・骨形成の低下」で、両側の 2 データすべてで再現した。**
   NABA core matrisome、collagen fibril organization、Hallmark EMT は 加齢 LepR+ / 加齢 CAR / Dox×4 / Dox×1 の 4 contrast すべてで FDR < 0.05 で ↓（NES −1.6〜−2.9）。
   ossification・osteoblast differentiation・bone mineralization も 4 本すべてで NES < −1.3 の ↓（FDR < 0.05 は 2〜3 本）。血管（Hallmark angiogenesis）は 加齢 LepR+ と Dox の両方で ↓。
   IFN-γ 応答 ↑ も両側で同じ向き（加齢 CAR +1.61、Dox ×4 +1.75 / ×1 +1.60 が FDR < 0.05、加齢 LepR+ は +1.24 で有意でない）。
   遺伝子では **Dmp1 ↓ が 3 contrast で共通**（加齢 LepR+ −6.5、加齢 CAR −0.83、Dox×4 −2.00。いずれも padj < 0.05）、Spp1 ↓ は 加齢 LepR+（−2.2, padj 0.051）と Dox（−1.58 / −0.99）。
2. **加齢で上がって、化療では上がらないもの = 炎症・SASP。ただしこれは CAR（GSE133922, 雄）だけで出る。**
   加齢 CAR: inflammatory response +2.27、TNF-α/NF-κB +1.84、IL6/JAK/STAT3 +1.92、SenMayo +2.07、myeloid activation +2.09（すべて FDR < 0.05）。遺伝子では Il1b +2.46、Serpine1 +1.49、Igfbp3 +1.64。
   Dox はどれも有意に上げない（SenMayo −1.38 / −1.37、IL6/JAK/STAT3・complement・coagulation はむしろ ↓）。
   **加齢 LepR+（GSE130885, 雌）では炎症は上がらない**（inflammatory response +1.27 で有意でない、手製炎症パネル −1.65 で ↓、Ccl2 −1.87）。加齢の炎症は 2 データで再現していない。
   加齢 CAR の炎症の形は LPS / poly I:C 刺激と似る（NES の Spearman 0.47 / 0.50）。Dox は LPS / poly I:C と似ない（0.05 / 0.24）。
3. **化療でのみ動くもの**: 増殖（E2F +3.3、G2M +2.7、Myc ↑。加齢では逆に ↓）、ER 応答（Hallmark estrogen early −1.6、Pgr −2.73、Greb1 −1.43）、
   TGF-β ↓、コレステロール恒常性 ↓、sprouting angiogenesis ↓、Reactome cellular senescence ↑（中身はヒストン・増殖系）。
   ER 標的は加齢側で比べられない（加齢 CAR は雄。加齢 LepR+ では Pgr が発現フィルタで落ち、Greb1 も動かない）。
4. **HSC 保持の niche 因子（Cxcl12, Kitl, Angpt1, Vcam1, Lepr）は加齢でも化療でも動かない**（|log2FC| < 1.3、padj ≥ 0.1）。
   動かしたのは LPS / poly I:C だけ（Kitl −1.87 / −2.22、Cxcl12 −0.72 / −0.53、Angpt1 −1.29 / −0.88）。Tnfsf11（RANKL）・Tnfrsf11b（OPG）・Adipoq も加齢・化療では動かない。
5. **一文**: 「加齢した LepR+/CAR 間質は、すでに ECM・骨形成プログラムが低く（CAR では炎症・SASP も高い）、化療は炎症を足さずに同じ ECM・骨形成・血管の軸をさらに下げ、ER 応答を消し、増殖を上げる。HSC 保持因子はどちらでも保たれる。したがって加齢 + 化療では『骨形成 niche の低下が上乗せされ、炎症は加齢由来のまま残る』と予想される」。
   ただし**交互作用は検定できない**（加齢マウスに化療を当てたトランスクリプトームは無い）。上の文は別実験の向きを並べた推論。
6. **Fig4B の補助パネルとしての可否**: Supplementary なら使える。主図は勧めない。
   使える点: 同じ細胞種（LepR+ / CAR）で、加齢 2 データと化療 2 contrast が ECM・骨形成 ↓ でそろう（`fig_niche.png` A の下半分）。「化療単独では炎症性 niche は再現されない」（`docs/prelim_chemo2.md` 結論 4）を、骨 bulk ではなくソートした間質で補強できる。
   弱い点: 加齢の炎症は 2 データで食い違う、性別がデータごとに違う、遺伝子単位の相関はほぼ 0（下記）、共通軸の骨形成 ↓ は「ソートに混ざる骨芽細胞の量」でも同じ形になる（GSE130885 の若齢群で Bglap が 428〜26,894 CPM とばらつく）。

## 1. データの設計（`out/prelim/niche/design_summary.csv`, `design_samples.csv`）

| | GSE130885（加齢） | GSE133922（加齢・炎症刺激） | GSE296758（化療） |
|---|---|---|---|
| 細胞 | Ter119−CD45− を FACS、100 細胞 bulk。6 集団（ALCAM+Sca1− 骨表面、**LepR+PDGFRα+ 中心骨髄**、LepR+PDGFRα− 中心骨髄、Sca1+PDGFRα+、CD31+ ×2） | FACS: CAR、PαS、類洞内皮、動脈内皮 | Lepr-Cre;tdTomato+ 骨髄 MSPC |
| 使った集団 | LEPR+_PDGFRA+_CBM（主）/ LEPR+_PDGFRA−_CBM（感度） | CAR | 全サンプル |
| 群・n | 若齢 8–12 週 3 / 老齢 24–25 か月 3 | 2 か月 4 / 2 年 4 / 2 か月 + LPS 4 / 2 か月 + poly I:C 4（2 週 4 もあるが不使用） | PBS×4 4 / Dox 10 mg/kg×1 3 / Dox 5 mg/kg×4 4。8–10 週で処置、最終投与 8 週後 |
| 性別（SOFT） | 記載なし | 記載なし | female（title） |
| 性別（発現から推定） | **雌**（Xist 49–380 CPM、Y 遺伝子 ≈ 0。1 例のみ Xist 8.8 で判定不能、使っていない集団） | **雄**（Ddx3y / Eif2s3y / Kdm5d / Uty 5–98 CPM。Xist の行が表に無い） | 雌（Xist+、Y 遺伝子の行が表に無い） |
| 系統 | C57BL/6 | C57BL/6JRj | Lepr-Cre;tdTomato |
| プラットフォーム | Smart-seq2 + Nextera XT, HiSeq 4000, mm9, featureCounts（RefSeq）raw counts | SMARTer Stranded Total RNA-Seq Pico, HiSeq 4000, mm10, htseq-count → DESeq2 normalized counts | SMART-Seq v4, NovaSeq 6000, GRCm39, featureCounts raw counts |
| 文献 | 未収載（SubSeries of GSE130899, Valletta / Nerlov） | PMID 31801092 | 未収載 |

- **GSE130885 の LepR+ の選び方**: 集団平均 CPM で LEPR+_PDGFRA+_CBM が Lepr 3,265、Cxcl12 37,232、Kitl 2,006、Adipoq 7,429 と CAR マーカーが最も高く（PDGFRA− は 1,131 / 4,780 / 392 / 1,456）、
  GSE296758 PBS MSPC との log CPM の Spearman（参照側上位 3,000 遺伝子）も 0.57 で 6 集団中最大（PDGFRA− は 0.30）。GSE133922 の CAR とも最大（0.39 vs 0.21）。→ 主解析は PDGFRA+（`gse130885_population_markers.csv`, `gse130885_population_similarity.csv`）。
- **LEPR+_PDGFRA− は感度分析として不適**: 老齢 3 例のうち 2 例で Lepr 4–8 CPM、Cxcl12 129–160 CPM しかなく（若齢は Lepr 490–5,239）、Cd79a（B 細胞）が 12–42 CPM ある。
  老齢側のソート純度の問題で、加齢の比較になっていない。この集団の DEG 1,414 個と、Dox との一致（NES Spearman 0.74、DEG 20 個中 19 個が同方向で Cd79a/b・Irf4・Alas2 ↑ を含む）は、**骨髄間質の減少と血球混入の増加という組成差どうしの一致**とみる（Dox×4 MSPC でも Cd79a +1.3、Alas2 +1.5、Hbb-bs +1.0 と血球混入が少し増えている）。表には残したが結論には使っていない。
- **GSE133922 の counts**: 公開されているのは DESeq2 の normalized counts だけだが、各列を「最小の非ゼロ値 = raw 1」として size factor を掛けると、79 列すべてで全値が整数に戻った（整数からのずれ最大 2e-9、値の最大公約数 1）。
  これを raw counts とみなし、DESeq2（pydeseq2）で size factor を推定し直して解析した。limma 等の代替は不要だった。
  列名 ADU / AGE / YOU と年齢の対応は SOFT に無い。PαS 2 か月が SOFT で 3 匹（replicate 3 欠番）、列 ADU_PAS が 3 本なので ADU = 2 か月、AGE = 2 年、YOU = 2 週とした（群サイズは 20 群すべてで SOFT と一致）。
- QC（`qc_pca.png`）: 両データとも PC1/PC2 で年齢がきれいに分かれる。CAR の LPS・poly I:C も別群に分かれ、poly I:C の 1 例はやや離れる。

## 2. 方法

- DE: `prelim/de.py` と同じ（filterByExpr 相当: CPM ≥ 10/中央値ライブラリ[M] をn_min サンプル以上、counts 合計 ≥ 15。pydeseq2 `~cond`、MLE、shrink なし）。DEG は |log2FC| ≥ 1 & padj < 0.05。比較に使う 2 群だけで fit。
- 化療側は既存 `out/prelim/chemo/de/B3_Dox4x_MSPC.csv`, `B3_Dox1x_MSPC.csv` を読むだけ（再計算していない）。
- GSEA: `prelim/chemo2/genesets.py` の collection（Hallmark 50 + niche 17 + 手製 5 = 72 set）、Wald 統計量の prerank、gseapy 1000 perm, seed 0, min 4, max 1000、FDR は 72 set で計算。設定は chemo2 と同じ。
- モジュール分類（`module_class.csv`）: 加齢の向き = 主 2 本（加齢 LepR+、加齢 CAR）のどちらかで FDR < 0.05、かつもう片方が逆向きに有意でない。化療の向き = Dox×4 か ×1 で同様。
- 遺伝子パネル（事前定義, `compare.py` の `GENE_PANEL`）: 依頼の niche 因子 11 個 + 炎症・脂肪・骨芽・ER 標的・老化・IFN。

| contrast | 検定遺伝子 | DEG ↑ / ↓ |
|---|---|---|
| 加齢 LepR+PDGFRα+（3/3） | 11,211 | 192 / 154 |
| 加齢 LepR+PDGFRα−（3/3, 不適） | 10,364 | 455 / 959 |
| 加齢 CAR（4/4） | 10,908 | 204 / 126 |
| LPS CAR（4/4） | 10,512 | 834 / 757 |
| poly I:C CAR（4/4） | 10,243 | 1,041 / 768 |
| Dox×4 MSPC（4/4, 既存） | 16,643 | 69 / 71 |
| Dox×1 MSPC（3/4, 既存） | 17,052 | 12 / 12 |

## 3. モジュール（`fig_niche.png` A, `heatmap_niche.png`, `heatmap_hallmark.png`, `nes_wide.csv`）

| 分類 | set | 加齢 LepR+ / CAR | Dox ×4 / ×1 |
|---|---|---|---|
| **加齢 ↓ 化療 ↓（同方向）** | NABA core matrisome | −1.58** / −2.19** | −2.50** / −2.14** |
| | GO collagen fibril org. | −1.94** / −2.34** | −2.80** / −2.00** |
| | H EMT | −1.66** / −2.36** | −2.90** / −2.23** |
| | GO bone mineralization | −1.99** / −1.72** | −1.95** / −1.35 |
| | GO ossification / osteoblast diff. | −1.92** / −1.53, −1.87** / −1.32 | −2.13** / −1.58**, −1.93** / −1.54** |
| | H angiogenesis | −1.76** / −1.22 | −2.16** / −2.02** |
| | H protein secretion | −1.66** / −1.52 | −2.47** / −2.23** |
| | GO vasculature dev.（加齢は割れる） | −1.52** / +1.26 | −1.88** / −1.97** |
| **加齢 ↑ 化療 ↑** | H IFN-γ response | +1.24 / +1.61** | +1.75** / +1.60** |
| | H allograft rejection | +1.28 / +1.93** | +1.99** / +2.01** |
| **加齢 ↑ のみ** | H inflammatory response | +1.27 / +2.27** | +1.38 / +1.28 |
| | H TNF-α/NF-κB | −0.98 / +1.84** | +1.05 / +1.39 |
| | SenMayo | +1.07 / +2.07** | −1.38 / −1.37 |
| | GO myeloid leukocyte activation | −1.03 / +2.09** | +1.16 / +1.19 |
| | H IFN-α response | +1.36 / +1.48** | +1.42 / +1.38 |
| **加齢 ↑ 化療 ↓（逆）** | H IL6/JAK/STAT3 | +1.48 / +1.92** | −1.46** / −1.21 |
| | H complement, H coagulation | +1.15 / +1.54**, +1.07 / +1.56** | −1.46** / −1.40, −2.23** / −1.86** |
| **加齢 ↓ 化療 ↑（逆）** | H G2M checkpoint, Myc v1 | −1.87** / −0.80, −1.63** / −1.62** | +2.71** / +2.64**, +2.38** / +2.18** |
| **化療のみ** | H E2F targets | −1.40 / −0.67 | +3.27** / +3.22** |
| | H estrogen early | −0.80 / +1.20 | −1.60** / −1.56** |
| | H TGF-β, sprouting angiogenesis, cholesterol homeostasis | NS | ↓ |
| | Panel adipocyte | −1.35 / +0.81 | +1.59** / +1.26 |
| **どちらも動かない** | Panel HSC retention, GO HSC homeostasis, Reactome SASP, p53, apoptosis, 破骨系 GO / WP, H adipogenesis, H Wnt | NS | NS |

** FDR < 0.05。分類の内訳: 同方向 11（↓ 9, ↑ 2）、逆方向 6、加齢のみ 8、化療のみ 21、どちらもなし 25（71 set。手製 senescence パネルはどの contrast でも 4 遺伝子未満で GSEA から外れた）。

- 加齢 2 データどうしの一致も弱い: NES の Spearman 0.46、両方で有意 7 set（同方向 5、逆 2）。共通は ECM・骨形成 ↓ と Myc ↓。炎症・SASP は CAR だけ。
- 加齢 vs 化療の NES Spearman: 加齢 LepR+ 0.24 / 0.19、加齢 CAR 0.16 / 0.17（Dox ×4 / ×1）。`docs/prelim_chemo2.md` の皮質骨の加齢 vs MSPC Dox（0.69）より低い。
  皮質骨 bulk で見えていた「加齢と化療が同じ向き」の大部分は、同じ細胞種にそろえると ECM・骨形成・血管の軸に絞られる。
- LPS / poly I:C は炎症・IFN・p53・apoptosis ↑ に加えて EMT ↑ で、ECM・骨形成の軸は加齢・化療と共通しない（ossification は LPS で ↓）。

## 4. 遺伝子（`fig_niche.png` B, `panel_genes.csv`, `deg_overlap.csv`, `deg_overlap_genes.csv`）

- niche 因子（log2FC, padj）: Cxcl12 / Kitl / Angpt1 / Vcam1 / Lepr はどれも加齢・化療で padj ≥ 0.1。**Spp1 だけが 加齢 LepR+（−2.18, 0.051）と Dox（−1.58, 2e-9 / −0.99, 6e-5）で ↓**。Tnfsf11・Tnfrsf11b・Adipoq は動かない。
  Il6 は 加齢 CAR で発現フィルタ落ち、加齢 LepR+ +1.39（0.85）、Dox +0.10。Ccl2 は 加齢で ↓ 傾向（−1.87 / −1.68, NS）、Dox ±0。
- 骨芽系: Dmp1 ↓（加齢 LepR+ −6.47, 加齢 CAR −0.83, Dox×4 −2.00）、Col1a1 ↓（−5.81 / −0.63 / Dox −0.69 NS）、Bglap ↓（加齢 CAR −2.38; 加齢 LepR+ −11.6 は Cook 外れ値で padj なし; Dox ±0）。
- 老化・炎症（加齢 CAR のみ）: Il1b +2.46、Serpine1 +1.49、Igfbp3 +1.64。Cdkn1a は どこでも動かず、Cdkn2a はすべてで発現フィルタ落ち。
- ER 標的（Dox のみ）: Pgr −2.73、Greb1 −1.43（×4）。Esr1 は動かない。
- IFN: Irf7 は 加齢 LepR+ で +2.74（0.017）、Dox では個々の ISG は動かない（モジュールの ↑ は弱い集積）。
- **DEG の重なり**（両方で検定された遺伝子を母集団）: 加齢 LepR+ ∩ Dox×4 = 6 個（期待 1.1、p = 0.0008）で、同方向は Dmp1 ↓・Postn ↓ の 2 個、逆が 4 個（Adra1b, Colgalt2, Megf10, Mmd2: 加齢 ↑ / Dox ↓）。
  加齢 CAR ∩ Dox×4 = 6 個（期待 1.4、p = 0.0024）で、同方向 5 個（Fabp4 ↑, Ptk2b ↑, A2m ↓, Adra2a ↓, Col24a1 ↓）、逆 1 個（Timp1: 加齢 ↑ / Dox ↓）。Dox×1 との重なりは 0。
- 遺伝子全体の Wald 統計量の Spearman は 加齢 vs Dox で −0.09〜+0.11 とほぼ 0（`spearman_gene_stat.csv`）。共通性はモジュール（ECM・骨形成）に集まり、遺伝子単位の相関としては出ない。

## 根拠が弱い点

- **別実験の比較**: 加齢と化療は別の研究室・別の時期・別のソート定義（PDGFRα+LepR+ 抗体ソート / CAR / Lepr-Cre;tdTomato+）・別のライブラリ法（100 細胞 Smart-seq2 / total RNA Pico / SMART-Seq v4）・別のゲノム版。交互作用は検定していない。
- **性別がそろわない**: SOFT には性別が無く、発現から推定すると GSE130885 は雌、GSE133922 は雄、GSE296758 は雌。**炎症・SASP の加齢 ↑ は雄の CAR だけで出て、雌の LepR+ では出ない**。性差なのかデータ差なのか区別できない。ER 応答の加齢側の比較もできない。
- **n**: 加齢 LepR+ 3/3、CAR 4/4、Dox ×1 は 3/4。100 細胞 Smart-seq2 は dropout が多く、log2FC が ±6〜12 になる遺伝子がある（Bglap −11.6, Adra1b +9.3）。
- **骨形成 ↓ の共通性は組成でも説明できる**: GSE130885 の若齢 LepR+PDGFRα+ では Bglap が 428〜26,894 CPM とばらつき、ソートに骨芽細胞が混ざる量が若齢側で多い可能性がある。Dox 側の骨形成 ↓ も「MSPC から骨芽系に向かう細胞の割合が減る」で同じ形になる。細胞内の転写変化か、集団の割合の変化かは bulk では分けられない。
- GSE296758 の Dox×4 vs PBS は flowcell と完全交絡（`docs/prelim_chemo.md`）。Dox 群は血球混入（Cd79a, Alas2, Hbb-bs）がわずかに多い。
- GSE133922 の counts 復元は「最小の非ゼロ値 = raw 1」という仮定に立つ。全値が整数に戻り最大公約数 1 なので妥当とみているが、著者の size factor が全 79 サンプルで推定されたことを前提にしている（DESeq2 で推定し直すので結果への影響は小さい）。列名と GSM の 1 対 1 対応は SOFT に無く、ADU / YOU の年齢対応は群サイズからの推定。
- 加齢側の若齢対照（8–12 週 / 2 か月）と Dox 側の採取時（処置 8–10 週 + 8 週 ≈ 16–18 週）で週齢が違う。
- 加齢 2 データどうしの一致が弱い（NES Spearman 0.46）ので、「加齢で動く軸」自体が ECM・骨形成 ↓ 以外は固まっていない。

## 生成物

- スクリプト（`prelim/niche/`）: `design.py`（SOFT・性別推定・集団選択・counts 復元）、`de.py`（加齢側 DE・PCA）、`compare.py`（GSEA・分類・遺伝子パネル・DEG 重なり・図）
- 出力（`out/prelim/niche/`）:
  - 図: `fig_niche.png|pdf|tiff`（170 × 113 mm。A: niche モジュール NES、B: niche 因子パネル log2FC）、`heatmap_niche.png|pdf`（niche 42 set × 7 contrast、感度列込み）、`heatmap_hallmark.png|pdf`、`qc_pca.png|pdf`
  - 表: `design_samples.csv`, `design_summary.csv`, `sex_markers.csv`, `gse130885_population_markers.csv`, `gse130885_population_similarity.csv`,
    `de/<label>.csv`, `de_summary.csv`, `gsea/<label>.csv`, `nes_long.csv`, `nes_wide.csv`, `module_class.csv`, `concordance_sets.csv`, `concordance_counts.csv`,
    `spearman_nes.csv`, `spearman_gene_stat.csv`, `panel_genes.csv`, `deg_overlap.csv`, `deg_overlap_genes.csv`, `qc_pca.csv`
  - counts: `counts/GSE130885_counts.tsv.gz`, `counts/GSE133922_counts.tsv.gz`（復元した整数 counts）

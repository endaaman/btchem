# in vitro 陽性対照: ヒト骨芽細胞・骨髄 MSC に DNA 傷害を与えると何が動くか（2026-10-13）

生体の骨（マウス Dox、ヒト生検）の RNA では「化療 → アポトーシス・老化・炎症」がはっきり出なかった（`docs/prelim_chemo2.md`, `docs/prelim_human.md`）。
同じ遺伝子・同じ gene set を、ヒトの骨系細胞に直接 DNA 傷害を与えた公開データで見て、**その軸が細胞レベルでは動くこと**（陽性対照）と、生体との差を確認した。
基準はこれまでと同じ（pydeseq2 0.5.4・MLE・`prelim/de.py` の `run()` と低発現フィルタ・DEG = |log2FC| ≥ 1 かつ padj < 0.05、prerank GSEA は Wald 統計量・gseapy 1000 perm・seed 0・min 4）。
再生成は `bash prelim/invitro/run_all.sh`。閾値・gene set は解析前に固定し、結果を見て動かしていない。

## 結論（先に）

- **ヒト骨芽細胞（NHO）と骨髄 MSC では、エトポシド（TOP2 阻害、Dox と同じ標的）でも放射線でも、p53・アポトーシス・老化・炎症がそろって強く上がる**。4 contrast すべてで
  Hallmark p53（NES 1.9–2.7）、Apoptosis（1.7–2.2）、TNF-α/NF-κB（1.8–2.8）、Inflammatory Response（1.6–2.4）、**SenMayo（1.8–2.7）** が FDR < 0.01。
  増殖は全滅（E2F / G2M NES −3.1〜−3.4、MKI67 −5〜−8、LMNB1 −5〜−6 log2FC）で、**8 日後の細胞は老化（増殖停止）状態**。CDKN1A（p21）は +2.1〜+3.4、GDF15 +2.6〜+6.2、MDM2 +0.7〜+1.8（すべて padj < 0.05）。
- **「化療 → アポトーシス因子」は、細胞では出る。ただし主役は p53 標的の一部で、全部ではない**。
  エトポシドでは BBC3（PUMA）+1.4/+0.3、PMAIP1（NOXA）+1.9/+1.2、FAS +1.1/+0.9、TNFRSF10B（DR5）+1.3/+1.1、CASP3 +0.6/+1.6（NHO / MSC、すべて padj < 0.05）。
  **BAX は動かないか下がる**（+0.1 / −0.8）。FASLG・TNF は発現なし（CPM 0）。放射線では PMAIP1 がむしろ下がる（−1.4 / −1.3）。
- **SASP（IL6, IL1B, CXCL8, CCL2）はエトポシドで大きく上がる**（IL6 +2.7/+1.3、IL1B +5.9/+5.4、CXCL8 +3.8/+6.0、CCL2 +3.1/+2.0）。放射線でも IL6・CCL2 は上がる。
  ただし絶対量は小さい（IL1B は対照 0.1–0.2 → 4–10 CPM）。**Reactome SASP は逆に負（NES −1.0〜−1.4、有意でない）**：この set の大半はヒストン遺伝子で、増殖停止で下がる分を拾っている。Reactome Cellular Senescence が有意に負（−1.7〜−1.9）なのも同じ理由（リーディングエッジは E2F1, CCNA2, CDK2, LMNB1 など増殖側）。老化の判定には SenMayo / 手製パネルを使う。
- **RANKL/OPG**: OPG（TNFRSF11B）は 4 contrast すべて上昇（+0.5〜+1.9）。RANKL（TNFSF11）は ほぼ発現なし（NHO エトポシドの +3.3 は 0 → 0.3 CPM の差で意味がない）。CSF1 は NHO で↑、MSC で↓。
  **in vitro の骨芽系細胞は DNA 傷害で破骨支持（RANKL/OPG 比）を上げる方向には動かない**（OPG 優位）。
- **骨形成は細胞で向きが割れる**: NHO は GO osteoblast diff. / ossification / mineralization が上がり（NES 1.5–1.8、FDR < 0.03。放射線で ALPL +4.6, COL1A1 +2.8）、
  MSC は RUNX2 −1.7/−2.0、COL1A1 −2.7/−1.9 で下がる（GO set は NES −1.2〜−1.5、FDR 0.06–0.2）。BGLAP・SP7 は対照で ~0 CPM で、出ている大きな log2FC はノイズ水準。
- **3D 骨（TE-Bone）× ドキソルビシン 4.16 µM 2 週（TempO-Seq パネル）**: p53 標的は上がる（GDF15 +3.3、TNFRSF10B +1.9、BAX +1.1、MDM2 +1.0、すべて padj < 0.05。Hallmark p53 NES 1.77, FDR 0.12）。
  炎症因子も上がる（IL1B +4.9、CXCL8 +3.5、MMP3 +6.7、IL6 +0.7）。一方で **基質・骨形成が大きく落ちる**（COL1A1 −4.0、SERPINE1 −3.7、NABA core matrisome NES −2.5、EMT −2.3、collagen fibril −2.2、FDR < 0.001）。
  ライブラリが Dox 群で 1/2〜1/5（15–57 万 vs 74–94 万リード）で、**細胞が減った（細胞毒性）状態を見ている**可能性が高い。CDKN1A はパネル上の発現が低くフィルタで落ちた。
- **生体マウス Dox（`out/prelim/chemo2`）との並置**（niche 関連 42 set、ヒト = NHO/MSC 4 contrast のうち 3 以上で同符号 FDR < 0.05、マウス = bulk 4 本の FDR < 0.05 の多数決）:
  - 同方向は IFN-α / IFN-γ response、adipogenesis、手製炎症パネルの 4 つだけ。
  - **ヒト in vitro だけで出るのが p53、Apoptosis、SenMayo、手製老化パネル、TNF-α/NF-κB、myeloid activation**（マウスはどれも FDR ≥ 0.06。p53 は NES +1.05〜+1.40 で向きは同じだが弱い）。
  - **逆方向は Inflammatory Response、complement、hypoxia、angiogenesis、EMT、Reactome cellular senescence の 6 つ**（マウス bone の Dox は炎症が下がる、LepR+ MSPC は EMT/血管が下がる）。
  - 増殖（E2F/G2M）はヒト細胞で↓、マウス LepR+ MSPC では↑（Dox 8 週後の再増殖）と時間相が違う。
  - **読み**: DNA 傷害そのものに対する細胞の応答（p53 → 老化 → SASP）は、ヒト骨系細胞で教科書どおり強く出る。生体の骨 RNA でそれが見えないのは、(1) 投与から採取まで 30 日〜8 週と長く、急性の p53 応答・死細胞が消えている、(2) 応答した細胞（老化した間質・骨芽細胞）が組織の一部で、骨髄・血液・筋の RNA に薄まる、(3) 生体では死細胞が除去され、残った細胞が再増殖する（マウス MSPC の E2F↑）、で説明できる。in vitro は 8 日・単一細胞種・高濃度で、応答が最も出やすい条件である。
- **Fig4 での使い道**: 補助パネル（または Supplementary）として「同じ gene set がヒト骨芽細胞・MSC では強く動く」＝**生体で出なかったのは方法の感度不足ではなく生物学的な違い（時間・希釈）**、を示すのに使える。
  `heatmap_vs_mouse` を絞った版（p53 / apoptosis / SenMayo / TNF-α / inflammatory / E2F / osteoblast diff. の 7 行 × ヒト 5 列 + マウス bulk 4 列）が候補。
  **ただし「乳がん化療（AC 等）で骨芽細胞がアポトーシスする」直接の根拠にはならない**（薬剤・用量・時間・細胞の状態が臨床と違う。下の限界）。

## データセットと設計（SOFT で確認）

| label | GSE | 細胞 | 処置 vs 対照 | n | 測定 |
|---|---|---|---|---|---|
| NHO_Eto | GSE284666 | ヒト初代骨芽細胞（NHO） | エトポシド 20 µM 8 日 vs DMSO 8 日（増殖中） | 4 / 4 | total RNA stranded、STAR → featureCounts（Ensembl v104）、1.1–1.9 億リード |
| NHO_IR | GSE284666 | 同上 | 放射線 10 Gy → 8 日 vs 同上 | 4 / 4 | 同上 |
| MSC_Eto | GSE284666 | 骨髄 MSC（BMMSC） | エトポシド 50 µM 8 日 vs DMSO | 4 / 4 | 同上 |
| MSC_IR | GSE284666 | 同上 | 放射線 15 Gy → 8 日 vs DMSO | 4 / 4 | 同上 |
| TEB_Dox | GSE135323 | TE-Bone（足場上のヒト骨芽細胞 3D 培養） | ドキソルビシン 4.16 µM、3 日投与 / 4 日休薬 × 2 週 vs 溶媒 | 3 / 3 | TempO-Seq S1500+（2,977 プローブ → 2,733 遺伝子に合算）、15–94 万リード |

- GSE284666: **ドナーの年齢・性別・入手先は SOFT に記載なし**。replicate 1–4 は同じ細胞（cell line 欄が "NHO" / "BMMSC"）の培養反復と読める＝**生物学的反復（ドナー）ではない**。対照は「DMSO 8 日の増殖中細胞」で、処置群は老化・増殖停止しているので、比較には「増殖 vs 停止」が含まれる。
  Ensembl → シンボルは `data/seq/GSE141595/human_ensembl2symbol.csv` で引いた（発現 ID 41,160 のうち 2,581 は引けず ID のまま。注目遺伝子はすべて引けている）。
- GSE135323: TE-Bone の Dox と対照は GSM4005429/31/33 と GSM4005423/25/27（SOFT の fastq 名で列を対応）。リンシチニブ群も samples.csv に入れたが使っていない。骨芽細胞の由来は SOFT に記載なし。低発現フィルタ後 1,337 遺伝子を検定。
- PCA（`qc/pca.png`）: NHO・MSC とも PC1（77–78%）で対照 vs 処置、PC2 でエトポシド vs 放射線がきれいに分かれる。TE-Bone は Dox 3 例が PC1 で対照から離れ、Dox 内のばらつきが大きい（ライブラリ 15 万の GSM4005429 が外れ寄り）。
- GSE17860（メルファラン 6 h、Agilent）は今回取得していない（6 時間で早すぎ、乳がんで使わない薬剤のため優先度を下げた）。

## DEG 数（padj < 0.05, |log2FC| ≥ 1、up / down）

| contrast | tested | up / down |
|---|---|---|
| NHO_Eto | 19,539 | 2,408 / 1,519 |
| NHO_IR | 18,275 | 2,132 / 1,717 |
| MSC_Eto | 19,676 | 2,623 / 1,843 |
| MSC_IR | 18,565 | 1,791 / 2,078 |
| TEB_Dox | 1,337（パネル） | 252 / 251 |

生体（マウス B_Dox など数十〜数百）とは桁が違う。他の閾値は `out/prelim/invitro/deg_counts.csv`。

## 注目遺伝子（log2FC、** padj < 0.05）

全体と群別 CPM は `out/prelim/invitro/focus_table.csv`、図は `focus_heatmap.png`。

| 遺伝子 | NHO Eto | NHO IR | MSC Eto | MSC IR | TE-Bone Dox | 注 |
|---|---|---|---|---|---|---|
| BAX | 0.1 | −0.1 | −0.8** | −0.4** | 1.1** | |
| BBC3 (PUMA) | 1.4** | 1.2** | 0.3** | 0.3 | −0.7 | |
| PMAIP1 (NOXA) | 1.9** | −1.4** | 1.2** | −1.3** | — | |
| FAS | 1.1** | 0.5** | 0.9** | 0.7** | −0.7 | |
| CASP3 | 0.6** | 0.2 | 1.6** | 1.0** | 0.3 | |
| TNFRSF10B (DR5) | 1.3** | 0.6** | 1.1** | 0.6** | 1.9** | |
| FASLG / TNF | — | — | — | — | — | 発現なし |
| CDKN1A | 3.4** | 2.7** | 2.1** | 2.1** | — | |
| CDKN2A | 1.1** | 0.2 | 0.2 | −0.2 | — | |
| SERPINE1 | 0.2 | 0.2 | 0.9** | 1.4** | −3.7** | |
| IGFBP3 | 1.5** | 1.7** | 0.3** | 0.0 | −1.3** | |
| IL6 | 2.7** | 3.6** | 1.3** | 2.5** | 0.7** | |
| IL1B | 5.9** | 2.4 | 5.4** | 0.6 | 4.9** | 対照 0.1–0.2 CPM |
| CXCL8 | 3.8** | −1.2 | 6.0** | 1.2** | 3.5** | |
| CCL2 | 3.1** | 3.9** | 2.0** | 1.8** | −0.9** | |
| MMP3 | 1.3** | −6.9** | 2.9** | 0.3 | 6.7** | |
| MMP13 | 7.8** | — | 0.2 | −0.3 | — | 0 → 0.3 CPM（ノイズ水準） |
| TNFSF11 (RANKL) | 3.3** | — | — | — | — | 0 → 0.3 CPM（ノイズ水準） |
| TNFRSF11B (OPG) | 1.0** | 1.7** | 0.5** | 1.9** | — | 110–440 CPM |
| CSF1 | 0.6** | 1.2** | −0.5** | −0.7** | — | |
| RUNX2 | −0.7** | 0.2 | −1.7** | −2.0** | — | |
| ALPL | 1.9** | 4.6** | 1.6** | 2.6** | — | |
| COL1A1 | 0.1 | 2.8** | −2.7** | −1.9** | −4.0** | |
| BGLAP / SP7 | (7.2** / 1.5) | — | (2.6** / −0.4) | (1.5 / −3.4**) | — | いずれも < 0.5 CPM |
| SPP1 | 2.4** | −1.3** | 3.9** | 0.6** | 4.8** | |
| ESR1 / GREB1 | 2.0** / 4.7** | 2.0** / 3.1** | 1.5** / 3.2** | 0.4 / 1.2 | — | いずれも ≤ 2.4 CPM。ER の機能を語れる量ではない |
| PGR | 0.5 | — | 1.2 | — | — | ほぼ発現なし |

— = 低発現フィルタで除外、または TempO-Seq パネルに無い。

## GSEA（`nes_wide.csv`, `heatmap_vs_mouse.png`, `heatmap_hallmark_iv.png`）

gene set は `prelim/chemo2/genesets.py` と同じ term（Hallmark 50 + MSigDB 17 + 手製 niche パネル 5）を MSigDB human 2025.1.Hs から引いた（マウスと set 名で 1:1 に並べるため。追加・削除なし）。
GMT は Broad から直接取得して `out/prelim/invitro/msigdb/` に置いた（h.all / c2.all / c5.go.bp）。

| set | NHO Eto | NHO IR | MSC Eto | MSC IR | TE-Bone Dox | マウス bulk Dox（B_Dox / Dox4x MSPC / Dox1x MSPC / Dox4x OB） |
|---|---|---|---|---|---|---|
| H p53 pathway | 2.69** | 2.03** | 2.39** | 1.91** | 1.77 (0.12) | 1.16 / 1.09 / 1.14 / 1.40（すべて FDR ≥ 0.12） |
| H apoptosis | 2.17** | 1.75** | 2.09** | 1.86** | 0.94 | −1.22 / −1.31 / −1.23 / 1.46 |
| SenMayo | 2.42** | 1.81** | 2.66** | 1.86** | 1.33 | −1.34 / −1.38 / −1.37 / 0.85 |
| H TNF-α/NF-κB | 2.78** | 1.96** | 2.35** | 1.84** | 1.34 | −2.20** / 1.05 / 1.39 / 1.75** |
| H inflammatory response | 2.30** | 1.76** | 2.35** | 1.64** | 1.25 | −1.81** / 1.38 / 1.28 / −0.99 |
| GO osteoblast diff. | 1.46** | 1.57** | −1.25 | −1.46 | 1.00 | |
| H E2F targets | −3.32** | −3.31** | −3.40** | −3.33** | 0.83 | |

** FDR < 0.05。TE-Bone は検定遺伝子 1,337 で set との重なりが小さく、FDR < 0.05 は 5 set（すべて基質・EMT・myogenesis・coagulation の低下）。
set 単位の向き分類は `concordance_vs_mouse.csv`（記述用。交互作用の検定ではない）。

## 限界

1. **薬剤**: エトポシドは乳がんの周術期レジメンではない（TOP2 阻害で Dox と機序が同じことが根拠）。放射線も化療ではない。Dox を使った TE-Bone はパネル測定で、全遺伝子の GSEA は参考値。
2. **用量・時間**: エトポシド 20–50 µM・8 日連続は老化誘導のための条件で、臨床の血中濃度・曝露時間より強い。TE-Bone の Dox 4.16 µM も細胞毒性域（ライブラリが 1/2〜1/5）。応答は「最大限に出る条件」の値。
3. **反復とドナー**: GSE284666 は 1 ドナー由来と思われる細胞の培養反復（n=4）で、ドナー間のばらつきを含まない。年齢・性別・閉経状態は不明（SOFT に記載なし、論文本文は未確認）。DEG 数が数千になるのはこのため（統計的な有意性を生体と比べない）。
4. **対照の状態**: 対照は増殖中の細胞で、処置群は増殖停止している。DEG・GSEA の相当部分は「増殖 vs 停止」（E2F/G2M、ヒストン、Reactome の senescence/SASP set の負）。
5. **ER 軸は見られない**: ESR1/PGR/GREB1 は培養細胞でほとんど発現せず（≤ 2.4 CPM）、エストロゲンも加えていない。2-hit 仮説の「閉経」側はこのデータでは検証できない。
6. **破骨支持**: RANKL がほぼ発現しない培養条件なので、RANKL/OPG の評価は OPG 側だけ。破骨細胞もいない。
7. **TE-Bone**: TempO-Seq S1500+ は毒性学用の 2,977 プローブで、骨形成マーカー（RUNX2, ALPL, BGLAP, SP7）と RANKL/OPG/CSF1 は載っていない。CDKN1A は載っているが発現が低く除外された。

## 出力

- スクリプト: `prelim/invitro/{prepare,ivcontrasts,de_iv,gsea_iv,figs_iv}.py`, `run_all.sh`
- counts: `data/seq/GSE284666/{counts.tsv.gz,samples.csv}`, `data/seq/GSE135323/{counts.tsv.gz,samples.csv}`
- `out/prelim/invitro/`: `de/<label>.csv`, `deg_counts.csv`, `focus_table.csv`, `volcano_all.{png,pdf,tiff}`（170 mm、A–E）, `volcano/<label>.{png,pdf}`,
  `focus_heatmap.{png,pdf,tiff}`（80 mm）, `heatmap_vs_mouse.{png,pdf,tiff}`（170 mm）, `heatmap_hallmark_iv.{png,pdf,tiff}`, `nes_long.csv`, `nes_wide.csv`, `concordance_vs_mouse.csv`, `gsea/<label>.csv`, `qc/pca.{png,pdf}`, `qc/group_cpm.csv`, `msigdb/`

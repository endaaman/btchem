# 追加公開データ調査 第三波（論文からの逆引き・GEO 以外のリポジトリ）

調査日: 2026-10-13。目的は「閉経後（加齢）の骨で、化療が骨転移 niche を作る」を dry で補強できる**データ**を探すこと。
第一波（`docs/datasets.md`）と第二波（`docs/datasets2.md`）は GEO / ArrayExpress / ProteomeXchange を総当たりし、
骨・骨髄 ×（加齢 / OVX / 化療）の RNA-seq はほぼ出尽くしていた。そこで今回は**論文から逆引き**し、Data availability と補足表に書かれた所在をたどった。
あわせて Zenodo / figshare / Dryad / Mendeley Data を直接検索した。
既に評価したデータ（第一波・第二波・`,docs/strategy_20261012.md` に出てくる accession 132 件）は再評価していない。

## 結論

- **理想の 2×2（加齢または OVX × 細胞傷害性化療、骨・骨髄）の新しいデータは、GEO 以外でも見つからなかった。**
  論文の Data availability（OA 全文 約 2,000 本）、Zenodo、figshare、Dryad、Mendeley Data のどれにも無い。
  加齢 × 化療を同じ実験で組んだものは、第二波の PXD078161（Izquierdo 2026）と、表現型だけの Yao 2020（OVX × Doxo）のままである。
  新しく見つかった「年齢 × 化療」の実験は busulfan の 1 本（PMID 39139448）だけで、雄、若齢 4–6 週 vs 成体 16–20 週、移植前処置で、データ寄託もない。
- **今回の収穫は「個体ごとの定量値」**。論文の補足 Excel に、トランスクリプトームではないが個体単位の数値が公開されているものが 3 本ある（すべて `data/lit/` に取得済み）。
  1. **Wasi 2024（PMID 39147353）`mmc2.xlsx`**: 成熟マウス（50 ± 4 週、C57BL/6J、雌雄）に Dox 2.5 / 5 mg/kg × 6 回。
     **骨細胞 TUNEL 陽性率が個体ごとに載っている**: 非投与 6.0 / 16.2 / 7.4 %（n=3）、低用量 30–58 %（n=5、平均 47 %）、高用量 77 / 84 %（n=2）。
     TRAP、骨標識、qPCR（Wnt1, Acp5, Ocn, Opg, Rankl, Sost, Rankl/Opg 比、2 週と 5 週）もある。若齢群がないので 2×2 にはならないが、
     **Fig4 の RNA で出なかった「化療による骨細胞アポトーシス」を、中年マウスの組織レベルで示す公開値**として使える。
  2. **Sci Rep 2025（PMID 40389522）`MOESM1_ESM.xlsx`**: 雌マウスに乳がんの実レジメン（AC、TAC）。µCT と骨の qPCR が個体ごとに載っている。
     19 週齢で TAC を始めた BALB/cJ を、開始 35 週後と 62 週後（約 20 か月齢）まで追った群がある（化療 → 加齢の順）。
     我々が個体値から計算すると、AC 2 週後の皮質骨で Cdkn1a ×2.1、Cdkn2a ×1.9、Mmp13 ×2.2、Ccl2 ×1.4（いずれも p < 0.05）。一方、Il6 ×0.73、Tnfsf11（RANKL）×0.98 で、
     **全骨の qPCR では SASP・RANKL は上がっていない**。TAC 7 週後は椎体の Tnfsf11 が ×0.57 に下がり、12 か月後は老化マーカーに差がない。
  3. **Nat Commun 2025（PMID 41469373）Source Data `MOESM4_ESM.xlsx`**: 雌 12 週、Doxo / PTX。SA-β-gal、p16、TRAP、**Rankl / Opg の qPCR**（骨常在細胞分画、精製骨髄脂肪細胞、CAR 細胞）が個体ごとに載っている。
     **PTX を投与した骨転移モデル（PyMT-BO1 を心腔内に接種）**の BLI、N.Oc/BS もある。若齢のみ。生の顕微鏡画像は Zenodo 17314594（4.8 GB、未取得）。
- **文献をつなぐと、OVX と化療が同じ細胞・同じ分子に収束する**: OVX による海綿骨の減少には、骨髄脂肪系前駆細胞（MALP、Adipoq+）の RANKL が必要である（Bone Res 2025、PMID 40102423）。
  化療による骨量減少は、老化した CAR 細胞・骨髄脂肪細胞の RANKL を介する（Nat Commun 2025）。データは無いが、
  **「閉経（1st hit）と化療（2nd hit）が骨髄の脂肪系 niche 細胞（LepR+ / CAR / MALP）の RANKL で重なる」**という仮説は、手元のデータで検証できる。
  手元の GSE296758（Dox × 4、LepR+）、GSE289491（Dox、CAR）、GSE133922（加齢 CAR）、GSE130885（加齢 LepR+）で Tnfsf11 / Tnfrsf11b を見ればよい（未実施）。
  ただし Sci Rep 2025 の全骨 qPCR では RANKL は動いておらず、細胞分画で見ないと出ない可能性がある。
- **加齢側の骨髄 niche scRNA が 2 本増えた**（第二波の案 B の加齢側を補う）:
  **GSE299773**（PNAS 2025、PMID 40063797、Lin−CD45−CD71− 骨髄間質、2 / 12 / 24 か月、各 1 サンプル。LepR+ の 5 亜集団と、加齢で増える炎症性の Sca1+Cxcl9+ LepR+ 細胞）と、
  GSE145477（骨髄間葉系、1 / 1.5 / 3 / 16 か月、MALP を定義したデータ）。
- **臨床（手順 4）では、個人単位のデータは公開されていない**。群平均として使えるのは次の 2 つ。
  - PMID 25923354（PLoS One 2015）の Table 2: 閉経前 29 例 / 閉経後 12 例で、術前化療 4 サイクルの前後の骨代謝マーカー。変化は閉経状態で大差ない（下表 §4）。
  - ClinicalTrials.gov NCT01298362: 閉経後の女性では、化療完了時の腰椎 BMD 変化は −0.19 %。
- **否定的・反対方向の知見も見つかった**（本文で触れるかの判断材料）:
  - Sci Rep 2025: 化療後の骨の老化マーカーは長期では消える。
  - PMID 29984313: 骨転移は若齢（6 週）マウスの方が 6 か月齢より多い。ただし 6 か月齢でも OVX を加えると骨破壊が増える。
  - PMID 38457059: 老化促進マウス（SAMP6）の骨髄 MSC の SASP は、乳がん細胞の生存を**抑える**。

### 使い方の提案

| 案 | 中身 | 弱点 |
|---|---|---|
| **E. 公開個体値の再プロット（Supplementary）** | Wasi 2024 の骨細胞 TUNEL（用量依存 10 → 47 → 81 %）と、Nat Commun 2025 の Rankl/Opg（骨常在細胞・CAR・骨髄脂肪細胞）を並べ、「RNA では見えない化療の骨細胞・niche への作用が、組織・細胞分画では出ている」と示す | 年齢 × 化療の 2×2 にはならない。Wasi は雌雄混合で、個体の性別が Excel に無い。図に再利用するときはライセンスの確認が要る（Sci Rep・Nat Commun は CC BY、Wasi は Elsevier の NIH 著者原稿） |
| **F. RANKL 収束仮説を手元データで見る** | GSE296758 / GSE289491（化療 LepR+ / CAR）と GSE133922 / GSE130885 / GSE299773（加齢 LepR+ / CAR）で Tnfsf11・Tnfrsf11b・Rankl/Opg 比を比べる | 別々の実験の比較。Sci Rep 2025 の全骨 qPCR では RANKL は動いていない |
| G. 本文・Discussion で引用 | 38161426（腫瘍のない宿主を AC で前処置すると、ER+ ZR-75-1 の肝・肺・**骨**転移が増える）、37699010（ドセタキセルが間質を傷害し、IL-6 / G-CSF で休眠がんが目覚める）、40102423 と 41469373（上の RANKL 収束） | データとしてはどれも骨・加齢の 2×2 ではない |

---

## 1. 指定 5 論文の公開定量データ（手順 3）

| PMID | 論文 | データの種類 | 所在（取得したもの） | 条件 | 本研究での評価 |
|---|---|---|---|---|---|
| **39147353** | Wasi 2024 *Bone* | 個体ごとの骨細胞 TUNEL（TUNEL+ / − 骨細胞数、%）、TRAP（骨内膜 % 面積）、ALP、骨標識（Ec.sL/BS）、µCT（pMOI、骨形成・吸収体積、3 部位）、qPCR（Wnt1, Acp5, Ocn, Tnfsf11b, Tnfsf11, Sost, Rankl/Opg 比。2 週と 5 週） | `data/lit/39147353/mmc2.xlsx`（ScienceDirect の `ars.els-cdn.com/content/image/1-s2.0-S8756328224002242-mmc2.xlsx`）。`mmc3.xlsx` は荷重の群別要約表、`mmc1.docx` は補足図表 | C57BL/6J 雌雄、50 ± 4 週。Dox 2.5 / 5 mg/kg × 6 回（2 週）± Yoda1 ± 脛骨荷重。高用量 Dox + 溶媒は雌 3 匹 | ○ **化療で骨細胞アポトーシスが増える個体値**（非投与 9.9 % → 低用量 47 % → 高用量 81 %。Kruskal–Wallis + Dunn で非投与 vs 高用量 p = 0.02）。若齢群なし。Excel に性別の列はない |
| **40389522** | Ali 2025 *Sci Rep* | 個体ごとの µCT（BV/TV、Tb.Th/N/Sp、皮質骨厚）と qPCR（Cdkn1a, Cdkn2a, Cd38, Il1a, Il1b, Il6, Hmgb1, Mmp13, Tnfsf11, Ccl2, Cxcl12, Spp1、Alpl, Bglap, Acp5, Ctsk）。椎体 L6 と皮質骨、腎・肝 | `data/lit/40389522/41598_2025_1717_MOESM1_ESM.xlsx`（Fig. 1–8 のシートごと）、`MOESM2_ESM.docx`（Supplementary Table S1 など）、本文 XML | 雌 C57BL/6 16 週: AC × 3 回 → 2 週後（n = 8–9）。TAC → 7 週後（n = 5–6）。雌 BALB/cJ 19 週: TAC → 開始 35 週後と 62 週後（約 20 か月齢、n = 5）。± D+Q / ピペルロングミン | ○ 乳がんの実レジメン × 雌で、**化療後に加齢まで追った唯一の個体値**。上の「結論」の値は、各行の前半を対照・後半を化療として我々が再計算したもの（列の分け方はシートの見出しと図の n から判断した） |
| **41469373** | Nat Commun 2025（Stewart 研） | Source Data（Figure 1–8、S1–S15 のシート）。SA-β-gal、p16 / tdTomato、TRAP、BFR、µCT、**Rankl / Opg / 比の qPCR**（骨常在細胞分画 n = 7–12、精製骨髄脂肪細胞 n = 3–4、CAR 細胞 n = 4–5）、骨転移モデル（PTX ± D+Q）の BLI・N.Oc/BS・Oc.S/BS | `data/lit/41469373/41467_2025_67793_MOESM4_ESM.xlsx`（57 MB）、`MOESM1_ESM.pdf`（補足図）。生画像は Zenodo 10.5281/zenodo.17314594（`RawData for Upload.zip` 4.8 GB、未取得）。scRNA は GSE289491（評価済み）と GSE230295（小児 ALL の骨髄で ✕） | 雌 12 週（INK-ATTAC、p16-CreERT2/tdTomato、ADQ-Cre/DTR、ADQ-Cre/Rankl fl/fl、WT）。Doxo、PTX | ○ **化療 → 脂肪系 niche 細胞の老化 → RANKL ↑ の個体値**と、**化療 × 骨転移モデル**の値。若齢のみ |
| **31932453** | Yao 2020 *Cancer Res* | 補足図のみ（数値表なし）。OVX × Doxo の µCT、PTX、放射線、OVX の SASP qPCR | `data/lit/31932453/Yao2020_supp_figs.pdf`（AACR の figshare 22424466）、本文 XML（PMC 著者原稿） | 第二波に記載 | ○（第二波と同じ）。数値は図から読むしかない |
| **42686970** | Izquierdo 2026 *GeroScience* | MOESM1–4 = 年齢・性別ごとの Spectronaut 比較表（MassIVE の `VG1–4_*_Candidates.xlsx` と同じ）。MOESM5 = Fig. S1（**加齢による変化と DX による変化の log2FC が正に相関**: 雄 Pearson r = 0.35、p = 9 × 10⁻¹²、雌 r = 0.25、p = 0.003） | `data/lit/42686970/11357_2026_2485_MOESM{1..5}_ESM.*`（Springer の static-content から直接取得。論文 HTML は bot 対策で取れない） | 第二波に記載 | ○ プロテオームは `docs/prelim_proteome.md` で解析済み。µCT・qPCR の値は補足に無い。Fig. S1 は「DX は骨で加齢に似た変化を起こす」の根拠として引用できる |

## 2. 論文からの逆引きで見つかったデータ（手順 1）

| PMID | データの種類 | 所在 | 条件（種・性別・年齢・処置） | 評価 |
|---|---|---|---|---|
| 40063797（PNAS 2025） | scRNA-seq、骨髄間質（Lin−CD45−CD71−） | **GSE299773**（= PRJNA1085191、RAW.tar）。SOFT を `data/seq/GSE299773/` に取得 | マウス、2 / 12 / 24 か月、各 1 サンプル。論文では雌雄を加齢させたとあるが、サンプルごとの性別は SOFT に記載なし | ○ 加齢側の LepR+ niche。加齢で周細動脈の LepR+ が炎症化（Sca1+Cxcl9+）。案 B（化療 LepR+ / CAR との比較）の加齢側に追加できる |
| 40102423（Bone Res 2025） | 文献（MALP の RANKL が OVX による骨量減少に必要）。新規 scRNA は無く、GSE145477 / GSE176171 / GSE253355 の再解析 | 「データは請求に応じて」 | マウス雌 3 か月、Adipoq-CreER × Rankl fl/fl、OVX | ○ 文献として。化療側（41469373）と同じ細胞・同じ分子に収束する根拠 |
| — （Zhong 2020、40102423 などが再利用） | scRNA-seq、骨内膜側の骨髄間葉系細胞 | **GSE145477**（RAW.tar）。SOFT を `data/seq/GSE145477/` に取得 | マウス B6SJL/F1、1 / 1.5 / 3 / 16 か月、各 1 サンプル。性別は記載なし | △ 加齢の MALP 参照。16 か月 1 点のみ |
| 37548937（2024） | bulk RNA-seq、大腿骨 | **GSE220679**（counts、`RNAseq_mouse.count.txt.gz`）。SOFT を `data/seq/GSE220679/` に取得 | マウス C57BL/6J **雄**、2 か月 vs 22 か月、n = 5 / 5 | △ 加齢の骨 bulk（雄）。手元の加齢皮質骨データと同じ役割 |
| 38904042（2024） | bulk RNA-seq、脛骨全体 | **GSE262441**（processed counts）。SOFT を `data/seq/GSE262441/` に取得 | マウス**雄**、若齢 / 老齢 × カロリー制限、n = 12 / 群 | △ 加齢の骨 bulk（雄）。n は多い |
| 37524694（Nat Commun 2023） | CyTOF（Mendeley 10.17632/hg7sd7hbk5.2）、scRNA（GSE237307）、CITE-seq（GSE237301） | GEO・Mendeley。SOFT（GSE237307）を `data/seq/GSE237307/` に取得 | マウス**雄**、老齢、INK-ATTAC、骨・骨髄の間葉系 | △ 加齢骨の老化細胞（p16+Ki67−BCL-2+）の同定。化療なし |
| 41592025（JCI 2026） | bulk RNA-seq（GSE283633: エトポシド処理 BMSC vs マクロファージ、若齢 / 老齢の単球。GSE294362）、CyTOF（Mendeley 10.17632/tk466mks85.1） | GEO・Mendeley | マウス老齢、p16+ 骨髄系細胞の除去 | △ 加齢骨髄の老化は骨髄系と間葉系で別物。化療は in vitro のエトポシドのみ |
| **38161426**（iScience 2024） | NanoString（PanCancer Immune / Pathways、肺）、図の値 | Zenodo 10.5281/zenodo.10045006 → `data/lit/38161426/nanostring_data.tar.gz`（RCC ファイル） | 雌 BALB/c 6–9 週、腫瘍なしで AC（Dox 2.7 + Cy 43 mg/kg）× 3–4 回 → 7–10 日後に腫瘍細胞を静注。NRG マウス + E2 ペレット + **ER+ ZR-75-1** | ○（文献）／△（データ）。**化療で前処置した宿主では、ER+ 細胞の肝・肺・骨転移が増える**（Fig. 1E、骨は 2 切片 / 匹、n = 3 / 群）。NanoString は肺のみ |
| 37699010（PLoS Biol 2023） | scRNA（GSE231350）、Source Data Excel、FlowRepository FR-FCM-Z6HQ / Z6J8 | GEO ほか | ドセタキセル、腫瘍間質オルガノイド、同系乳がんモデル（乳腺） | △ 骨ではない。「タキサンが間質を傷害 → IL-6 / G-CSF → 休眠がんが目覚める」の機序の引用元 |
| 36038666（Leukemia 2022） | 培養骨髄 MSC の RNA-seq（tMN n = 3 vs 年齢を合わせた健常）、サイトカイン | 寄託なし。補足 docx（`data/lit/36038666/`、66 MB） | **ヒト高齢**（中央値 69–73 歳）。殺細胞性治療を受けた後の治療関連骨髄腫瘍。骨髄腫の化療前後ペア（Fig. 8） | △ 高齢ヒトの骨髄間質が化療後に老化・SASP を示す。ただし培養細胞で、若年の比較がない。乳がんは含まない |
| 39139448（Front Cell Dev Biol 2024） | 図の値のみ（血管漏出、µCT、骨代謝） | 寄託なし | マウス**雄**、若齢 4–6 週 vs 成体 16–20 週 × busulfan 40 / 80 mg/kg（移植前処置） | △ **今回唯一の新しい「年齢 × 化療 × 骨髄」**。骨髄血管の変化は若齢では回復し、成体では 6 週後も残る。雄で、閉経後に相当する年齢ではない |
| 41502782（J Bone Oncol 2026） | 図の値（µCT、骨髄細胞集団のフロー、骨芽細胞・破骨細胞数） | 補足は抗体表などで、数値なし | 雌 BALB/c と BALB/c ヌード、7 週。Dox 4 / 6 mg/kg を週 1 回 × 4 | △ 雌、乳がんモデルと同じ系統での Dox の骨・骨髄への影響。若齢のみ |
| 42122233（Cancers 2026） | 図の値（µCT、血清・骨髄の PINP / OC / TRACP-5b） | 請求に応じて | 雌 BALB/c 4 週、Dox 6 mg/kg × 4、回復期あり | △ 幼若期の Dox で骨量減少が持続 |
| 40057216（Bone 2025） | 5-FU 後の骨髄間葉系 scRNA（MALP が増え、筋線維芽細胞化） | 「請求に応じて」（寄託なし） | マウス、5-FU 単回 | ✕ データ非公開。文献としては「化療で MALP が増える」の根拠 |
| 29984313（Tomography 2018） | 図の値 | 寄託なし | 雌、若齢 6 週 vs 6 か月 ± OVX、大腿動脈注入の骨転移モデル | △ 反対方向: 骨転移は若齢の方が多い。6 か月齢でも OVX で骨破壊が増える |
| 38457059（Stem Cell Rev Rep 2024） | 図の値 | — | SAMP6（骨の老化促進）の骨髄 MSC の培養上清 | △ 反対方向: 老化 MSC の SASP は乳がん細胞の生存を抑える |
| 41619722（Cell Rep Med 2026） | ヒト骨転移 62 例の scRNA（13 原発、正常骨髄・原発巣のペアあり） | NODE OEP005136。臨床情報は `mmc3.xlsx` | ヒト脊椎骨転移 | △ 骨転移 niche の参照。治療歴で層別できるかは臨床表次第（未確認） |

## 3. GEO 以外のリポジトリの直接検索（手順 2）

| リポジトリ | ID | 内容 | 評価 |
|---|---|---|---|
| Mendeley Data | 10.17632/xs76zjzvc5 | 4T1 担がん雌 BALB/c + DOX で、骨髄の脂肪蓄積と骨髄リピドミクス（フォルモノネチン介入）。`original data.7z` 4.2 GB（未取得） | △ 化療で骨髄脂肪が増える（若齢・担がん）。リピドミクスの生データ |
| Mendeley Data | 10.17632/ht5846dgdk ほか | 老齢マウスの CCL3 と骨髄脂肪化・骨量減少の元データ（PMID 36378535） | △ 加齢のみ |
| Zenodo | 17314594 | Nat Commun 2025 の生の顕微鏡画像と定量表（Source Data と同じ内容の xlsx） | ○ 上の §1 と同じ |
| Zenodo | 10045006 | iScience 2024 の NanoString（肺） | △ 取得済み |
| Zenodo | 7319842 | 「Osteocytes regulate senescence of bone and bone marrow」（eLife 2022）のデータ | △ 加齢のみ |
| Zenodo | 15587768 / 20429119 / 12684327 | 加齢マウス骨髄の CITE-seq、ヒト加齢 CD146+ MSC の生データ（E-MTAB-16622）、老化した血管 niche と MSC | △ 加齢のみ（造血側が中心） |
| figshare | 22424466 | Yao 2020 の補足図 PDF | ○ 取得済み（§1） |
| figshare | 10325978 | 化療（若齢マウス）による骨・筋の異常をビスホスホネートで改善（Front Endocrinol 2019） | △ 若齢、雄中心 |
| figshare | 1398219 / 1398224 | PMID 25923354 の図（閉経前 / 閉経後の骨代謝マーカー） | ○ 値は本文 Table 2（§4） |
| Dryad | 10.25338/B8693K（= Zenodo 8310144） | 人工骨髄内の骨肉腫細胞に Dox を投与した bulk RNA-seq | ✕ 腫瘍細胞側 |
| Dataverse（Borealis） | 10.5683/SP3/MTHK1O | ラット椎体転移モデル × ドセタキセルの時期（PMID 40245075） | △ 化療 × 骨転移・骨質。加齢・閉経なし |

**Zenodo / figshare / Dryad / Mendeley に、加齢または OVX × 化療の骨・骨髄データは無かった。** ヒットの大半は腫瘍細胞（骨肉腫・白血病）、ドラッグデリバリー材料、加齢単独、OVX × 生薬だった。

## 4. 乳がん臨床: 化療前後の骨代謝マーカー（手順 4）

| 出典 | 対象 | 項目 | 所在 | 評価 |
|---|---|---|---|---|
| **PMID 25923354**（PLoS One 2015） | 中国、術前化療（CEF/EC 29 例、EC-T 12 例）、**閉経前 29 / 閉経後 12**、健常対照 30 / 30 | 血清 ICTP、β-CrossLaps、BAP、PINP、E2、FSH、25(OH)D、PTH。診断時（C0）と各サイクル 24 時間後（C1–C4） | 本文 Table 2・3（群平均 ± SD）、Table 5（12 例の個人値: E2、各マーカーの変化率、6 か月後の BMD T スコア）。補足ファイルは無い | ○ C0 → C4 の変化（Table 2 から計算）: 閉経前 ICTP +18 %、β-CTX +36 %、BAP −10 %、PINP −27 %。閉経後 ICTP +29 %、β-CTX +32 %、BAP −18 %、PINP −24 %。**吸収↑・形成↓は閉経状態によらない**。閉経前は E2 が 149 → 7 pg/mL に落ちる（化療による卵巣機能低下） |
| ClinicalTrials.gov **NCT01298362**（PMID 27212697） | ギリシャ、ER+ の閉経後早期乳がん、化療（アントラサイクリン / タキサン）→ AI 群と AI 単独群、290 例 | DXA の BMD | 結果が登録済み（群平均のみ） | △ **閉経後では化療完了時の腰椎 BMD 変化は −0.19 %**、AI 12 か月で −0.59 %。閉経後の骨量への化療そのものの影響は小さい |
| PMID 26998426（J Bone Oncol 2016、OPTION 試験 NCT00427245） | 閉経前、化療 ± ゴセレリン | 骨型 ALP、尿 NTX（36 か月まで） | 本文の図表のみ | △ 化療で卵巣機能が落ちると骨代謝回転が上がり、ゴセレリンで防げる（化療 → 医原性閉経の経路） |
| PMID 36899284（Arch Osteoporos 2023） | 閉経後 109 例、化療 | BMD、FRAX | 請求に応じて | △ 第二波の 37410202 と同じコホートとみられる（未確認） |
| PMID 27146497（Support Care Cancer 2016） | 早期乳がん、アントラサイクリン + タキサン（高用量プレドニゾロン併用） | BMD（0 → 4 か月） | 本文のみ | △ |
| PMID 37691079（Eur Radiol 2024） | 骨髄毒性のある化療を受けた患者 19 例（乳がんに限らない）と、年齢・性別をそろえた対照 38 例 | 椎体骨髄の MR プロトン密度脂肪率（PDFF） | 本文のみ | △ **ヒトでも化療後に骨髄脂肪率が +8.7 % 上がる**（p = 0.01）。年齢で層別した解析はない |

**乳がん患者で、化療前後の骨代謝マーカーや骨髄の測定値を個人単位で公開しているデータ（Project Data Sphere、Vivli、YODA、試験の公開データ）は見つからなかった。**
ClinicalTrials.gov で「breast cancer × chemotherapy × 骨代謝マーカー / BMD / DTC」を検索した（下記）。結果が登録されている試験は NCT01298362 と NCT00295867（ゾレドロン酸 × 骨髄微小転移）だけだった。

---

## 方法

- **PubMed**（E-utilities `esearch db=pubmed`、2015–2026、`[tiab]`）: 下記 21 クエリで 1,285 件。書誌と OA 状況は Europe PMC（`searchPOST`、`EXT_ID`）で取り、PubMed → GEO は `elink` で引いた。
- **Europe PMC**（`/search`、`resultType=core`、`SRC:MED`、2015–2026）: 下記 20 クエリ。全文にもマッチするので 1 クエリあたり上位 2,000 件（関連度順）に絞った（計約 35,000 件）。
  そこから、タイトルに骨・骨髄の語があり、抄録で化療 × 加齢 / OVX、化療 × オミクス、または化療 × niche / 骨細胞 / 脂肪 / 老化 / 転移 / 休眠を満たすもの 501 件を選んだ。
- **逆引き（被引用）**: Europe PMC の `/MED/<PMID>/citations` で、Yao 2020、Luo 2016、Ottewell 2015、Wasi 2024、Sci Rep 2025、Nat Commun 2025、Chandra 2022、Izquierdo 2026、Lawson 2015、PMID 30119678（骨髄間質の老化が ER+ 休眠細胞を起こす）、Demaria 2017 を引用している論文 計 1,811 件。
- **全文**: 上記のうち OA の論文の `fullTextXML`（計 1,988 本）から、Data availability / Accession / Source data の節、補足ファイル名、本文中の accession（GSE、PRJNA、SRP、E-MTAB、PXD、MSV、IPX、ST、MTBLS、HRA、CRA、OMIX、phs、EGAS）、リポジトリの URL / DOI（Zenodo、figshare、Dryad、Mendeley、OSF、Synapse、NGDC、GitHub）を正規表現で抜いた。
  既評価の 132 accession を除き、骨・骨髄 ×（化療 または 加齢 / OVX）の論文を目視で選別した。OA でない論文は抄録で判断した。
- **補足ファイルの取得**: Springer / Nature は `static-content.springer.com/esm/art%3A<DOI>/MediaObjects/<file>`（論文 HTML は bot 対策で取れないが、ESM は直接取れる）。Elsevier は Crossref で PII を引いて `ars.els-cdn.com/content/image/1-s2.0-<PII>-mmcN.<ext>`。AACR は figshare。
  PMC の `/bin/` と Europe PMC の `supplementaryFiles` は、この環境では bot 対策の HTML が返り、取れなかった。
- **リポジトリ**: Zenodo（`/api/records?q=`）、figshare（`/v2/articles/search`、POST）、Dryad（`/api/v2/search`）、Mendeley Data（`/api/research-data/search?query=`。`search=` や `q=` は無関係の結果を返すので注意）を、下記 27 語（Mendeley は 23 語）で検索した。
- **臨床試験**: ClinicalTrials.gov API v2（`query.cond=breast cancer`、`query.term` に下記 4 通り）。
- 評価: ◎ = Fig4/5 に直接使える、○ = 補助パネル・Supplementary・検証・本文の引用に使える、△ = 参照程度、✕ = 目的に合わない。

### 検索語

- PubMed（各語は `[tiab]`）
  - "chemotherapy-induced bone loss" / "cancer treatment-induced bone loss"
  - 化療語（chemotherapy, doxorubicin, cyclophosphamide, paclitaxel, docetaxel, cisplatin, fluorouracil, methotrexate, 5-FU, busulfan, carboplatin, adriamycin, taxane）を、次とそれぞれ組み合わせた:
    "marrow adipo*" / "marrow fat" / "marrow adiposity" / "bone marrow adipose"、senescence / senescent / SASP × bone / osteocyte* / osteoblast* / "bone marrow stromal" / skeletal、
    "pre-metastatic niche" / "premetastatic niche" / "metastatic niche" × bone / "bone marrow"、dormancy / dormant / "disseminated tumor cells" × bone / "bone marrow"、osteocyte(s)、
    RNA-seq / transcriptom* / single-cell / proteom* / microarray × bone / "bone marrow" / osteoblast* / osteocyte* / stromal、
    "bone marrow niche" / "bone marrow microenvironment" / "bone marrow stromal" / LepR / "CXCL12-abundant" / "sinusoidal endothelial" × single-cell / RNA-seq、
    "bone loss" / osteoclast* / osteoblast* / "mesenchymal stromal" / "bone mineral density" / "bone mass"
  - 年齢・閉経: (aged / old / aging / ageing / elderly / middle-aged / older) × 化療語 × ("bone loss" / "bone marrow stroma*" / osteocyte* / osteoblast* / "bone microarchitecture" / "marrow adipo*" / "bone marrow niche") × (mice / rat)、
    (ovariectomy / ovariectomized / ovariectomised / "estrogen deficiency") × 化療語、(aged / old / aging) × 化療語 × ("bone marrow niche" / "bone marrow stromal" / "mesenchymal stromal" / endosteal)、
    ("aged bone marrow" / "aging bone marrow" / "old bone marrow" / "aged mice" / "old mice" / "skeletal aging") × (niche / stromal / LepR / osteocyte* / "marrow adipocyte*") × オミクス語、
    (ovariectomy / "estrogen deficiency" / menopaus* / postmenopausal) × ("bone marrow" / osteocyte* / "cortical bone") × オミクス語
  - 骨転移: (estrogen / ovariectomy / menopaus*) × ("bone metastasis" / "bone metastases") × (niche / dormancy / "bone marrow")、(aged / "old mice" / aging) × "bone metastas*"、
    ("marrow adipocyte*" / "marrow adiposity" / "marrow fat") × ("bone metastasis" / "breast cancer" / "prostate cancer" / dormancy / metasta*)
  - 臨床: "breast cancer" × (chemotherapy / neoadjuvant / adjuvant) × ("bone turnover" / CTX / P1NP / PINP / sclerostin / RANKL / osteoprotegerin / "bone mineral density" / "bone marker*") × (premenopausal / postmenopausal / menopaus*)、
    "breast cancer" × "bone marrow" × ("disseminated tumor cells" / DTC) × オミクス語、"breast cancer" × ("bone marrow aspirate*" / "bone marrow biops*" / "bone marrow plasma") × chemotherapy
- Europe PMC: 上とほぼ同じ組み合わせ（`[tiab]` なし、本文を含む）。キー: ctibl, chemo_marrow_adipo, chemo_bone_senescence, pmn_bone_chemo, dormancy_bm_chemo, aged_bm_niche, ovx_bm_niche, ovx_osteocyte_seq, chemo_osteocyte, dox_bone_seq, cyclophos_bone, taxane_stroma, chemo_niche_seq, aged_chemo_bone, ovx_chemo, bc_chemo_btm, bc_dtc_marrow_seq, senescence_bone_aging_seq, estrogen_niche_dtc, chemo_bm_adipocyte_tumor
- Zenodo / figshare / Dryad / Mendeley: "doxorubicin bone" / "doxorubicin bone marrow" / "doxorubicin osteocyte" / "chemotherapy bone loss" / "chemotherapy bone marrow stroma" / "chemotherapy bone marrow niche" / "cyclophosphamide bone" /
  "paclitaxel bone marrow" / "docetaxel bone" / "cisplatin bone loss" / "chemotherapy osteocyte" / "chemotherapy senescence bone" / "aged bone marrow niche" / "aging bone marrow stromal single-cell" /
  "ovariectomy bone marrow single-cell" / "ovariectomy osteocyte RNA-seq" / "estrogen deficiency bone marrow" / "bone marrow adiposity chemotherapy" / "breast cancer bone metastasis niche" /
  "breast cancer dormancy bone marrow" / "premetastatic niche bone" / "5-fluorouracil bone marrow niche" / "bone marrow adipocyte aging" / "senescence bone aging mice" /
  "breast cancer chemotherapy bone turnover markers" / "micro-CT doxorubicin" / "chemotherapy induced bone loss mice"（Mendeley は "chemotherapy bone marrow fat" / "methotrexate bone" / "breast cancer chemotherapy bone mineral density" を追加）
- ClinicalTrials.gov: chemotherapy AND (bone turnover OR CTX OR P1NP OR NTX OR sclerostin OR RANKL OR osteoprotegerin)、chemotherapy AND bone mineral density AND (premenopausal OR postmenopausal)、"chemotherapy induced bone loss"、chemotherapy AND bone marrow aspirate AND disseminated tumor cells

### 今回も見つからなかったもの（第二波の「無い」リストへの追加）

- **加齢 / OVX × 化療の骨・骨髄データを Zenodo / figshare / Dryad / Mendeley Data に置いたもの。**
- **骨細胞アポトーシス率を年齢 × 化療で比べた個体値**（Wasi 2024 は成熟マウスのみ、PMID 21938727 は幼若ラットのみ）。
- **化療 × 骨髄 niche（LepR+ / CAR / MALP）の scRNA で、12 か月齢以上のマウスを使ったもの**（今回見た化療側の研究、5-FU × MALP（40057216）、Doxo × CAR（41469373）、busulfan（39139448、scRNA なし）は、すべて若齢〜成体）。
- **乳がん患者の、化療前後の骨代謝マーカー・骨髄の個人単位の公開データ。**
- 白血病以外の固形がん患者の、化療前後の骨髄間質のオミクス（36038666 は骨髄腫・治療関連骨髄腫瘍で、培養細胞）。

## ダウンロードしたもの

- `data/lit/39147353/`（Wasi 2024）: `mmc2.xlsx`（個体値: µCT、ALP、骨標識、TRAP、TUNEL、qPCR）、`mmc3.xlsx`（荷重の要約表）、`mmc1.docx`（補足図表）、`PMC11475016.efetch.xml`（本文）。
- `data/lit/40389522/`（Sci Rep 2025）: `41598_2025_1717_MOESM1_ESM.xlsx`（Fig. 1–8 の個体値）、`MOESM2_ESM.docx`、`PMC12089391.xml`。
- `data/lit/41469373/`（Nat Commun 2025）: `41467_2025_67793_MOESM4_ESM.xlsx`（Source Data）、`MOESM1_ESM.pdf`（補足図）、`PMC12848019.xml`。
- `data/lit/31932453/`（Yao 2020）: `Yao2020_supp_figs.pdf`（figshare 22424466）、`PMC7056549.efetch.xml`。
- `data/lit/42686970/`（Izquierdo 2026）: `11357_2026_2485_MOESM1–4_ESM.xlsx`（年齢・性別ごとのプロテオーム比較表）、`MOESM5_ESM.docx`（Fig. S1）。
- `data/lit/36038666/`（Leukemia 2022）: `41375_2022_1686_MOESM1_ESM.docx`（補足の方法・表・図）。
- `data/lit/38161426/`（iScience 2024）: `nanostring_data.tar.gz`（Zenodo 10045006、肺の NanoString RCC）。
- `data/seq/GSE299773/`、`GSE145477/`、`GSE220679/`、`GSE262441/`、`GSE237307/`: SOFT のみ。
- 取得していないもの: Zenodo 17314594 の生画像（4.8 GB）、Mendeley xs76zjzvc5（4.2 GB）、GSE299773 / GSE145477 の RAW.tar。

# 追加公開データ調査 第二波（加齢 / 閉経 × 化療、骨・骨髄 niche）

調査日: 2026-10-13。目的は打合せ（`data/20261013_meeting.md`）で決まった「背景（加齢 / 閉経）が合わさると化療で骨 niche が動く」を補強できるデータを探すこと。
`docs/datasets.md`（第一波）、`,docs/strategy_20261012.md`、`docs/prelim_chemo.md` で評価済みのデータは再評価していない
（GSE222752 / GSE238244 / GSE205237 / GSE296757 / GSE296758 / GSE296752 / GSE318370 / GSE289491 / GSE276391 / GSE108891 / GSE59894 / GSE208289 / GSE27421 / GSE324491 / NIHMS1776347 ほか）。

## 結論

- **理想の 2×2（加齢または OVX × 細胞傷害性化療を同じ実験の中で比べた、骨・骨髄の RNA-seq / scRNA-seq）は無い。**
  GEO・ArrayExpress・BioProject のどこにも無かった（検索語は末尾に記載）。ヒトでも「高齢と若年の化療前後の骨・骨髄間質」のデータは、乳がんでは無い。
- **トランスクリプトーム以外では、理想の 2×2 が 1 本だけある: PXD078161（MassIVE MSV000101777）。**
  若齢（3 か月）と老齢（20 か月）、雌と雄のマウスに、対照 / ドキソルビシン（DX）/ セノリティクス D+Q（DQ）/ DX+DQ を与え、**脛骨皮質骨（骨髄を除去）のプロテオーム**を測ったもの。
  論文は Izquierdo et al., *Geroscience* 2026（PMID 42686970、CC BY-NC-ND）。著者の結論は「化療による骨の加齢様変化（海綿骨の減少、骨髄脂肪化、老化・炎症遺伝子の上昇）は**雌でより強く、年齢によって異なる**」。
  ただし以下の制約があり、2×2 の交互作用をそのまま検定することはできない。年齢・性別ごとに「DX vs 対照」の変化を並べて比べる使い方になる。
  - 各群 n = 1–3（Sample_Key による。若齢雌は対照 2 / DX 2、老齢雌は対照 2 / DX 2 で、search ファイル名 `…_no51ctrl_…` から老齢雌の対照 1 匹（CB-51）が除外された可能性がある）。
  - 4 つの年齢・性別群は**別々のバッチ**で測定・解析されている（注入量が 500 ng と 1 µg で違い、若齢雌は技術的な重複測定。Spectronaut の解析ファイルも群ごとに別）。
  - 定量できたタンパク質は群あたり 125–400 程度（コラーゲン・ECM が中心）。DX vs 対照で q < 0.05 になったタンパク質は、若齢雌 12 / 老齢雌 2 / 若齢雄 6 / 老齢雄 14（`*_Candidates.xlsx` の `All Proteins` シートから数えた）。
  - DX の用量と時点は MassIVE の記載にない。論文本文は Springer の bot 対策でこの環境から取得できなかった（OA なのでブラウザからは読める）。
- **OVX × ドキソルビシンの 2×2 は表現型（µCT・子宮重量・qPCR）なら有る**: Yao et al., *Cancer Res* 2020（PMID 31932453、PMC7056549）。
  6 週齢の雌 FVB/NJ マウスで Sham / OVX × Veh / Doxo（5 mg/kg を週 1 回 × 4）。Doxo の骨量減少は OVX より大きい。OVX で骨量減少が止まった後（OVX 7 週後）に Doxo を足すと、さらに海綿骨・皮質骨が減る。骨の SASP（IL6, CCL3, MCSF, MMP12 の qPCR）は Doxo で上がり、OVX では上がらない（Supplementary Fig. S5）。
  ただし若齢マウスで、トランスクリプトームではない。

### 次善の組み合わせ案

| 案 | 化療側 | 加齢 / 閉経側 | 並べる単位 | 弱点 |
|---|---|---|---|---|
| **A. 皮質骨・同一実験** | PXD078161 の DX vs 対照（若齢雌・老齢雌） | 同じ PXD078161 の中で、年齢による違い | 同じタンパク質の log2FC を若齢と老齢で並べる（ECM・コラーゲン・解糖系・炎症） | n が極小。バッチが年齢と交絡。RNA ではないので Fig4A（RNA）との突き合わせは pathway 単位になる |
| **B. LepR+ / CAR 細胞（骨髄 niche）** | GSE296758（**雌** LepR+ MSPC、Dox 5 mg/kg × 4 → 8 週後、counts、n=4/4）＋ GSE289491（雌 CAR scRNA、Dox 96 h、n=1）＋ GSE296759（同じ細胞の ATAC、n=3/3） | **GSE130885**（LepR+PDGFRα± など、8–12 週 vs 24–25 か月、featureCounts、n=3）／ **GSE133922**（CAR 細胞、2 か月 vs 2 年、normalized counts、n=4。LPS・poly I:C の炎症刺激群もある）／ GSE227144（LepR+ MSC、3 か月 vs 18 か月、FPKM、n=3） | 同じ細胞種で DEG / Hallmark の重なりを見る（炎症・SASP、脂肪化、骨芽分化の低下、Cxcl12 / Kitl などの niche 因子、ER 標的） | 加齢側の性別は SOFT に記載なし。別データどうしの比較になる |
| **C. ヒト骨髄・年齢で層別** | GSE212615（AML 患者のコア骨髄生検、寛解導入化療（アントラサイクリン＋シタラビン）の前と 2 週後のペア、11 例、raw counts） | 同じデータの年齢（24–68 歳）で層別する（52 歳以上 7 例 / 46 歳以下 4 例） | 化療前後の変化量が高齢層で大きいかを見る（粗い年齢 × 化療） | AML でレジメンが違う。2 週後は骨髄が低形成になり、細胞組成の変化が支配的になる。n が小さい |
| D. 表現型の文献値 | Yao 2020（OVX × Doxo）、Izquierdo 2026（年齢 × Dox）、Sci Rep 2025（雌・TAC、終了 12 か月後まで追跡） | 同左 | 本文・Discussion で「加齢 / 閉経下でも化療で骨量が減り、SASP が出る」の根拠として引用する | 図に載せる定量データとしては使いにくい |

推奨: **A を Fig4B の差し替え候補**（同じ実験の中で年齢と化療の両方を持つ唯一のデータ）とし、**B を補助**（RNA で、雌・臨床に近い反復投与で、骨髄 niche のデータ）に使う。
C は Supplementary の候補。D は本文の引用に回す。

---

## 方法

- GEO: E-utilities（`esearch db=gds` → `esummary`）。年齢語 12 語 × 薬剤語 15 語 × 組織語の総当たり（計 180 検索）、薬剤語 12 語 × niche 語 19 語、年齢 / OVX 語 8 語 × niche 語 14 語で集め、
  タイトルと summary で一次選別した。候補は FTP の `GSExxxxxx_family.soft.gz` を取得し、Series の summary / overall_design / supplementary と全 Sample の characteristics を集計した。
  下表の種・性別・年齢・組織・処置・n・ファイル名は SOFT の記載から取った。**記載がない項目は「記載なし」とし、推測では埋めていない。**
- PubMed（`esearch db=pubmed` → `efetch`、`elink` で GEO と照合）。本文は NCBI PMC（`efetch db=pmc`）と Europe PMC で読んだ。
- ProteomeXchange（PROXI API `proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets?search=…`）、PRIDE API、MassIVE FTP（FTPS）。
- Metabolomics Workbench（REST `study/study_title/<語>/summary`）。
- ArrayExpress / BioStudies（`www.ebi.ac.uk/biostudies/api/v1/search`）。
- 評価: ◎ = Fig4/5 に直接使える、○ = 補助パネル・Supplementary・検証に使える、△ = 参照程度、✕ = 目的に合わない。

---

## 1. 加齢 / OVX × 化療を同じ実験の中で比べたもの（最優先）

| accession | 種・性別・年齢・組織 | 群・処置 | n | 処理済みファイル | PMID | 評価 |
|---|---|---|---|---|---|---|
| **PXD078161**（MassIVE MSV000101777） | マウス**雌雄**、**3 か月 / 20 か月**、脛骨の皮質骨（骨端を除き、骨髄を生理食塩水で洗い流したもの。`Methods.docx` の記載） | 対照 / DX / DQ / DX+DQ（年齢・性別ごと） | Sample_Key: 若齢雌 2 / 2 / 3 / 3、老齢雌 2 / 2 / 3 / 3、若齢雄 5 / 3 / 6 / 6（CB-41〜49・87〜89 と CB-50・81〜86・90 の 2 ブロックの合計。後者は Sample_Key では 3 か月と記載され、raw の VG3_yMales に CB-85 が含まれる）。老齢雄は Sample_Key に記載なし（raw / search には VG4_aMales として CB-74〜80 がある） | `quant_stats/VG{1-4}_*_Candidates.xlsx`（Spectronaut の群間比較表。全組み合わせの log2 比・p・q）、`metadata/Sample_Key.xlsx`、`other/Methods.docx`、raw（wiff）、search（sne、0.3–68 GB） | 42686970 | ◎ **唯一の「年齢 × 化療 × 骨」の同一実験**。プロテオーム。制約は「結論」に記載 |
| Yao 2020（データ寄託なし） | マウス雌 FVB/NJ 6 週（OVX 実験）、C57BL/6 16 週（SASP の qPCR）、雄も一部 | Sham / OVX × Veh / Doxo 5 mg/kg 週 1 × 4。OVX 7 週後に Doxo を足す実験。PTX・放射線の群もある | 群あたり 4–9 | 図の値のみ（µCT、子宮重量、qPCR） | 31932453 | ○ **OVX × Doxo の 2×2（表現型）**。Doxo の骨量減少は OVX より大きく、OVX 後にも上乗せされる。SASP は Doxo で上がり OVX では上がらない |
| Izquierdo 2026（PXD078161 の論文） | 同上 | 同上。µCT、骨髄脂肪化、老化・炎症・骨形成遺伝子の qPCR | 記載は本文（未取得） | 図の値 | 42686970 | ○ 表現型と qPCR は本文の図。CC BY-NC-ND なので、ken がブラウザから PDF を取れば数値を確認できる |
| Wasi 2024（寄託なし） | マウス C57BL/6J 雌雄、成熟（約 50 週） | Dox 2.5 / 5 mg/kg × 6 回（2 週）± Yoda1・荷重。年齢をそろえた無処置対照（雌 4） | 3–8（雌雄を合算） | 図の値（µCT、骨細胞 TUNEL%、TRAP、ALP） | 39147353（PMC11475016） | △ 成熟マウスで Dox により**骨細胞アポトーシス（TUNEL）が増える**。ただし若齢群がないので 2×2 にはならない |
| PMID 41711309 | ラット雌、OVX ＋ 脛骨インプラント | OVX 後にシスプラチン 5 / 2.5 mg/kg ± タモキシフェン。非 OVX の対照あり | 計 100 | 図の値（組織形態計測、TRAP / OCN / RUNX2 の免疫染色） | 41711309 | △ OVX × シスプラチン。インプラント周囲の骨で、シスプラチン非投与の OVX 群との比較が主（OVX 有無 × 化療の完全な 2×2 ではない） |
| PMID 19890044 | マウス、中年・老齢 | シクロホスファミド ± 性ステロイド除去（去勢） | — | 図の値 | 19890044 | ✕ 骨髄の造血・リンパ球回復が対象で、骨・niche ではない |
| PMID 32601899 | マウス、若齢 / 老齢 | 5-FU、day 4 / day 10。VSEL・HSC・MSC の数 | — | 図の値（フロー） | 32601899 | △ **若齢 / 老齢 × 5-FU の骨髄**。MSC（CD44+）の増殖が老齢で高い。性別は抄録に記載なし。トランスクリプトームなし |

**GEO / ArrayExpress / BioProject に、年齢または OVX と化療を同じ実験で組んだ骨・骨髄のトランスクリプトームは無かった。**
最も近い GSE302366（LepR+ 骨髄間質の scRNA、若齢 / 老齢 × 骨髄除去「手術」による損傷）は、損傷が化療ではなく手術。

---

## 2. 化療 × 骨髄 niche（第一波以外）

| accession | 種・性別・年齢・組織 | 群・処置 | n | 処理済みファイル | PMID | 評価 |
|---|---|---|---|---|---|---|
| **GSE296759**（SuperSeries GSE296761、GSE296758 の兄弟） | マウス**雌**（サンプル名に female と記載）、8–10 週で処置、Lepr-Cre;tdTomato の骨髄 MSPC | PBS × 4 / Dox 5 mg/kg × 4（静注、1 週間隔）、8 週後 | 3 / 3 | ATAC-seq の narrowPeak（サンプルごと） | 未収載 | ○ GSE296758（RNA）と同じ細胞・同じ処置のクロマチン。「トランスクリプトーム以外」の節にも記載 |
| GSE123078（SubSeries of GSE108892） | マウス（性別記載なし）、7 週（Lepr）/ 9 週（VE-cad）、骨髄の Lepr+ と VE-cad+ 細胞 | PBS / 5-FU 150 mg/kg、4 日目 | 細胞種と処置ごとに 1 | サンプルごとの raw counts txt | 30971824 | △ 第一波の GSE108891 の検証実験。n=1 |
| GSE104887 | マウス（性別・年齢記載なし）、N-cadherin+ の骨・骨髄間質、Col2.3+ 骨芽細胞、Cxcl12+、Lepr+、マクロファージ | Ncad 細胞: 無処置 4 / 5-FU 3（骨・骨髄それぞれ） | 3–6 | FPKM xlsx のみ | 30650358 | △ 5-FU（CMF / FEC の成分）による骨表面 niche の変化。counts なし |
| GSE135866 | マウス C57BL/6（性別・年齢記載なし）、骨髄 MSC（初代培養） | エタノール（mock）/ ENU 100 mg/kg | 3 / 3 | `GSE135866_MSC_ENU_tpmCounts.txt.gz` | 未収載 | ✕ ENU はアルキル化変異原で、乳がんのレジメンではない。培養細胞 |
| GSE291050 / GSE160091 | ヒト前立腺間質細胞（PSC27）をマウスに異種移植 | ミトキサントロン ± 薬剤 | 3 | DESeq2 表 / count | 41792133 / 34263179 | ✕ 骨ではない |
| GSE216842 / GSE190278 | ヒト初代間質細胞（線維芽細胞・HUVEC・前立腺間質） | ブレオマイシン（in vitro） | 3–6 | DESeq2 表 | 41493135 | ✕ in vitro、骨ではない |
| GSE190148 | マウス HSC | 5-FU 後の H3K27ac CUT&Tag | 8 | — | — | ✕ 造血系 |
| PMC12771323（寄託なし） | マウス雌 BALB/c と BALB/c nude、7 週 | 生理食塩水 / Dox 4 or 6 mg/kg 週 1 × 4 | — | 図の値（µCT、骨髄細胞群のフロー） | — | △ 若齢雌の Dox による骨と骨髄細胞群の変化。骨転移モデル用の基礎データ |

---

## 3. 加齢 / OVX × 骨髄 niche（化療側と同じ細胞種で並べる用）

| accession | 種・性別・年齢・組織 | 群 | n | 処理済みファイル | PMID | 評価 |
|---|---|---|---|---|---|---|
| **GSE130885**（SubSeries of GSE130899） | マウス C57BL/6（性別記載なし）、骨髄間質を FACS（100 細胞の bulk）: ALCAM+Sca1-、LepR+PDGFRα-、LepR+PDGFRα+、Sca1+PDGFRα+、血管内皮 | 若齢（8–12 週）/ 老齢（24–25 か月） | 細胞種ごとに 3 / 3 | `GSE130885_Valletta_et_al_RNAseq_stromal_featurecounts.txt.gz`（**raw counts**） | 未収載 | ○ **LepR+ 細胞の加齢を raw counts で持つ**。GSE296758（雌 LepR+ MSPC × Dox）と同じ細胞種で並べられる |
| **GSE133922** | マウス C57BL/6JRj（性別記載なし）、CAR 細胞・PDGFRα+Sca1+・類洞内皮・動脈内皮 | 生後早期 / 2 か月 / 2 年。さらに 2 か月で LPS・poly I:C | 細胞種・年齢ごとに 4 | `GSE133922_ALLsamples.normalizedCounts.txt.gz`、FPKM | 31801092 | ○ **CAR 細胞の加齢と炎症刺激を同じ系で持つ**。「加齢 niche の炎症」を CAR で示せる |
| GSE227144（SuperSeries GSE227148） | マウス（性別記載なし）、骨髄 MSC | 3 か月 / 18 か月 | 3 / 3 | FPKM のみ | 37037837 | △ counts なし。内皮は GSE227145（同じデザイン） |
| GSE169162（GSE168613 / GSE168586 / GSE168661 / GSE189217） | マウス、骨髄の間質・内皮（scRNA、Smart-seq と 10x）、アレイ | 若齢 / 老齢 | 多数 | h5 ほか | 36650381 | ○ タイトルは「IL-1 による間質 niche の炎症が造血の加齢を駆動する」。加齢 niche の炎症の参照として Fig5 の議論に使える |
| GSE316278 | マウス、骨髄 HSC と CD45-Lin-Ter119- 間質の scRNA | 若齢雄 / 若齢雌（8 週）/ 老齢 vehicle / 老齢セノリティクス 753b（100 週） | 各 1 ライブラリ | 10x | 42523366 | △ n=1。老齢側の性別記載なし |
| GSE299773 | マウス、Lin-CD45-CD71- 骨髄間質の scRNA | 2 / 12 / 24 か月 | 各 1 | 10x | 40063797 | △ n=1。12 か月（閉経相当の年齢）を含むのは利点 |
| GSE315430 | マウス、骨髄間質の scRNA | 成体 / 老齢 | 各 1 | 10x | — | △ Adipo-CAR の加齢変化。n=1 |
| GSE302366 | マウス、LepR+ 骨髄間質の scRNA（BD Rhapsody） | 対照 / 損傷 1 回 / 損傷 2 回（15 週）/ 老齢（18 か月）/ 老齢 + 損傷 | 各 1 | barcodes など | 42618571 | △ 「老齢 × 損傷」はあるが、損傷は骨髄除去の手術で化療ではない |
| GSE173371 | マウス、LepR+ 骨格幹・前駆細胞 | 12 週 / 52 週 / 52 週 + サリチル酸 | 3 / 3 / 3 | 正規化 counts xlsx | 未収載 | △ 中年（52 週）との比較 |
| GSE305558 | **ヒト**、整形外科手術時の骨髄から内皮と MSC を分取した scRNA | 若年（28–48 歳）5 / 高齢（58–79 歳）7、男女を含む | 12 | h5 | 41954292 | ○ ヒトの加齢骨髄 niche（内皮・MSC）。化療データはないが、ヒトでの「加齢 niche の炎症・線維化」の参照 |
| GSE326234 | **ヒト閉経後女性**、骨髄の CAR 細胞 scRNA | 無治療 2 / テリパラチド 2（1 ライブラリにまとめて登録） | 1 | mtx と cell metadata | — | △ 閉経後女性の CAR 細胞の参照として使える |
| GSE272266 / GSE272265 | マウス LSK（造血） | 雌雄、OVX + E2 / DHT / vehicle | 2 | — | 40178977 | ✕ 造血細胞。niche ではない |

**OVX で骨髄 niche（LepR+ / CAR / 内皮 / 骨髄脂肪細胞）を in vivo で分取したトランスクリプトームは見つからなかった**（OVX の骨髄は GSE238244 の脛骨丸ごと、GSE194075 の培養間質、GSE343407 のラット scRNA（第一波で評価済み）だけ）。

---

## 4. ヒト: 化療前後の骨髄 × 年齢

| accession | 種・性別・年齢・組織 | 群 | n | 処理済みファイル | PMID | 評価 |
|---|---|---|---|---|---|---|
| **GSE212615**（SuperSeries GSE213210） | ヒト AML、**FFPE のコア骨髄生検**（固相の間質を含む）。性別・年齢はサンプルごとに記載: 24 / 26 / 31 / 46 / 52 / 55 / 58 / 65 / 66 / 66 / 68 歳。女性 4–5 例（R5 は 1st が F、2nd が M と SOFT 内で食い違う） | 診断時 / 寛解導入化療（アントラサイクリン＋シタラビン）の 2 週後（ペア）。奏効 6 / 非奏効 5 | 11 ペア | サンプルごとの raw counts txt（Ensembl ID、`GSE212615_RAW.tar`） | 36354085 | ○ ヒトで同じ患者の化療前後の骨髄（間質を含む）を、年齢つきで持つ唯一のデータ。年齢で層別すれば粗い年齢 × 化療になる。制約は「次善の組み合わせ案」C に記載 |
| GSE213056 | 同上の Agilent アレイ版 | 同上 | 11 ペア | RAW | 36354085 | △ RNA-seq 版（GSE212615）があるので補助 |
| GSE279904 | ヒト AML、寛解導入後の骨髄吸引液の scRNA。年齢・性別・レジメンを記載 | CR 12 / 不応 10（化療後のみ） | 22 | h5 | 40691440 | △ 化療前のデータがない。リンパ球が対象 |
| GSE318370（第一波で評価済み） | ヒト ALL の骨髄培養間質、診断時 / 寛解導入 4 週後 | — | 8 / 5 | counts | 未収載 | （再掲）年齢は SOFT に記載なし |

**化療を受けた乳がん（固形がん）患者の骨・骨髄間質を、化療前後で、年齢つきで測ったトランスクリプトームは無かった**（第一波と同じ結論）。

---

## § トランスクリプトーム以外

### (1) プロテオーム / リン酸化プロテオーム（ProteomeXchange / PRIDE / MassIVE）

| accession | 内容 | n | ファイル | PMID | 評価 |
|---|---|---|---|---|---|
| **PXD078161**（MassIVE MSV000101777） | 若齢 / 老齢 × 雌雄 × 対照 / DX / DQ / DX+DQ、脛骨皮質骨（骨髄除去）、DIA（ZenoTOF 7600+） | 1–3 / 群 | §1 を参照。`data/seq/PXD078161/` に取得済み | 42686970 | ◎ §1 を参照 |
| PXD007048 | ヒト骨髄の細胞種別プロテオーム、加齢（同名の GEO GSE115348 / GSE114961 がある） | — | PRIDE | — | △ ヒト加齢骨髄の参照。化療なし |
| PXD054205 | 老齢マウス骨のインスリンシグナル（リン酸化プロテオーム） | — | PRIDE | — | △ 加齢骨。化療なし |
| PXD068037 | ApoE4 雌マウスの骨細胞機能（Buck Institute、同じ骨プロテオミクス系） | — | MassIVE | — | ✕ 目的外。PXD078161 と同じ測定系の参照 |
| PXD067728 / PXD082189 / PXD035745 ほか | OVX の骨・BMSC に生薬などを投与した実験 | — | — | — | ✕ 化療を含まない。OVX 単独でも介入が主目的 |

ProteomeXchange で「doxorubicin」「cyclophosphamide」「cisplatin」「paclitaxel」「fluorouracil」「chemotherapy」をそれぞれ検索し、骨・骨髄・間質を含むものを拾った。
**PXD078161 以外に、化療 × 骨・骨髄のプロテオームは無かった**（ヒットは腫瘍細胞・卵巣がん・骨肉腫・筋）。

### (2) ATAC-seq / ChIP / CUT&Tag / メチル化

| accession | 内容 | n | ファイル | PMID | 評価 |
|---|---|---|---|---|---|
| **GSE296759** | 雌 LepR+ 骨髄 MSPC、PBS × 4 / Dox 5 mg/kg × 4、8 週後の ATAC-seq | 3 / 3 | narrowPeak | 未収載 | ○ GSE296758（同じ細胞の RNA）と組み合わせ、Dox が niche 細胞の骨芽系遺伝子座を閉じるかを見られる |
| GSE240286（SuperSeries GSE240292） | 若齢（3 か月）/ 中年（12 か月）マウスの SSPC（継代 1 の培養）の ATAC-seq | 3 / 3 | bigWig | 未収載 | △ 加齢側の ATAC。培養細胞なので GSE296759（新鮮分取）とは条件が違う |
| GSE264569 | NF-κB 活性化した骨髄環境での HSPC の ATAC | 3 | raw counts | 40875292 | ✕ 造血系 |
| GSE190148 | 5-FU 後の HSC の H3K27ac CUT&Tag | 8 | — | — | ✕ 造血系 |

**化療 × 骨・骨細胞のメチル化データ、OVX × 骨・骨髄 niche の ATAC / メチル化データは見つからなかった。**

### (3) 代謝・サイトカインパネル

- Metabolomics Workbench を doxorubicin / cyclophosphamide / chemotherapy / bone marrow / bone / ovariectomy / osteoporosis で検索した。**化療 × 骨・骨髄、加齢 / OVX × 化療の代謝データは無かった**（ヒットは腫瘍細胞株、脳、マクロファージ、皮質骨の肥満モデル ST004373 など）。
- 動物のサイトカインパネル: Yao 2020 の骨 SASP は qPCR（Luminex ではない）。
- ヒト骨髄のサイトカイン: PMID 7593450 は、乳がん手術時の腸骨骨髄の培養上清で、閉経前 12 / 閉経後 18 / エストロゲン補充中止直後 5 の IL-1α・TNF-α・IL-6・PGE2・GM-CSF を測ったもの（化療前）。閉経・エストロゲン中止で骨髄の炎症性サイトカインが上がるという 1st hit 側の古典的根拠として引用できる。

### (4) 論文の定量データ（PubMed で探したもの。表の場所を記録）

| PMID | 種・条件 | 加齢 / OVX × 化療 | 定量項目と所在 | 評価 |
|---|---|---|---|---|
| **42686970**（Izquierdo 2026 Geroscience） | マウス雌雄、3 か月 / 20 か月 × DX ± D+Q | **年齢 × 化療（同一実験）** | µCT、骨髄脂肪、qPCR（老化・炎症・Alpl / Runx2 / Dmp1）、皮質骨プロテオーム。本文の図（PDF は未取得、CC BY-NC-ND） | ◎ |
| **31932453**（Yao 2020 Cancer Res、PMC7056549） | マウス雌 6 週、Sham / OVX × Veh / Doxo。OVX 後に Doxo を追加 | **OVX × 化療（同一実験）** | Fig. 1（OVX × Doxo の BV/TV、皮質骨厚、子宮重量）、Fig. 2F–H（OVX 後の Doxo 追加）、Fig. 4A と Supplementary Fig. S5（SASP の qPCR、OVX では上がらない）、Supplementary Fig. S8（16 週齢の OVX） | ○ |
| 39147353（Wasi 2024 Bone、PMC11475016） | マウス雌雄、約 50 週、Dox | 加齢（成熟）下の化療のみ | 骨細胞の TUNEL 陽性率、TRAP、ALP、µCT（Fig. 4–6、Supplementary に雌雄別の値） | △ 骨細胞アポトーシスの定量がある数少ない例 |
| 40389522（Sci Rep 2025、PMC12089391） | マウス**雌** C57BL/6 16 週・BALB/cJ 19 週、**AC / TAC（ドキソルビシン＋シクロホスファミド＋ドセタキセル）** ± D+Q / ピペルロングミン | 化療後 6 / 12 か月まで追跡（化療 → 加齢の順） | µCT、皮質骨の qPCR（老化マーカー）。「全データは本文と Supplementary」と記載 | ○ **乳がんの実レジメン × 雌**。化療直後は老化マーカーが上がり、骨量の減少は 12 か月後まで残るが、そのころには老化マーカーは検出されない。セノリティクスでは防げない |
| 41469373（Nat Commun 2025、PMC12848019。GSE289491 の論文） | マウス雌 12 週、Doxo / PTX | 化療のみ | SA-β-gal、p16 / tdTomato、TRAP、BFR、RANKL / OPG の qPCR。Source Data ファイルあり | ○ 化療による老化は**骨髄の脂肪系細胞（CAR・骨髄脂肪細胞）に限られ、骨細胞では起きない**。骨細胞濃縮分画では Rankl / Opg 比が変わらない（Fig. S10） |
| 21938727（J Cell Physiol 2012） | ラット幼若、MTX | 化療のみ | 骨細胞アポトーシス 4.3 倍、TRAP+ 破骨細胞 1.8 倍、IL-6 / IL-11 | △ MTX、幼若 |
| 32601899（Stem Cell Rev Rep 2020） | マウス若齢 / 老齢 × 5-FU | 年齢 × 化療 | 骨髄の VSEL / HSC / MSC 数、CFU（フロー） | △ 骨髄の幹細胞数のみ |
| 41711309（J Periodontol 2026） | ラット雌 OVX、シスプラチン ± タモキシフェン | OVX × 化療（インプラント周囲） | 組織形態計測、TRAP / OCN / RUNX2 | △ |
| **25923354**（PLoS One 2015） | ヒト乳がん、閉経前 29 / 閉経後 12、術前化療 4 サイクル | **閉経状態 × 化療（ヒト）** | 血清 PINP / BAP / ICTP / β-CrossLaps、E2、FSH、25(OH)D を各サイクル後に測定 | ○ 閉経状態に関係なく、化療で吸収マーカー↑・形成マーカー↓ |
| 37642751（Support Care Cancer 2023） | ヒト乳がん 118 例、月経状態別、EC / TC / TEC / EC-T(H) | 閉経状態 × 化療（ヒト） | BMD と骨代謝指標の変化 | ○ BMD 低下は主に閉経前群（化療による卵巣機能低下とみられる） |
| 37410202（J Bone Miner Metab 2023） | ヒト閉経後乳がん 109 例、AC → T（デキサメタゾン併用） | 閉経後 × 化療（単群） | BMD、CTX、PINP、炎症性サイトカイン、酸化ストレス | ○ 閉経後で化療前後の炎症性サイトカインを測った例 |
| 25575458（Mol Cell Endocrinol 2015） | ヒト閉経前乳がん 40 例、化療誘発卵巣機能不全 | 化療 → 閉経 | OPG、BMD、FSH | △ |

**骨細胞アポトーシス率（TUNEL / cleaved caspase-3）や空の骨小腔率を、加齢 / OVX × 化療の 2×2 で定量した論文は見つからなかった。**
骨細胞 TUNEL は Wasi 2024（成熟マウス × Dox、若齢群なし）と PMID 21938727（幼若ラット × MTX）だけ。

---

## 検索した語と、ヒットしなかった組み合わせ

### 検索語

- GEO（`db=gds`、`gse[ETYP]`）
  - 年齢・閉経: aged / "old mice" / "aging" / "young and old" / "old and young" / ovariectomy / ovariectomized / "estrogen deficiency" / "middle-aged" / geriatric / elderly / senescence / menopause / postmenopausal / "18 month" / "20 month" / "24 month"
  - 薬剤: doxorubicin / adriamycin / chemotherapy / "chemotherapy-induced" / cyclophosphamide / fluorouracil / 5-FU / paclitaxel / docetaxel / cisplatin / carboplatin / methotrexate / epirubicin / "cytotoxic" / busulfan /（参考として irradiation）
  - 組織: "bone marrow" / bone / osteoblast / osteocyte / stromal / niche / skeletal / "bone marrow stromal" / "bone marrow niche" / "bone marrow microenvironment" / "mesenchymal stromal" / LepR / "leptin receptor" / CXCL12 / "CAR cells" / "bone marrow endothelial" / sinusoidal / "bone marrow adipocyte" / "marrow adiposity" / "bone marrow adipose" / "bone loss" / "bone marrow macrophage" / "bone marrow mesenchymal" / single-cell / female / femur / tibia
  - アッセイ: ATAC-seq / methylation / ChIP-seq / bisulfite / CUT&Tag
  - ヒト: "before and after" / "pre- and post" / paired / "post-treatment" × biopsy / trephine / mesenchymal / stromal
- PubMed: (aged OR old OR aging OR elderly) × (各薬剤) × (bone loss / osteocyte / osteoblast / "bone marrow stromal" / "bone marrow niche" / "marrow adipo*") × young、
  ("aged mice" OR "old mice" OR "aged female" OR "middle-aged mice") × (各薬剤) × (bone / bone marrow / osteocyte / skeletal)、"chemotherapy-induced bone loss"、
  (osteocyte apoptosis / osteocyte senescence / "empty lacunae") × (各薬剤)、(ovariectomy / ovariectomized / ovariectomised) × (各薬剤) × (bone loss / "bone metastasis" / osteoclast / "bone marrow")、
  乳がん × 化療 × 閉経前・閉経後 × (RANKL / OPG / bone turnover markers / CTX / P1NP / sclerostin / IL-6 / cytokine)、chemotherapy-induced menopause / amenorrhea / ovarian failure × bone
- ProteomeXchange（PROXI の search）: doxorubicin / chemotherapy / cyclophosphamide / cisplatin / paclitaxel / fluorouracil / senolytic / "bone marrow" / bone / osteocyte / osteoporosis / ovariectomy / ovariectomized / "skeletal aging"
- Metabolomics Workbench（タイトル）: doxorubicin / cyclophosphamide / chemotherapy / bone marrow / bone / ovariectom / osteoporosis
- BioStudies: "aged doxorubicin bone" / "aged chemotherapy bone marrow stroma" / "ovariectomy doxorubicin" / "ovariectomy chemotherapy bone" / "5-fluorouracil aged bone marrow niche" / "cyclophosphamide bone marrow stromal" / "doxorubicin bone marrow mesenchymal"

### 見つからなかった組み合わせ（「無い」ことの根拠）

- **加齢 × 化療の骨・骨髄トランスクリプトーム**（bulk / sc / sn、マウス・ラット・ヒト）。年齢 × 化療を同じ実験で組んだものはプロテオームの PXD078161 だけ。
- **OVX × 化療の骨・骨髄のオミクス**（トランスクリプトーム・プロテオーム・エピゲノム・代謝のどれも）。OVX × Doxo は表現型と qPCR（Yao 2020）のみ。
- **OVX × 骨髄 niche（LepR+ / CAR / 内皮 / 骨髄脂肪細胞）を in vivo で分取したトランスクリプトーム**。
- **化療 × 骨髄 niche で、高齢マウス（12 か月以上）を使ったもの**（化療側はすべて 7–16 週の若齢）。
- **ヒトの固形がん（乳がんを含む）患者の骨・骨髄間質の、化療前後ペアのトランスクリプトーム**。ヒトの化療前後ペアは AML（GSE212615 / GSE213056）と ALL（GSE318370）だけ。
- **骨細胞アポトーシス率（TUNEL / caspase-3）や空の骨小腔率を、加齢 / OVX × 化療の 2×2 で定量した論文**。
- **化療 × 骨・骨髄の代謝データ、メチル化データ**。
- **化療 × 骨・骨髄のプロテオーム**（PXD078161 以外）。

## ダウンロードしたもの

- `data/seq/PXD078161/`: `Sample_Key.xlsx`、`Methods.docx`、`params.xml`、`VG1–VG4_*_Candidates.xlsx`（Spectronaut の群間比較表）、`PXD078161_proteomexchange.xml`。raw（wiff）と search（sne、0.3–68 GB）は取得していない。
- `data/seq/GSE212615/`: SOFT、`GSE212615_RAW.tar`（5.6 MB）と、それを展開した `raw/`（22 サンプルの raw counts）。
- `data/seq/GSE130885/`: SOFT、`GSE130885_Valletta_et_al_RNAseq_stromal_featurecounts.txt.gz`。
- `data/seq/GSE133922/`: SOFT、`GSE133922_ALLsamples.normalizedCounts.txt.gz`。
- `data/seq/GSE296759/`、`data/seq/GSE227144/`、`data/seq/GSE213056/`: SOFT のみ。

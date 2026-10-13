# 追加公開データ調査 第三波（in vitro の骨系細胞 × 化療 / ヒト化療後骨髄 / メタデータ直接検索）

調査日: 2026-10-13。第一波（`docs/datasets.md`）・第二波（`docs/datasets2.md`）で GEO のタイトル・summary の総当たりはほぼ尽きたので、
今回は (1) in vitro / ex vivo の骨系細胞 × 化療薬、(2) AML 以外のヒト化療後骨髄、(3) SRA / BioSample / ARCHS4 / recount3 のサンプル属性からの直接検索、の 3 方向に変えた。
第一波・第二波で評価済みのもの（GSE324491 / GSE296757–9 / GSE318370 / GSE212615 / GSE289491 / GSE108891 / GSE104887 ほか）は再評価していない。

## 結論

- **骨細胞（MLO-Y4 / Ocy454 / IDG-SW3 / 初代骨細胞）× 化療薬のトランスクリプトームは無い。** GEO・BioSample・SRA・ARCHS4（最新版、ヒト 109.9 万 / マウス 122.5 万サンプル）・recount3（マウスは走査未完）のどれにも無かった。
  骨細胞株の公開データは PTH・機械刺激・低酸素・転移（GSE276373、GSE244453）のみ。
- **ヒトの骨系細胞 × 化療薬（または同じ機序の DNA 傷害薬）は 3 本見つかった。** どれもドナーの年齢・性別は SOFT に記載がない。
  - **◎ GSE284666（SenCat、Mol Cell 2026）**: ヒト初代骨芽細胞（NHO）と骨髄 MSC に **エトポシド**（TOP2 阻害。ドキソルビシンと同じ標的）または放射線で老化を誘導（8 日）。各 n=4、サンプルごとの featureCounts。
    「DNA 傷害が骨芽細胞・骨髄 MSC のどの軸（p53 / p21 / SASP / 骨芽分化 / niche 因子）を動かすか」を、ヒトの骨系細胞で、counts から DE できる唯一のデータ。
  - **◎ GSE135323（Toxicol Sci 2020）**: 組織工学で作った 3D の骨（足場上のヒト骨芽細胞、TE-Bone）に **ドキソルビシン 4.16 µM** を 2 週間（3 日投与 / 4 日休薬 × 2）。対照 3 / Dox 3。
    **ヒト骨芽細胞 × ドキソルビシンそのもの**の唯一のデータ。ただし TempO-Seq の S1500+ パネル（2,977 プローブ）で全トランスクリプトームではない。
  - ○ GSE17860（PLoS One 2012）: ヒト初代骨芽細胞（PromoCell）に **メルファラン**（アルキル化薬。CTX と同系）50 µg/ml を 6 時間、n=4/4（Agilent、共通リファレンス）。化療を受けた骨髄間質の培養上清を当てた間接群（n=3 ペア）もある。
- **マウスの骨系細胞 × 化療薬（in vitro）**: ○ GSE283633（骨髄 BMSC、エトポシドで老化、n=4/4、counts）、○ GSE301693（骨髄間質細胞株 ST2、シスプラチン 2 µM 24 h、n=3/3、counts）、
  △ GSE290573（RAW264.7＝破骨前駆の株、ドキソルビシン 100 nM 24 h、n=5/5、FPKM）。破骨細胞 × 化療薬は他に無い（GSE243008 はタイトルに methotrexate とあるが MTX 群は寄託されていない）。
- **ドナーの年齢・性別が分かるヒト MSC × 化療薬は、メチル化の 1 本だけ**: △ GSE227160（女性ドナーの骨髄 MSC にドキソルビシンで老化誘導、EPIC アレイ、n=4、年齢記載なし）。
  RNA では、年齢・性別つきのヒト骨髄 MSC は**化療なしの加齢データ**として ○ GSE39540（61 ドナー、女性 55 歳以上 28 / 45 歳未満 11、U133A 2.0）と GSE35959（高齢女性・骨粗鬆症ドナー）がある。
  「エトポシドで骨髄 MSC に出る変化（GSE284666）」と「閉経後女性ドナーの骨髄 MSC で上がっている遺伝子（GSE39540）」の重なりを見る、という別データ同士の比較ならできる。
- **ヒトの化療後骨髄（AML 以外）で、骨髄の間質を含み、年齢・性別が分かるものは無い。** 見つかったのは
  △ GSE220946（小児神経芽腫の骨髄 scRNA 15 例。治療あり 11 / なし 4 と記載があるが、レジメン・年齢・性別・採取時期の記載なし、転移の有無と交絡）だけ。
  乳がんの術前化療中の骨髄、リンパ腫・固形がんの自家移植前の骨髄、化療誘発性骨減少のヒト骨生検のトランスクリプトームはいずれも無かった。
- **メタデータの直接検索（狙い 3）で新たに出たのは GSE284666（BioSample の属性から）と SRA のみの SRP655200（△）だけ。** ARCHS4 最新版のサンプル属性に化療薬名を持つシリーズはヒト 189 本・マウス 15 本で、骨系は GSE135323 と既知のもの（GSE104887 / GSE108891 = 5-FU の in vivo）しかなかった。
  GEO のタイトル・summary 検索で拾えていなかったデータは、ほぼ無いと考えてよい。

### 使い方の案（生体データとの突き合わせ）

| 案 | in vitro 側 | 生体側（既存） | 並べる単位 | 弱点 |
|---|---|---|---|---|
| **A. ヒト骨系細胞の DNA 傷害応答** | GSE284666 の NHO・BM-MSC「エトポシド vs 対照」（n=4/4）＋ GSE135323 の TE-Bone「Dox vs 対照」（n=3/3、パネル内） | GSE205237（雄の皮質骨、Dox）・GSE296757/8（雌の骨芽細胞 / LepR+、Dox × 4）・PXD078161（年齢 × Dox のプロテオーム） | 同じ遺伝子セット（p53 標的、アポトーシス、SASP、骨芽分化、Cxcl12/Kitl）の log2FC を並べる。in vitro で強く出る軸が生体でどこまで残るかを見る | 細胞種・用量・時間がそろわない。エトポシドは乳がんのレジメンではない（機序が Dox と同じ TOP2 阻害であることを根拠にする） |
| **B. ヒト骨髄 MSC の「化療」と「閉経後」** | GSE284666 の BM-MSC エトポシド | GSE39540 の女性ドナー 55 歳以上 vs 45 歳未満（化療なし） | DEG / Hallmark の重なり（老化・炎症・脂肪化・骨芽分化） | 別データ・別プラットフォーム。化療と年齢の交互作用は検定できない |
| C. アルキル化薬 | GSE17860 の HOB メルファラン 6 h（＋間質経由の間接群） | GSE324491 の CPA（プロドラッグなので in vitro では効かない）と対比 | 早期応答（TGF-β、CXCL12、VCAM1） | 6 時間で早すぎる。メルファランは乳がんで使わない |

推奨: **A を Fig4B の補助（Supplementary）候補**にする。「生体の骨で Dox による DEG がほとんど出ない（Fig4B）のに、同じ機序の薬をヒト骨系細胞に直接当てるとアポトーシス・老化の軸がはっきり動く」という対比になるので、生体で出なかった理由（用量・時間・細胞組成の希釈）を論じる材料になる。
ただし in vitro の陽性対照として扱い、生体の結果を言い換える材料にはしない。

---

## 方法

- GEO（`esearch db=gds` → `esummary`）: 細胞語 18 語（osteoblast / osteocyte / MLO-Y4 / Ocy454 / IDG-SW3 / MC3T3-E1 / hFOB / "mesenchymal stem cells" / "mesenchymal stromal cells" / "bone marrow stromal" / BMSC / osteoclast / RAW264.7 / "bone explant" / calvaria など）× 薬剤語 16 語の総当たり（288 検索）と、
  老化・ドナー年齢・骨器官培養などの組み合わせ 23 本、ヒト化療後骨髄の組み合わせ 10 本。ヒットはタイトルではなく **summary 本文で「骨系細胞語」と「化療薬語」が共起するか**を機械判定して一次選別した。
- PubMed → GEO（`esearch db=pubmed` → `elink dbfrom=pubmed db=gds`）: 論文側の検索で GEO に寄託のあるものを拾った（骨系細胞 × 化療 × トランスクリプトーム 1,461 報 → 204 シリーズ、ヒト骨髄 × 化療 × 小児・固形がん・リンパ腫 1,444 報 → 216 シリーズ）。
- BioSample（`esearch db=biosample` → `efetch`）: 5,169 件を取得し、**属性値**に薬剤名と骨系細胞語が両方入るものを BioProject 単位に集計した。
- SRA（`esearch db=sra` → `esummary`）: `"transcriptomic"[Source]` で絞って Study 単位に集計した。
- ARCHS4: 最新版 `human_gene_v2.latest.h5`（2026-07、62 GB）と `mouse_gene_v2.latest.h5` の `meta/samples` だけを HTTP Range で読み（title / source_name / characteristics）、薬剤語を含むサンプルを GSE 単位に集計した。2020 年版の `*_gsm_meta.rda` も同様に走査した。
- recount3: `sra.recount_project.MD.gz` で全プロジェクト（ヒト 8,802 / マウス 10,116）を列挙し、プロジェクトごとの `sra.sra.<project>.MD.gz` の `sample_attributes`・`sample_title`・`experiment_title` を走査した。**マウスは 9,008 プロジェクトまでで打ち切った（未完）。** 再開するときは、作業ディレクトリ `/tmp/claude-1000/-home-ken-src-github-com-endaaman-btchem/99618536-1210-4a9b-a834-3e1075540012/scratchpad/w3/recount3/` で `mouse_projects.txt` のうち `md_mouse/<project>.gz` が無いものを `https://duffel.rail.bio/recount3/mouse/data_sources/sra/metadata/<project の末尾 2 文字>/<project>/sra.sra.<project>.MD.gz` から取り直し、`python3 -I /tmp/claude-1000/-home-ken-src-github-com-endaaman-btchem/99618536-1210-4a9b-a834-3e1075540012/scratchpad/bin/r3scan.py md_mouse` を回す（作業ディレクトリは一時領域なので、消えていたらプロジェクト一覧の取得からやり直しになる）。
- ArrayExpress（BioStudies API の `arrayexpress/search`、AND 検索）。
- 候補は FTP の `GSExxxxxx_family.soft.gz` を取得し、種・細胞・処置・n・ファイルは SOFT から取った。**SOFT に記載がない項目は「記載なし」とし、推測で埋めていない。** ドナー情報は論文本文（Europe PMC）で確認できたものだけ書いた。
- 評価: ◎ = 生体データとの突き合わせに直接使える、○ = 補助・Supplementary に使える、△ = 参照程度、✕ = 目的に合わない。

---

## 1. in vitro / ex vivo の骨系細胞 × 化療薬

| ID | 種・細胞・ドナー年齢性別 | 処置・用量・時間 | n | 処理済みファイル | PMID | 評価 |
|---|---|---|---|---|---|---|
| **GSE284666**（SenCat。GSE285826 の一部） | ヒト初代骨芽細胞（NHO）と骨髄 MSC（BMMSC）。**ドナーの年齢・性別・入手先は SOFT に記載なし**（論文は OA でなく本文未確認） | 対照（DMSO 8 日）/ **エトポシド 20 µM（NHO）・50 µM（BMMSC）8 日** / 放射線 10 Gy（NHO）・15 Gy（BMMSC）8 日 | 細胞ごとに 4 / 4 / 4 | サンプルごとの `*_geneCOUNT.txt.gz`（STAR → featureCounts、Ensembl ID） | 42276073 | ◎ ヒトの骨芽細胞と骨髄 MSC で、TOP2 阻害薬による老化の応答を counts から DE できる。同じシリーズに他の初代細胞種（ドキソルビシンで老化させたのは冠動脈内皮と血管平滑筋のみ）もあるので、骨系に特有の応答かどうかも見られる |
| **GSE135323** | ヒト骨芽細胞を足場上で 3D 培養した TE-Bone（骨芽細胞の由来・ドナーは SOFT に記載なし）。TE-Tumor（骨芽細胞＋Ewing 肉腫）もある | TE-Bone: **ドキソルビシン 4.16 µM、3 日投与 / 4 日休薬 × 2 週** vs 溶媒対照（他にリンシチニブ）。TE-Tumor にはシスプラチン・MTX・ビンクリスチン・デキサメタゾン・MAP も | TE-Bone: 対照 3 / Dox 3 | `GSE135323_TCTC_Mapped_Reads.csv.gz`（TempO-Seq S1500+、2,977 プローブ × 132 サンプル。列名は fastq 名で、SOFT の `Sample_description` と対応） | 31626302 | ◎（ただし注意つき）ヒト骨芽細胞 × ドキソルビシンの唯一のデータ。毒性学用の S1500+ パネルなので p53 標的・ストレス応答は入っているが、全遺伝子の GSEA はできない。TE-Bone の Dox は GSM4005429/31/33、対照は GSM4005423/25/27 |
| GSE17860（SubSeries GSE17839 / GSE17842 / GSE17843） | ヒト初代骨芽細胞（HOB、PromoCell から購入と論文に記載。ドナー情報なし） | GSE17839: **メルファラン 50 µg/ml 6 h** vs 溶媒（共通リファレンス）/ GSE17842: メルファラン処理した骨髄間質の培養上清 vs 対照上清（色素交換）/ GSE17843: TGF-β1 | 4 / 4、3 ペア、3 ペア | `GSE17839_filtered_LogRatio.txt.gz` ほか（Agilent 4x44K） | 22363485 | ○ アルキル化薬の直接作用と、傷んだ間質を介した間接作用を分けて見られる。論文の結論は TGF-β1 活性（Smad2 リン酸化）の上昇と、造血前駆・CXCR4+ 細胞の支持能の低下。エトポシド（VP-16）は表現型のみでアレイなし |
| GSE283633 | マウス（系統・性別・週齢は SOFT に記載なし）、骨から得た BMSC（初代培養） | 対照 / **エトポシドで老化誘導**（用量・日数は SOFT に記載なし）。同じシリーズに骨髄マクロファージ（同処置）と単球（若齢 vs 老齢） | BMSC 4 / 4 | サンプルごとの `*_featurecount.tabular.gz` | 未収載 | ○ マウス骨髄 BMSC の DNA 傷害性老化。GSE296758（雌 LepR+、生体 Dox）と並べられる |
| GSE301693 | マウス骨髄間質細胞株 ST2 | 無処置 / **シスプラチン 2 µM 24 h** / シスプラチン＋アロエエモジン | 3 / 3 / 3 | `GSE301693_gene.count.matrix.annotation.xls.gz` | 未収載 | ○ 白金製剤（TNBC のカルボプラチンに近い）による骨髄間質の応答。著者はフェロトーシスと骨芽分化の抑制を主張 |
| GSE290573 | マウス RAW264.7（破骨前駆として使われるマクロファージ株） | 培地 / **ドキソルビシン 100 nM 24 h** | 5 / 5 | `GSE290573_summary_FPKM.txt.gz` | 42155834 | △ 目的は心毒性。破骨分化（RANKL）はしていない。Dox による SASP の参照程度 |
| GSE227160 | ヒト骨髄 MSC（**女性**）・脂肪 MSC（女性）・皮膚線維芽細胞（男女）。年齢記載なし | 若い継代 / 中間継代 / 複製老化 / 放射線老化 / **ドキソルビシン老化** | BM-MSC: 各 2–4（Dox 4） | `GSE227160_Matrix_signal_intensities_GEO.csv.gz`（EPIC メチル化アレイ） | 36980268 | △ RNA ではない。女性骨髄 MSC の Dox 老化はこれだけ |
| SRP655200（SRA のみ。GEO なし） | ヒト骨髄 MSC と MSC 由来脂肪細胞、健常ドナー 4 名（Lonza、**21–47 歳、中央値 27 歳**、性別は "not applicable"） | B-ALL 株 Nalm-6 と共培養し、無処置 / **デキサメタゾン＋ダウノルビシン＋ビンクリスチン 72 h** | ドナー 4 × 2 条件 × 2 細胞 | 処理済みなし（FASTQ のみ） | 未収載 | △ 白血病細胞との共培養でステロイドを含むカクテル。ドナーが若い |
| GSE216842 | ヒト初代間質（乳腺・肺・臍帯・HUVEC）。骨髄由来ではない | ブレオマイシン | 3–6 | DESeq2 表 | — | ✕ 骨系ではない（第二波で既出） |
| GSE248107 | ヒト胸部大動脈由来 MSC | 放射線 | 4 / 3 | `GSE248107_raw_counts.csv.gz` | 39011489 | ✕ 化療薬ではなく、骨髄由来でもない |
| GSE243008 | マウス骨髄由来細胞 → 破骨分化 | タイトルに methotrexate とあるが、寄託されたのはアミノ酸欠乏・SHMT2 阻害・PHGDH KO のみ | — | counts | 未収載 | ✕ MTX 群なし |
| GSE128728 | 二重ヒットリンパ腫 PDX のマウス骨髄マクロファージ | Dox / CTX / アレムツズマブ | 5–11 | counts | 未収載 | ✕ 腫瘍を持つマウスのマクロファージ（in vivo）。間質ではない |
| GSE326160 / GSE133098 | ヒト臍帯 MSC | 複製老化・KDM 操作（薬剤なし） | 1–2 | counts | — / 31704649 | ✕ |

**骨細胞（MLO-Y4 / Ocy454 / IDG-SW3）× 化療薬は無い。** MLO-Y4 の公開データ（21 シリーズ）は低酸素・機械刺激・磁場・レーザー・転移との共培養、Ocy454 / IDG-SW3（22 シリーズ）は PTH・ビタミン D・機械刺激・分化のみ。
**MC3T3-E1 × 化療薬**も無い（α 線照射の miRNA のみ）。**hFOB × 化療薬**も無い。**破骨細胞 × 化療薬**は RAW264.7 の Dox（GSE290573）だけ。
**骨器官培養（calvaria・骨片の ex vivo）× 化療薬**も無い。

### ドナーの年齢・性別が分かるヒト骨髄 MSC（化療なし。上の in vitro 化療データと並べる側）

| ID | 内容 | n | 処理済みファイル | PMID | 評価 |
|---|---|---|---|---|---|
| GSE39540 | ヒト骨髄 MSC、61 ドナー（17–89 歳）。性別・年齢・採取部位（寛骨臼 / 腸骨稜）をサンプルごとに記載。**女性 55 歳以上 28 / 45–54 歳 7 / 45 歳未満 11**、男性 15 | 61 | RAW（U133A 2.0） | 22927939 | ○ 閉経後女性ドナーの MSC の加齢シグネチャを作れる。案 B の生体側 |
| GSE35959（SuperSeries） | ヒト骨髄 MSC。中年（42–67 歳、女性 4・男性 1）/ 高齢（79–89 歳、女性 3・男性 1）/ **高齢女性の原発性骨粗鬆症 5 例（79–94 歳）**/ 長期培養で老化（5） | 19 | RAW（U133 Plus 2.0） | 23028809 | ○ 骨粗鬆症の高齢女性ドナーがある |
| GSE39035 | ヒト骨髄 MSC、若年 4（21–24 歳）/ 高齢 4（75–82 歳）× 継代 4・8。性別記載なし | 16 | RAW（Agilent） | 23271708 | △ |

---

## 2. ヒトの化療後骨髄（AML 以外）

| ID | 内容 | n | 処理済みファイル | PMID | 評価 |
|---|---|---|---|---|---|
| GSE220946 | 小児神経芽腫の骨髄組織（手術室で採取）の scRNA。`treatment: yes/no` と `metastatic` の記載のみ。**年齢・性別・レジメン・採取時期は SOFT に記載なし** | 治療あり 11（転移 5 / 非転移 6）/ なし 4（転移 2 / 非転移 2） | サンプルごとの `*.count.csv.gz` | 38358826 | △ 化療を受けた骨髄の間質細胞が含まれていれば比較できるが、小児で、治療の内容も時期も分からない |
| GSE266330 | 骨転移 34 例（乳がん 10・腎がん 21 ほか）と健常骨髄 10 の scRNA | 63 | barcodes / features / matrix | 40412393 | △ 骨転移そのもの。治療歴は SOFT に記載なし |
| GSE162554 | 小児がん患者の化療中 / 後のペア（32 ペア）のメチル化。骨髄は 3 ペアだけで、残りは末梢血 | 64 | β 値 RData | 35354948 | ✕ |
| GSE25624 | 神経芽腫の骨髄（転移 9 / 限局 23 / 健常若年 8） | 40 | Agilent | 21994039 | ✕ 治療前の診断時で、化療の情報なし |
| GSE51907 | 造血器腫瘍患者のシクロホスファミド投与前後（6 / 24 h） | 41 | NimbleGen | 24466173 | ✕ 末梢血単核球 |
| GSE150624 | 転移性乳がん患者の CTC を移植したマウスの骨髄常在がん細胞 | 12 | RAW | 32575420 | ✕ がん細胞側 |
| GSE57570 | 化療前後の血清 miRNA | 25 | log2 | 26024523 | ✕ |

**乳がんの術前・術後化療中の骨髄（間質・全骨髄・骨生検）、リンパ腫・固形がんの自家移植前採取骨髄の間質、化療誘発性骨減少のヒト骨生検のトランスクリプトームは無い。**
ヒトで化療前後の骨髄間質を比べられるのは、引き続き AML（GSE212615 / GSE213056）と ALL（GSE318370）だけ。

---

## 3. メタデータの直接検索（SRA / BioSample / ARCHS4 / recount3）

| 経路 | 対象 | 骨系 × 化療薬で拾えたもの | うち GEO のタイトル・summary 検索で未把握だったもの |
|---|---|---|---|
| BioSample の属性 | 検索語でヒットした 5,169 件。属性値に薬剤名が入るのは 55 件 | GSE284666（NHO / BMMSC × エトポシド）、GSE135323、GSE296757/8（既知）、GSE303814（ラット MSC × Dox 老化 → ナノ製剤の評価が主目的）、GSE333489（ヒト MSC × エトポシド、遺伝子編集の毒性評価） | **GSE284666**（タイトルは「SenCat」で骨の語を含まない） |
| SRA（`"transcriptomic"[Source]`） | 骨系 × 薬剤 89 run / 骨・骨髄・間質 × 化療 1,618 run | 上記に加えて SRP655200（ヒト BM-MSC × ダウノルビシン等、Nalm-6 共培養） | **SRP655200**（GEO なし） |
| ARCHS4 最新版（2026-07） | ヒト 1,098,771 / マウス 1,225,044 サンプル。属性に化療薬名を持つシリーズはヒト 189 / マウス 15 | ヒト: GSE135323 のみ（他は U2OS・SW1353 などの肉腫株）。マウス: GSE104887 / GSE108891（5-FU の in vivo、既知）と造血系のみ | なし（GSE284666 は ARCHS4 に未収載） |
| ARCHS4 2020 年版 | ヒト 287,071 / マウス 371,797 サンプル | ヒト: 骨系なし（ES 由来の神経構築物 × 5-FU の GSE63935 のみ）。マウス: 5-FU の造血系と GSE104887 / GSE108891 | なし |
| recount3（**マウスは未完**） | ヒト: 全 8,802 プロジェクト（SRA 由来、2019 年頃まで）を走査済み。マウス: 10,116 中 9,008 プロジェクトで取得を打ち切り（残り約 1,100 は未走査） | ヒト: 骨系なし（ES 由来の神経構築物 × 5-FU の SRP050892 と ENCODE のみ）。マウス（取得済み分）: SRP093722（5-FU 後の GMP、造血系）と SRP119850（= GSE104887、既知）のみ | なし |
| ArrayExpress | `osteoblast AND doxorubicin` など 18 通り | E-MTAB 独自のものは無い（GEO の再掲と E-MTAB-3732 などの大規模集成のみ） | なし |

---

## 検索語と、見つからなかった組み合わせ

### 検索語

- GEO（`db=gds`、`gse[ETYP]`）
  - 細胞: osteoblast(s) / osteocyte / MLO-Y4 / MLOY4 / Ocy454 / IDG-SW3 / MC3T3-E1 / hFOB / "human fetal osteoblast" / "primary human osteoblasts" / "human osteoblasts" / "mesenchymal stem cells" / "mesenchymal stromal cells" / "bone marrow stromal" / "bone marrow mesenchymal" / BMSC / osteoclast(s) / osteoclastogenesis / RAW264.7 / "bone explant" / "organ culture" / "ex vivo" / calvaria / calvarial
  - 薬剤: doxorubicin / epirubicin / daunorubicin / idarubicin / mitoxantrone / anthracycline / cyclophosphamide / 4-hydroperoxycyclophosphamide / mafosfamide / ifosfamide / "4-HC" / fluorouracil / 5-FU / paclitaxel / taxol / docetaxel / vincristine / carboplatin / cisplatin / oxaliplatin / etoposide / methotrexate / cytarabine / chemotherapy / "DNA damage" / genotoxic / "therapy-induced senescence" / senolytic
  - ドナー: "donor age" / "elderly donors" / "aged donors" / "old donors" / postmenopausal
  - ヒト骨髄: "bone marrow" × (neoadjuvant / "breast cancer patients" / lymphoma / Hodgkin / neuroblastoma / Ewing / rhabdomyosarcoma / "solid tumor" / "stem cell mobilization" / apheresis / "autologous transplant" / "before and after" / paired / R-CHOP / ABVD / "bone loss" / osteoporosis) × (stroma / mesenchymal / aspirate / biopsy / trephine / single-cell)
- PubMed → GEO（elink）: 上の薬剤 × 細胞 × (transcriptom* / RNA-seq / microarray / "gene expression profil*" / single-cell)、("bone marrow" / "bone biopsy" / trephine / "marrow stromal") × chemotherapy × (lymphoma / neuroblastoma / "solid tumor" / "breast cancer" / sarcoma / pediatric / autologous / "bone loss")、骨系細胞 × 薬剤 × (aged / elderly / "donor age" / postmenopausal / senescence)
- BioSample: 薬剤 13 語 × (osteoblast / osteocyte / MLO-Y4 / Ocy454 / IDG-SW3 / MC3T3-E1 / hFOB / osteoclast)、薬剤 × ("mesenchymal stem" / "mesenchymal stromal" / BMSC / "bone marrow stromal" / "stromal cells")
- SRA: 同上の薬剤 × 細胞、化療 × ("bone marrow" / bone / mesenchymal / stromal / osteo*)、いずれも `"transcriptomic"[Source]`
- ARCHS4 / recount3: サンプル属性・タイトルに対する正規表現（薬剤: doxorubicin / adriamycin / doxo / epirubicin / cyclophosphamide / mafosfamide / 4-HC / fluorouracil / 5-FU / paclitaxel / taxol / docetaxel / carboplatin / cisplatin / CDDP / etoposide、細胞: 上の細胞語＋ bone / marrow / stroma / skeletal / NHO / ST2）
- ArrayExpress: osteoblast / osteocyte / mesenchymal / osteoclast / stromal × 各薬剤 / chemotherapy、"bone AND explant AND doxorubicin"、"marrow AND chemotherapy AND lymphoma"

### 見つからなかった組み合わせ

- **骨細胞（株・初代・骨細胞濃縮骨片）× 化療薬**のトランスクリプトーム（GEO・BioSample・SRA・ARCHS4・recount3（マウスは 9 割まで）・ArrayExpress のどれにも無い）。
- **骨芽細胞株（MC3T3-E1 / hFOB）× 化療薬**。ヒト骨芽細胞は初代（NHO・HOB）と 3D の TE-Bone のみ。
- **シクロホスファミドの活性体（4-HC / マフォスファミド）× 骨系細胞**（in vitro で CTX を効かせるには活性体が要るが、使ったデータは無い。GSE324491 は CPA そのもの）。
- **タキサン（パクリタキセル / ドセタキセル）× 骨系細胞・骨髄 MSC**。
- **5-FU × 骨系細胞の in vitro**（in vivo の骨髄 niche は GSE104887 / GSE108891、ヒト ES 由来骨芽系は GSE324491 のみ）。
- **破骨細胞（RANKL で分化させたもの）× 化療薬**。
- **骨器官培養（calvaria / 骨片）× 化療薬**。
- **年齢・性別の記載があるヒト骨芽細胞・骨髄 MSC に化療薬を当てた RNA データ**（年齢つきは化療なしの加齢データだけ。化療ありは GSE227160 のメチル化で女性のみ記載）。
- **ヒトの化療後骨髄（AML 以外）で、間質を含み、年齢・性別があるもの**。乳がん患者の化療中の骨髄、リンパ腫・固形がんの自家移植前骨髄の間質、化療誘発性骨減少の骨生検は、いずれも無い。

## ダウンロードしたもの

- `data/seq/GSE284666/`: `GSE284666_family.soft.gz`、`raw/` に NHO と BMMSC の 24 サンプル分の `*_geneCOUNT.txt.gz`（対照・エトポシド・放射線 × 各 4。他の細胞種と `GSE284666_RAW.tar`（42 MB）は取得していない）。
- `data/seq/GSE135323/`: `GSE135323_family.soft.gz`、`GSE135323_TCTC_Mapped_Reads.csv.gz`。
- ○ 以下は SOFT のみ確認（作業ディレクトリ側。`data/` には置いていない）。

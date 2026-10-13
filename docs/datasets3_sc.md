# 追加公開データ調査 第三波（単一細胞・単一核アトラスと細胞アトラスのポータル）

調査日: 2026-10-13。第一波（`docs/datasets.md`）と第二波（`docs/datasets2.md`）で骨・骨髄 ×（加齢 / OVX / 化療）のバルク RNA-seq はほぼ調べ尽くしたので、
今回は単一細胞・単一核のデータとアトラス系のポータルから探した。目的は「閉経後（加齢）の骨で化療が骨転移 niche を作る」を dry で補強すること。
第一波・第二波で評価済みのもの（GSE108891 / GSE123078 / GSE289491 / GSE296752 / GSE296758 / GSE169162 / GSE316278 / GSE299773 / GSE315430 / GSE302366 / GSE305558 / GSE326234 / GSE343407 / GSE279904 / GSE325591 ほか）は再評価していない。

## 結論

- **(3) 加齢 × 化療、OVX × 化療を同じ実験で組んだ単一細胞データは無い。** CELLxGENE・Single Cell Portal・HCA・GEO のどこにも無かった（第二波の bulk の結論と同じ）。
  いちばん近い設計は GSE276784 で、骨髄間質（MSC 系）の scRNA を「老齢（24 か月）」と「若齢 + 炎症刺激 poly(I:C)」で並べている。ただし刺激は化療ではなく、老齢と刺激の組み合わせ群も無い。
- **(1) 雌の加齢骨髄間質の単一細胞データは 2 本見つかった。これが今回の主な収穫。**
  - **GSE210584**（◎候補）: **雌** C57BL/6J、若齢成体 5 匹と**中年** 9 匹を **1 匹ずつ別ライブラリ**で 10x scRNA にかけたもの。骨髄の造血細胞と niche 細胞（MSC・内皮）を分取している。
    論文（Young et al., *Blood* 2024, PMID 38598841）の結論は「中年で MSC の Kitl / Igf1 / アディポネクチンが下がり、HSC の加齢を駆動する。IL-1β / IL-6 / TNF などの炎症性サイトカインとは独立」。
    個体が反復になる（n = 5 vs 9）ので、擬似バルクにすれば雌の中年 MSC の DEG がちゃんと検定できる。**中年はマウスの生殖老化の時期に重なるので、「閉経前後」に最も近い単一細胞データ**。
    中年の月齢は GEO に記載がない（論文本文は非 OA で未確認。要確認）。
  - **GSE275179**（○）: **全サンプル雌**。大腿骨と頭蓋骨の骨髄で、若齢（10–14 週）と「geriatric」（95 週以上）の間質（造血系除去）・HSPC・全細胞を BD Rhapsody で測ったもの（Koh et al., *Nature* 2024, PMID 39537918）。
    群ごとに 1–3 ライブラリで、大腿骨の間質は若齢 1 / 老齢 1。雌の老齢骨髄間質の参照には使えるが、統計検定には足りない。
- **(1) 雌雄両方で性差を見られるデータ**: GSE267695（若齢の雌 1 と雄 1、骨髄間質 scRNA。雄の方が LepR+ MSC が多く、Kdm5c–Cxcl12 に性差がある。PMID 39836478）。加齢の軸は無い。
  Tabula Muris Senis の骨髄（CELLxGENE）は雌（3 / 18 / 21 か月）を含むが**間質が 1 細胞も入っていない**（洗い流した骨髄の造血細胞だけ）。
  **Tabula Muris Senis の bulk（GSE132040）には「Bone」（大腿骨と脛骨）があり、雌が 1–21 か月の 8 時点 × 各 n = 2**（雄は 1–27 か月 × n = 4）。単一細胞ではないが、雌の骨で 12–15 か月の閉経相当期をまたぐ時系列が取れるのはこれだけ。
- **ヒトの雌の加齢骨髄間質**: CELLxGENE の「Bone Marrow Reference」（Chu et al., *Nat Commun* 2026）が複数の公開データ（GSE190965 / GSE169396 / GSE196678 ほか）を統合しており、
  女性ドナーの間質細胞が閉経前（19 / 21 / 45 歳）と閉経後年齢（53 / 58 / 61 / 72 / 74 / 77 歳）の両方にある（下表）。元データごとに単離法が違うので、元データ内で比べる（GSE190965 の中なら女性 19・21 歳 vs 53・58 歳）のが安全。
- **(2) 化療後の骨髄 niche の単一細胞データで、雌・年齢が分かる新しいもの**は見つからなかった。既知の GSE289491（雌の CAR 細胞、Dox 96 時間、n = 1）と GSE108891 / GSE123078（5-FU、性別記載なし）のまま。
  ヒトでは、乳がん骨転移を含む骨転移の scRNA（GSE266330 の乳がん 12 例、HTAPP の骨生検 5 例）に前治療歴があるが、化療の前後比較ではなく、治療歴は論文の補足表にしか無い。

### 使い道の提案

| 案 | 中身 | 何を補強するか | 弱点 |
|---|---|---|---|
| **A. 雌の中年 MSC シグネチャ** | GSE210584 の MSC を個体ごとに擬似バルク化 → 中年 vs 若齢の DEG（n = 9 vs 5） | 「閉経前後の雌の骨髄間質で何が変わるか」を雌・個体反復つきで示す。第二波の案 B（GSE296758 の雌 LepR+ MSPC × Dox）と同じ細胞種で、DEG / Hallmark の重なりを見る | 化療は別データ。MSC の細胞数は造血系より少ない |
| B. 雌の骨の時系列（bulk） | GSE132040 の Bone、雌 1–21 か月 × n = 2 | Fig4A（OVX）の遺伝子が雌の自然加齢でも 12 か月以降に動くかを時系列で示す | bulk。n = 2。骨髄を含むかどうかは論文に記載なし（「femurs & tibiae」とだけある） |
| C. 「老齢の間質 ≒ 傷害を受けた間質」 | GSE276784 の老齢 iMSC-L と poly(I:C) 後 iMSC-L（Smart-seq、counts 行列あり） | 加齢と急性傷害で同じ炎症性間質（Sca-1 高、IFN 応答）ができる、という論法の先行例。化療で同じ状態ができるかを GSE289491 / GSE296758 と突き合わせる | 刺激は化療ではない。性別は雌雄混合・プール |
| D. ヒト女性の閉経前後の骨髄間質 | CELLxGENE Bone Marrow Reference（GSE190965 の女性 2 vs 2 を中心に） | ヒトでも閉経後年齢で間質の炎症・脂肪化が上がるかの参照 | n が極小。元データ間のバッチ |

案 A を「4A（閉経・加齢）の雌の単一細胞版」として最優先で試す価値がある。4B（化療）の置き換えになるデータは今回も見つからなかった。

---

## 方法

- **CELLxGENE Discover**: Curation API `https://api.cellxgene.cziscience.com/curation/v1/datasets`（全 2,240 データセット）と `/collections` を取得し、
  tissue のラベルが bone / marrow / skelet / femur / osteo / vertebra / calvar / iliac に当たる 66 データセットを列挙した。そのうち間質系の細胞種（stromal / osteo / mesench / fibro / chondro / adipo / endotheli / pericyte / smooth / reticular）を含むものについて、
  公開 h5ad の `obs` を **HTTP の範囲読み（h5py + fsspec）でダウンロードせずに読み**、donor・sex・development_stage・disease・cell_type を集計した。
- **Single Cell Portal（Broad）**: `/single_cell/api/v1/search?type=study&terms=…`（bone marrow stroma / bone marrow niche / osteocyte / osteoblast / bone marrow aging / aged bone marrow / ovariectomy / chemotherapy bone marrow / 5-FU / doxorubicin / skeletal stem / bone / LepR）と、候補の `/site/studies/<SCP>` の記載。
- **HCA Data Portal（Azul）**: `/index/projects?filters={"organ":{"is":[…]}}`（bone marrow / bone tissue / bone element / skeletal element / femur / rib / skull / hindlimb / musculoskeletal system）。
- **GEO**: E-utilities（`esearch db=gds` → `esummary`）。単一細胞語（"single-cell" / "single cell" / scRNA-seq / snRNA-seq / "single-nucleus" / "single nucleus"）× 組織語 17 × 条件語 19 の総当たり（323 検索、220 件）と、
  追加の組み合わせ（骨髄 × 化療薬 15、骨 × OVX・エストロゲン 7、骨髄間質 × 性差 7、骨組織 × 加齢 8、ヒト化療後の骨髄間質、骨髄 × 加齢 × 化療、骨 × 加齢 × 化療、OVX × 化療、乳がん × 骨髄 / 骨転移 × 単一細胞）。
  候補は FTP の `GSExxxxxx_family.soft.gz` を取って、Series の summary / overall_design / supplementary と全 Sample の characteristics を集計した。
- 論文: PubMed の抄録と、Europe PMC / NCBI PMC の全文（OA のもの）で、性別・月齢・処置を確認した。
- **表の種・性別・年齢・処置・n は SOFT・API のメタデータ・h5ad の obs、または論文本文から取った。記載がない項目は「記載なし」とし、推測では埋めていない。**
- 評価: ◎ = Fig4/5 に直接使える、○ = 補助パネル・Supplementary・検証に使える、△ = 参照程度、✕ = 目的に合わない。

---

## 1. 雌の加齢（または OVX）骨髄間質・骨芽細胞・骨細胞（目的 (1)）

| ID | 種・性別・年齢 | 組織・細胞 | 群・処置 | n | 所在・形式 | PMID | 評価 |
|---|---|---|---|---|---|---|---|
| **GSE210584** | マウス C57BL/6J **雌**（性別は論文抄録。SOFT には記載なし）。若齢成体 / 中年（SOFT の summary に 12–14 か月と記載） | 脊椎・胸骨・上腕骨・脛骨・大腿骨・腸骨稜を砕いた骨髄。造血細胞と niche 細胞（酵素消化後に分取。MSC・内皮を含む） | 若齢 / 中年（無処置） | **個体ごとに 1 ライブラリ: 若齢 5 / 中年 9** | GEO、サンプルごとの 10x `barcodes / features / matrix.mtx`（CellRanger 6 の filtered） | 38598841 | **◎候補**。雌・個体反復つき・中年（閉経相当期）。MSC の擬似バルクで雌の中年 MSC の DEG が取れる。中年は 12–14 か月（閉経相当期） |
| **GSE275179** | マウス C57BL/6 **雌**（全 19 サンプル）。若齢 10–14 週 / geriatric 95 週以上（論文 Methods。scRNA は young と geriatric） | 大腿骨と頭蓋骨の骨髄。造血系を MACS で除いた間質、KIT+ HSPC、全細胞（BD Rhapsody） | 若齢 / 老齢 × 大腿骨 / 頭蓋骨 | ライブラリ単位で 1–3（大腿骨の間質は若齢 1 / 老齢 1、大腿骨の全細胞は若齢 1 / 老齢 2） | GEO、サンプルごとの `barcodes / features / matrix.mtx`、`GSE275179_BD_CLS{1-3}(_v2).txt.gz` | 39537918 | ○ 雌の老齢骨髄間質の参照。統計には n が足りない |
| **GSE132040**（Tabula Muris Senis bulk） | マウス C57BL/6JN、**雌 1 / 3 / 6 / 9 / 12 / 15 / 18 / 21 か月 × n = 2**、雄 1–27 か月 × n = 3–4 | 「Bone」（論文では femurs & tibiae）と「Marrow」ほか計 17 臓器 | 無処置の加齢時系列 | Bone: 雌 16 / 雄 39。Marrow: 雌 15 / 雄 39 | GEO、全臓器の counts（`GSE132040_190214_…csv.gz`、40 MB）と `GSE132040_MACA_Bulk_metadata.csv.gz` | 32669715 | ○ **bulk** だが、雌の骨で閉経相当期（12–15 か月）をまたぐ時系列はこれだけ。Fig4A の OVX 遺伝子の動きを雌の自然加齢で確認できる |
| Tabula Muris Senis 骨髄（CELLxGENE。10x と Smart-seq2） | マウス、雌 3 / 18 / 21 か月（10x は各 2 匹）と雄 1–30 か月 | 洗い流した骨髄 | 加齢 | 10x 40,220 細胞（雌 11,707）、Smart-seq2 14,517 細胞 | CELLxGENE h5ad（0.57 / 0.34 GB）。figshare / AWS にもある | 32669715 | ✕（目的 (1) について）**間質・骨芽細胞が 0 細胞**。細胞種は顆粒球・赤芽球・B 細胞・マクロファージなど造血系だけ。雌の加齢骨髄マクロファージの参照なら △ |
| GSE267695 | マウス C57BL/6J、若齢の雌 1 / 雄 1 | 骨髄間質（CD45- Lin- Ter119-）scRNA | 性別のみ | 1 / 1 | GEO RAW | 39836478 | △ 若齢での間質の性差（雄で LepR+ MSC が多い、Kdm5c–Cxcl12）。加齢の軸がない |
| GSE276784 | マウス C57BL/6、8–12 週と 24 か月。性別は「雌雄プール」50 / 雄 37 / 雌 2（SOFT） | 中心骨髄と骨内膜の間質（CD45- Ter119-）の 10x scRNA、炎症性 MSC 様細胞（iMSC-L）の Smart-seq | 無処置 / poly(I:C) 10 mg/kg 48 時間ごと × 3 日。老齢は無処置のみ | 10x はライブラリ単位 1–3、Smart-seq は細胞単位（老齢 iMSC-L は 1 匹由来） | GEO、`GSE276784_SMARTseq_counts_matrix.csv.gz`、`_cell_groups.csv.gz`、RAW | 41060335 | ○ 「老齢の間質」と「急性炎症後の間質」が同じ iMSC-L 状態（Sca-1 高、IFN 応答）になる。案 C の根拠。化療ではない |
| GSE129726（SubSeries of GSE130299） | マウス、10–12 週 / 20 か月超（性別記載なし） | 骨髄内皮（CD45- CD31+）10x | 若齢 / 老齢 | 3 / 2 ライブラリ | GEO RAW | 31685996 | △ 加齢の骨髄内皮 |
| GSE138689 | マウス C57BL/6、8–10 週と 12 か月（性別記載なし） | Lepr-Cre 標識の間質（骨髄と骨片）10x | 成体 / 加齢（12 か月）/ ロシグリタゾン / 放射線 5 Gy / 骨折 | 各 1–2 | GEO、`GSE138689_expression.txt.gz`、RAW。SCP1939 に同じ系の行列 | 34957577 | △ LepR+ の 12 か月が 1 ライブラリ。ストレス群は放射線で化療ではない |
| GSE161946 / GSE166441 | マウス**雄**、2 か月 / 24 か月（GSE166441 は 2 / 12 / 24 か月） | 骨格幹細胞（SSC）の Smart-seq2 670 細胞 / 骨折仮骨の SSC bulk | 加齢 | — / 3 ずつ | csv / TPM xlsx | 34381212 | △ 雄のみ |
| GSE145477（SCP1017） | マウス**雄** Col2/Td、1 / 1.5 / 3 / 16 か月 | 骨内膜側の Td+ 間葉系（MALP を定義したデータ） | 加齢 | 各 1 ライブラリ（2–3 匹プール） | GEO RAW、SCP の行列 | 32286228 | △ 雄のみ |
| GSE237307 | マウス**雄**、24 か月（INK-ATTAC） | 骨と骨髄の間葉系（CD45- Lin-） | 老齢のみ | 2 レーン | GEO RAW | 37524694 | ✕ 若齢群なし |
| GSE240291（SubSeries of GSE240292） | マウス（性別記載なし）、若齢 / 中年 | 骨と骨髄（造血・内皮と間質を 1:1 で混ぜたもの） | 若齢 WT / 中年 WT / 中年 LeprCre;Ncstn cKO | 各 1（5 匹プール） | GEO RAW | 未収載 | △ 中年の SSPC で Notch 上昇 |
| GSE285020 | マウス**雄**、2 / 12 / 20 か月 | 皮質骨（骨芽〜骨細胞）bulk、8 か月の scRNA 1 本 | 加齢、運動、Trpv4 作動薬 | 3 ずつ | GEO RAW | 未収載 | △ 雄、bulk |
| GSE197594 | マウス**雌** 12 週 | CD63+ MSC、CD11b+ 骨髄球、LSK の scRNA | WT / Il1rn-KO | 各 1 | GEO RAW | 36596811 | ✕ 加齢なし（IL-1 過剰モデル） |

### ヒト

| ID | 性別・年齢 | 組織・細胞 | 群 | n | 所在・形式 | PMID | 評価 |
|---|---|---|---|---|---|---|---|
| **CELLxGENE「Bone Marrow Reference」**（コレクション 0391c84c、Chu et al. 2026） | 女性: 19 / 21 / 45 / 53 / 58 / 61 / 72 / 74 / 77 歳。男性: 22–66 歳。不明 13 | 骨髄の統合アトラス（117,308 細胞、間質 32,576 細胞）。元データは GSE190965 / GSE169396 / GSE196678 / GSE156640 / GSE147390 / GSE202813 / GSE224152 | 健常・変形性関節症・リンパ腫 | 女性の間質細胞数（obs を集計）: 19 歳 34 / 21 歳 5,067 / 45 歳 733 / 53 歳 521 / 58 歳 1,151 / 61 歳 415 / 72 歳 589 / 74 歳 242 / 77 歳 28 | CELLxGENE h5ad 0.65 GB（obs に `gse` 列あり） | doi 10.1038/s41467-026-69863-6 | ○ ヒト女性の閉経前後の骨髄間質。元データ間で単離法が違う（CD271+ 分取、大腿骨頭の酵素消化など）ので、元データの中で比べる |
| GSE190965（上の主な元データ） | 19 / 21 / 22 / 25 / 32 歳（若年）と 52 / 53 / 58 / 61 歳（高齢）。性別は CELLxGENE の obs で女性 19・21・53・58 歳 | 骨髄の CD45- CD235a-（4）と CD45- CD235a- CD271+（5）scRNA | 若年 / 高齢 | 9 | GEO、`spliced / unspliced.mtx`（生リードはスウェーデンの法律で非公開） | 36876630 | ○ 女性 2 vs 2 で閉経前後の比較ができる。n が極小 |
| GSE308956（CELLxGENE 2317da21） | 52–84 歳、男女 | 大腿骨頭の骨髄 niche（間質・内皮・線維芽細胞を含む）scRNA | CHIP 6 / 非 CHIP 10 | 16 | CELLxGENE h5ad 3 GB。GEO は処理済みのみ | 未収載（doi 10.1038/s41590-026-02668-3） | △ 全員高齢で若年対照がない |
| GSE253355 | 記載なし（SOFT） | 人工股関節手術時の大腿骨頭骨髄。間葉系濃縮・HSPC・全細胞 | 健常 | 12 | `GSE253355_MSC_Subset_Seurat.rds.gz` ほか | 38714197 | △ ヒト骨髄 niche の地図。年齢・性別は論文補足を要確認 |
| GSE169396（HCA 68df3629） | 女性 45・61 歳、男性 31・66 歳 | 大腿骨頭 scRNA | 正常 / 骨粗鬆症 / 変形性関節症 | 4 | GEO RAW | 34111027 | △ 上の統合アトラスに含まれる |
| CELLxGENE 93eebe82（Triana 2021） | 若年 25・26・34 歳 / 高齢 59・65・69 歳、男女 | 骨髄（造血中心） | 若年 / 高齢 | 3 / 3 | CELLxGENE h5ad | doi 10.1038/s41590-021-01059-0 | ✕ 間葉系細胞が計 18 個しかない |

---

## 2. 化療後の骨髄 niche の単一細胞（目的 (2)）

| ID | 種・性別・年齢 | 組織・細胞 | 群・処置 | n | 所在・形式 | PMID | 評価 |
|---|---|---|---|---|---|---|---|
| （既出）GSE289491 | マウス雌（第二波） | CAR 細胞 scRNA | Dox 96 時間 | 1 | — | — | （再掲）新しい化療 sc は見つからなかった |
| （既出）GSE108891 / GSE123078（SCP1248 に統合） | マウス、性別記載なし、7–9 週 | Lepr+ / VE-cad+ / Col2.3+ niche | PBS / 5-FU | 1 | — | 30971824 | （再掲）SCP1248（Dolgalev & Tikhonova 2021 の niche 統合データ）は Tikhonova・Baryawno・Baccin ほかを 1 つにまとめた正規化行列で、5-FU 群もこの中にある |
| **GSE266330** | ヒト、**乳がん骨転移 12 例**（BC_01–12）ほか 6 がん種の骨転移と、健常骨髄 10（5 例 × 2） | 骨転移組織の scRNA（cell hashing で 3–5 例を多重化） | がん種別。前治療歴は 42 例分が論文 Table S1 にある（本文の記載） | 63 ライブラリ | GEO RAW | 40412393（*Cell Genom* 2025、OA） | △ 前治療ありの乳がん骨転移 niche。化療の前後比較ではない。著者は「治療と免疫アーキタイプに相関なし」としている |
| HTAPP 乳がん転移アトラス（CELLxGENE a96133de、Klughammer et al. *Nat Med* 2024） | ヒト女性、骨生検は 5 例（36 / 38 / 43 / 57 / 60 歳、すべて ER+） | 転移生検の sc（新鮮）/ sn（凍結） | 前治療歴は obs に無い（論文補足） | 骨 5 例。間質細胞は 1–2,768 個と少ない例が多い | CELLxGENE h5ad、生データは dbGaP phs002371 | 10.1038/s41591-024-03215-z | △ 閉経前後の年齢が混ざるが n が小さく、間質も少ない |
| GSE202813 | ヒト、腎がん骨転移と良性の骨髄（間質濃縮） | 骨髄 scRNA | 前治療歴は 1 例だけ「パゾパニブ・ソラフェニブ・Ra-223・アキシチニブ」、ほかはなし | 33 | GEO RAW | 38281962 | ✕ 化療ではない |

---

## 3. 加齢 × 化療、OVX × 化療の単一細胞（目的 (3)、最優先）

**該当なし。** 近いものを並べる。

| ID | 内容 | 足りないもの |
|---|---|---|
| GSE276784 | 老齢間質（24 か月）と、若齢 + poly(I:C) の間質。老齢と若齢刺激後で同じ炎症性 MSC 様細胞ができる | 刺激が化療ではない。老齢 × 刺激の群がない |
| GSE302366（第二波で評価済み） | LepR+ 間質の scRNA、老齢 × 骨髄除去による損傷 | 損傷が手術。各 1 ライブラリ |
| GSE138689 | LepR+ 間質、12 か月 / 放射線 | 放射線は化療ではなく、加齢と組み合わせていない |
| GSE316278（第二波で評価済み） | 若齢の雄 / 雌、老齢 vehicle / 老齢 senolytic（753b）の間質 scRNA | 化療なし。老齢側は雌雄プール |
| GSE331488 / GSE331489（SuperSeries GSE333859） | 21 組織の sci-ATAC（EdU 標識）、1 / 5 / 21 か月、雌雄（多くは雌雄プール） | 単一細胞クロマチン。化療なし。骨髄の間質がどれだけ入っているかは未確認 |

---

## 当たった先と検索条件

| 当たった先 | 条件 | 結果 |
|---|---|---|
| CELLxGENE Discover（Curation API） | tissue に bone / marrow / skelet / femur / osteo / vertebra / calvar / iliac を含む 66 データセット。間質系の細胞種を含むものは obs を集計 | マウスの骨・骨髄は Tabula Muris Senis（造血のみ）と CNS 境界の免疫アトラスだけ。間質を含むのはヒトの統合アトラス（Chu 2026）、大腿骨頭 CHIP、HTAPP。化療の条件を持つものはない |
| Single Cell Portal | 上記 13 語 | SCP361 / 362（Baryawno、8–10 週、AML）、SCP1017（Zhong、雄）、SCP1248（niche 統合）、SCP1939（LepR+）、SCP978（ImmGen、7 週の雌雄各 1）。ovariectomy / osteocyte は 0 件、doxorubicin は腎糸球体のみ |
| HCA Data Portal（Azul） | organ = bone marrow / bone tissue / bone element / skeletal element / femur / rib / skull / hindlimb / musculoskeletal system | 骨髄 9 件は造血・免疫中心。骨組織 6 件のうち大腿骨頭（4 例）以外は胎児・腫瘍・Tabula Muris（3 か月）。加齢 × 骨間質の project はない |
| GEO（E-utilities） | 単一細胞語 × 組織語 × 条件語（方法を参照）。計 400 弱の検索 | 本書の表に載せたもの。化療 × 骨髄 niche の sc は既出の GSE108891 / GSE123078 / GSE289491 / GSE296752 から増えなかった |
| PubMed / Europe PMC / NCBI PMC | 候補論文の性別・月齢・処置 | GSE210584 の雌は抄録で確認。GSE275179 の月齢と雌は本文 Methods。GSE276784 の老齢 iMSC-L は雄 1 匹、pIC は雌雄プール（本文の図の説明）。TMS bulk の Bone は femurs & tibiae（本文） |

### 見つからなかった組み合わせ

- **加齢 × 化療、OVX × 化療を同じ実験で組んだ骨・骨髄の単一細胞データ**（マウス・ラット・ヒト、sc / sn とも）。
- **OVX の雌の骨髄間質・骨芽細胞・骨細胞の単一細胞データ**（OVX の sc は第一波のラット GSE343407 と、造血細胞の GSE272266 だけ。GSE249471 はヒト BMSC の培養）。
- **化療後の骨髄 niche の単一細胞データで、12 か月以上のマウスを使ったもの**、および**雌雄両方を化療したもの**。
- **雌の骨細胞（皮質骨）の加齢 snRNA**（GSE325591 は第一波で評価済みで、加齢の比較はない）。
- **ヒト固形がん（乳がんを含む）患者の骨髄間質の化療前後の単一細胞データ**。ヒトの化療前後ペアは白血病の bulk（GSE212615 / GSE318370）だけで、sc は寛解導入後のリンパ球（GSE279904）だけ。
- アトラス系（Tabula Muris Senis、Tabula Sapiens、HCA、Calico の加齢アトラス）には化療の条件が無い。Tabula Sapiens の骨髄は 7 例で、間質は内皮と線維芽細胞のみ。

## 保存したもの（`data/seq/`）

- `GSE132040/`: SOFT と `GSE132040_MACA_Bulk_metadata.csv.gz`（臓器・月齢・性別の表）。counts（40 MB）は取得していない。
- `GSE210584/`、`GSE275179/`、`GSE276784/`、`GSE190965/`、`GSE266330/`、`GSE267695/`: SOFT のみ。
- h5ad はダウンロードせず、`obs` だけを範囲読みで集計した（読み取りスクリプトは scratchpad に置いた。repo には入れていない）。

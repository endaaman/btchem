---
title: "探して無かったもの — 加齢・閉経 × 化療の骨データは存在しない"
date: 2026-10-13
author: 遠田
---

10/13 の打合せで、「閉経後（加齢）の動物に化療（ドキソルビシンなど）を当てて、骨や骨髄で変化が出ている公開データ」を骨髄まで広げて探し直すことになりました。
3 回に分けて、考えられる限りの場所を探しました。この章では、何を探し、どこを当たり、何が無かったかを、根拠（検索語と件数）つきで報告します。
「無い」ことは論文の Limitation に書くことになるので、その英文案も末尾に付けました。

# 結論

- **加齢または卵巣摘出（OVX）と、細胞傷害性の化療を、同じ実験の中で組み合わせた骨・骨髄のトランスクリプトーム（RNA のデータ）は、公開されているものが 1 本もありませんでした。** bulk・単一細胞・単一核のいずれも、マウス・ラット・ヒトのいずれもです。
- RNA 以外も含めると、**「年齢 × 化療 × 骨」を同じ実験で組んだオミクスはプロテオームの 1 本（PXD078161）だけ**でした。これは解析しましたが、測定のばらつきが大きく、年齢による違いを読むことができませんでした（下で詳しく説明します）。
- **OVX × 化療**は、骨量（µCT）と qPCR の論文が 1 本（Yao 2020）あるだけで、オミクスは RNA・タンパク質・エピゲノム・代謝のどれも見つかりませんでした。
- **ヒトでは、化療を受けた乳がん患者の骨・骨髄間質を化療の前後で測ったデータは無く**、化療の前後をペアで比べた骨髄データは白血病（AML・ALL）のものだけでした。
- 化療側の骨・骨髄データは、**すべて 7–16 週の若いマウス**で行われていました。12 か月以上のマウスに化療を当てた骨髄 niche のデータはありません。

これは本研究の**新規性の裏付け**（誰もこの組み合わせを調べていない）であると同時に、**dry 解析で「2 つが重なると何が起きるか」を直接示すことはできない**という限界でもあります。
前章の比較がすべて「別々の実験の向きを並べたもの」になっているのはこのためです。

# 探したもの（理想のデータ）

理想は次のようなデータです。

1. **動物の 2 × 2**: 若齢 / 老齢（または sham / OVX）× 対照 / 化療 の 4 群を同じ実験で作り、骨（骨細胞を含む皮質骨）や骨髄の間質を RNA-seq で測ったもの。雌で、乳がんの周術期レジメン（アントラサイクリン、シクロホスファミド、タキサン）に近い薬が望ましい。
2. **ヒト**: 乳がん患者の骨・骨髄を、化療の前後で（できれば閉経前と閉経後に分けて）測ったもの。
3. 1・2 が無い場合の次善: 年齢・性別の分かる化療後の骨・骨髄データ、化療を受けた骨転移モデルの骨側のデータ、骨細胞アポトーシスを年齢 × 化療で定量した論文の数値など。

# 探索 1（10/12）— GEO を中心にした第一波

**当たった先**: GEO（米国 NCBI の遺伝子発現データベース。公開 RNA データの大半がここに置かれます）、ArrayExpress（欧州の同種のデータベース）、NCBI BioProject（GEO を通さず生データだけ置いたもの）、MSigDB（遺伝子セット）。
候補ごとに GEO のメタデータ（`family.soft`）を取得し、種・性別・組織・処置・n を全サンプル分読みました。

**結果**: 閉経・OVX の骨、化療の骨、ステロイドの骨、乳がん骨転移モデルの骨を合わせて 55 件ほどを評価しました（前章の補助・不採用の多くがここで見つかったものです）。次のものは見つかりませんでした。

- 乳がん骨転移モデル × 化療の、骨側のトランスクリプトーム（近いのは GSE174638 で、骨転移モデルにドキソルビシン + ゾレドロン酸を与えていますが、測ったのは原発腫瘍の免疫パネルだけ）
- 乳がん骨転移モデル × OVX の、骨側のトランスクリプトーム
- タキサン（パクリタキセル・ドセタキセル）× 骨・骨細胞・骨髄間質
- シクロホスファミド × 骨組織（骨髄全体のデータだけ）
- AC / TC レジメン（併用）× 骨
- 化療を受けた乳がん患者の骨・骨髄間質（あるのは骨髄の播種腫瘍細胞そのものか、原発腫瘍のデータ）

# 探索 2（10/13 午前）— 加齢・OVX × 化療の総当たり

打合せを受けて、「加齢 / OVX と化療の組み合わせ」に絞って総当たりしました。

**GEO**: 年齢の語 12 語 × 薬剤の語 15 語の総当たり（組織の語はまとめて OR で掛け、計 180 検索）、薬剤の語 12 語 × niche の語 19 語、年齢・OVX の語 8 語 × niche の語 14 語で集め、タイトルと要約で選別したあと、候補のメタデータを全部読みました。

- 年齢・閉経の語: aged / "old mice" / aging / "young and old" / ovariectomy / ovariectomized / "estrogen deficiency" / "middle-aged" / geriatric / elderly / senescence / menopause / postmenopausal / "18 month" / "20 month" / "24 month" など
- 薬剤の語: doxorubicin / adriamycin / chemotherapy / "chemotherapy-induced" / cyclophosphamide / fluorouracil / 5-FU / paclitaxel / docetaxel / cisplatin / carboplatin / methotrexate / epirubicin / cytotoxic / busulfan
- 組織の語: "bone marrow" / bone / osteoblast / osteocyte / stromal / niche / LepR / CXCL12 / "CAR cells" / sinusoidal / "bone marrow adipocyte" / "bone loss" / "bone marrow mesenchymal" / femur / tibia など

**ほかに当たった先**: PubMed（論文の検索と、論文から GEO への逆引き）、ProteomeXchange・PRIDE・MassIVE（タンパク質のデータベース）、Metabolomics Workbench（代謝物のデータベース）、ArrayExpress / BioStudies。
ProteomeXchange では doxorubicin / chemotherapy / cyclophosphamide / cisplatin / paclitaxel / fluorouracil / senolytic / bone / osteocyte / osteoporosis / ovariectomy などを検索しました。

**結果**:

- **加齢 × 化療を同じ実験で組んだものは、プロテオームの PXD078161 だけ**でした。
- OVX × ドキソルビシンを 2 × 2 で組んだ実験は、表現型（µCT・子宮重量・qPCR）の論文が 1 本（Yao 2020 *Cancer Res*、PMID 31932453）だけで、データの寄託はありません。この論文では、ドキソルビシンの骨量減少は OVX より大きく、OVX で減った後にさらに上乗せされ、骨の SASP（qPCR）はドキソルビシンでは上がるが OVX では上がりませんでした。ただし 6 週齢の若いマウスです。
- 最も近い RNA のデータは GSE302366（LepR+ 骨髄間質の単一細胞、老齢 × 損傷）でしたが、損傷は骨髄を掻き出す手術で、化療ではありません。
- ProteomeXchange で化療 × 骨・骨髄に当たったのは PXD078161 だけで、他のヒットは腫瘍細胞・卵巣がん・骨肉腫・筋でした。Metabolomics Workbench にも化療 × 骨・骨髄の代謝データはありませんでした。

# 探索 3（10/13 午後）— 3 方向から並行して

GEO のタイトル・要約の総当たりはほぼ出尽くしたので、探し方そのものを 3 通りに変えました。

## 3a. 単一細胞のアトラスとポータル

**当たった先と件数**:

- **CELLxGENE Discover**（単一細胞データの最大の集積サイト）: 全 2,240 データセットのうち、組織が骨・骨髄に当たる 66 データセットを列挙し、間質の細胞を含むものは中身のメタデータ（ドナー・性別・年齢・病態・細胞種）をダウンロードせずに読んで集計しました。
- **Single Cell Portal**（Broad 研究所）: bone marrow stroma / bone marrow niche / osteocyte / osteoblast / bone marrow aging / ovariectomy / chemotherapy bone marrow / 5-FU / doxorubicin / skeletal stem / LepR など 13 語。ovariectomy と osteocyte は 0 件、doxorubicin は腎臓のデータだけでした。
- **HCA Data Portal**（Human Cell Atlas）: 骨髄・骨組織・大腿骨・頭蓋骨など 9 つの臓器区分。骨髄 9 件は造血・免疫が中心で、加齢 × 骨の間質の project はありませんでした。
- **GEO**: 単一細胞の語（6 語を OR でまとめたもの）× 組織の語 17 語 × 条件の語 19 語の総当たり（323 検索）と追加の組み合わせ、計 400 弱の検索。

**結果**: 加齢 × 化療、OVX × 化療の単一細胞データは**どこにもありません**でした。アトラス系（Tabula Muris Senis、Tabula Sapiens、HCA、加齢アトラス）には化療の条件がありません。
化療後の骨髄 niche の単一細胞は、既知のもの（GSE289491、GSE108891 など、すべて若齢）から増えませんでした。
収穫は雌の加齢骨髄間質の GSE210584（前章で説明）です。

## 3b. 培養細胞と、メタデータの直接検索

**ねらい**: GEO のタイトルや要約に「骨」や「化療」の語が無くても、サンプルの属性欄（細胞名・処置）に書いてあれば拾えるように、サンプル単位の情報を直接調べました。

**当たった先と件数**:

- **GEO**: 細胞の語 18 語（osteoblast / osteocyte / MLO-Y4 / Ocy454 / IDG-SW3 / MC3T3-E1 / hFOB / BMSC / osteoclast / RAW264.7 / calvaria など）× 薬剤の語 16 語の総当たり（288 検索）。
- **PubMed → GEO**: 骨系細胞 × 化療 × トランスクリプトームの論文 1,461 報から 204 シリーズ、ヒト骨髄 × 化療 × 小児・固形がん・リンパ腫の論文 1,444 報から 216 シリーズを逆引き。
- **BioSample**（サンプルの属性のデータベース）: 5,169 件を取得し、属性値に薬剤名と骨系細胞の語が両方入るものを集計。
- **SRA**（生データのデータベース）: RNA のデータに絞って、骨系 × 薬剤 89 run、骨・骨髄・間質 × 化療 1,618 run を集計。
- **ARCHS4**（GEO のヒト・マウス RNA-seq をすべて再処理した集積）: 最新版のヒト 1,098,771 サンプル・マウス 1,225,044 サンプルの属性を全部走査。化療薬の名前を持つシリーズはヒト 189・マウス 15 で、骨系は既知のものだけでした。
- **recount3**（同様の再処理集積）: ヒトは全 8,802 プロジェクトを走査済み。**マウスは 10,116 プロジェクト中 9,008 で打ち切り**（残り約 1,100 は未走査）。

**結果**:

- **骨細胞（細胞株・初代・骨片）× 化療薬のトランスクリプトームはありませんでした。** 骨細胞株の公開データ（MLO-Y4 21 シリーズ、Ocy454 / IDG-SW3 22 シリーズ）は、PTH・機械刺激・低酸素・転移などだけでした。
- ヒトの骨系細胞 × 化療薬（または同じ仕組みの DNA 傷害薬）は 3 本見つかり、うち GSE284666（エトポシド）と GSE135323（ドキソルビシン）を陽性対照として解析しました（前章）。
- **ドナーの年齢・性別が分かるヒトの骨芽細胞・骨髄 MSC に化療薬を当てた RNA データはありません**（年齢つきは化療なしの加齢データだけ）。
- メタデータの直接検索で新たに出たのは、GEO のタイトルに「骨」の語を含まない GSE284666 と、GEO を通さない SRP655200 の 2 本だけでした。**GEO のタイトル・要約の検索で拾えていなかったデータは、ほぼ無いと考えてよい**、というのがこの探索の結論です。

## 3c. 論文から逆引き・GEO 以外のリポジトリ

**ねらい**: データが GEO 以外（論文の補足ファイル、Zenodo などの汎用リポジトリ）に置かれている可能性をつぶしました。

**当たった先と件数**:

- **PubMed**（2015–2026、21 クエリ、1,285 件）と **Europe PMC**（20 クエリ、各上位 2,000 件、計約 35,000 件から 501 件を選別）。
- **被引用からの逆引き**: Yao 2020、Wasi 2024、Izquierdo 2026、Ottewell 2015 など鍵になる 11 論文を引用している論文 1,811 件。
- **全文の機械読み**: 上の論文のうちオープンアクセスの全文 **1,988 本**から、Data availability の節、補足ファイル名、本文中の登録番号（GSE、PRJNA、SRP、E-MTAB、PXD、MSV など 14 種類）、リポジトリの URL（Zenodo、figshare、Dryad、Mendeley、OSF、GitHub など）を抜き出し、既に評価した 132 件を除いて目視で選別しました。
- **汎用リポジトリ**: Zenodo、figshare、Dryad、Mendeley Data を 27 語（Mendeley は 23 語）で検索。例: "doxorubicin bone" / "doxorubicin osteocyte" / "chemotherapy bone loss" / "chemotherapy bone marrow niche" / "aged bone marrow niche" / "ovariectomy bone marrow single-cell" / "bone marrow adiposity chemotherapy" / "breast cancer bone metastasis niche" など。
- **臨床試験**: ClinicalTrials.gov で乳がん × 化療 × 骨代謝マーカー / 骨密度 / 骨髄の播種腫瘍細胞を検索。

**結果**:

- **加齢または OVX × 化療の骨・骨髄データは、論文の Data availability にも、Zenodo / figshare / Dryad / Mendeley にもありませんでした。** 汎用リポジトリのヒットの大半は、腫瘍細胞（骨肉腫・白血病）、薬物送達の材料、加齢だけ、OVX × 生薬でした。
- 新しく見つかった「年齢 × 化療」の実験は busulfan の 1 本（PMID 39139448）だけで、雄、若齢 4–6 週 vs 成体 16–20 週、移植前処置の薬で、データの寄託もありません。
- 収穫は、論文の補足 Excel に**個体ごとの数値**が載っている 3 本（Wasi 2024 の骨細胞 TUNEL、Ali 2025 *Sci Rep* の化療後の長期追跡、*Nat Commun* 2025 の RANKL / OPG）で、前章で説明しました。
- **乳がん患者の化療前後の骨代謝マーカーや骨髄の測定値を、個人単位で公開しているデータ**（Project Data Sphere、Vivli、YODA、試験の公開データ）はありませんでした。結果が登録されている試験は 2 本（NCT01298362、NCT00295867）で、群の平均だけです。
  参考になる群平均として、PMID 25923354（閉経前 29 例 / 閉経後 12 例、術前化療 4 サイクル）では、骨吸収マーカーの上昇と骨形成マーカーの低下が**閉経状態によらず**同程度でした（例: β-CTX 閉経前 +36%、閉経後 +32%）。

# 見つからなかった組み合わせ（まとめ）

3 回の探索を通じて、以下はどこにもありませんでした。

1. **加齢 × 化療の骨・骨髄のトランスクリプトーム**（bulk / 単一細胞 / 単一核、マウス・ラット・ヒト）
2. **OVX × 化療の骨・骨髄のオミクス**（RNA・タンパク質・エピゲノム・代謝のどれも）
3. **OVX の雌の骨髄間質・骨芽細胞・骨細胞の単一細胞データ**、および OVX で骨髄 niche の細胞（LepR+ / CAR / 内皮 / 骨髄脂肪細胞）を生体から分取したデータ
4. **化療 × 骨髄 niche で、12 か月以上のマウスを使ったもの**（化療側はすべて 7–16 週）、および雌雄両方に化療を当てたもの
5. **乳がん骨転移モデル × 化療**、**乳がん骨転移モデル × OVX** の骨側のデータ
6. **タキサン × 骨・骨髄間質**、**シクロホスファミドの活性体 × 骨系細胞**、**AC / TC 併用 × 骨**
7. **骨細胞（細胞株・初代・骨片）× 化療薬**、**骨器官培養 × 化療薬**、**破骨細胞 × 化療薬**（RAW264.7 の 1 本を除く）
8. **骨細胞のアポトーシス率（TUNEL・カスパーゼ 3）や空の骨小腔の率を、年齢 × 化療の 2 × 2 で定量した論文**（骨細胞 TUNEL は成熟マウスだけの Wasi 2024 と、幼若ラットだけの 1 本）
9. **化療を受けた乳がん（固形がん）患者の骨・骨髄間質の、化療前後のトランスクリプトーム**、および**個人単位の骨代謝マーカー**
10. **化療 × 骨・骨髄のプロテオーム**（PXD078161 を除く）、**代謝データ**、**DNA メチル化データ**

**注意（探し残し）**: recount3 のマウスは約 1 割（約 1,100 プロジェクト）が未走査です。ただし、同じ GEO の RNA-seq を網羅している ARCHS4 の最新版は全サンプルを走査済みなので、ここに未知のデータが残っている可能性は低いと考えています。
オープンアクセスでない論文は抄録で判断したので、本文の補足にだけデータの所在が書かれている論文は拾えていない可能性があります。

# 唯一の候補 PXD078161 が使えなかった理由

## どんなデータか

Izquierdo et al., *GeroScience* 2026（PMID 42686970）のデータで、タンパク質のデータベース（ProteomeXchange PXD078161 / MassIVE MSV000101777）に置かれています。
若齢（3 か月）と老齢（20 か月）の雌雄のマウスに、対照 / ドキソルビシン（DX）/ 老化細胞除去薬 D+Q / DX + D+Q を与え、脛骨の皮質骨（骨髄を生理食塩水で洗い流し、脱灰して抽出）のタンパク質を質量分析で測ったものです。
著者の結論は「化療による骨の加齢様の変化（海綿骨の減少、骨髄の脂肪化、老化・炎症遺伝子の上昇）は雌でより強く、年齢と性別で異なる」です。設計としては、まさに探していた 2 × 2 です。

## 何をしたか

著者が置いた群間比較の表（Spectronaut の出力）を使い、年齢・性別ごとに「DX vs 対照」のタンパク質の変化（log2 比）を並べ、若齢と老齢で比べました。遺伝子セットの解析（GSEA）も、RNA と同じセットで回しました。

![PXD078161 の解析（`out/prelim/proteome/fig_proteome_dx_age.png`）。A = DX の効果の若齢雌（横軸）と老齢雌（縦軸）の比較、B = 年齢・性別ごとの遺伝子セットの NES（右は RNA の参照）](prelim/proteome/fig_proteome_dx_age.png){width=100%}

## 使えなかった理由

1. **匹数が少なすぎる**: 雌は若齢・老齢とも対照 2 匹 / DX 2 匹。若齢雌は同じ試料を 2 回測った技術的な繰り返しで、老齢雄は匹数の記載がありません。
2. **年齢と測定のバッチが完全に重なっている**: 4 つの年齢・性別の群は、それぞれ別の日に、別の注入量（500 ng と 1 µg）で測られ、解析ファイルも群ごとに別です。年齢による違いと、測定の回の違いを分けられません。
3. **効果の再現性がゼロ**: 共通に測れた 116 タンパク質で、DX の効果は若齢雌と老齢雌で相関しませんでした（Spearman ρ = -0.03）。ところが**同じ年齢の雌と雄でも ρ = 0.06** で、どの 2 群の組み合わせでも相関の絶対値は 0.26 以下でした。
   つまり、「年齢で DX の効果が変わった」のではなく、「DX の効果そのものが群をまたいで再現しない」状態で、年齢の違いを読む土台がありません。
4. **測れているタンパク質が少なく、偏っている**: 群あたり 125–397 種類で、骨の基質（コラーゲンなど）と血漿のタンパク質が中心です。炎症、IFN、エネルギー代謝の遺伝子セットは、どの群でもタンパク質が 5 個未満で検定できませんでした。本研究の仮説で見たい軸（炎症・老化）がほぼ測れていません。
5. **統計の作り方**: 著者の q 値は、ペプチド × サンプルの組を独立の観測として数えた t 検定で、個体の数を反映しません（2 匹 vs 2 匹でも自由度が 1,940 になるタンパク質がある）。そのため q 値ではなく変化の大きさで読みましたが、上の 3 のとおり再現しませんでした。
6. **投与の条件が分からない**: DX の用量・回数・採取時期は、データベースのメタデータにも付属の方法の文書にも書かれておらず、論文本文はこの環境から取得できませんでした（ken がブラウザで PDF を取れば確認できます）。

**唯一そろっていたのは「DX でコラーゲン（基質）が下がる方向」**で、4 群とも同じ向きでしたが、どれも有意ではありません。これは RNA の結果（化療の LepR+ 間質で基質が低下）と向きが合う、という程度の補強にしかなりません。
したがって **Figure 4 には載せず**、本文で引くなら論文の表現型（µCT・骨髄脂肪化・qPCR）を文献として引用する方が強い、と判断しました。
なお論文の補足図には「加齢による変化と DX による変化が正に相関する」（雄 r = 0.35、雌 r = 0.25）という結果があり、「DX は骨で加齢に似た変化を起こす」の文献的な根拠としては使えます。

# 「無い」ことの意味と、次の一手

- **新規性**: 「閉経後（加齢）の骨 × 周術期化療」を同じ枠で調べたデータは公開されておらず、本研究の臨床所見（閉経後でだけ化療と骨転移が関連する）は、基礎側でも誰も検証していない組み合わせです。
- **限界**: dry 解析は、加齢と化療を別々の実験で見て向きを並べることしかできず、「重なると何が起きるか」（交互作用）は示せません。
- **ウェットで埋めるなら**: 打合せのとおり、まずは現状の dry 解析で投稿し、査読で追加実験を求められた時点で検討する方針でよいと考えます。その場合の最小の設計は「若齢 / 中年以上（または sham / OVX）の雌 × 対照 / AC 相当の反復投与、骨（皮質骨）と骨髄 LepR+ 間質の RNA-seq、各群 4 匹以上」で、今回の解析の枠組み（同じ遺伝子セット・同じ閾値）がそのまま使えます。

# Limitation 文案（英文）

> Sixth, no public dataset combined aging or ovariectomy with cytotoxic chemotherapy within the same experiment for bone or bone marrow. We systematically searched GEO, ArrayExpress, BioProject, SRA and BioSample, sample-level metadata of uniformly reprocessed repositories (ARCHS4 and recount3), single-cell resources (CELLxGENE Discover, Single Cell Portal and the Human Cell Atlas), proteomic and metabolomic repositories (ProteomeXchange, PRIDE, MassIVE and Metabolomics Workbench), general-purpose repositories (Zenodo, figshare, Dryad and Mendeley Data), and the data availability statements of approximately 2,000 open-access articles. The only dataset with an aging-by-chemotherapy design was a cortical bone proteome of young and aged mice treated with doxorubicin (PXD078161); however, with two animals per group, age groups measured in separate batches and no reproducibility of the doxorubicin effect even between sexes of the same age, it did not allow the effect of age to be evaluated. All public transcriptomic datasets of chemotherapy-exposed bone or bone marrow were derived from young adult animals (7–16 weeks of age), and no transcriptome of bone or bone marrow stroma from patients with breast cancer before and after chemotherapy was available. Consequently, the effects of aging or menopause and of chemotherapy were analysed in independent datasets and compared only in terms of the direction of pathway-level changes; their interaction could not be tested, and our transcriptomic findings should be regarded as hypothesis-generating.

（v1 の Limitation 案の First–Fifth に続く 6 番目として書きました。番号は本番の構成に合わせて変えます。）

# 清水先生に確認したいこと

1. **「無い」ことを論文のどこに書くか**: Limitation だけでなく、Introduction（「この組み合わせは調べられていない」＝研究の動機）にも 1 文入れるのがよいと考えますが、いかがでしょうか。
2. **PXD078161（Izquierdo 2026）の扱い**: データとしては使わず、論文の表現型（化療による骨の加齢様変化が雌でより強い）を Discussion で引用する方針でよいか。論文本文（用量・時点）を確認したいので、PDF を入手できれば共有をお願いします。
3. **追加実験の見通し**: 査読で求められた場合に備え、上の最小設計（年齢または OVX × AC 相当の反復投与、雌、骨と LepR+ 間質の RNA-seq）を共同研究先（中島先生など）に打診できそうか、感触を伺いたいです。

# 文献調査 2 — 加齢/閉経の骨（A）と化療の骨（B）を骨転移 niche に繋ぐ文献

作成 2026-10-13。第一波（`docs/literature.md`、75 本）の続き。**再調査はしていない**ので、
重複する文献（Ottewell, Yao, Luo, Demaria, Shaked 系など）は PMID だけ再掲し、要点は第一波を参照。

対象の問い: dry 解析で出た **A（加齢/閉経の骨 = 炎症・老化・p53 ↑）** と
**B（化療の骨 = 骨形成/ECM/血管 ↓、間質 ER 標的の消失、IFN ↑、in vitro では p53/老化 ↑）** が、
骨転移 niche（DTC の休眠解除・増殖・破骨活性化）にどう繋がるかを指摘した文献を探す。

**検証手順**: NCBI E-utilities（esearch/efetch）で PMID・書誌・抄録を実取得。要点は抄録の記述のみから書き、
本文を読んでいない推測は書かない。抄録に無いことは「言えない」として (d) に回した。

---

## (a) A と B を繋ぐ最も使える 9 本

1. **Singh N, Kusumbe AP et al. *JCI Insight* 2019（PMID 31292293）** — 加齢した骨の secretome が HSC・MSC・DTC の増殖を誘発し、
   その原因は **骨特異的な pericyte の減少（= 静止シグナルの喪失）**。逆に放射線・化療は PDGF-B を介した type H 血管のリモデリングと
   pericyte 拡大を起こし、**静止促進 secretome を作って治療抵抗性にする**。
   → **A（加齢で niche の静止維持能が落ちる）と B（化療が骨の血管 niche を作り替える）を同じ臓器・同じ軸で扱った唯一の論文**。本論文の「化療で angiogenesis ↓」の受け皿。
2. **Luo X, Stewart SA et al. *Cell Rep* 2016（PMID 26725121）** — 老化骨芽細胞を誘導するだけで乳癌の骨コロニー形成が増える。老化骨芽細胞は石灰化せず、
   **SASP の IL-6 を介して局所の破骨細胞形成を増やす**。IL-6 中和で両方が止まる。
   → 「宿主側の老化だけで骨 niche が整う」先行例。A の炎症・老化が出口（破骨）に繋がる筋の中核。
3. **Tivari S et al. *Cell Commun Signal* 2018（PMID 30119678）** — 骨髄ストローマに酸化・低酸素・**エストロゲン除去（フルベストラント）**で
   secretory senescence を誘導すると、休眠 ER+ MCF-7 が再増殖する。rhIL-6 / IL-8 / TGF-β1 単独でも再増殖を誘導し、
   休眠細胞は ERα 低下・N-cadherin/SLUG 上昇（間葉様）を示す。
   → **A の炎症性 SASP（IL6）と「間質の ER シグナル喪失」の両方が、ER+ の休眠解除に繋がることを in vitro で示した唯一の文献**。本論文の仮説に最も近い。
4. **Price TT, Sipkins DA et al. *Sci Transl Med* 2016（PMID 27225183）** — 休眠乳癌 DTC は骨髄の **E-selectin / SDF-1 に富む perisinusoidal 血管周囲領域**に局在する。
   E-selectin 阻害は骨への entry を、CXCR4 阻害は tether を切り休眠細胞を循環に動員する。患者骨髄でも SDF-1+ 血管に隣接して微小転移が見られる。
   → 「休眠を決めているのは血管 niche」という場所の同定。B の血管 niche 低下を dormancy に結ぶ前提。
5. **Carlson P, Ghajar CM et al. *Nat Cell Biol* 2019（PMID 30664790）** — DTC は血管周囲 niche（PVN）に居ることで、細胞周期の状態と無関係に**化療から保護される**。
   内皮由来 vWF / VCAM-1 と integrin の相互作用を阻害すると DTC が化療に感受性になり、骨転移が予防される。
   → **化療と血管 niche と DTC を同時に扱う**。B が niche を変えること自体が DTC の運命を左右する、の直接の足場（ただし「化療が PVN を壊す」とは言っていない。下記 (d)）。
6. **Kaur J et al. *Cancer Res* 2024（PMID 39312185）** — 乳癌骨転移巣の **osteocyte が早期老化し、破骨誘導性の SASP** を獲得する（scRNA-seq、
   患者試料でも p16Ink4a・telomere dysfunction）。senolytics で骨吸収が止まり骨量が保たれる。
   → 第一波の Farr 2017 / Yao 2020 と繋げると「老化 osteocyte → 破骨 → 骨破壊」が転移の文脈で閉じる。
7. **Farr JN, Khosla S et al. *Nat Med* 2017（PMID 28825716）** — 加齢マウスで老化細胞を除去（INK-ATTAC / senolytics / JAK 阻害）すると骨量・骨質が改善。
   機序は**骨吸収の低下**で、老化細胞の馴化培地は骨芽細胞の石灰化を阻害し破骨細胞形成を促す。
   → A の「老化 ↑」が骨で何を起こすかの標準引用（SASP → 骨形成 ↓ ＋ 破骨 ↑ ＝ 我々の B の骨形成 ↓ と同じ出力）。
8. **Bado IL, Zhang XH-F et al. *Dev Cell* 2021（PMID 33878299）** — **骨の osteogenic niche は ER+ 骨微小転移の ER 発現・活性を一過性かつ可逆的に下げ**（EZH2 による
   エピゲノム再プログラム）、内分泌療法抵抗性と basal/stem 様形質を生む。niche から離れると部分的に戻る。
   → 「骨 niche では ER シグナルが落ちている」。我々の B（化療後の間質で Pgr/Greb1 が消える）と同じ方向の現象を、**腫瘍細胞側**で示した文献。
9. **Albrengues J, Egeblad M et al. *Science* 2018（PMID 30262472）** — 持続炎症が NET を誘導し、NE/MMP9 による laminin の切断 → integrin α3β1 → **休眠細胞の覚醒**。
   → 「A の炎症が休眠を解く」一般機序（臓器は肺）。第一波の He 2025（PMID 40614736, 化療 → NET → 覚醒）と対にすると、A と B が同じ NET 軸に合流する言い方ができる。

（補助）**Bakir B et al. *Int J Mol Sci* 2025（PMID 41465319）** — 骨での dormancy escape の総説。破骨細胞による吸収が TGF-β/IGF-1 を遊離し、
**骨芽細胞由来の静止シグナルの喪失が endosteal niche を壊す**、骨髄脂肪細胞が代謝的に支える、と整理されている。本論文の Discussion の枠組みをそのまま借りられる。

---

## (b) 問いごとの表

### 1. 老化/炎症した骨髄 niche × 骨芽細胞・血管 niche の喪失 → 休眠 DTC の再活性化

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 本論文での使い方 |
|---|---|---|---|---|
| 27225183 | Sci Transl Med 2016 | マウス異種移植 + 患者骨髄 IHC | 休眠乳癌 DTC は E-selectin/SDF-1 に富む perisinusoidal niche に居る。E-selectin が骨への entry、SDF-1/CXCR4 が tether を担い、CXCR4 阻害で休眠微小転移が動員される | 「休眠の場所 = 血管 niche」。B の angiogenesis ↓ を論じる前提 |
| 30664790 | Nat Cell Biol 2019 | マウス、乳癌 DTC | DTC は血管周囲 niche に居ることで細胞周期に依らず化療から守られる。integrin/vWF/VCAM-1 を切ると化療感受性になり骨転移が予防される | 化療 × 血管 niche × DTC を同時に扱う唯一級の論文。「化療期は niche 介入の window」 |
| 31292293 | JCI Insight 2019 | マウス骨、加齢 + 放射線/化療 | 加齢骨の secretome が HSC/MSC/DTC を増殖させる。原因は pericyte 減少。治療（放射線・化療）は type H 血管を PDGF-B 経由で作り替え pericyte を増やし、静止促進 secretome で抵抗性を作る | **A と B を同じ軸で扱える中核**。ただし向きは我々と逆（下記 6） |
| 27761369 | J Bone Oncol 2016 | 総説 | 骨の血管微小環境は DTC の休眠を支え、かつ再活性化も指示する。血管 niche を標的にすれば転移再発を遅らせうる | Discussion の導入（血管 niche が両方向に働く） |
| 27074508 | Nature 2016 | マウス骨、加齢 | 加齢で HSC niche を作る血管（type H、動脈、PDGFRβ+ 血管周囲細胞）が強く減る。内皮 Notch 活性化で回復する | A 側の「加齢で骨の血管 niche が減る」直接証拠。我々のマウス加齢 angiogenesis ↓ と整合 |
| 29514796 | Cancer Res 2018 | 前立腺癌 + マウス大腿骨 | **分化した**骨芽細胞の分泌因子（GDF10/TGF-β2）が TGFβRIII-p38MAPK-RB を介して DTC を静止させる。未分化骨芽細胞では起きない | **骨芽細胞 niche の喪失 = 休眠維持シグナルの喪失**。B の ossification ↓ を dormancy に接続する最良の分子根拠 |
| 22137794 | Cancer Cell 2011 | 乳癌骨転移休眠モデル | VCAM-1（NF-κB 依存）が α4β1+ 破骨前駆細胞を呼び、indolent な微小転移を overt な骨溶解病変に転換させる。抗 VCAM-1/α4 で阻止 | 「炎症 → 破骨前駆の動員 → 休眠離脱」の分子経路（A の出口） |
| 30262472 | Science 2018 | マウス、肺 | 持続炎症が NET を誘導し、NE/MMP9 の laminin 切断 → α3β1 → 休眠細胞の覚醒 | A の炎症が休眠を解く一般機序（臓器は肺、骨では未検証） |
| 41465319 | Int J Mol Sci 2025 | 総説 | 骨での dormancy escape: 破骨による吸収が TGF-β/IGF-1 を遊離、**骨芽細胞由来の静止シグナルの喪失**、骨髄脂肪細胞の代謝的支援、免疫監視の破綻 | Discussion の枠組み・最新総説として |
| 28956197 | Angiogenesis 2017 | 総説（ラット/マウス） | 化療・放射線は骨髄 sinusoidal endothelium を傷害する。骨髄微小血管は造血・骨形成・骨転移と相互作用する | 「化療が骨髄血管を壊す」の総説的根拠（B の angiogenesis ↓ の裏付け） |
| 23728425 | Nat Cell Biol 2013 | — | （第一波 1-7）血管周囲 niche の TSP-1 が休眠を維持し、内皮の活性化・sprouting で破れる | 第一波から再掲 |
| 26632274 | Nat Commun 2015 | — | （第一波 1-5）破骨細胞によるリモデリングが休眠骨髄腫細胞を再活性化 | 第一波から再掲 |

### 2. 間質の ER シグナル喪失と ER+ DTC の制御

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 本論文での使い方 |
|---|---|---|---|---|
| 33878299 | Dev Cell 2021 | ER+ 乳癌骨微小転移（マウス・バーコード） | 骨の osteogenic niche が ER の発現・活性を一過性/可逆的に下げ（EZH2 依存のエピゲノム再プログラム）、内分泌抵抗性と basal/stem 様形質を生む。niche を離れると戻る | **「骨 niche にいる ER+ 細胞は ER-low になる」**。我々の「化療後の間質で ER 標的が消える」と並べ、骨の ER シグナルが可塑的であることを示す |
| 28368409 | Oncogene 2017 | 総説 | ER は乳癌の driver であると同時に骨の発生・恒常性でも中心的。ER を標的にする治療は必然的に骨芽細胞・破骨細胞を含む微小環境にも作用する | **「抗 ER 治療 = 微小環境治療でもある」**を 1 文で引ける。間質 ER を論じる正当化 |
| 30119678 | Cell Commun Signal 2018 | ヒト MCF-7 休眠モデル + マウス間質 | **フルベストラントによるエストロゲン除去**を含むストローマ傷害で secretory senescence を誘導すると、休眠 ER+ 細胞が再増殖。IL-6/IL-8/TGF-β1 単独でも覚醒し、休眠細胞では ERα が低い | 間質の ER 喪失 → SASP → ER+ DTC 覚醒、という一本道をそのまま引用できる |
| 24687923 / 25603921 | Clin Cancer Res 2014 / Int J Cancer 2015 | OVX マウス | （第一波 4-1, 1-4）OVX のみで骨内の休眠 DTC が増殖し骨転移化、ZOL/OPG-Fc で阻止。骨外腫瘍には無関係 | 1st hit 側の中核（第一波から再掲） |
| 24451543 | Nature 2014 | マウス HSC | エストロゲンが ERα を介して HSC の自己複製を増やす（雌・妊娠時）。骨髄の幹細胞がエストロゲンに直接応答する | 「骨髄の間質・幹細胞は ER 応答性である」の一般的裏付け（DTC ではない） |

### 3. IFN 応答が加齢骨・化療後骨の両方で上がる意味

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 本論文での使い方 |
|---|---|---|---|---|
| 22820642 | Nat Med 2012 | 乳癌骨転移 | 腫瘍細胞側の **Irf7 経路の沈黙が免疫逃避を介して骨転移を促す**（IFN は宿主の NK/T による監視を支える） | IFN ↑ を「抑制側」として読む主柱 |
| 26198985 | Cancer Immunol Res 2015 | Ifnar1−/− マウス、複数の乳癌モデル | **宿主**の I 型 IFN シグナルの喪失が骨転移を加速し、NK 機能を損なう | 宿主側 IFN の喪失 = 骨転移増。我々の IFN ↑ を「防御反応」と読む根拠 |
| 32314873 | EMBO Rep 2020 | 前立腺癌骨転移（scRNA + 患者） | 骨で増殖する癌細胞では腫瘍内在性 I 型 IFN が失われ、免疫原性と治療応答が落ち、骨細胞の活性化を促す。HDAC 阻害で回復 | 「骨で休眠が破れるとき IFN は失われている」= IFN ↑ は休眠側 |
| 30546090 | Oncogene 2019 | ER− 乳癌、Dox/MTX 後 | 化療を生き延びた細胞で **IRF7/IFN-β/IFNAR 軸が持続活性化し、免疫学的 dormancy を作る**。IFNAR 阻害で dormancy が解除される。患者血中 IFN-β でも確認 | **化療 → IFN ↑ → 休眠維持**。我々の「化療で IFN ↑」の最も直接的な先行例（ただし腫瘍細胞側） |
| 11117749 | Nature 2000 | マウス骨/免疫 | IFN-γ が TRAF6 を分解して RANKL-RANK シグナルを遮断し、破骨細胞形成を強く抑制する | **IFN ↑ は破骨を抑える方向**。我々で RANKL/破骨が動かないことと整合 |
| 11961557 | Nature 2002 | マウス骨 | RANKL が c-Fos 依存に IFN-β を誘導し、IFN-β が c-Fos を抑えて破骨細胞分化を抑制する（自己制御ループ） | 同上。「IFN ↑ ＝ 骨破壊 ↑」と書いてはいけない根拠 |
| 37512126 | Medicina 2023 | マウス局所照射 + MLO-Y4 骨細胞 | 照射が **STING-TBK1 経路**を介して骨微小環境を破壊し細胞老化を起こす。IR による骨微小環境の破壊は骨量減少を招き、**腫瘍転移に好都合な環境を作る** | ゲノム傷害 → cGAS/STING → IFN → 骨 niche 破壊、を骨で示した唯一の近傍文献（化療ではなく放射線） |

### 4. 老化骨芽細胞/骨細胞の SASP が破骨を活性化 / 化療誘発老化が同じ経路に乗るか

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 本論文での使い方 |
|---|---|---|---|---|
| 28825716 | Nat Med 2017 | 加齢マウス | 老化細胞の除去（遺伝学的・senolytics・JAK 阻害）で骨量・骨質が改善。機序は骨吸収の低下。老化細胞の馴化培地は骨芽細胞の石灰化を阻害し破骨細胞形成を促す | A の老化 ↑ が「骨形成 ↓ ＋ 破骨 ↑」を生む標準引用 |
| 26725121 | Cell Rep 2016 | マウス骨 + 乳癌 | 老化骨芽細胞の誘導だけで乳癌の骨コロニー形成が増える。老化骨芽細胞は石灰化せず、SASP の **IL-6** で局所破骨細胞形成を増やす。IL-6 中和で阻止 | 老化 → IL-6 → 破骨 → 転移、を一本で持つ最重要文献（第一波 2-2 の再掲） |
| 39312185 | Cancer Res 2024 | マウス・患者の骨転移 | 骨転移巣の osteocyte が早期老化し破骨誘導性 SASP を獲得。senolytics で骨吸収抑制・骨量保持 | 老化 osteocyte が破骨の引き金であることを**転移の文脈**で示す最新例 |
| 36305580 | Elife 2022 | DMP1+ osteocyte 部分除去マウス | osteocyte を減らすと骨粗鬆・サルコペニア・寿命短縮。骨形成 ↓ / 破骨 ↑ に傾き、骨髄は骨髄系へ偏り、骨系・骨髄系の両方で **SASP を獲得**する | 「osteocyte の喪失そのものが骨と骨髄に SASP を広げる」= 2nd hit の下流 |
| 31932453 | Cancer Res 2020 | 化療マウス | （第一波 2-1）化療後の骨減少は老化/SASP 依存で、エストロゲン非依存。p38MAPK-MK2 阻害・senolytics で保護 | 化療誘発老化が A と同じ経路に乗ることの直接証拠 |
| 30913313 | J Bone Miner Res 2019 | 若齢マウス OVX/ORCH + ヒト | **若年成体の性ホルモン欠乏（OVX/ORCH 2 か月）は老化マーカーも SASP も誘導しなかった**。エストロゲン欠乏と細胞老化は骨粗鬆症の**独立した**機序 | **「閉経が老化細胞の蓄積を速める」とは書けない**ことの根拠。A の炎症と老化を分けて書く |
| 32794138 | Curr Osteoporos Rep 2020 | 総説 | osteocyte の細胞老化の総説 | 用語・概念の導入 |
| 35974106 | Nat Commun 2022 | ヒト・マウス骨/骨髄 | SenMayo 遺伝子セットの作成と検証（加齢ヒト骨生検で濃縮、senolytics 後に低下、scRNA で老化細胞を同定） | 我々が SenMayo を使う方法論の引用 |

### 5. 加齢/閉経 niche × 化療を同時に扱った総説・展望

| PMID | 雑誌・年 | 系 | 主張（抄録より） | 本論文での使い方 |
|---|---|---|---|---|
| 31292293 | JCI Insight 2019 | 原著 | 加齢（pericyte 減少 → 静止喪失）と治療（type H 血管リモデリング → 静止促進）を**同じ論文で**扱う | (a) 1 を参照。本論文に最も近い枠組み |
| 33324652 | Front Cell Dev Biol 2020 | 総説 | ストレス・加齢・疾患における骨の血管新生と血管 niche のリモデリング | 「骨の血管 niche は加齢とストレスで変わる」の総説枠 |
| 33875860 | Nat Rev Clin Oncol 2021 | 総説 | （第一波 1-12）骨転移微小環境を標的とする治療 | 第一波から再掲 |
| 35477650 / 26029361 | Semin Oncol Nurs 2022 / BoneKEy 2015 | 総説 | （第一波 2-9, 2-8）CTIBL を直接作用（骨細胞毒性）と間接作用（卵巣機能・ステロイド・不活動）に分けて整理 | 「閉経後女性の化療性骨減少」の総説枠。第一波から再掲 |
| 37999904 / 31645711 | Clin Exp Metastasis 2024 / Nat Rev Cancer 2019 | 総説 | （第一波 3-7, 3-6）治療に対する宿主反応が転移を助ける | 第一波から再掲 |
| 36568212 | Front Oncol 2022 | 肺癌、シスプラチン | **化療で休眠化した癌細胞由来のエクソソーム（IGF-2/IGFBP2）が骨髄間質の解糖を活性化し、骨髄の pre-metastatic niche 形成を促す** | 「化療 → 骨髄 niche 形成」を骨髄で示した数少ない例（腫瘍由来 EV 経由、乳癌ではない） |
| 41465319 | Int J Mol Sci 2025 | 総説 | 骨の dormancy escape（上記 1） | Discussion の骨格 |

→ **「加齢/閉経 niche × 化療」を骨転移の文脈で正面から扱った総説は見つからなかった**（第一波の「2-hit の先行例なし」と同じ結論）。

### 6. 逆の報告・仮説に不利な報告

| PMID | 雑誌・年 | 主張（抄録より） | 扱い |
|---|---|---|---|
| 30546090 | Oncogene 2019 | 化療後に持続する IRF7/IFN-β/IFNAR 軸は**免疫学的 dormancy を維持**し、IFNAR 阻害で dormancy が解除される | 我々の「化療で IFN ↑」は**転移を抑える**方向かもしれない。両論併記が必須 |
| 26198985 | Cancer Immunol Res 2015 | 宿主 I 型 IFN の**喪失**が骨転移を加速する | 同上 |
| 11117749 / 11961557 | Nature 2000 / 2002 | IFN-γ も IFN-β も破骨細胞分化を抑制する | IFN ↑ を破骨活性化の根拠にはできない。我々で破骨系が動かないことと一致 |
| 31292293 | JCI Insight 2019 | 放射線・化療は骨の pericyte を増やし**静止促進** secretome を誘導、加齢では逆に pericyte が減って増殖性になる | **向きが我々と逆**（化療が niche を「守る」側）。我々の angiogenesis ↓ と矛盾しうる |
| 30913313 | J Bone Miner Res 2019 | 若年成体の OVX/ORCH では老化マーカーも SASP も増えない | 「閉経が老化細胞を増やす」は言えない |
| 38457059 | Stem Cell Rev Rep 2024 | 老化促進マウス（SAMP6）の骨髄 MSC は早期老化し、その **SASP は乳癌細胞の viability を下げ、dormancy を誘導する**側に働く | 老化間質 = 転移促進、とは限らない |
| 35650435 | Nature 2022 | 加齢した肺間質は休眠からの脱出を促すが、**加齢した皮膚は増殖を抑えて播種を促す** = 加齢の効果は臓器依存 | 「加齢 = 一律に転移促進」とは書けない |
| 40389522 | Sci Rep 2025 | （第一波 2-5）AC-T 相当で senolytic は骨減少を防げず、後期には老化マーカーも検出されない | 第一波から再掲（老化仮説の限界） |
| 31806543 | Lancet Oncol 2020 | （第一波 6-6）D-CARE 陰性 | 第一波から再掲 |

---

## (c) Discussion の段落案（英文、約 200 語）

> Our transcriptomic survey argues that menopause/ageing and chemotherapy do not converge on a single bone
> program: ageing and postmenopausal bone is dominated by inflammation and senescence (IL6, TNF, SenMayo, p53),
> whereas doxorubicin-exposed murine bone stroma loses its osteogenic, matrisome and angiogenic programs and its
> oestrogen-responsive signature without any inflammatory induction. Each change has, independently, been linked
> to the release of disseminated tumour cells (DTCs) from dormancy. Senescent osteoblasts alone are sufficient to
> increase breast cancer colonisation of bone through IL-6-driven osteoclastogenesis (PMID 26725121), senescent
> osteocytes acquire a pro-osteoclastogenic secretome (PMID 39312185), and stromal secretory senescence — including
> that induced by oestrogen withdrawal — reawakens dormant ER+ cells (PMID 30119678). Conversely, dormancy in bone
> is actively imposed by the very compartments that chemotherapy erodes: differentiated osteoblasts secrete GDF10
> and TGF-β2 that enforce quiescence (PMID 29514796), and dormant DTCs reside in E-selectin/SDF-1-rich
> perisinusoidal vessels (PMID 27225183) whose age- and therapy-dependent remodelling alters DTC behaviour
> (PMID 27074508, 31292293, 30664790). A plausible, but as yet untested, synthesis is that the first insult supplies
> an osteoclastogenic, inflammatory signal while the second removes the osteoblastic and vascular signals that keep
> DTCs quiescent. **No study has examined both insults together in bone, and we therefore present this as a
> hypothesis rather than an interaction we have demonstrated; our data show two separate shifts, not their product.**

（語数 ≒ 215。最後の 1 文が交互作用の直接証拠が無いことの明示。）

---

## (d) 繋がらなかった点

1. **「加齢/閉経 niche × 化療」を骨（骨転移）で同時に扱った文献は、原著・総説とも見つからなかった。**
   第一波の結論（2-hit の先行例なし）と同じ。本論文の novelty であると同時に、Discussion は仮説の形でしか書けない。
2. **「化療が骨芽細胞 niche・血管 niche を壊すことで DTC が目覚める」を直接示した論文は無い。**
   揃っているのは (i) 骨芽細胞/血管 niche が休眠を維持する（PMID 29514796, 27225183, 23728425）、
   (ii) 化療が骨髄血管・骨芽細胞を傷める（PMID 28956197、第一波 39147353）—— の 2 つで、
   **(i)+(ii) を繋いだ実験は存在しない**。Carlson 2019（PMID 30664790）は近いが、主張は逆向き
   （niche が DTC を化療から守る）であって「化療が niche を壊す」ではない。
3. **IFN の向きが読めない。** 我々のデータでは加齢骨・化療骨の両方で IFN 応答が上がるが、文献上は
   I 型 IFN は (a) 宿主側で転移抑制（PMID 26198985, 22820642）、(b) 腫瘍側で dormancy 維持（PMID 30546090, 32314873）、
   (c) 破骨細胞分化の抑制（PMID 11117749, 11961557）と、**いずれも「転移を抑える」側**に働く。
   「化療と加齢に共通する唯一の軸が IFN」という我々の所見は、そのままでは 2-hit を支持しない。
   Discussion では「IFN は宿主の防御反応で、他の 2 軸（炎症・niche 喪失）を相殺しうる」と書くのが正確。
4. **間質（骨側）の ER が DTC を直接制御するという報告は見つからなかった。** 見つかったのは
   (i) 骨 niche が**癌細胞側**の ER を下げる（PMID 33878299）、(ii) ER 標的薬は必然的に骨微小環境にも効く（PMID 28368409）、
   (iii) エストロゲン除去を含むストローマ傷害が ER+ 休眠細胞を起こす（PMID 30119678）まで。
   「骨髄間質の ERα が DTC 休眠を制御する」を示した conditional KO 等の実験は**該当なし**。
   我々の「化療後に間質の Pgr/Greb1 が消える」は、先行文献の無い所見として書くことになる。
5. **化療 → cGAS–STING → IFN を骨髄で示した文献は見つからなかった。** 最も近いのは局所照射の
   PMID 37512126（STING-TBK1 → 骨微小環境の破壊・老化・転移に好都合な環境）で、薬剤ではなく放射線。
6. **「高齢の骨髄が転移を抑える」という骨限定の逆報告は見つからなかった**（肺・皮膚の臓器依存性は PMID 35650435、
   老化 MSC が dormancy を誘導する例は PMID 38457059 まで）。したがって「加齢骨は一律に転移促進」とは書けないが、
   骨で明確に反証されてもいない。
7. 検索で確認できなかった固有名: Kolb/Bussard の **dormancy** 論文（骨転移・骨芽細胞の論文は多数あるが
   休眠を主題にしたものは PubMed で特定できず、`Kolb AD AND dormancy bone` はヒット 2 件・いずれも別主題）。
   Price 2016 は PMID 27225183 として確認済み（第一波の「未確認」欄は本稿で解決）。

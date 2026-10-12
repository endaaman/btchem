# 文献調査 — 周術期化学療法と骨転移（2-hit 仮説）

作成 2026-10-12。対象論文: ER+ 早期乳がん 785 例（BMA 不使用）で、周術期化学療法が**骨転移のみ**
HR 2.05、**閉経後でのみ**術後化療が独立不良因子 HR 1.72。仮説は「閉経（エストロゲン欠乏）による
骨の炎症・破骨細胞誘導（1st hit）＋ 化療・ステロイドによる骨細胞アポトーシス/老化（2nd hit）
→ 骨の pre-metastatic niche」。

**検証手順**: すべて NCBI E-utilities（esearch/esummary/efetch）で PMID・書誌・抄録を実取得して確認。
PMID が確認できなかったものは末尾「未確認／見つからなかったもの」に分離した（本文中には載せない）。

---

## 1. 骨の pre-metastatic niche と DTC の dormancy / 再活性化

| # | 文献 | PMID / DOI | 要点 | 使い所 |
|---|---|---|---|---|
| 1-1 | Kaplan RN et al. *Nature* 2005 | PMID 16341007 / 10.1038/nature04186 | VEGFR1+ の骨髄由来造血前駆細胞が転移に先立って標的臓器に動員され、pre-metastatic niche を形成する。腫瘍細胞到着前に「土壌」が変化するという概念の原典。 | Intro: niche 概念の導入（腫瘍細胞より先に土壌が変わる）。 |
| 1-2 | Peinado H et al. *Nat Rev Cancer* 2017 | PMID 28303905 / 10.1038/nrc.2017.6 | pre-metastatic niche の総説。全身性因子（分泌因子・EV・骨髄動員）が臓器特異的に niche を準備する枠組みを整理。 | Intro: niche の一般論を 1 文で引く総説枠。 |
| 1-3 | Braun S et al. *NEJM* 2005 | PMID 16120859 / 10.1056/NEJMoa050434 | 乳がん 4,703 例のプール解析。診断時の骨髄 DTC が独立予後因子で、骨髄は早期から DTC が存在する器官。 | Intro: 「骨には早期から DTC が居る」＝ 2nd hit で目覚める前提の臨床的根拠。 |
| 1-4 | Ottewell PD et al. *Int J Cancer* 2015 | PMID 25603921 / 10.1002/ijc.29439 | OVX で骨の dormant MDA-MB-231 が増殖して骨転移化する。これは破骨細胞依存で、OPG-Fc（RANKL 阻害）で抑制される。 | Discussion の中核: 「骨代謝回転の上昇が dormant DTC を起こす／RANKL 阻害で止まる」。 |
| 1-5 | Lawson MA et al. *Nat Commun* 2015 | PMID 26632274 / 10.1038/ncomms9983 | 骨髄の endosteal niche において、破骨細胞による骨リモデリングが dormant 骨髄腫細胞を再活性化する（抗吸収薬で dormancy 維持）。 | Discussion: 破骨細胞活性 = dormancy 解除スイッチ、という細胞レベルの機序。 |
| 1-6 | Wang H et al. *Cancer Cell* 2015 | PMID 25600338 / 10.1016/j.ccell.2014.11.017 | osteogenic niche（骨芽細胞との異種接着・Cx43/N-cadherin）が播種直後の乳癌細胞の初期増殖を促す。 | Intro/Discussion: 骨の常在細胞そのものが niche の能動的構成要素であること。 |
| 1-7 | Ghajar CM et al. *Nat Cell Biol* 2013 | PMID 23728425 / 10.1038/ncb2767 | 血管周囲 niche（内皮由来 TSP-1）が乳癌 DTC の dormancy を維持し、内皮の活性化で dormancy が破れる。 | Discussion: dormancy は niche 側の状態で決まる（宿主側の変化が引き金）。 |
| 1-8 | Shiozawa Y et al. *J Clin Invest* 2011 | PMID 21436587 / 10.1172/JCI43414 | 前立腺癌細胞が HSC niche（CXCL12/CXCR4 軸）を奪って骨髄に足場を作る。niche の競合という考え方。 | Intro: CXCL12 軸を骨 niche の分子として導入。 |
| 1-9 | Riquelme MA et al. *Front Endocrinol* 2020 | PMID 33162934 / 10.3389/fendo.2020.567844 | 骨細胞（osteocyte）と骨転移の総説。osteocyte が RANKL・sclerostin・IL-6・CCL5 等を介して腫瘍と双方向に作用し、骨転移を修飾する。 | Fig4 の「皮質骨＝骨細胞寄り」を骨転移に接続する総説。 |
| 1-10 | Pin F, Bonewald LF et al. *Curr Osteoporos Rep* 2021 | PMID 34773212 / 10.1007/s11914-021-00712-9 | 「Osteocytes and Cancer」総説。癌・癌治療が osteocyte を障害し、osteocyte 由来因子が腫瘍増殖と骨破壊を支える。 | Intro: osteocyte を 2-hit の収束点として据える根拠。 |
| 1-11 | Yi Y et al. *Oncogene* 2026 | PMID 42426271 / 10.1038/s41388-026-03887-6 | 乳癌骨転移における osteocyte を統合的制御因子として位置づけ直す最新総説。 | Intro の最新 citation（骨細胞中心の見方は arguable でなく今の潮流だと示す）。 |
| 1-12 | Hofbauer LC et al. *Nat Rev Clin Oncol* 2021 | PMID 33875860 / 10.1038/s41571-021-00499-9 | 骨転移微小環境を標的とする治療の総説（RANKL/破骨細胞・免疫・niche）。 | Discussion: 「臨床で介入できるのは実質 RANKL/破骨細胞だけ」の根拠付け。 |
| 1-13 | Nakashima T et al. *Nat Med* 2011 | PMID 21909105 / 10.1038/nm.2452 | osteocyte が骨の主要な RANKL 供給源であることを条件付き KO で証明。 | 機序の要: 骨細胞障害 → RANKL → 破骨細胞 という連結。 |
| 1-14 | Xiong J et al. *Nat Med* 2011 | PMID 21909103 / 10.1038/nm.2448 | 上記の独立した確認（matrix-embedded cells が骨吸収を制御）。 | 1-13 と併記して「osteocyte-RANKL」を堅くする。 |
| 1-15 | Kennedy OD et al. *Bone* 2012 | PMID 22342796 / 10.1016/j.bone.2012.01.025 | 骨の微小損傷モデルで、アポトーシスした osteocyte の近傍の生存 osteocyte が RANKL を産生し局所の骨吸収を起こす（集団ごとに役割が分かれる）。 | 2nd hit の機序: 骨細胞アポトーシスが RANKL 産生と局所吸収に直結する。 |

## 2. 化学療法が骨に与える影響（CTIBL・骨細胞アポトーシス・老化/SASP）

| # | 文献 | PMID / DOI | 要点 | 使い所 |
|---|---|---|---|---|
| 2-1 | Yao Z, Stewart SA et al. *Cancer Res* 2020 | PMID 31932453 / 10.1158/0008-5472.CAN-19-2348 | 「Therapy-Induced Senescence Drives Bone Loss」。化療後の骨減少はエストロゲン非依存的な老化細胞・SASP による。老化細胞除去または p38MAPK-MK2 阻害で骨が保たれる。 | **Fig4B の解釈の柱**。化療の骨作用が「卵巣機能低下の二次効果」だけでないことの直接証拠。 |
| 2-2 | Luo X, Stewart SA et al. *Cell Rep* 2016 | PMID 26725121 / 10.1016/j.celrep.2015.12.016 | 老化（senescent）骨芽細胞を誘導すると、腫瘍由来シグナルが無くても乳癌の骨コロニー形成が増える。機序は SASP 因子による局所破骨細胞誘導。 | **2-hit 仮説に最も近い先行例**。「宿主側ストローマの変化だけで niche が整う」。 |
| 2-3 | Wasi M et al. *Bone* 2024 | PMID 39147353 / 10.1016/j.bone.2024.117235 | 成熟マウスで doxorubicin が osteocyte アポトーシスを増やし、リモデリングを乱して骨構造を悪化させる。機械刺激/Yoda1 では救済できないほど強い。 | **Dox → 骨細胞アポトーシス**の直接証拠（Fig4B の in vivo 裏付け）。 |
| 2-4 | Izquierdo DJ et al. *GeroScience* 2026 | PMID 42686970 / 10.1007/s11357-026-02485-4 | Dox で海綿骨減少・骨髄脂肪増加・老化/炎症遺伝子上昇。効果は性・年齢依存で**雌でより顕著**、senolytic (D+Q) で部分回復。 | 「化療の骨障害は雌／高齢でより強い」= 閉経後群に偏る臨床所見と整合。 |
| 2-5 | Ali MM et al. *Sci Rep* 2025 | PMID 40389522 / 10.1038/s41598-025-01717-5 | AC-T 相当レジメンで急速な骨量減少が起き 12 か月後も持続するが、senolytic 2 種で予防できず、後期には老化マーカーも検出されない。 | **反証として必ず引く**。老化細胞仮説の限界＋「骨障害は長期持続する」点は仮説に有利。 |
| 2-6 | Shapiro CL et al. *J Clin Oncol* 2001 | PMID 11454877 / 10.1200/JCO.2001.19.14.3306 | 化療後に卵巣機能不全となった女性は 6 か月で腰椎 BMD −4.0%、12 か月までさらに −3.7%、骨代謝マーカー上昇。卵巣機能が保たれた群では減少なし。 | 臨床の 1st+2nd hit の重なり: 化療が「医原性閉経」を作り骨代謝回転を上げる。 |
| 2-7 | Waqas K et al. *J Bone Oncol* 2021 | PMID 33948427 / 10.1016/j.jbo.2021.100355 | 早期乳癌の CTIBL 管理の最新ガイダンス（閉経前後別の骨折リスクと抗吸収薬の使い方）。 | Discussion: 「骨を守る介入は既に標準化されている」の足場。 |
| 2-8 | Hadji P. *BoneKEy Rep* 2015 | PMID 26029361 / 10.1038/bonekey.2015.60 | 乳癌の治療関連骨減少（CTIBL）総説。化療・卵巣抑制・AI それぞれの骨への寄与を整理。 | Intro で CTIBL を一語で導入する総説枠。 |
| 2-9 | Rayner-Myers SD et al. *Semin Oncol Nurs* 2022 | PMID 35477650 / 10.1016/j.soncn.2022.151280 | 補助化療の骨減少を「直接作用（骨細胞への細胞毒性）」と「間接作用（卵巣機能・ステロイド・不活動）」に分けて統合レビュー。 | 直接/間接の切り分けを提示するときの引用。 |
| 2-10 | Oostra DR et al. *Mol Cell Endocrinol* 2015 | PMID 25575458 / 10.1016/j.mce.2014.12.028 | 補助化療後の骨減少と血中 OPG の関連を前向きに検討（RANKL/OPG 系が臨床で動く）。 | 「RANKL/OPG 軸がヒトでも化療後に動く」橋渡し。 |
| 2-11 | Chandra A et al. *Aging Cell* 2022 | PMID 35363946 / 10.1111/acel.13602 | 放射線後の骨減少・骨髄脂肪増加は p21+ 老化細胞の除去で予防でき、p16+ では予防できない。 | 「ゲノム傷害性治療 → 老化細胞 → 骨減少」の一般性（Cdkn1a に注目する根拠）。 |
| 2-12 | Chandra A et al. *J Bone Miner Res* 2020 | PMID 32023351 / 10.1002/jbmr.3978 | 局所放射線後の骨減少が senolytic（D+Q）で軽減。老化細胞負荷が骨減少の因果に乗る。 | 2-11 と対で引く。 |
| 2-13 | Demaria M et al. *Cancer Discov* 2017 | PMID 27979832 / 10.1158/2159-8290.CD-16-0241 | 化療誘発の老化細胞と SASP が有害事象を増やし、さらに**再発・転移を促進**する。老化細胞除去で改善。 | 2nd hit（老化）を「骨減少」から「転移促進」へ接続する要文献。 |

## 3. 化学療法の「宿主反応」による転移促進

| # | 文献 | PMID / DOI | 要点 | 使い所 |
|---|---|---|---|---|
| 3-1 | Karagiannis GS et al. *Sci Transl Med* 2017 | PMID 28679654 / 10.1126/scitranslmed.aan0026 | 術前化療（AC→パクリタキセル）後の残存腫瘍で TMEM 密度・MENA^INV が増え、血管内侵入と遠隔転移が増える。TIE2 阻害で可逆。 | **Intro の主柱**: 「化療は腫瘍を縮めながら転移リスクを上げ得る」臨床付き証拠。 |
| 3-2 | Keklikoglou I et al. *Nat Cell Biol* 2019 | PMID 30598531 / 10.1038/s41556-018-0256-3 | タキサン・アントラサイクリンが annexin A6 に富む腫瘍由来 EV を誘導 → 肺内皮 NF-κB 活性化、CCL2、Ly6C+CCR2+ 単球動員で pre-metastatic niche を形成。患者 EV でも ANXA6 検出。 | **「化療が niche を作る」**の分子的中核。骨で同じことが起きうるという類推。 |
| 3-3 | Chang YS et al. *PNAS* 2017 | PMID 28784776 / 10.1073/pnas.1700455114 | 非腫瘍宿主細胞のストレス応答遺伝子 Atf3 が、パクリタキセルによる転移増悪に必要。PTX は CTC を増やし肺の「土壌」を腫瘍に好適化（炎症性単球増加）。 | 「2nd hit は宿主側で起こる」ことを遺伝学的に示す例。 |
| 3-4 | Volk-Draper L et al. *Cancer Res* 2014 | PMID 25274031 / 10.1158/0008-5472.CAN-14-0067 | パクリタキセルが TLR4 依存に炎症を誘導し、局所だけでなく血液・脾・**骨髄**・肺で炎症メディエーターを上げ、TLR4+ 腫瘍の転移を増やす。 | 化療の全身炎症が**骨髄を含む**点が本論文に直結。 |
| 3-5 | Daenen LG et al. *Cancer Res* 2011 | PMID 21975929 / 10.1158/0008-5472.CAN-11-0627 | シスプラチン/パクリタキセルで**前処理した宿主**（薬は既に消失）に腫瘍細胞を入れると肺転移が増える。内皮 VEGFR-1 上昇が機序で、抗 VEGFR-1 で回復。 | 「薬が抜けた後も宿主が転移許容的に残る」= 周術期化療 → 数年後の骨転移という時間差の説明。 |
| 3-6 | Shaked Y. *Nat Rev Cancer* 2019 | PMID 31645711 / 10.1038/s41568-019-0209-6 | 「The pro-tumorigenic host response to cancer therapies」総説。治療が骨髄由来細胞動員・内皮活性化・サイトカインを介して転移を助ける機序を体系化。 | Intro の総説枠（この分野があることを 1 文で示す）。 |
| 3-7 | Deo A, Shaked Y et al. *Clin Exp Metastasis* 2024 | PMID 37999904 / 10.1007/s10585-023-10243-5 | 上記の最新版総説（宿主反応と耐性・転移・臨床含意）。 | 3-6 の更新 citation。 |
| 3-8 | Shaked Y et al. *Cancer Cell* 2008 | PMID 18772115 / 10.1016/j.ccr.2008.08.001 | 化療が急性に内皮前駆細胞を動員する（宿主反応の古典）。 | 骨髄由来細胞動員の原典的 citation。 |
| 3-9 | He D et al. *Cancer Cell* 2025 | PMID 40614736 / 10.1016/j.ccell.2025.06.007 | Dox・シスプラチンが dormant 乳癌 DTC を目覚めさせ肺転移を増やす。DormTracer で直接証明、機序は NET 形成。 | **最新・最強の「化療が dormancy を破る」証拠**（ただし臓器は肺）。 |
| 3-10 | Shi S et al. *Cancer Res* 2026 | PMID 42696464 / 10.1158/0008-5472.CAN-26-0151 | 化療が老化肝細胞を介した NET で肝の dormant 膵癌 DTC を再活性化。 | 「化療 → 宿主細胞老化 → DTC 覚醒」の臓器横断パターン（骨で未検証と言える根拠）。 |
| 3-11 | Vorontsova A et al. *J Clin Med* 2020 | PMID 33276524 / 10.3390/jcm9123912 | 化療後の腫瘍微小環境における骨髄由来細胞の二面性の総説。 | 骨髄由来細胞動員の整理に。 |

## 4. 閉経・エストロゲン欠乏と骨転移

| # | 文献 | PMID / DOI | 要点 | 使い所 |
|---|---|---|---|---|
| 4-1 | Ottewell PD et al. *Clin Cancer Res* 2014 | PMID 24687923 / 10.1158/1078-0432.CCR-13-1246 | OVX（閉経後モデル）で骨吸収が上がり、骨内の DTC が増殖して腫瘍を形成（OVX 83% vs sham 17%）。骨外の腫瘍には影響なし。ZOL は OVX 群の腫瘍形成を完全に阻止、sham 群では無効。 | **論文の中核 citation**。「閉経後の骨でのみ」「骨だけ」が揃っており、本研究の臨床所見と鏡像。 |
| 4-2 | Canuas-Landero VG, Ottewell PD et al. *Front Endocrinol* 2021 | PMID 34733240 / 10.3389/fendo.2021.749428 | エストラジオールを pre/peri/post 相当濃度に作り分けると、ZOL の骨転移抑制は低エストラジオール下でのみ得られる。軟部転移では効果が逆方向にも出る。 | 「閉経状態が骨 niche の薬剤応答を決める」の定量的裏付け。 |
| 4-3 | George CN, Ottewell PD et al. *J Bone Oncol* 2020 | PMID 32995253 / 10.1016/j.jbo.2020.100317 | エストロゲンと ZOL が骨・免疫環境をどう変えるかの総説。閉経前後で抗腫瘍効果が違う理由の仮説を整理。 | Discussion: 閉経依存性の機序候補を列挙するときの総説。 |
| 4-4 | Holen I et al. *Clin Exp Metastasis* 2016 | PMID 26585891 / 10.1007/s10585-015-9770-x | ER+ 乳癌の骨転移モデルで、骨微小環境を標的にすると転移が抑制される。 | ER+ 文脈（本論文は ER+ 限定）での骨 niche 介入の妥当性。 |
| 4-5 | Cenci S, Pacifici R et al. *J Clin Invest* 2000 | PMID 11086024 / 10.1172/JCI11066 | エストロゲン欠乏は T 細胞の TNF-α 産生を増やして骨減少を起こす（免疫介在性）。 | **1st hit = 骨の炎症**の原典。 |
| 4-6 | Roggia C et al. *PNAS* 2001 | PMID 11717453 / 10.1073/pnas.251534698 | 骨髄の TNF 産生 T 細胞の増加がエストロゲン欠乏性骨減少の鍵。 | 4-5 と対。骨髄が炎症の場であること。 |
| 4-7 | Eghbali-Fatourechi G et al. *J Clin Invest* 2003 | PMID 12697741 / 10.1172/JCI17215 | 閉経直後の女性で骨髄細胞表面 RANKL が上昇し、骨吸収上昇と相関。エストロゲン/ラロキシフェンで低下。 | **ヒトでの 1st hit の直接証拠**（RANKL 上昇）。 |
| 4-8 | Yu M et al. *J Clin Invest* 2021 | PMID 33586672 / 10.1172/JCI143137 | OVX による骨減少は腸管由来の TNF+ T・Th17 細胞の骨髄移行を要する。 | 1st hit の炎症を最新の機序で補強。 |
| 4-9 | Tomkinson A et al. *J Clin Endocrinol Metab* 1997 | PMID 9284757 / 10.1210/jcem.82.9.4200 | ヒト骨で、エストロゲン離脱（GnRH 作動薬）により osteocyte アポトーシスが増える。 | **1st hit も骨細胞アポトーシスを起こす** → 2 つの hit が同じ細胞に収束する論拠。 |
| 4-10 | Emerton KB et al. *Bone* 2010 | PMID 19925896 / 10.1016/j.bone.2009.11.006 | マウス OVX 後に osteocyte アポトーシスが骨吸収に先行して起こり、空の lacunae と吸収活性化が続く。 | 4-9 のマウス版。Fig4A の読みの足場。 |
| 4-11 | Allocca G, Ottewell PD et al. *J Bone Oncol* 2019 | PMID 31236323 / 10.1016/j.jbo.2019.100244 | 乳癌骨転移 niche と HSC niche の in vivo での重複を検討。 | niche 競合の議論に。 |

## 5. ステロイド（支持療法のデキサメタゾン等）と骨細胞・転移

| # | 文献 | PMID / DOI | 要点 | 使い所 |
|---|---|---|---|---|
| 5-1 | Weinstein RS et al. *J Clin Invest* 1998 | PMID 9664068 / 10.1172/JCI2799 | グルココルチコイドが骨芽細胞形成を抑制し、骨芽細胞・**骨細胞のアポトーシス**を誘導する（GIOP の中心機序）。 | **2nd hit のステロイド側原典**。Fig4C の解釈の柱。 |
| 5-2 | Weinstein RS et al. *J Clin Endocrinol Metab* 2000 | PMID 10946902 / 10.1210/jcem.85.8.6714 | ヒト大腿骨頭で、ステロイド性骨壊死に osteocyte アポトーシスが伴う。 | ヒトでも同じことが起きる（臨床的妥当性）。 |
| 5-3 | Weinstein RS et al. *Endocrinology* 2011 | PMID 21771887 / 10.1210/en.2011-0170 | OPG（RANKL 阻害）投与がグルココルチコイド誘発 osteocyte アポトーシスを予防する。 | **RANKL 阻害がステロイド性の骨細胞障害も止める** → 「化療前にデノスマブ」の議論に直結。 |
| 5-4 | Jia J et al. *FASEB J* 2011 | PMID 21705669 / 10.1096/fj.11-182519 | グルココルチコイドの用量が osteocyte の運命（オートファジー生存 vs アポトーシス）を決める。 | 用量依存性（支持療法の反復低用量 vs 高用量）の議論に。 |
| 5-5 | Xia X et al. *J Bone Miner Res* 2010 | PMID 20564240 / 10.1002/jbmr.160 | グルココルチコイドが骨細胞にオートファジーを誘導する。 | 5-4 と対。Fig4C の DEG にオートファジー系が出たときの解釈枠。 |
| 5-6 | Obradović MMS et al. *Nature* 2019 | PMID 30867597 / 10.1038/s41586-019-1019-4 | 乳癌進行に伴うストレスホルモン上昇が遠隔転移巣でグルココルチコイド受容体を活性化し、コロニー形成を増やし生存を短縮（ROR1 経由）。**臨床でのステロイド使用に注意を促す**。 | **Discussion の重要 citation**: 支持療法ステロイドが転移側にも働き得る。 |
| 5-7 | Weinstein RS. *Endocrinol Metab Clin North Am* 2012 | PMID 22877431 / 10.1016/j.ecl.2012.04.004 | GIOP と骨壊死の総説。 | 総説枠。 |
| 5-8 | Cheng CH et al. *Int J Mol Sci* 2022 | PMID 35163300 / 10.3390/ijms23031376 | エストロゲン欠乏とグルココルチコイド過剰が骨代謝回転に与える影響を**並べて**整理した総説。 | **1st hit と 2nd hit を同じ枠で語れる**既存総説（2-hit の下敷き）。 |

## 6. 臨床: 周術期化療・骨減少の疫学と BMA の adjuvant 効果

| # | 文献 | PMID / DOI | 要点 | 使い所 |
|---|---|---|---|---|
| 6-1 | EBCTCG. *Lancet* 2015 | PMID 26211824 / 10.1016/S0140-6736(15)60908-4 | 18,766 例の IPD メタ解析。全体では骨再発 RR 0.83 のみ明確。**閉経後 11,767 例では**再発 0.86、遠隔再発 0.82、**骨再発 0.72**、乳癌死 0.82 と有意。閉経前では効果なし。 | **最重要の臨床 citation**。「骨 niche は閉経後でのみ薬理学的に抑えられる」= 本論文の閉経後特異性と一致。 |
| 6-2 | Coleman R et al. *Lancet Oncol* 2014（AZURE 最終） | PMID 25035292 / 10.1016/S1470-2045(14)70302-X | 全体では DFS 改善なし（HR 0.94）。骨転移は減少（初発事象 HR 0.78）。**閉経後 5 年超の 1,041 例では IDFS HR 0.77**、それ以外では 1.03。 | 閉経状態による効果差の単一試験での証拠。ONJ 1.7% も記載。 |
| 6-3 | Gnant M et al. *Lancet* 2015（ABCSG-18） | PMID 26040499 / 10.1016/S0140-6736(15)60995-3 | 閉経後 HR+ 早期乳癌で AI 併用のデノスマブ 60mg/6M が臨床骨折を有意に減らす。 | 「閉経後では安全に RANKL を抑えられる」実務的根拠。 |
| 6-4 | Gnant M et al. *Lancet Oncol* 2019（ABCSG-18 DFS） | PMID 30795951 / 10.1016/S1470-2045(18)30862-3 | 中央値 73 か月で DFS が改善（HR 0.82, 95%CI 0.69-0.98, p=0.026。ただし多重性未補正の記述的解析）。 | **「化療前にデノスマブ」を支持する側**の最良データ（ただし postmenopausal・AI 併用集団）。 |
| 6-5 | Gnant M et al. *NEJM Evid* 2022 | PMID 38319865 / 10.1056/EVIDoa2200162 | ABCSG-18 の長期成績。 | 6-4 の追補。 |
| 6-6 | Coleman R et al. *Lancet Oncol* 2020（D-CARE） | PMID 31806543 / 10.1016/S1470-2045(19)30687-4 | 高リスク早期乳癌 4,509 例、標準（neo)adjuvant 治療＋デノスマブ。**骨転移無発生生存は改善せず（HR 0.97, 0.82-1.14）**。 | **最大の反証**。「化療と同時にデノスマブ」を薄い根拠で推すと必ず突かれる。 |
| 6-7 | Coleman R et al. *Adv Ther* 2021（D-CARE 骨健康） | PMID 34185259 / 10.1007/s12325-021-01812-9 | D-CARE で骨関連アウトカム（骨折・BMD）は改善した（腫瘍アウトカムは不変）。 | 「骨は守れるが転移は止まらなかった」という切り分け。 |
| 6-8 | Gnant M et al. *Ann Oncol* 2015（ABCSG-12） | PMID 25403582 / 10.1093/annonc/mdu544 | 閉経前＋卵巣機能抑制下では ZOL 併用で DFS 改善傾向（薬理的閉経では骨 niche が「閉経後」的に振る舞う）。 | 「閉経＝ホルモン環境の問題」であり年齢ではないことの臨床的傍証。 |
| 6-9 | Hadji P et al. *Ann Oncol* 2016 | PMID 26681681 / 10.1093/annonc/mdv617 | 欧州パネルの adjuvant BP コンセンサス（閉経後に推奨）。 | 現行の推奨を 1 文で引く。 |
| 6-10 | Shapiro CL. *Clin Breast Cancer* 2021 | PMID 34045175 / 10.1016/j.clbc.2021.04.009 | 乳癌の BMA 総説（適応・用量・期間）。 | 実務的推奨の citation。 |
| 6-11 | Shapiro CL et al. *Eur J Cancer* 2011（CALGB 79809） | PMID 21324674 / 10.1016/j.ejca.2010.11.024 | 化療誘発卵巣機能不全の閉経前女性で ZOL が BMD を保持。 | 「化療期の骨保護は実行可能」の RCT。 |
| 6-12 | Lin Q et al. *J Clin Invest* 2026 | PMID 42138086 / 10.1172/JCI199285 | ER+ 乳癌骨転移で約 40% がデノスマブ下でも骨病変が進行。CRKL/circCCDC50/NFATc1 による **RANKL 非依存**の破骨細胞活性化が耐性機序。 | **デノスマブ万能論への歯止め**（Discussion の limitation / 将来展望）。 |
| 6-13 | Pan H et al. *NEJM* 2017 | PMID 29117498 / 10.1056/NEJMoa1701830 | ER+ 乳癌は内分泌療法 5 年終了後も 20 年にわたり再発し続ける（骨が主な部位）。 | 「周術期の niche 変化が数年〜十数年後に表現される」時間スケールの正当化。 |
| 6-14 | Fornetti J et al. *J Bone Miner Res* 2018 | PMID 30476357 / 10.1002/jbmr.3618 | 骨転移の生物学・臨床の総説（vicious cycle、骨代謝回転と転移）。 | Intro の背景総説。 |

## 7. 「2-hit」的な骨転移仮説の先行例

厳密に「骨転移の 2-hit モデル」を名乗る文献は、今回の検索では**見つからなかった**
（検索語: two-hit model bone metastasis niche priming / sequential insults 等 → ヒット 0）。
最も近い先行例は以下の 3 本で、いずれも「hit は 1 つ」の枠で書かれている。

- **Luo 2016 (PMID 26725121)** — 宿主ストローマ（老化骨芽細胞）の変化**単独**で骨コロニー形成が増える。
  ＝ 本論文の 2nd hit だけを実験的に切り出した形。
- **Ottewell 2014 / 2015 (PMID 24687923 / 25603921)** — OVX（＝1st hit）**単独**で dormant DTC が
  破骨細胞依存に覚醒する。
- **Cheng 2022 (PMID 35163300)** — エストロゲン欠乏とグルココルチコイド過剰を同じ枠で比較した総説
  （骨代謝の話に限られ、転移には触れていない）。

→ **「1st hit（閉経）× 2nd hit（化療＋ステロイド）の重ね合わせ」として明示的に定式化した例は無い**
ので、本論文の novelty 主張はここに置ける。ただし「2-hit」という語は Knudson の腫瘍抑制遺伝子仮説を
想起させるので、`two-hit` ではなく **"sequential two-insult model" / "priming-and-trigger model"** の
ような言い方を併記するのが安全。

---

## 支持する証拠 / 反する証拠

### 支持（strong）
1. **閉経後の骨でのみ DTC が覚醒し、抗吸収薬が効く** — Ottewell 2014/2015（前臨床）と EBCTCG 2015・
   AZURE（臨床）が独立に同じ方向。本論文の「閉経後でのみ HR 1.72」と鏡像をなす。
2. **骨特異性** — OVX は骨内腫瘍だけを増やし、骨外には影響しない（PMID 24687923）。本論文の
   「骨転移のみ HR 2.05、肺・肝は無関連」と一致。臓器選択性が仮説の最大の強み。
3. **化療は骨の細胞を直接傷める** — Dox で osteocyte アポトーシス（PMID 39147353）、
   エストロゲン非依存の老化/SASP による骨減少（PMID 31932453）、雌でより顕著（PMID 42686970）。
4. **1st hit と 2nd hit が同じ細胞（osteocyte）に収束する** — エストロゲン離脱でも
   （PMID 9284757, 19925896）、ステロイドでも（PMID 9664068）骨細胞アポトーシスが起き、
   骨細胞は RANKL の主要供給源（PMID 21909105, 21909103）でアポトーシス周囲は RANKL を出す（PMID 22342796）。
5. **宿主側の変化だけで niche は整う** — 老化骨芽細胞のみで骨コロニー形成増（PMID 26725121）。
6. **化療の宿主効果は薬剤消失後も残り、dormancy を破る** — PMID 21975929、40614736、42696464、27979832。
7. **時間スケールが合う** — ER+ は 20 年再発し続ける（PMID 29117498）、骨髄 DTC は早期から存在
   （PMID 16120859）。周術期の数か月の niche 変化が後年の骨転移として現れる筋は通る。

### 反する証拠・弱点
1. **D-CARE 陰性（PMID 31806543）** — 化療と同時のデノスマブは骨転移無発生生存を改善しなかった
   （HR 0.97）。「化療前/同時に RANKL を抑えれば 2nd hit を止められる」の最大の反証。
   骨健康は改善している（PMID 34185259）ので、「骨量を守ること ≠ niche を無効化すること」。
2. **老化細胞仮説の再現性** — senolytic で化療性骨減少を防げず、後期には老化マーカーも検出されない
   （PMID 40389522）。SASP を 2nd hit の唯一の機序に置くのは危険。
3. **骨での直接証拠が無い** — 「化療が骨転移を増やす」ことを示した前臨床モデルは検索で見つからず
   （肺・肝では PMID 40614736, 42696464 があるが骨は空白）。乳癌骨転移モデル × 化療、
   × 閉経の RNA-seq も存在しない（strategy メモの記載と一致）。これは novelty でもあり limitation でもある。
4. **データセットのミスマッチ** — Fig4 の 3 データセットは腫瘍を持たない正常マウス、系統・性・年齢・
   投与法がばらばら（B は雄単回 Dox、C は DEX 5 日、A は雌 OVX）。ヒトの周術期レジメン
   （AC-T ＋ 制吐のデキサメタゾン＋ AI）と直接対応しない。
5. **交絡（臨床側）** — 化療を受ける患者は元々 high risk（T/N/grade/Ki-67）。多変量で調整しても
   残余交絡と適応バイアスは残る。化療後の内分泌療法・AI 曝露、骨折・骨減少、BMD、骨スキャンの
   追跡頻度（surveillance bias: 骨痛のある人がよく撮られる）も骨転移の検出に効く。
6. **RANKL 非依存の経路** — デノスマブでも 40% が進行（PMID 42138086）。RANKL 阻害を
   「唯一打てる手」と書くなら、効かない経路があることも書く。
7. **ステロイドの寄与を分離できない** — 臨床データでは制吐ステロイドの累積量が取れていない
   （おそらく）。Fig4C を根拠に「ステロイドも 2nd hit」と言うのは機序仮説にとどまる。

---

## Introduction の段落構成案（英語箇条書き）

- **P1 — Bone is the dominant and earliest metastatic site in ER+ breast cancer.**
  Bone-first relapse dominates ER+ disease and risk persists for ≥20 years after endocrine therapy
  (PMID 29117498); DTCs are already present in marrow at diagnosis and are prognostic (PMID 16120859);
  the bone/osteogenic niche actively supports early colonization (PMID 25600338, 30476357).
- **P2 — Metastatic outgrowth in bone is a niche-controlled, not a cell-autonomous, event.**
  Pre-metastatic niche concept (PMID 16341007, 28303905); dormancy is maintained/broken by the niche
  (PMID 23728425); osteoclastic remodeling reactivates dormant cells (PMID 26632274);
  osteocytes are the central regulators and the main source of RANKL (PMID 21909105, 21909103,
  33162934, 34773212, 42426271).
- **P3 — Perioperative chemotherapy has well-documented pro-metastatic host effects.**
  TMEM/MENA^INV after neoadjuvant chemo (PMID 28679654); chemo-elicited EVs build a lung
  pre-metastatic niche (PMID 30598531); host ATF3 and TLR4 inflammation (PMID 28784776, 25274031);
  host remains permissive after drug clearance (PMID 21975929); chemo awakens dormant DTCs via NETs
  (PMID 40614736); systematized as the pro-tumorigenic host response (PMID 31645711, 37999904).
  **But: all of this has been studied in lung and liver — the bone compartment is unexamined.**
- **P4 — Chemotherapy and its supportive steroids specifically damage bone-resident cells.**
  CTIBL (PMID 26029361, 33948427, 35477650); estrogen-independent senescence/SASP-driven bone loss
  (PMID 31932453), with a female/age bias (PMID 42686970); Dox induces osteocyte apoptosis
  (PMID 39147353); glucocorticoids induce osteoblast/osteocyte apoptosis (PMID 9664068, 10946902).
- **P5 — Menopause itself primes the bone compartment.**
  Estrogen loss raises marrow RANKL in women (PMID 12697741), drives TNF+ T/Th17-mediated bone loss
  (PMID 11086024, 11717453, 33586672), and induces osteocyte apoptosis (PMID 9284757, 19925896);
  OVX alone awakens dormant DTCs in bone only, and only there is ZOL/OPG-Fc effective
  (PMID 24687923, 25603921, 34733240) — mirrored clinically by the menopause-restricted benefit of
  adjuvant bisphosphonates (PMID 26211824, 25035292).
- **P6 — Gap and aim.**
  No study has asked whether perioperative chemotherapy changes the *bone* metastatic risk, nor whether
  that effect depends on menopausal status. We report a 785-patient ER+ early breast cancer cohort
  without BMA exposure, show a bone-restricted association (HR 2.05) that is confined to postmenopausal
  women (HR 1.72), and then interrogate public bone transcriptomes of the three relevant exposures
  (estrogen loss, doxorubicin, dexamethasone) to propose a sequential two-insult model of bone niche priming.

## Discussion の段落構成案（英語箇条書き）

- **D1 — Summary of findings**, with the organ specificity as the key internal control: chemotherapy
  associated with bone but not lung/liver metastasis; interaction with menopausal status.
- **D2 — Why bone, and why only after menopause: the first insult.**
  Estrogen loss creates an inflamed, high-turnover, RANKL-rich marrow (PMID 12697741, 11086024, 33586672)
  with osteocyte apoptosis (PMID 9284757, 19925896); in mice this alone converts dormant DTCs into
  overt bone metastases (PMID 24687923, 25603921). Our postmenopausal-restricted HR is the clinical
  counterpart.
- **D3 — The second insult: cytotoxics plus steroids acting on osteocytes.**
  Our Fig4B/4C transcriptomes (describe actual DEGs/pathways; report negatives honestly);
  anchor in PMID 31932453 (senescence/SASP, estrogen-independent), 39147353 (osteocyte apoptosis),
  42686970 (female bias), 9664068/21771887 (steroid-induced osteocyte apoptosis, preventable by OPG).
- **D4 — Convergence (Fig5): both insults converge on the osteocyte → RANKL → osteoclast axis,**
  which is the one node where dormant DTCs are known to be released (PMID 26632274, 25603921) and the
  one node that is druggable in the clinic (PMID 33875860). Nearest precedent for a host-only,
  tumor-independent niche is senescent-osteoblast-driven colonization (PMID 26725121).
- **D5 — Clinical implication, stated cautiously.**
  Adjuvant anti-resorptives reduce bone recurrence only in postmenopausal women
  (PMID 26211824, 25035292, 25403582), and adjuvant denosumab improved DFS in ABCSG-18
  (PMID 30795951, descriptive analysis). **However D-CARE was negative (PMID 31806543)** while still
  improving bone health (PMID 34185259). So our data motivate a *hypothesis* — that the window of
  vulnerability is the perioperative chemotherapy period in postmenopausal women — not a treatment
  recommendation. Any trial would have to differ from D-CARE in population (postmenopausal ER+ only),
  timing (start before/with the first cycle), and endpoint (bone-first recurrence).
- **D6 — Steroids as a modifiable co-factor**: glucocorticoid receptor activation at metastatic sites
  worsens outcome (PMID 30867597); steroid-induced osteocyte apoptosis is preventable by RANKL blockade
  (PMID 21771887). Steroid-sparing antiemetics are a testable, low-cost intervention.
- **D7 — Limitations**: retrospective single-cohort design, residual confounding/indication bias,
  no steroid dose data, surveillance bias for bone imaging; public datasets are tumor-free mice with
  heterogeneous sex/strain/dosing and no chemo × bone-metastasis or menopause × bone-metastasis dataset
  exists; senescence may not be the operative mechanism (PMID 40389522); RANKL-independent
  osteoclastogenesis exists (PMID 42138086).

---

## 査読者に突かれそうな点と応答案

| # | 想定される指摘 | 応答案 |
|---|---|---|
| R1 | 化療を受けたのは元々高リスク患者で、交絡（適応バイアス）では。 | 多変量で T/N/grade/Ki-67/内分泌療法を調整済み。決定的なのは**臓器特異性**: 同じ交絡なら肺・肝の転移も増えるはずだが HR は無関連。「高リスクだから」では骨だけ増える理由を説明できない。 |
| R2 | 骨転移の検出バイアス（骨痛で撮る、化療後は受診が多い）。 | 画像のフォロー間隔・初発転移部位の定義を Methods に明記し、可能なら診断までの画像回数で感度分析。first-site 解析と any-site 解析を併記。 |
| R3 | D-CARE が陰性なのに「化療前にデノスマブ」と言うのは矛盾。 | 矛盾しない。D-CARE は閉経前後混在の高リスク集団で、エンドポイントは骨転移無発生生存の複合。EBCTCG と AZURE が示すのは**閉経後に限って**骨再発が減るということ。我々の主張は「仮説：対象を閉経後 ER+ に絞り、投与を化療開始前からにした試験が必要」であり、推奨ではない。D-CARE で骨健康は改善した（骨量保護と niche 無効化は別）という切り分けも書く。 |
| R4 | マウスの正常骨の transcriptome から転移の話をするのは飛躍。 | 同意し、Fig4/5 は「仮説生成 (hypothesis-generating)」と明記。化療 × 骨転移モデル、閉経 × 骨転移モデルの公開データは存在しないことを Limitation に書き、必要な検証実験（OVX + AC-T + DTC トレーシング）を Future directions に提案する。 |
| R5 | 老化/SASP が機序というが senolytic で防げない報告がある。 | PMID 40389522 を自分から引き、「老化は 2nd hit の唯一の経路ではない」「アポトーシス（PMID 39147353, 9664068）とそれに続く RANKL 上昇が共通の下流である」という書き方にする。我々のデータも unbiased に出し、老化マーカーが出なければ出ないと書く。 |
| R6 | 閉経後の患者は AI を使っており、骨への影響は AI によるのでは。 | 内分泌療法の種類（TAM vs AI）を共変量に入れ、AI 曝露で層別した感度分析を提示。AI だけでは「化療あり群でのみ骨転移が増える」理由にならない点も述べる。 |
| R7 | 制吐ステロイドの寄与を主張する根拠が弱い。 | 臨床データに累積ステロイド量が無いことを明記し、ステロイドは「仮説上の co-factor」とトーンを落とす。根拠は機序側（PMID 9664068, 21771887, 30867597）に限る。 |
| R8 | 「2-hit」は Knudson の用語で誤解を招く。 | 本文では "sequential two-insult (priming-and-trigger) model" と定義してから、以降略語を使う。Fig5 のラベルも合わせる。 |
| R9 | 化療誘発卵巣機能不全で閉経前患者も「閉経後」になるはず。なぜ閉経前で効果が無い？ | 重要な指摘。化療性卵巣機能不全は急速な骨減少を起こす（PMID 11454877）一方、ABCSG-12 では卵巣抑制下の閉経前患者で ZOL が効いた（PMID 25403582）。閉経前群では**曝露時点では**エストロゲンが保たれており 1st hit が先行していないこと、サンプル数と追跡期間の制約を述べ、化療性無月経の情報が取れる症例での探索的解析を提示する。 |
| R10 | ER+ だけで BMA 未使用の 785 例という設定は選択的ではないか。 | むしろ強み: BMA は骨転移リスクを直接修飾するので、未使用コホートは交絡のない natural history を見られる。ER+ 限定は骨優位の転移パターンを均質化する。この点を Methods/Strengths に明示。 |

---

## 未確認／見つからなかったもの（本文では使わない）

今回の検索で PMID を確認できなかった、または該当文献が見つからなかったもの。引用するなら要再検索。

- Hiratsuka らの「肺内皮での MMP9 誘導と premetastatic phase」— **未確認**（検索語を変えて再試行要）。
- Price ら「E-selectin / SDF-1 と骨髄血管周囲 niche での乳癌 dormancy」— **未確認**。
- Kennedy OD らの 2014 年の「アポトーシス骨細胞の近傍が RANKL を出す」単独論文 — **未確認**
  （代わりに同グループの 2012 Bone, PMID 22342796 を使用）。
- NLRP3 インフラマソームと乳癌骨転移 — **該当なし**（骨転移に絞ると見つからなかった）。
- 補体（C5a 等）と**骨の** pre-metastatic niche — **該当なし**（肺の総説 PMID 30061895 のみ）。
- 化療による骨髄 CAR 細胞（CXCL12-abundant reticular cell）の障害 — **該当なし**。
  CXCL12 軸自体は PMID 21436587 で引ける。
- シクロホスファミド／タキサン単剤による骨芽細胞・骨細胞アポトーシスの個別論文 — **該当なし**
  （複合レジメンの PMID 40389522、Dox の PMID 39147353 で代替）。
- GSE205237（Dox 単回 30 日後 皮質骨）の原著論文 — **未特定**。GEO の citation 欄を直接確認すること。
- S100A8/A9 と**骨**髄の pre-metastatic niche — **骨では該当なし**（イメージングの PMID 28744322 は
  一般的な PMN）。本論文で S100a8/9 が DEG に出た場合は、骨での先行文献が無いと書くのが正確。
- 「骨転移の two-hit モデル」を明示した先行文献 — **該当なし**（上記 7 節参照）。

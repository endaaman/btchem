# 要約モジュール図（下書き, 2026-10-13）

図: `out/prelim/summary/summary_modules.{png,pdf}`（170 × 80 mm）。値: `out/prelim/summary/summary_modules_values.csv`
（図から外した候補 4 行も `shown=False` で入れてある）。スクリプト: `prelim/summary/summary_modules.py`。

表現: FDR < 0.05 だけ色付き（赤 ↑ / 青 ↓、濃さ = NES、数値は有意セルだけ小さく）。FDR 0.05–0.25 は薄い枠線（枠の色で向き）。FDR ≥ 0.25 は空白。

## legend の最初の文（案）

> Senescence-associated programs (SenMayo, p53) are elevated in postmenopausal and older women's bone and in human osteoblasts exposed to a topoisomerase II inhibitor, but not in mouse bone or bone-marrow stroma after doxorubicin in vivo, whereas the bone-formation/ECM/angiogenesis axis moves in opposite directions between human biopsies (up) and mouse aged or doxorubicin-exposed stroma (down).

（一文にすると長い＝主張が一本にまとまっていない証拠。下の自己評価を参照。）

## 行（8）の選び方

依頼の 10 候補から 2 つ外した。外した行の値も csv にある。

| 行 | set | 選んだ / 外した理由 |
|---|---|---|
| Inflammation | Hallmark inflammatory response | 1st hit の主軸（加齢で ↑、Dox 骨で ↓）。 |
| IFN-α response | Hallmark IFN-α | 依頼どおり。ただし **IFN-γ に替えると化療 3 列すべてが有意 ↑ になる**（MSPC +1.75 / Dox 骨 +1.91 / NHO +2.25。IFN-α は MSPC で +1.42, FDR 0.055）。化療に共通する唯一の行として見せるなら γ の方が良い（下の提案）。 |
| Senescence (SenMayo) | SAUL_SEN_MAYO | 清水仮説の「老化」。Reactome SASP / cellular senescence は中身がヒストン・増殖系なので使わない（`prelim_invitro.md`）。 |
| p53 pathway | Hallmark p53 | 「アポトーシス」の代表。Hallmark apoptosis とほぼ同じ形（閉経後 +2.36、NHO +2.17）だが、apoptosis は雌 CAR で +1.68（FDR 0.028）が出る以外は差が無く、TOP2 阻害の直接の応答として p53 を残した。 |
| Osteoclast support | GO positive reg. of osteoclast differentiation | 2-hit の出口（破骨）。ほぼ全部空白だが、「出ない」こと自体を見せるために残した。 |
| Bone formation | GO ossification | niche 軸。 |
| Extracellular matrix | NABA core matrisome | niche 軸。 |
| Angiogenesis | Hallmark angiogenesis | niche 軸（血管）。 |
| ~~TNF-α / NF-κB~~ | Hallmark | Inflammation と同じ形（違いは雌 CAR の FDR 0.056 くらい）。冗長なので外した。 |
| ~~Estrogen response~~ | Hallmark estrogen early | 雄の列では意味が無い。しかも閉経後で **↑**（+2.16）と直感に反し、閉経後 vs 閉経前のバッチ交絡を疑わせる。要約図には向かない（このこと自体は Limitation で使える）。 |

## 列（7）の選び方と出どころ

| 列 | データ | 出どころ |
|---|---|---|
| Human ♀ bone biopsy, postmeno. vs premeno. | GSE141595 | **再計算**（下記） |
| Human ♀ bone biopsy, old vs young | GSE72815 | **再計算** |
| Mouse ♂ cortical bone, 30 vs 2 mo | NIHMS1776347 | `out/prelim/chemo2/nes_wide.csv` `Age_2.5y` |
| Mouse ♀ CAR cells, mid-age (12–14 mo) vs young | GSE210584 | `out/prelim/aging_f/nes_wide_all.csv` `sc210584_CAR` |
| Mouse ♀ LepR+ MSPC, Dox ×4 (8 wk) | GSE296758 | `chemo2/nes_wide.csv` `B3_Dox4x_MSPC` |
| Mouse ♂ cortical bone, Dox ×1 (30 d), 組成補正版 | GSE205237 | `chemo2/nes_wide.csv` `B_Dox_compadj` |
| Human osteoblasts, etoposide (in vitro) | GSE284666 NHO | `out/prelim/invitro/nes_wide.csv` `NHO_Eto` |

**ヒト 2 列は GSEA を回し直した。** 既存の `out/prelim/human/enrich/gsea_*.csv` は Hallmark 50 + SenMayo しか無く、
ossification / osteoclast diff. / NABA が空欄になるため。既存の DE（`out/prelim/human/de/*.csv`）の Wald 統計量を、
in vitro と同じ set 集合（`prelim/invitro/gsea_iv.py` の `collection()`：Hallmark 50 + C2/C5 17 + 手製 5、MSigDB 2025.1.Hs）・
同じ設定（1000 perm, seed 0, min 4, max 1000）で prerank した（`out/prelim/summary/gsea_human/`）。
set 集合が変わると FDR も変わるので、ヒト列は全行をこの再計算版から取った。既存版との差は小さい
（閉経後 TNF-α +2.35 → +2.33、高齢 SenMayo +1.65/0.007 → +1.65/0.011、高齢 p53 +1.41/FDR 0.050 → +1.41/0.060。
高齢 p53 はこの差で FDR < 0.05 から外れて枠線になっている）。

## 自己評価（厳しめ）

**主張は伝わらない。少なくとも「2-hit」の主張は図から読めない。** 読めるのは別の 2 つのこと。

1. **列を分けているのは「insult の種類」ではなく「種・試料」。** 左右の括弧（1st / 2nd insult）の中で色がそろっていない。
   一番目立つパターンは下 3 行（骨形成・ECM・血管）で、**ヒト 3 列が赤、マウスの加齢骨と MSPC Dox が青**。
   つまり読者が最初に受け取るのは「ヒトとマウスで niche 軸の向きが逆」で、これは仮説の支持ではなく、
   ヒト生検の組成（血管・間質の割合）か種差の問題（`prelim_human.md` の ρ = −0.30 と同じ話）。
   Dox 骨（組成補正版）も赤で、マウス内でも化療 2 列が割れている。
2. **老化（SenMayo, p53）は「ヒト ↑、マウス in vivo は出ない」。** これは図から一目で読める唯一のきれいな対比で、
   in vitro の列が「ヒト細胞では DNA 傷害で教科書どおり出る」ことを示す。ただしこれは 2-hit の図ではなく、
   「in vivo マウスで老化が出なかったのは感度ではなく時間・希釈」という限界の説明になる。
3. **Osteoclast support はほぼ白。** 2-hit の出口（破骨）の行が空なので、仮説の結末を図が支えない。
4. 列 4（雌 CAR）は 8 行すべて枠線か空白で、情報が無い。列 6 は組成補正版でも炎症 ↓・ECM ↑ と加齢の逆。
5. 見た目の改善（空白・枠線で非有意を沈める）は効いていて、40×14 版よりは「どこが有意か」は一目で分かる。
   問題はデータの側で、有意セルの並びが 1st / 2nd の括弧と対応していない。

### 次に試すなら（提案、未実施）

- **IFN-α を IFN-γ に替える**: 化療 3 列 + 加齢マウス + 閉経後がそろって赤になり、「化療と加齢に共通する唯一の軸 = IFN」が 1 行で言える（`prelim_chemo2.md` 結論 2 と同じ主張）。
- **列を「ヒト / マウス」で先に分け、その中で insult を並べる**: 種で割れる事実を隠さずに、種の中での加齢 vs 化療を比べる形にする。マウス内では「加齢と MSPC Dox が niche 軸 ↓ でそろう、炎症は加齢だけ」が読める（`prelim_niche.md` 結論 1–2）。ヒト側は化療 in vivo が無いので、in vitro 列は別枠（参考）にする。
- 行は 6 本（Inflammation / IFN-γ / SenMayo / Bone formation / ECM / Angiogenesis）まで絞れる。Osteoclast support は「出ない」ことを本文で書き、図からは外してよい。

# PXD078161 皮質骨プロテオーム: 加齢背景で DX の骨への効果は変わるか（予備解析）

解析日: 2026-10-13。スクリプト `prelim/proteome/proteome.py`（`uv run python -I prelim/proteome/proteome.py`）。出力 `out/prelim/proteome/`。
データ: PXD078161 / MassIVE MSV000101777（Izquierdo et al., *Geroscience* 2026, PMID 42686970）。`data/seq/PXD078161/` の `*_Candidates.xlsx`（Spectronaut の群間比較表）だけを使った。

## 結論（先に）

1. **このデータでは「加齢背景で DX の骨への効果が変わる」とは言えない。言えないのは「変わらない」からではなく、測定のばらつきが大きく、効果どうしを比べられないから。**
   - 共通に定量された 116 タンパク質で、DX の log2FC は若齢雌と老齢雌でまったく相関しない（Spearman ρ = −0.03, P = 0.75）。
   - ただし**同じ年齢の雌と雄（若齢雌 vs 若齢雄）でも ρ = 0.06** で、どの 2 群の組み合わせでも |ρ| ≤ 0.26（`pairwise_correlation.csv`）。DX 効果のタンパク質単位の再現性が群間でゼロなので、「老齢でだけ動く」「若齢でだけ動く」というリストは、年齢差ではなくノイズ・バッチ差の可能性を否定できない。
   - 若齢 vs 老齢の雄は ρ = −0.26（P = 0.004）と弱い逆相関だが、老齢雄は Sample_Key に記載がなく n も不明（下記）。解釈しない。
2. **何が動くか: 一貫して見えるのは「DX でコラーゲン（ECM）が下がる」方向だけ。それも年齢差を示すものではない。**
   - log2FC での prerank GSEA は 4 群とも FDR < 0.05 の set がほぼ無い（若齢雌の Hallmark MYC targets v1 ↓ NES −2.0 のみ。中身は CCT7/HSPD1/NPM1 などのシャペロン・核小体タンパク質）。
   - NABA collagens は 4 群とも NES −0.7〜−1.0、GO collagen fibril organization も −0.8〜−1.3、GO bone mineralization も −0.8〜−1.5 と**向きは 4 群でそろって ↓**（いずれも FDR ≥ 0.5）。
   - 感度分析（順位を t = log2FC / SE にする）では、若齢雌で collagen fibril org.（−2.09）・ossification（−2.02）・NABA collagens（−1.90）・osteoblast diff.（−1.87）・vasculature dev.（−1.83）が FDR < 0.05、若齢雄でも NABA collagens・collagen fibril org.・proteoglycans が FDR < 0.05。**老齢雌・老齢雄ではこれらは FDR < 0.05 にならない**。
     ただしこれは「若齢でだけコラーゲンが下がる」ではない: log2FC そのものは老齢雌でも同程度に ↓（NABA collagens NES −1.0 vs 若齢雌 −0.9）で、差は SE（＝ Spectronaut の ratio 数）から来ている。Col1a1 は若齢雌で 1,941 ratios（ペプチド × サンプル対を独立とみなした数）あり、SE が極端に小さい。**t 順位の差は統計の作り方の差で、生物学的な年齢差とは読めない**。
   - タンパク質単位で著者基準（q < 0.05 かつ |log2FC| ≥ 0.58）を満たすもの: 若齢雌 10（Col1a1 −0.64, Col11a2 −0.83, Col11a1 −0.68, Col5a1 −1.08, Aspn −0.89, Prelp +1.03 と、ケラチン 4 個・Hba-a2）、老齢雌 2（Col2a1 −1.14, Ighg2c −1.98）、若齢雄 2（Myh1 +1.86, Hbbt1）、老齢雄 9（Spp1 −1.52, Fmod −0.99, Prelp −1.02, Serpinh1 −0.82, Atp5b −0.90, Anxa2 −0.62 と Mup1/Mup14/Krt1）。**どの 2 群にも共通する有意タンパク質は無い**。
3. **niche 関連の軸のうち、炎症・老化・代謝はこのデータでは見られない。**
   - Hallmark inflammatory response / TNFα-NF-κB / IL6-JAK-STAT3 / IFN-γ / OXPHOS は 4 群とも定量タンパク質が 5 未満で検定できない（皮質骨を脱灰してグアニジン抽出した ECM 中心のプロテオームなので、細胞内・サイトカインはほぼ測れていない）。
   - SenMayo は 5–13 タンパク質（ほとんどが分泌 ECM・凝固系）で、老齢雌 +1.2 / 若齢雄 +1.0 / 老齢雄 −1.1（いずれも FDR ≥ 0.7）。老化の指標にはならない。
   - Hallmark glycolysis は 6–17 タンパク質で若齢雌 −0.8 / 老齢雌 +0.8 / 若齢雄 +1.0（FDR ≥ 0.88）。個別には Pkm（+0.53 / +0.62）と Eno1（+0.53 / +0.79）が若齢雌・老齢雌とも ↑ 方向だが、q ≥ 0.26。論文抄録の「DX で解糖系酵素が上がる」と向きは合うが、この表からは有意とは言えない。
   - 補体（H complement）は若齢雌 +0.9 / 老齢雌 +1.4 / 若齢雄 +1.0 / 老齢雄 −0.9、凝固 +0.6 / +1.0 / −0.8 / +0.5。血漿タンパク質の残り（骨髄洗い流しの残差）の影響を受けやすく、組成差と区別できない。
4. **RNA との突き合わせ（パスウェイ単位）**: 一致とも不一致とも言えない。
   - RNA 加齢（2.5y vs 2m 皮質骨 ♂）と同じ向きの割合は、老齢雌 9/22 set（41%）、若齢雌 13/22（59%）、若齢雄 20/37（54%）、老齢雄 11/17（65%）。どの群もほぼコイン投げ。
   - RNA 側で FDR < 0.05、プロテオーム側で FDR < 0.25 の set は 0–2 個しかない（`rna_concordance_summary.csv`）。
   - 唯一「向きが RNA と合う」と言えるのは、上記の**コラーゲン・骨形成・ECM の ↓**: RNA の加齢（collagen fibril −2.9, ossification −2.4）と Dox ×4 の LepR+ MSPC（−2.8, −2.1）でも ↓。ただし RNA の単回 Dox 骨（GSE205237）では逆に ↑（+1.9, +1.6）で、`docs/prelim_chemo2.md` の結論（化療は骨形成・ECM niche を削る方向だがデータで割れる）を超える情報は無い。
5. **Fig4 には載せないことを推奨**。載せるなら Supplementary で、「年齢・性別ごとに別バッチで測った皮質骨プロテオームでも、DX でコラーゲン / 骨形成系タンパク質が全群で ↓ 方向（有意ではない）」という一文の根拠として、`fig_proteome_dx_age` の B（ヒートマップ）の ECM・骨形成の行だけ。
   「加齢背景で化療の効果が変わる」の根拠には使えない（A の散布図は無相関で、しかも同年齢の雌雄でも無相関なので、年齢効果を読む土台が無い）。本文で引くなら論文の表現型（µCT・骨髄脂肪化・qPCR）を文献として引く方が強い。

## データの構造

| 群 | ファイル | 注入量 | n（Sample_Key: 対照 / DX / DQ / DX+DQ） | 定量タンパク質（DX vs 対照） | q < 0.05 | q < 0.05 & \|log2FC\| ≥ 0.58 | \|log2FC\| ≥ 0.58 |
|---|---|---|---|---|---|---|---|
| 若齢雌 3 か月 | VG1 | 500 ng、技術的重複 ×2 | 2 / 2 / 3 / 3 | 197 | 12 | 10 | 70 |
| 老齢雌 20 か月 | VG2 | 1 µg | 2 / 2 / 3 / 3 | 187 | 2 | 2 | 79 |
| 若齢雄 3 か月 | VG3 | 1 µg | 5 / 3 / 6 / 6（2 ブロックの合計） | 397 | 6 | 2 | 77 |
| 老齢雄 20 か月 | VG4 | 1 µg | **Sample_Key に記載なし**（raw には CB-74〜80 の 7 匹） | 125 | 14 | 9 | 28 |

- `All Proteins` シートに 6 比較（DX/Ctrl, DQ/Ctrl, DXDQ/Ctrl, DQ/DX, DXDQ/DQ, DXDQ/DX）が縦に積まれている。列は `AVG Log2 Ratio`, `Pvalue`, `Qvalue`, `Standard Error`, `Degrees Of Freedom`, `# of Ratios`, `Genes`, `UniProtIds` など。
- **VG2（老齢雌）と VG4（老齢雄）は比較の向きが `Ctrl / DX`** なので、`Condition Numerator / Denominator` 列で判定して符号を反転し DX/Ctrl にそろえた。反転しないと老齢群の向きが全部逆になる。
- 1 タンパク質 1 行（Genes の先頭シンボル、重複は `# of Ratios` 最大を残す。若齢雄で 2 件）。若齢雌の Valid = False は 4 行（除外していない）。
- 統計: Methods は「Storey 法、paired t-test、群ごとの多重補正、|log2FC| ≥ 0.58 かつ q ≤ 0.05」。表の `Degrees Of Freedom` = `# of Ratios` − 1 で、**ペプチド（precursor）× サンプル対を独立な観測として t 検定している**。個体 n = 2 vs 2 でも Col1a1 は df 1,940 になる。q 値は個体レベルの不確実性を反映しないので、ここでは効果量（log2FC）を主に読んだ。

### 測定条件（`Methods.docx`）

- 脛骨: 筋を除き、"diaphysis removed"（原文のまま。どの部位を残したかは本文で要確認）、骨髄を生理食塩水で洗い流す → 遠心で残りの骨髄を除去 → コラゲナーゼ I/II 5 分 → 水中超音波 3 回 → 1.2 M HCl で一晩脱灰 → 凍結粉砕 → 6 M グアニジン 72 h 抽出 → S-Trap、Lys-C / トリプシン。
- DIA（ZenoSWATH）、ZenoTOF 7600、120 分グラジエント、1 µg または 500 ng 注入。Spectronaut 16、骨専用スペクトルライブラリ（71,206 ペプチド、150 測定から作成）。
- **DX の用量・投与回数・採取時点は Methods.docx にも MassIVE にも無い**。

### 論文本文

- PMID 42686970 は PMC に入っていない（Europe PMC: `inPMC = N`, `isOpenAccess = N`, PMCID なし。E-utilities の elink も PMC へのリンク無し）。
- Springer（`link.springer.com/article/10.1007/s11357-026-02485-4`）は 3 KB のボット対策ページしか返らず、本文は取れなかった。**用量・時点は未確認**。ken がブラウザで PDF を取れば確認できる。
- 抄録（Europe PMC）の主結果: DX で海綿骨減少・骨髄脂肪化・老化 / 炎症遺伝子の上昇（雌でより強く、年齢・性別で異なる）。DQ で Alpl / Runx2 / Dmp1 が部分回復、骨髄脂肪化が減少（特に老齢雌）。皮質骨プロテオームでは「DX で ECM リモデリングと解糖系へのシフト（炎症性コラーゲンアイソフォームと解糖系酵素の増加）」、DQ で部分的に逆転。

## 方法

- タンパク質単位: 各群の DX/Ctrl log2FC を共通タンパク質で並べ、Pearson / Spearman。カテゴリは 2 通り: 著者基準（q < 0.05 かつ |log2FC| ≥ 0.58）と効果量のみ（|log2FC| ≥ 0.58、探索的）。
- GSEA: gseapy prerank、log2FC 順位、1000 perm、seed 0、min_size 5、max_size 1000。gene set は `prelim/chemo2/genesets.py` の `collection()`（Hallmark 50 + niche MSigDB 17 + 手製パネル）に、ECM 中心のプロテオーム向けとして NABA の下位分類 6 set（collagens / ECM glycoproteins / proteoglycans / ECM regulators / ECM-affiliated / secreted factors）を**事前に**足した。マウス GMT なのでシンボルを大文字化して照合。FDR はコレクション全体で計算。感度分析として t = log2FC / SE の順位でも回した（NES の Spearman は log2FC 版と 0.57–0.83）。
- RNA との突き合わせ: `out/prelim/chemo2/nes_wide.csv` の Age_2.5y（加齢 2.5y vs 2m 皮質骨 ♂）、B3_Dox4x_MSPC、B_Dox と、同じ term の NES の符号を比べた。NABA 下位分類は RNA 側で回していないので比較対象外。

## 結果の補足

### 若齢雌 vs 老齢雌（116 共通タンパク質、効果量のみの分類）

- 両方で同じ向き 8: Lox ↓、Sfrp2 ↑、C1qa / C1qc ↑、Myh1 ↑（筋混入）、Krt2 / Krt5 / Krt10 ↑（ケラチン＝取り扱い由来の混入と考えるのが自然）。
- 逆向き 8: Hapln1（若 +2.6 / 老 −1.8）、Ighg2c、Thbs2、Ngp、ケラチン 4 個。
- 老齢雌だけ 32: ↓ Col2a1 −1.14、Sparc −0.71、Ibsp −0.75、Mfge8、Apcs ／ ↑ Adipoq +1.49、Anxa5 +1.50、Spon1、Spp2、Mgp +1.05、Fmod、Calr、Eno1 +0.79、Pkm +0.62、凝固因子（F2, F10, Proc, Proz）。q < 0.05 は Col2a1 だけ。
- 若齢雌だけ 13: ↓ Col1a1、Col5a1、Col11a1、Col11a2、Aspn、Mup3 ／ ↑ Col6a1、Col6a2、Dpt、Prelp、Clec3b、Gpx3。
- 老齢雌の Adipoq ↑ は論文の「DX で骨髄脂肪化」と向きは合う（骨髄の残り）。ただし q = 0.44、Hallmark adipogenesis も +1.3（FDR 0.59）で、裏付けは弱い。
- 一覧: `compare_YF_AF.csv`（雄は `compare_YM_AM.csv`）。

### 根拠の弱い点（まとめ）

1. 群あたり n = 2（雌）。老齢雄は n 不明。若齢雌は技術的重複。
2. 4 群が別バッチ・別注入量・別 Spectronaut 解析で、**年齢はバッチと完全に交絡**。年齢間でタンパク質の絶対量は比べられず、各群内の DX/Ctrl 比どうしを並べることしかできない。
3. 定量 125–397 タンパク質で ECM・血漿タンパク質が中心。炎症・老化・代謝の軸はほぼ測れていない。
4. q 値は ratio 単位の t 検定で個体 n を反映しない。t 順位の GSEA はこの偏りをそのまま持ち込む。
5. ケラチン・Ig・ヘモグロビン・筋（Myh1/2/4, Acta1）・Mup（肝・尿由来）など、組成・混入で動くものが上位に多い。
6. DX の用量・時点が未確認。

## 出力（`out/prelim/proteome/`）

- `fig_proteome_dx_age.{png,pdf,tiff}`: 170 mm 幅。A = 若齢雌 vs 老齢雌の DX 効果の散布図、B = 年齢・性別 × パスウェイの NES ヒートマップ（右に RNA の参照列）。
- `figS_proteome_males.{png,pdf,tiff}`: 80 mm 幅。若齢雄 vs 老齢雄の散布図。
- `summary_groups.csv`（群ごとの n・定量数・有意数）、`dx_vs_ctrl_long.csv` / `dx_vs_ctrl_wide.csv`（向きをそろえた DX/Ctrl）、`pairwise_correlation.csv`、`compare_YF_AF.csv`、`compare_YM_AM.csv`
- `gsea_long.csv`、`gsea_wide_log2fc.csv`、`gsea_wide_t.csv`、`genesets_meta.csv`
- `rna_concordance_long.csv`、`rna_concordance_summary.csv`

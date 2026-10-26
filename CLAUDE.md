# btchem — 乳がん周術期化学療法と骨転移（清水と共同、MedComm 2本目）

- uv を使う。実行は `uv run python <script>`。
- `data/` → `~/Sync/Projects/Shimizu/btchem`（git 外）。公開データは `data/seq/<GSE>/`、清水からの資料は `data/` 直下。
- 出力は `out/`（git 外）。予備解析は `out/prelim/`。
- 図の規定（170/80 mm 幅、8–10 pt、パネルラベル 12 pt 太字、600 dpi、png/pdf/tiff）は `prelim/common.py`。
- 手本は `,/参考1`（1本目 osteoporosis-bone-meta-fig）。KM 等は同じ構成にする。
- 戦略・データセット評価は `,docs/strategy_20261012.md`。解析メモは `docs/`。
- 清水の指示（`data/2本目MEDCOMM_遠田作業依頼資料20260930.docx`）に従う。解析は unbiased に回し、出なかった結果もそのまま書く。
- 清水の依頼資料にある BioRender のアカウント情報は転記しない。
- 清水向け資料（`docs/report/`）は清水の目線で書く。版の名前（v1/v2/v3）、資料の作り方（並べ替えた・落とさずに）、経緯（以前は〜・訂正）は書かず、
  現在の理解だけを書く。清水の依頼資料の文面を原文で引用しない（突き返しに読める）。清水が見たのは 10/12 の資料（v1）までで、v2（10/16）は見せていない。
- 出力の世代管理: 図・デザイン・資料の出力先は `out/<OUT_VER>/{fig,design,report}`。世代は repo 直下の `OUT_VER`（例 `v2`）で指定し、
  `fig/common.py`・`design/render.sh`・`docs/report/build.sh` がそれを読む。一区切りで git tag `vN` を打ち、次の世代に進むときは `OUT_VER` を上げる。
  `out/` 直下にはバージョンの無いもの（`prelim/`）だけ置く。`data/` はそのまま（git 外、上書きしない）。
  v1 = 2026-10-12 の第一波（`out/v1`、tag `v1`）。v2 = 2026-10-13〜 の第二波（tag `v2`）。
  v3 = 2026-10-26〜 清水の 10/14 依頼資料（`data/20261014_遠田先生解析依頼資料.docx`、`out/v3` にも写し）を背骨にした組み直し（Fig4A 案1–4 / 4B / 4C 案1–2）。

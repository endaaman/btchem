# btchem — 乳がん周術期化学療法と骨転移（清水と共同、MedComm 2本目）

- uv を使う。実行は `uv run python <script>`。
- `data/` → `~/Sync/Projects/Shimizu/btchem`（git 外）。公開データは `data/seq/<GSE>/`、清水からの資料は `data/` 直下。
- 出力は `out/`（git 外）。予備解析は `out/prelim/`。
- 図の規定（170/80 mm 幅、8–10 pt、パネルラベル 12 pt 太字、600 dpi、png/pdf/tiff）は `prelim/common.py`。
- 手本は `,/参考1`（1本目 osteoporosis-bone-meta-fig）。KM 等は同じ構成にする。
- 戦略・データセット評価は `,docs/strategy_20261012.md`。解析メモは `docs/`。
- 清水の指示（`data/2本目MEDCOMM_遠田作業依頼資料20260930.docx`）に従う。解析は unbiased に回し、出なかった結果もそのまま書く。
- 清水の依頼資料にある BioRender のアカウント情報は転記しない。
- 出力の世代管理: 一区切りごとに `out/{fig,design,report}` を `out/vN/` にコピーして保存し、同じ名前の git tag（`vN`）を打つ。
  `git checkout vN` して `fig/all.sh`・`design/render.sh`・`docs/report/build.sh` を回せば `out/vN` を再現できる。
  `data/` はそのまま置いておく（git 外、上書きしない）。fig と report は必ずセットで撮る。v1 = 2026-10-12 の第一波。

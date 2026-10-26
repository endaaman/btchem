#!/usr/bin/env bash
# design/*.drawio を out/<OUT_VER>/design/ に png (-s 4) と pdf で書き出す。
# 使い方: design/render.sh [name ...]   （省略時は全部）
# drawio 単位: 2 px = 1 pt。-s 4 で 8 px/pt = 576 dpi 相当（-b 20 の余白込み）。
set -euo pipefail
cd "$(dirname "$0")/.."
ver=$(cat OUT_VER)
out="out/$ver/design"
mkdir -p "$out"
names=("$@")
[ ${#names[@]} -eq 0 ] && names=(fig1 fig5_revised fig5_shimizu ga)
# 出力はフォルダ指定で一時ディレクトリに出してから移す。Arch の drawio-desktop は `electron app.asar "$@"` で起動するため
# app.asar も入力に数えられ、-o にファイル名を渡すと "output must be a folder when exporting multiple files" で落ちる。
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
for n in "${names[@]}"; do
  drawio --no-sandbox --disable-gpu -x -f png -s 4 -b 20 -o "$tmp/" "design/$n.drawio" >/dev/null 2>&1 || true
  drawio --no-sandbox --disable-gpu -x -f pdf --crop -b 20 -o "$tmp/" "design/$n.drawio" >/dev/null 2>&1 || true
  mv "$tmp/$n.png" "$out/$n.png"
  mv "$tmp/$n.pdf" "$out/$n.pdf"
  echo "$out/$n.{png,pdf}"
done

#!/bin/sh
# docs/report/*.md -> out/<OUT_VER>/report/*.pdf（pandoc + lualatex + ltjsarticle、日本語、余白 15 mm）
# 実行: sh docs/report/build.sh [name ...]   （repo root から。引数なしで全部）
set -e
cd "$(dirname "$0")/../.."
ver=$(cat OUT_VER)
out="out/$ver/report"
mkdir -p "$out"
if [ $# -gt 0 ]; then files=""; for n in "$@"; do files="$files docs/report/$n.md"; done; else files=$(ls docs/report/*.md); fi
for md in $files; do
  name=$(basename "$md" .md)
  pandoc "$md" -o "$out/$name.pdf" \
    --pdf-engine=lualatex \
    --resource-path=.:docs/report:out/$ver:out \
    -V documentclass=ltjsarticle -V geometry:margin=15mm -V fontsize=10pt \
    -V mainfont='Noto Sans CJK JP' -V monofont='Noto Sans Mono CJK JP' \
    -V linkcolor=blue -V colorlinks=true \
    -H docs/report/header.tex
  echo "built $out/$name.pdf"
done

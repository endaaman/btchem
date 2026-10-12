"""docs/report/parts/*.md（図ごとの原稿）を 1 本に結合して docs/report/shimizu_20261012.md を作る。
各 part の YAML は捨て、title を章見出し（#）にし、本文の見出しを 1 段下げる。
「この資料について」節は落とし、「清水先生に確認したいこと」節は末尾に章ごとにまとめる。
実行: uv run python docs/report/merge.py && sh docs/report/build.sh
"""
import re
from pathlib import Path

PARTS = ['fig1_design', 'fig2_fig3_km', 'fig4_dry', 'fig5_ga', 'literature']
SRC = Path('docs/report/parts')
OUT = Path('docs/report/shimizu_20261012.md')
DROP = ('この資料について',)
ASK = '清水先生に確認したいこと'

HEAD = """---
title: "乳がん周術期化学療法と骨転移 — 図の作成状況と公開データ解析の報告"
date: 2026-10-12
author: 遠田
---

# 要点

- **臨床（Fig. 1–3）**: 1 本目のコホート（785 例）に周術期化療の有無が入っていたので、骨転移の KM（Fig. 3）は実データで作成済み。
  化療 (+) vs (−) の Cox 単変量 HR は全体 1.89（1.20–2.99）、閉経前 1.24（0.59–2.61, n.s.）、閉経後 2.46（1.37–4.40）。
  化療の記載が無い 43 例は除外（解析 742 例）。肺・肝（Fig. 2）は拡張データが届いてから。
- **公開データ解析（Fig. 4, S1–S3）**: マウス骨（加齢 / Dox / DEX）とヒト骨生検を unbiased に解析した。
  加齢骨で炎症＋アポトーシス、DEX で炎症＋OPG 低下（RANKL/OPG 比の上昇）が出た。
  一方、**化療（Dox）とステロイドでアポトーシス因子の上昇は、どの公開データでも転写レベルでは確認できなかった**。
  Fig. 4A の加齢データは雄マウスのみで、Limitation に書く。
- **概念図（Fig. 5, GA）**: データに合わせた修正版（1st hit = 炎症＋老化、2nd hit = RANKL/OPG シフト＋骨形成抑制）を推奨。
  "two-hit" は Knudson を想起させるので "sequential two-insult" 等を提案。
- 確認したい事項は末尾にまとめた。

"""


def split_sections(text):
    """'# ' 見出し単位に分割 → [(title, body)]。先頭の見出し前は title=None。"""
    out, title, buf = [], None, []
    for line in text.splitlines():
        m = re.match(r'^# (.+)$', line)
        if m:
            out.append((title, '\n'.join(buf)))
            title, buf = m.group(1).strip(), []
        else:
            buf.append(line)
    out.append((title, '\n'.join(buf)))
    return out


def demote(body):
    return re.sub(r'^(#+) ', lambda m: '#' + m.group(1) + ' ', body, flags=re.M)


chapters, asks = [], []
for name in PARTS:
    text = (SRC / f'{name}.md').read_text()
    m = re.match(r'^---\n(.*?)\n---\n', text, flags=re.S)
    meta = dict(re.findall(r'^(\w+):\s*"?(.*?)"?\s*$', m.group(1), flags=re.M))
    text = text[m.end():]
    secs = split_sections(text)
    chap = [f"# {meta['title']}\n"]
    for title, body in secs:
        if title is None:
            if body.strip():
                chap.append(body)
            continue
        if title.startswith(DROP):
            continue
        if title == ASK:
            asks.append((meta['title'], body.strip()))
            continue
        chap.append(f'## {title}\n' + demote(body))
    chapters.append('\n'.join(chap))

tail = ['# 清水先生に確認したいこと（まとめ）\n']
for chap_title, body in asks:
    tail.append(f'## {chap_title}\n\n{demote(body)}\n')

OUT.write_text(HEAD + '\n\n'.join(chapters) + '\n\n' + '\n'.join(tail))
print(f'wrote {OUT} ({OUT.stat().st_size} bytes)')

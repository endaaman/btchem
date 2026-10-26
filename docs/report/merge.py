"""docs/report/<parts>/*.md（章ごとの原稿）を 1 本に結合して docs/report/<out>.md を作る。
各 part の YAML は捨て、title を章見出し（#）にし、本文の見出しを 1 段下げる。
「この資料について」節は落とし、「清水先生に確認したいこと」節は末尾に章ごとにまとめる。
実行: uv run python docs/report/merge.py v1|v2 && sh docs/report/build.sh
"""
import re
import sys
from pathlib import Path

GEN = {
    'v1': dict(parts_dir='parts', parts=['fig1_design', 'fig2_fig3_km', 'fig4_dry', 'fig5_ga', 'literature'],
               out='shimizu_20261012', head='head_v1.md'),
    'v2': dict(parts_dir='parts2',
               parts=['clinical', 'datasets', 'analysis', 'not_found', 'fig_options', 'hypothesis_options'],
               out='shimizu_20261016', head='head_v2.md'),
    'v3': dict(parts_dir='parts3',
               parts=['guide', 'fig4a', 'fig4b', 'fig4c', 'clinical', 'crosscut', 'fig_options', 'fig5_ga',
                      'hypothesis', 'literature', 'not_found', 'datasets_list', 'v1_dry', 'methods'],
               out='shimizu_20261026', head='head_v3.md'),
}
DROP = ('この資料について',)
ASK = '清水先生に確認したいこと'


def split_sections(text):
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


def main(gen):
    cfg = GEN[gen]
    base = Path('docs/report')
    src = base / cfg['parts_dir']
    head = (base / cfg['head']).read_text()
    chapters, asks = [], []
    for name in cfg['parts']:
        p = src / f'{name}.md'
        if not p.exists():
            print(f'skip (missing): {p}')
            continue
        text = p.read_text()
        m = re.match(r'^---\n(.*?)\n---\n', text, flags=re.S)
        meta = dict(re.findall(r'^(\w+):\s*"?(.*?)"?\s*$', m.group(1), flags=re.M))
        text = text[m.end():]
        chap = [f"# {meta['title']}\n"]
        for title, body in split_sections(text):
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
    tail = [f'# {ASK}（まとめ）\n']
    for chap_title, body in asks:
        tail.append(f'## {chap_title}\n\n{demote(body)}\n')
    out = base / f"{cfg['out']}.md"
    out.write_text(head + '\n\n'.join(chapters) + '\n\n' + '\n'.join(tail))
    print(f'wrote {out} ({out.stat().st_size} bytes)')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'v2')

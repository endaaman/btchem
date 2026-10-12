"""ヒト contrast の DE。prelim/de.py の run()（pydeseq2・低発現フィルタ）をそのまま使う。

  uv run python -I prelim/human/de_h.py [label ...]
出力: out/prelim/human/de/<label>.csv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import de  # noqa: E402
from hcontrasts import HC, OUT  # noqa: E402


def main(labels):
    (OUT / 'de').mkdir(parents=True, exist_ok=True)
    for lb in labels:
        r = de.run(HC[lb])
        r.to_csv(OUT / 'de' / f'{lb}.csv', index=False)


if __name__ == '__main__':
    main(sys.argv[1:] or list(HC))

"""in vitro contrast の DE。prelim/de.py の run()（pydeseq2・低発現フィルタ・MLE）をそのまま使う。

  uv run python -I prelim/invitro/de_iv.py [label ...]
出力: out/prelim/invitro/de/<label>.csv
TE-Bone（TempO-Seq）も同じ関数に渡す（プローブを遺伝子に合算した counts。全遺伝子でなく 2,733 遺伝子のパネル）。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import de  # noqa: E402
from ivcontrasts import IC, OUT  # noqa: E402


def main(labels):
    (OUT / 'de').mkdir(parents=True, exist_ok=True)
    for lb in labels:
        r = de.run(IC[lb])
        r.to_csv(OUT / 'de' / f'{lb}.csv', index=False)


if __name__ == '__main__':
    main(sys.argv[1:] or list(IC))

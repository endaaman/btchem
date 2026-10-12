"""prelim/genes.py の注目パネルをヒトシンボルに変換（大文字化 + 名前が違うものだけ手当て）。

マウス → ヒトの直交オルソログは大文字化で近似する。パネル内で名前が変わるのは Fasl→FASLG, Trp53→TP53 の 2 つ。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import genes  # noqa: E402

RENAME = {'Fasl': 'FASLG', 'Trp53': 'TP53'}


def to_h(g):
    return RENAME.get(g, g.upper())


INFLAM = [to_h(g) for g in genes.INFLAM]
APOP = [to_h(g) for g in genes.APOP]
ALL_GENES = list(dict.fromkeys(INFLAM + APOP))
M2H = {g: to_h(g) for g in genes.ALL_GENES}

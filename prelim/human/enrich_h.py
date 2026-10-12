"""ヒト contrast の enrichment。prelim/enrich.py の ora() / gsea()（閾値・背景・GSEA 設定）をそのまま使い、
ライブラリをヒト版に差し替える。prerank GSEA には SenMayo（MSigDB C2 SAUL_SEN_MAYO, 2023.2.Hs）を Hallmark に足して回す。

  uv run python -I prelim/human/enrich_h.py [label ...]
前提: out/prelim/human/genesets/c2.cgp.v2023.2.Hs.symbols.gmt（data.broadinstitute.org から取得済み）
出力: out/prelim/human/enrich/ora_<label>.csv, gsea_<label>.csv, summary_top.csv
"""
import sys
from pathlib import Path

import pandas as pd
import gseapy as gp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import enrich  # noqa: E402
from hcontrasts import HC, OUT  # noqa: E402

GMT = OUT / 'genesets' / 'c2.cgp.v2023.2.Hs.symbols.gmt'
HALLMARK = 'MSigDB_Hallmark_2020'
_cache = {}


def senmayo():
    for line in open(GMT):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'SAUL_SEN_MAYO':
            return f[2:]
    raise KeyError('SAUL_SEN_MAYO')


def lib(name):
    if name not in _cache:
        g = gp.get_library(name=name, organism='Human')
        g = {k: [x.upper() for x in v] for k, v in g.items()}
        if name == HALLMARK:  # ORA / GSEA とも Hallmark 50 + SenMayo 1 の 51 セットとして扱う
            g['SenMayo (SAUL_SEN_MAYO)'] = senmayo()
        _cache[name] = g
    return _cache[name]


enrich.LIBS = ['GO_Biological_Process_2023', 'KEGG_2021_Human', HALLMARK]
enrich.DE = OUT / 'de'
enrich.OUT = OUT / 'enrich'
enrich.lib = lib


def main(labels):
    enrich.main(labels)


if __name__ == '__main__':
    main(sys.argv[1:] or list(HC))

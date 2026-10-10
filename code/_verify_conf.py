# -*- coding: utf-8 -*-
"""_verify_conf.py -- verify the four conference papers that Crossref does not
index, using DBLP and (for PMLR) the proceedings site."""
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "ref_verify", "raw")
os.makedirs(OUT, exist_ok=True)
UA = "Mozilla/5.0 (compatible; PaperA-ref-audit/1.0; mailto:8073592@qq.com)"

QUERIES = {
    "Cuturi2013": "Sinkhorn distances lightspeed computation of optimal transport",
    "Genevay2019": "Sample complexity of Sinkhorn divergences",
    "Janati2020": "Debiased Sinkhorn barycenters",
    "Feydy2019": "Interpolating between optimal transport and MMD using Sinkhorn divergences",
}


def get(url, timeout=40):
    req = urllib.request.Request(url, headers={"User-Agent": UA,
                                               "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as f:
        return f.read().decode("utf-8", "replace")


for key, q in QUERIES.items():
    print("=" * 90)
    print(key, "|", q)
    url = ("https://dblp.org/search/publ/api?q=%s&format=json&h=5"
           % urllib.parse.quote(q))
    try:
        data = get(url)
        io.open(os.path.join(OUT, key + ".dblp.json"), "w",
                encoding="utf-8").write(data)
        d = json.loads(data)
        hits = (d.get("result", {}).get("hits", {}) or {}).get("hit", []) or []
        if isinstance(hits, dict):
            hits = [hits]
        for h in hits:
            i = h.get("info", {})
            a = i.get("authors", {}).get("author", [])
            if isinstance(a, dict):
                a = [a]
            print("   T:", (i.get("title") or "")[:80])
            print("   A:", ", ".join(x.get("text", "") for x in a))
            print("   V:", i.get("volume"), "| P:", i.get("pages"),
                  "| Y:", i.get("year"), "| venue:", (i.get("venue") or "")[:50])
            print("   ee:", (i.get("ee") or "")[:70])
            print("   -")
    except Exception as e:        # noqa
        print("   DBLP ERROR:", e)
    time.sleep(1.5)

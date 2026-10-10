# -*- coding: utf-8 -*-
"""
_verify_openalex.py -- second independent bibliographic source (OpenAlex)
for every entry of references.bib.  Writes ref_verify/openalex_check.csv
"""
import io
import json
import os
import re
import csv
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.join(HERE, "paper_a", "manuscript")
OUT = os.path.join(HERE, "ref_verify")
RAW = os.path.join(OUT, "raw")
os.makedirs(RAW, exist_ok=True)
UA = "PaperA-ref-audit/1.0 (mailto:8073592@qq.com)"

import importlib.util
spec = importlib.util.spec_from_file_location("vr", os.path.join(HERE, "_verify_refs.py"))
vr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vr)

bib = vr.parse_bib(os.path.join(MS, "references.bib"))


def oa_search(title, author_last):
    q = "%s %s" % (title, author_last or "")
    url = ("https://api.openalex.org/works?search=%s&per-page=5&mailto=8073592@qq.com"
           % urllib.parse.quote(q))
    try:
        d = json.load(urllib.request.urlopen(
            urllib.request.Request(url, headers={"User-Agent": UA}), timeout=60))
    except Exception as e:        # noqa
        return None, str(e)
    want = set(vr.norm_title(title).split())
    best, bs = None, 0.0
    for w in d.get("results", []):
        t = w.get("title") or ""
        ws = set(vr.norm_title(t).split())
        s = len(want & ws) / float(max(1, len(want | ws)))
        if s > bs:
            bs, best = s, w
    return (best, bs) if best is not None and bs >= 0.5 else (None, bs)


rows = []
for e in bib:
    key = e["key"]
    bt = vr.clean_tex(e.get("title", ""))
    au = vr.split_authors(e.get("author", ""))
    last = au[0][0] if au else ""
    best, score = oa_search(bt, last)
    rec = {"key": key, "bib_title": bt, "score": score}
    if best is None:
        rec.update({"oa_title": "", "oa_authors": "", "oa_venue": "",
                    "oa_year": "", "oa_volume": "", "oa_pages": "",
                    "oa_doi": "", "oa_type": "", "oa_url": ""})
        rec["score"] = -1.0 if isinstance(score, str) else float(score)
    else:
        b = best.get("biblio") or {}
        src = ((best.get("primary_location") or {}).get("source") or {})
        pages = ""
        if b.get("first_page"):
            pages = "%s-%s" % (b.get("first_page"), b.get("last_page"))
        authors = [a["author"]["display_name"] for a in
                   (best.get("authorships") or [])]
        rec.update({
            "oa_title": best.get("title") or "",
            "oa_authors": "; ".join(authors),
            "oa_venue": src.get("display_name") or "",
            "oa_year": best.get("publication_year") or "",
            "oa_volume": b.get("volume") or "",
            "oa_pages": pages,
            "oa_doi": best.get("doi") or "",
            "oa_type": best.get("type") or "",
            "oa_url": (best.get("primary_location") or {}).get(
                "landing_page_url") or "",
        })
        io.open(os.path.join(RAW, key + ".openalex.json"), "w",
                encoding="utf-8").write(json.dumps(best, ensure_ascii=False)[:300000])
    # comparisons
    rec["title_ok"] = vr.norm_title(bt) == vr.norm_title(rec["oa_title"])
    bl = set(x[0].lower() for x in au)
    rl = set()
    for a in (rec["oa_authors"].split("; ") if rec["oa_authors"] else []):
        w = a.split()
        if w:
            rl.add(w[-1].lower().replace("í", "i").replace("é", "e")
                   .replace("ñ", "n").replace("ü", "u").replace("ö", "o"))
    bl2 = set(x.replace("í", "i").replace("é", "e").replace("ñ", "n")
              .replace("ü", "u").replace("ö", "o") for x in bl)
    rec["authors_ok"] = bool(bl2) and (bl2 <= rl or rl <= bl2)
    rec["year_ok"] = (str(e.get("year", "")) == str(rec["oa_year"]))
    pg = lambda s: re.sub(r"[^0-9]", "", (s or ""))
    rec["pages_ok"] = (not rec["oa_pages"]) or \
        (pg(rec["oa_pages"]) == pg(e.get("pages", "")) or
         pg(rec["oa_pages"]).startswith(pg(e.get("pages", ""))[:4]))
    rec["vol_ok"] = (not rec["oa_volume"]) or \
        (str(rec["oa_volume"]) == str(e.get("volume", "")))
    rows.append(rec)
    print("  %-24s score=%.2f T=%s A=%s Y=%s V=%s P=%s | %s"
          % (key, rec["score"], rec["title_ok"], rec["authors_ok"],
             rec["year_ok"], rec["vol_ok"], rec["pages_ok"],
             (rec["oa_venue"] or "-")[:38]))
    time.sleep(0.4)

cols = ["key", "score", "bib_title", "oa_title", "oa_authors", "oa_venue",
        "oa_year", "oa_volume", "oa_pages", "oa_doi", "oa_type", "oa_url",
        "title_ok", "authors_ok", "year_ok", "vol_ok", "pages_ok"]
with io.open(os.path.join(OUT, "openalex_check.csv"), "w", encoding="utf-8",
             newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow(r)
print("wrote ref_verify/openalex_check.csv")

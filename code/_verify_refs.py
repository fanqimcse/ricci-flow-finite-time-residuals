# -*- coding: utf-8 -*-
"""
_verify_refs.py -- three-way verification of the bibliography of Paper A.

Phase A  parse references.bib
Phase B  cross-check the \\cite keys used in manuscript.tex
Phase C  download the authoritative record for every entry
         (arXiv API for 10.48550 DOIs, Crossref for everything else)
         and save the raw response to ref_verify/raw/
Phase D  compare title / authors / journal / volume / pages / year
         and write ref_verify/bibliography_check.csv

Run:  python _verify_refs.py [--fetch] [--pdfs]
"""
import io
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.join(HERE, "paper_a", "manuscript")
OUT = os.path.join(HERE, "ref_verify")
RAW = os.path.join(OUT, "raw")
PDFS = os.path.join(OUT, "pdfs")
os.makedirs(RAW, exist_ok=True)
os.makedirs(PDFS, exist_ok=True)

UA = ("PaperA-reference-audit/1.0 (mailto:8073592@qq.com) "
      "python-urllib")
MAILTO = "8073592@qq.com"


def get(url, timeout=45, retries=3):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    last = None
    for k in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return f.read()
        except Exception as e:          # noqa
            last = e
            time.sleep(2 + 3 * k)
    raise RuntimeError("fetch failed: %s (%s)" % (url, last))


# --------------------------------------------------------------- Phase A ----
def parse_bib(path):
    txt = io.open(path, encoding="utf-8").read()
    # strip comment lines
    txt = re.sub(r"(?m)^\s*%.*$", "", txt)
    entries = []
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,]+)\s*,", txt):
        kind, key = m.group(1).lower(), m.group(2).strip()
        # find the matching closing brace
        i = m.end()
        depth = 1
        while i < len(txt) and depth:
            if txt[i] == "{":
                depth += 1
            elif txt[i] == "}":
                depth -= 1
            i += 1
        body = txt[m.end():i - 1]
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*", body):
            name = fm.group(1).lower()
            j = fm.end()
            if j < len(body) and body[j] == "{":
                d = 1
                k = j + 1
                while k < len(body) and d:
                    if body[k] == "{":
                        d += 1
                    elif body[k] == "}":
                        d -= 1
                    k += 1
                val = body[j + 1:k - 1]
            else:
                k = body.find(",", j)
                val = body[j:k if k > 0 else len(body)]
            fields[name] = " ".join(val.split())
        entries.append({"kind": kind, "key": key, **fields})
    return entries


def clean_tex(s):
    if s is None:
        return ""
    s = s.replace("\\&", "&")
    s = re.sub(r"\\url\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\[\"'`^~]\s*\{?([A-Za-z])\}?", r"\1", s)
    s = re.sub(r"\\([a-zA-Z]+)\s*", r"\1", s)
    s = s.replace("{", "").replace("}", "").replace("\\", "")
    return " ".join(s.split())


def split_authors(s):
    s = clean_tex(s)
    parts = re.split(r"\s+and\s+", s)
    out = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        if "," in p:
            last, first = p.split(",", 1)
            out.append((last.strip(), first.strip()))
        else:
            w = p.split()
            out.append((w[-1], " ".join(w[:-1])))
    return out


def norm_title(s):
    s = clean_tex(s).lower()
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    return " ".join(s.split())


# --------------------------------------------------------------- Phase B ----
def cite_order(path):
    txt = io.open(path, encoding="utf-8").read()
    txt = re.sub(r"(?m)^\s*%.*$", "", txt)
    order, seen = [], set()
    for m in re.finditer(r"\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]*)\}", txt):
        for k in m.group(1).split(","):
            k = k.strip()
            if k and k not in seen:
                seen.add(k)
                order.append(k)
    return order, seen


# --------------------------------------------------------------- Phase C ----
def arxiv_id(doi):
    m = re.match(r"10\.48550/arXiv\.(.+)$", doi or "")
    return m.group(1) if m else None


ARXIV_NS = "{http://www.w3.org/2005/Atom}"


def fetch_arxiv(aid, cache=None):
    url = ("http://export.arxiv.org/api/query?id_list=%s&max_results=1"
           % urllib.parse.quote(aid))
    if cache and os.path.exists(cache) and os.path.getsize(cache) > 500:
        data = io.open(cache, encoding="utf-8", errors="replace").read()
    else:
        data = get(url).decode("utf-8", "replace")
    import xml.etree.ElementTree as ET
    ns = {"a": "http://www.w3.org/2005/Atom",
          "arxiv": "http://arxiv.org/schemas/atom"}
    root = ET.fromstring(data)
    entry = root.find("a:entry", ns)
    if entry is None:
        return {"source": "arxiv-api(empty)", "url": url, "title": "",
                "authors": [], "journal": "", "year": "", "doi": "",
                "volume": "", "pages": "", "raw": data}
    def t(name):
        el = entry.find("a:" + name, ns)
        return " ".join((el.text or "").split()) if el is not None else ""
    authors = [" ".join((e.text or "").split())
               for e in entry.findall("a:author/a:name", ns)]
    jref = t("journal_ref")
    adoi = ""
    el = entry.find("arxiv:doi", ns)
    if el is not None and el.text:
        adoi = el.text.strip()
    return {
        "source": "arxiv-api",
        "url": "https://arxiv.org/abs/%s" % aid,
        "title": t("title"),
        "authors": authors,
        "journal": jref,
        "year": t("published")[:4],
        "doi": adoi,
        "volume": "", "pages": "",
        "raw": data,
    }


def fetch_crossref(doi, cache=None):
    url = "https://api.crossref.org/works/%s?mailto=%s" % (
        urllib.parse.quote(doi, safe="/"), MAILTO)
    data = None
    if cache and os.path.exists(cache) and os.path.getsize(cache) > 500:
        try:
            head = io.open(cache, encoding="utf-8", errors="replace").read(2000)
            if '"message"' in head or head.lstrip().startswith("{"):
                data = io.open(cache, encoding="utf-8", errors="replace").read()
        except Exception:      # noqa
            data = None
    if data is None:
        data = get(url).decode("utf-8", "replace")
    j = json.loads(data)["message"]
    au = []
    for a in j.get("author", []) or []:
        au.append("%s %s" % (a.get("given", ""), a.get("family", "")))
    issued = (j.get("issued", {}).get("date-parts") or [[None]])[0][0]
    return {
        "source": "crossref-doi",
        "url": url,
        "title": " ".join((j.get("title") or [""])[0].split()),
        "authors": [x.strip() for x in au],
        "journal": " ".join((j.get("container-title") or [""])[0].split()),
        "year": str(issued) if issued else "",
        "doi": j.get("DOI", ""),
        "volume": j.get("volume", "") or "",
        "pages": j.get("page", "") or "",
        "publisher": j.get("publisher", ""),
        "raw": data,
    }


def fetch_crossref_query(title, author, cache=None):
    q = urllib.parse.urlencode({
        "query.bibliographic": "%s %s" % (title, author),
        "rows": "5", "mailto": MAILTO,
        "select": "DOI,title,author,container-title,issued,volume,page,publisher,type",
    })
    url = "https://api.crossref.org/works?%s" % q
    data = None
    if cache and os.path.exists(cache) and os.path.getsize(cache) > 500:
        try:
            head = io.open(cache, encoding="utf-8", errors="replace").read(2000)
            if '"message"' in head or head.lstrip().startswith("{"):
                data = io.open(cache, encoding="utf-8", errors="replace").read()
        except Exception:      # noqa
            data = None
    if data is None:
        data = get(url).decode("utf-8", "replace")
    items = json.loads(data)["message"].get("items", [])
    best, bestscore = None, 0.0
    want = set(norm_title(title).split())
    for it in items:
        t = " ".join((it.get("title") or [""])[0].split())
        ws = set(norm_title(t).split())
        score = len(want & ws) / float(max(1, len(want | ws)))
        if score > bestscore:
            bestscore, best = score, it
    if best is None or bestscore < 0.45:
        return None, bestscore, data
    au = ["%s %s" % (a.get("given", ""), a.get("family", ""))
          for a in best.get("author", []) or []]
    issued = (best.get("issued", {}).get("date-parts") or [[None]])[0][0]
    return {
        "source": "crossref-query(score=%.2f)" % bestscore,
        "url": url,
        "title": " ".join((best.get("title") or [""])[0].split()),
        "authors": [x.strip() for x in au],
        "journal": " ".join((best.get("container-title") or [""])[0].split()),
        "year": str(issued) if issued else "",
        "doi": best.get("DOI", ""),
        "volume": best.get("volume", "") or "",
        "pages": best.get("page", "") or "",
        "publisher": best.get("publisher", ""),
        "raw": data,
    }, bestscore, data


def lastnames(val):
    """last names, lowercased, from either a bibtex author string or a list."""
    if isinstance(val, str):
        return [a[0].lower() for a in split_authors(val)]
    out = []
    for a in val or []:
        w = a.split()
        out.append((w[-1] if w else a).lower())
    return out


def main(do_fetch=True, do_pdfs=False):
    bib = parse_bib(os.path.join(MS, "references.bib"))
    print("parsed %d bib entries" % len(bib))
    order, seen = cite_order(os.path.join(MS, "manuscript.tex"))
    keys = [e["key"] for e in bib]
    dup = [k for k in set(keys) if keys.count(k) > 1]
    uncited = [k for k in keys if k not in seen]
    missing = [k for k in order if k not in keys]
    print("cited %d; duplicates %s; uncited %s; cited-but-missing %s"
          % (len(order), dup, uncited, missing))
    first_cite = {k: i + 1 for i, k in enumerate(order)}
    order_ok = [first_cite[k] for k in keys if k in first_cite]
    print("first-citation order monotonic:",
          order_ok == sorted(order_ok))

    rows = []
    for e in bib:
        key = e["key"]
        rec = {"key": key, "kind": e["kind"],
               "bib_title": clean_tex(e.get("title", "")),
               "bib_authors": clean_tex(e.get("author", "")),
               "bib_journal": clean_tex(e.get("journal", "")
                                        or e.get("booktitle", "")),
               "bib_year": e.get("year", ""),
               "bib_volume": e.get("volume", ""),
               "bib_pages": e.get("pages", ""),
               "bib_doi": e.get("doi", ""),
               "cited": key in seen,
               "first_cite": first_cite.get(key, "")}
        if do_fetch:
            rawfile = os.path.join(RAW, key + ".txt")
            got = None
            try:
                aid = arxiv_id(e.get("doi", ""))
                if aid:
                    got = fetch_arxiv(aid, cache=rawfile)
                elif e.get("doi"):
                    got = fetch_crossref(e["doi"], cache=rawfile)
                else:
                    got, sc, _ = fetch_crossref_query(
                        rec["bib_title"],
                        lastnames(e.get("author", ""))[0]
                        if lastnames(e.get("author", "")) else "",
                        cache=rawfile)
                    if got is None:
                        got = {"source": "NOT-FOUND", "url": "", "title": "",
                               "authors": [], "journal": "", "year": "",
                               "doi": "", "volume": "", "pages": "",
                               "raw": ""}
            except Exception as ex:      # noqa
                got = {"source": "FETCH-ERROR: %s" % ex, "url": "", "title": "",
                       "authors": [], "journal": "", "year": "", "doi": "",
                       "volume": "", "pages": "", "raw": ""}
            io.open(rawfile, "w", encoding="utf-8").write(
                (got.get("raw") or "")[:400000])
            for f in ("source", "url", "title", "journal", "year", "doi",
                      "volume", "pages"):
                rec["ref_" + f] = got.get(f, "")
            rec["ref_authors"] = "; ".join(got.get("authors", []) or [])
            rec["title_match"] = (norm_title(rec["bib_title"])
                                  == norm_title(rec["ref_title"]))
            bl = set(lastnames(e.get("author", "")))
            rl = set()
            for a in got.get("authors", []) or []:
                w = a.split()
                if w:
                    rl.add(w[-1].lower())
            rec["author_match"] = (bool(bl) and bl == rl)
            rec["authors_missing"] = "; ".join(sorted(bl - rl))
            rec["authors_extra"] = "; ".join(sorted(rl - bl))
            rec["year_match"] = (str(rec["bib_year"])
                                 == str(rec["ref_year"])[:4])
            rec["journal_match"] = (
                not rec["ref_journal"]
                or norm_title(rec["bib_journal"]).find(
                    norm_title(rec["ref_journal"])[:18]) >= 0
                or norm_title(rec["ref_journal"]).find(
                    norm_title(rec["bib_journal"])[:18]) >= 0)
            rec["vol_match"] = (not rec["ref_volume"]
                                or str(rec["bib_volume"])
                                == str(rec["ref_volume"]))
            pg = lambda s: re.sub(r"[^0-9A-Za-z]", "", (s or "")).lower()
            rec["pages_match"] = (not rec["ref_pages"]
                                  or pg(rec["bib_pages"]) == pg(rec["ref_pages"]))
            if do_pdfs and arxiv_id(e.get("doi", "")):
                aid = arxiv_id(e["doi"])
                try:
                    pdf = get("https://arxiv.org/pdf/%s" % aid, timeout=90)
                    io.open(os.path.join(PDFS, key + ".pdf"), "wb").write(pdf)
                    rec["pdf"] = "yes"
                except Exception as ex:      # noqa
                    rec["pdf"] = "err:%s" % str(ex)[:40]
        rows.append(rec)
        print("  %-24s %-22s title=%s authors=%s year=%s"
              % (key, rec.get("ref_source", "")[:22],
                 rec.get("title_match"), rec.get("author_match"),
                 rec.get("year_match")))

    import csv
    cols = ["key", "cited", "first_cite", "bib_authors", "bib_title",
            "bib_journal", "bib_year", "bib_volume", "bib_pages", "bib_doi",
            "ref_source", "ref_title", "ref_authors", "ref_journal",
            "ref_year", "ref_volume", "ref_pages", "ref_doi", "ref_url",
            "title_match", "author_match", "authors_missing", "authors_extra",
            "year_match", "journal_match", "vol_match", "pages_match", "pdf"]
    with io.open(os.path.join(OUT, "bibliography_check.csv"), "w",
                 encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("wrote ref_verify/bibliography_check.csv")


if __name__ == "__main__":
    main(do_fetch="--fetch" in sys.argv or True,
         do_pdfs="--pdfs" in sys.argv)

# -*- coding: utf-8 -*-
"""Fetch every network in the analysis set from the Network Data Repository.

Usage
-----
    python download_networks.py [--out DIR] [--limit N] [--sleep SEC]

Reads ``data/dataset_sources.csv`` (or ``dataset_sources.csv`` next to this
file), which lists for each network its landing page on
http://networkrepository.com and the verified direct download URL, downloads
the archive and unpacks it into ``<out>/<dataset_id>/``.

The data are distributed by the Network Data Repository for research use; see
http://networkrepository.com for its terms and cite Rossi & Ahmed (AAAI 2015)
when reusing them.
"""
import argparse, csv, io, os, ssl, sys, time, zipfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(HERE, "dataset_sources.csv")
if not os.path.exists(DEFAULT_CSV):
    DEFAULT_CSV = os.path.join(HERE, "..", "data", "dataset_sources.csv")

UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124.0 Safari/537.36")
HEADERS = {"User-Agent": UA,
           "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
           "Accept-Language": "en-US,en;q=0.9"}


def fetch(url, timeout=60):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(url, headers=HEADERS)
    return urllib.request.urlopen(req, timeout=timeout, context=ctx)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--out", default=os.path.join(HERE, "networks"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--sleep", type=float, default=0.2)
    a = ap.parse_args()

    rows = list(csv.DictReader(io.open(a.csv, encoding="utf-8")))
    if a.limit:
        rows = rows[: a.limit]
    os.makedirs(a.out, exist_ok=True)

    ok = skip = fail = 0
    for i, r in enumerate(rows, 1):
        ds = r["dataset_id"]
        url = (r.get("download_url_zip") or "").strip()
        dest = os.path.join(a.out, ds)
        if os.path.isdir(dest) and os.listdir(dest):
            skip += 1
            continue
        if not url:
            print("[%3d] %-52s no URL" % (i, ds)); fail += 1; continue
        try:
            data = fetch(url).read()
            zf = zipfile.ZipFile(io.BytesIO(data))
            os.makedirs(dest, exist_ok=True)
            zf.extractall(dest)
            ok += 1
            print("[%3d] %-52s ok (%d bytes)" % (i, ds, len(data)))
        except Exception as e:
            fail += 1
            print("[%3d] %-52s FAIL %s" % (i, ds, e))
        time.sleep(a.sleep)
    print("\ndone: %d downloaded, %d already present, %d failed -> %s"
          % (ok, skip, fail, a.out))
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

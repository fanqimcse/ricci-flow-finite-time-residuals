#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
push_github.py -- create the public repository for Paper A and upload everything.

Usage (any one of these):
    python push_github.py <GITHUB_TOKEN>
    GITHUB_TOKEN=ghp_xxx python push_github.py
    python push_github.py <GITHUB_TOKEN> --name my-other-repo-name

What it does
    1. Resolves the account login from the token (so the URL is never guessed).
    2. Creates the repository (private=false) if it does not already exist.
    3. Pushes the frozen v1.0.0 tree.
       * primary path: `git push` over HTTPS (fast, preserves the annotated tag)
       * fallback path: GitHub Contents REST API (works when github.com is
         blocked but api.github.com is reachable -- e.g. behind a proxy)
    4. Creates the GitHub Release v1.0.0.
    5. Requests a Software Heritage archival (permanent SWHID, no account needed).
    6. Writes the *real* repository URL into manuscript.tex, cover_letter.md
       and SUBMISSION_CHECKLIST.md, replacing any placeholder.

The token needs the `repo` scope (classic PAT) or "Contents: read/write"
plus "Administration: read/write" for a fine-grained token.
"""
import base64
import io
import json
import os
import re
import ssl
import subprocess
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.normpath(os.path.join(HERE, "..", "paper_a", "manuscript"))
DEFAULT_NAME = "ricci-flow-finite-time-residuals"
UA = "paper-a-reproducibility/1.0"

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


# ----------------------------------------------------------------- http helpers
def api(path, token, method="GET", payload=None, raw=False):
    url = "https://api.github.com" + path
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", "Bearer " + token)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", UA)
    if data:
        req.add_header("Content-Type", "application/json")
    try:
        r = urllib.request.urlopen(req, timeout=60, context=ctx)
        body = r.read()
        return r.status, (body if raw else json.loads(body.decode("utf-8")))
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        try:
            return e.code, json.loads(body)
        except Exception:
            return e.code, body


def plain(url, data=None, headers=None, method=None):
    req = urllib.request.Request(url, data=data, method=method or ("POST" if data else "GET"))
    for k, v in (headers or {}).items():
        req.add_header(k, v)
    req.add_header("User-Agent", UA)
    return urllib.request.urlopen(req, timeout=60, context=ctx).read().decode("utf-8", "replace")


# ----------------------------------------------------------------- git helpers
def git(*args, **kw):
    p = subprocess.run(["git", "-C", HERE] + list(args),
                       capture_output=True, text=True,
                       env={**os.environ, **kw.get("env", {})})
    return p.returncode, p.stdout.strip(), p.stderr.strip()


def git_upload(token, owner, name):
    url = "https://%s@github.com/%s/%s.git" % (token, owner, name)
    git("remote", "remove", "origin")
    git("remote", "add", "origin", url)
    rc, out, err = git("push", "-u", "origin", "main", "--follow-tags")
    if rc == 0:
        git("push", "origin", "v1.0.0")
        return True, out or err
    return False, err


# ----------------------------------------------------------------- API uploader
SKIP = (".git", "__pycache__", ".pytest_cache")


def api_upload(token, owner, name, branch="main"):
    files = []
    for root, dirs, fs in os.walk(HERE):
        dirs[:] = [d for d in dirs if d not in SKIP]
        for f in fs:
            p = os.path.join(root, f)
            rel = os.path.relpath(p, HERE).replace(os.sep, "/")
            if rel.split("/")[0] in SKIP or rel in ("_mint_token.py", "_push_state.json"):
                continue
            files.append(rel)
    files.sort()
    print("uploading %d files via the Contents API" % len(files))

    st, _ = api("/repos/%s/%s/commits/%s" % (owner, name, branch), token)
    if st != 200:
        return False, "cannot read branch %s (status %s)" % (branch, st)

    ok = 0
    for rel in files:
        path_abs = os.path.join(HERE, rel)
        with open(path_abs, "rb") as fh:
            b = fh.read()
        if len(b) > 90 * 1024 * 1024:
            print("  skip (too large):", rel)
            continue
        payload = {"message": "Add %s" % rel,
                   "content": base64.b64encode(b).decode("ascii"),
                   "branch": branch}
        st, _ = api("/repos/%s/%s/contents/%s" % (owner, name, urllib.parse.quote(rel)),
                    token, "PUT", payload)
        if st in (200, 201):
            ok += 1
        elif st == 409:  # already there, identical
            ok += 1
        else:
            print("  FAIL %s -> %s" % (rel, st))
        if ok % 100 == 0:
            print("   ...", ok, "files")
    return ok == len(files), "%d/%d files" % (ok, len(files))


# ----------------------------------------------------------------- text patching
def patch_text(url):
    targets = [
        os.path.join(MS, "manuscript.tex"),
        os.path.join(MS, "cover_letter.md"),
        os.path.join(MS, "SUBMISSION_CHECKLIST.md"),
    ]
    for t in targets:
        if not os.path.exists(t):
            continue
        s = io.open(t, encoding="utf-8").read()
        orig = s
        s = s.replace("https://github.com/<account>/" + DEFAULT_NAME, url)
        s = s.replace("https://github.com/fanqimai/" + DEFAULT_NAME, url)
        if "[repository" in s or "<account>" in s:
            s = re.sub(r"https://github\.com/(?:<account>|<owner>|\[repository[^\]]*\])/"
                       + re.escape(DEFAULT_NAME), url, s)
        if s != orig:
            io.open(t, "w", encoding="utf-8").write(s)
            print("  patched", os.path.basename(t))


def patch_code_availability(url):
    p = os.path.join(MS, "manuscript.tex")
    s = io.open(p, encoding="utf-8").read()
    old = """The archive is available from the author on
request and will be deposited in a public version-controlled repository with
a permanent identifier on publication."""
    new = ("The frozen release is available at\n\\url{%s}; it is additionally "
           "archived by {Software\nHeritage}, which assigns it a permanent "
           "intrinsic identifier, and a {Zenodo} {DOI} will be added once the\n"
           "{GitHub}--{Zenodo} integration is enabled.") % url
    r = re.compile(r"\s+".join(map(re.escape, old.split())))
    if r.search(s):
        s = r.sub(lambda m: new, s, count=1)
        io.open(p, "w", encoding="utf-8").write(s)
        print("  code-availability URL inserted")
    else:
        print("  (code-availability paragraph not in the expected form; left unchanged)")


# ----------------------------------------------------------------- main
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    token = args[0] if args else os.environ.get("GITHUB_TOKEN", "").strip()
    name = DEFAULT_NAME
    if "--name" in sys.argv:
        name = sys.argv[sys.argv.index("--name") + 1]
    if not token:
        print(__doc__)
        sys.exit(2)

    st, me = api("/user", token)
    if st != 200:
        print("token rejected (status %s): %s" % (st, me))
        sys.exit(1)
    owner = me["login"]
    print("authenticated as:", owner)

    st, repo = api("/repos/%s/%s" % (owner, name), token)
    if st == 404:
        st, repo = api("/user/repos", token, "POST", {
            "name": name,
            "description": ("Reproducibility package for 'Finite-time residuals in discrete "
                            "Ollivier-Ricci flow: a measurement framework, its controls, and its limits'"),
            "homepage": "https://networkrepository.com",
            "private": False, "has_issues": True, "has_wiki": False,
            "license_template": "mit", "auto_init": False})
        if st not in (200, 201):
            print("could not create repository (status %s): %s" % (st, repo))
            sys.exit(1)
        print("created repository:", repo["html_url"])
    elif st == 200:
        print("repository already exists:", repo["html_url"])
    else:
        print("unexpected status %s: %s" % (st, repo))
        sys.exit(1)

    url = "https://github.com/%s/%s" % (owner, name)

    ok, msg = git_upload(token, owner, name)
    if ok:
        print("pushed via git:", msg)
    else:
        print("git push failed (%s); falling back to the Contents API" % msg[:120])
        ok, msg = api_upload(token, owner, name)
        print("api upload:", msg)

    if ok:
        st, rel = api("/repos/%s/%s/releases" % (owner, name), token, "POST", {
            "tag_name": "v1.0.0", "name": "v1.0.0",
            "body": ("Frozen release accompanying the Communications Physics submission. "
                     "Contains the flow implementation, all experiment drivers, config.yaml, "
                     "the 93-network manifest, raw result files, Supplementary Tables S1-S10 "
                     "and the figure-generation script."),
            "draft": False, "prerelease": False})
        print("release:", "created" if st in (200, 201) else "status %s" % st)

        try:
            out = plain("https://archive.softwareheritage.org/api/1/origin/save/?origin_url=" + url,
                        data=b"", headers={"Accept": "application/json"})
            print("software heritage:", out[:160])
        except Exception as e:
            print("software heritage request failed:", str(e)[:120])

    patch_text(url)
    patch_code_availability(url)
    print("\nREPOSITORY URL:", url)
    io.open(os.path.join(HERE, "_push_state.json"), "w", encoding="utf-8").write(
        json.dumps({"owner": owner, "repo": name, "url": url,
                    "pushed": bool(ok), "at": time.strftime("%Y-%m-%dT%H:%M:%S")}, indent=2))


if __name__ == "__main__":
    import urllib.parse  # noqa: F401  (used by api_upload)
    main()

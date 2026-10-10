# -*- coding: utf-8 -*-
"""
_update_release_r6.py -- stage the round-6 artefacts into github_release/
so that the public repository (and hence Zenodo) can be updated to v1.0.1.

It only COPIES files and prints what changed; it never commits or pushes.
"""
import hashlib
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REL = os.path.join(HERE, "github_release")
MS = os.path.join(HERE, "paper_a", "manuscript")

# (source, destination inside github_release)
MANUSCRIPT_FILES = [
    "manuscript.tex",
    "references.bib",
    "SUBMISSION_CHECKLIST.md",
    "RESPONSE_TO_REVIEW_ROUND6.md",
    "Supplementary_Table_S1_dataset_sources.csv",
    "Supplementary_Table_S11_long_horizon_three_networks.csv",
    "Supplementary_Table_S12_weight_semantics.csv",
    "Supplementary_Table_S13a_checkpoints.csv",
    "Supplementary_Table_S13b_final_probes.csv",
    "Supplementary_Table_S13c_initial_weight_perturbation.csv",
    "Supplementary_Table_S14_reference_audit_three_way.csv",
    "REFERENCE_AUDIT_ROUND6.md",
]
CODE_FILES = [
    "exp_round6.py",
    "make_supplementary_r6.py",
    "_make_s1_weights.py",
    "_r6_rawstats.py",
    "_r6_sem_summary.py",
    "_fill_response_r6.py",
    "_fill_response_r6b.py",
    "_update_release_r6.py",
    "_verify_refs.py",
    "_verify_openalex.py",
    "_verify_conf.py",
    "_find_arxiv.py",
    "_make_ref_audit_report.py",
]
CODE_EXTRA = [
    ("paper_a", "manuscript", "_patch_r6.py"),
    ("paper_a", "manuscript", "_patch_r6b.py"),
    ("paper_a", "manuscript", "_patch_r6c.py"),
    ("paper_a", "manuscript", "_patch_r6d.py"),
    ("paper_a", "manuscript", "_patch_r6e.py"),
    ("paper_a", "manuscript", "_patch_r6f.py"),
    ("paper_a", "manuscript", "_patch_r6g.py"),
    ("paper_a", "manuscript", "_patch_r6h.py"),
    ("paper_a", "manuscript", "_patch_refs_r6.py"),
]


def sha(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()[:16]


def copy(src, dst, log):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst) and sha(src) == sha(dst):
        log.append(("=", os.path.relpath(dst, REL)))
        return
    shutil.copy2(src, dst)
    log.append(("+" if not os.path.exists(dst) else "~",
                os.path.relpath(dst, REL)))


def main():
    log = []
    for f in MANUSCRIPT_FILES:
        s = os.path.join(MS, f)
        if os.path.exists(s):
            copy(s, os.path.join(REL, "manuscript", f), log)
        else:
            log.append(("!", "missing " + f))
    for f in CODE_FILES:
        s = os.path.join(HERE, f)
        if os.path.exists(s):
            copy(s, os.path.join(REL, "code", f), log)
        else:
            log.append(("!", "missing " + f))
    for parts in CODE_EXTRA:
        s = os.path.join(HERE, *parts)
        if os.path.exists(s):
            copy(s, os.path.join(REL, "code", parts[-1]), log)
        else:
            log.append(("!", "missing " + parts[-1]))

    # round-6 raw experiment outputs
    for d in ("exp_round6", "exp_round6_focus", "exp_round6_rest"):
        src = os.path.join(HERE, d)
        if not os.path.isdir(src):
            continue
        for root, _, files in os.walk(src):
            if "__pycache__" in root:
                continue
            for fn in files:
                if fn.endswith((".csv", ".json")):
                    rel = os.path.relpath(os.path.join(root, fn), src)
                    copy(os.path.join(root, fn),
                         os.path.join(REL, "results", d, rel), log)

    changed = [x for x in log if x[0] in "+~"]
    print("staged: %d added/updated, %d unchanged, %d missing"
          % (len(changed), sum(1 for x in log if x[0] == "="),
             sum(1 for x in log if x[0] == "!")))
    for tag, p in log:
        if tag in "+~!":
            print("  %s %s" % (tag, p))


if __name__ == "__main__":
    main()

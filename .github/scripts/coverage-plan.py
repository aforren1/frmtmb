#!/usr/bin/env python3
"""Which coverage jobs a commit needs (test-coverage.yaml, job `plan`).

Writes two lines for $GITHUB_OUTPUT:
  core=true|false       run the core job
  matrix=[{...}, ...]   the extension jobs to run, as matrix `include`

A job is skipped only when no file it reads changed: core's coverage
reads everything outside extensions/ that is not documentation or
another workflow; an extension's reads that, its own directory and the
directories of the siblings it installs. When the changed files cannot
be listed (a new branch, a force push, a shallow history), every job
runs. Codecov carries a skipped flag's last report forward
(.github/codecov.yml), so a skipped job does not lower the totals.

Environment: EVENT (github.event_name), BEFORE (github.event.before),
BASE (the pull request's base sha), HEAD (github.sha). CHANGED_FILES,
newline separated, replaces the git call, for a dry run.
"""
import json
import os
import subprocess
import sys

# `siblings` are the frmtmb extensions the package names in its
# DESCRIPTION (Imports or Suggests), in install order; `absent` the
# packages a test may skip as not installed (the skip gate fails on any
# other). The comments in test-coverage.yaml say why each is there.
EXTENSIONS = [
    {"pkg": "frmtmb.coupling", "siblings": ""},
    {"pkg": "frmtmb.eam", "siblings": "frmtmb.sample"},
    {"pkg": "frmtmb.latent", "siblings": ""},
    {"pkg": "frmtmb.learn", "siblings": "frmtmb.eam"},
    {"pkg": "frmtmb.ode", "siblings": ""},
    {"pkg": "frmtmb.sample",
     "siblings": "frmtmb.latent frmtmb.ode frmtmb.eam frmtmb.survival",
     "absent": "RTMBode"},
    {"pkg": "frmtmb.survival", "siblings": ""},
    {"pkg": "frmtmb.spatial", "siblings": ""},
]


def neutral(path):
    """A file no coverage job reads: covr runs the tests only."""
    if path.startswith(("dev/", "docs/", "pkgdown/")):
        return True
    # a testthat snapshot (tests/testthat/_snaps/*.md) is read by a test
    tests = path.startswith("tests/") or (
        path.startswith("extensions/") and "/tests/" in path)
    if path.endswith(".md") and not tests:
        return True
    if path in ("codemeta.json", "_pkgdown.yml"):
        return True
    if path.startswith("extensions/") and path.endswith("/_pkgdown.yml"):
        return True
    # another workflow; this one and the scripts it runs are read
    if path.startswith(".github/workflows/"):
        return path != ".github/workflows/test-coverage.yaml"
    return False


def changed_files():
    listed = os.environ.get("CHANGED_FILES")
    if listed is not None:
        return [p for p in listed.splitlines() if p.strip()]
    event = os.environ.get("EVENT", "")
    head = os.environ.get("HEAD", "")
    base = os.environ.get("BASE" if event == "pull_request" else "BEFORE",
                          "")
    if not base or not head or set(base) == {"0"}:
        return None
    # --no-renames: with rename detection git prints only the new path of
    # a moved file, and the job that read the old one would be skipped
    r = subprocess.run(["git", "diff", "--name-only", "--no-renames", base,
                        head], capture_output=True, text=True)
    if r.returncode != 0:
        return None
    return [p for p in r.stdout.splitlines() if p.strip()]


def plan(files):
    if files is None:
        return True, list(EXTENSIONS), "changed files unknown: run all"
    core = False
    touched = set()
    for p in files:
        if neutral(p):
            continue
        if p.startswith("extensions/"):
            parts = p.split("/")
            if len(parts) > 1:
                touched.add(parts[1])
            continue
        core = True
    run = []
    for e in EXTENSIONS:
        reads = {e["pkg"]} | set(e["siblings"].split())
        if core or reads & touched:
            run.append(e)
    why = "%d files; core %s; extensions touched: %s" % (
        len(files), "changed" if core else "unchanged",
        " ".join(sorted(touched)) or "none")
    return core, run, why


def main():
    core, run, why = plan(changed_files())
    sys.stderr.write("coverage plan: %s\n" % why)
    sys.stderr.write("core job: %s; extension jobs: %s\n" % (
        core, " ".join(e["pkg"] for e in run) or "none"))
    print("core=%s" % ("true" if core else "false"))
    print("matrix=%s" % json.dumps(run, separators=(",", ":")))


if __name__ == "__main__":
    main()

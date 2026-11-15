"""Review of lane hygiene, item 3177: dry run of coverage-plan.py.

Runs the lane's script itself (its changed_files() with the git call it
makes) against a bare scratch clone, for the last N first-parent
commits of main, as a push of that one commit would (BEFORE = first
parent). Then compares the plan with the one --no-renames gives, and
runs the constructed rename cases.

  uv run --no-project python dev/reviews/hygiene-review-plan.py <bare-git-dir> <N>
"""
import importlib.util
import os
import subprocess
import sys

gitdir, n = sys.argv[1], int(sys.argv[2])
spec = importlib.util.spec_from_file_location(
    "plan", ".github/scripts/coverage-plan.py")
plan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plan)
os.environ["GIT_DIR"] = gitdir
os.environ["GIT_OPTIONAL_LOCKS"] = "0"


def git(*a):
    return subprocess.run(["git", *a], capture_output=True, text=True,
                          check=True).stdout


commits = git("rev-list", "--first-parent", "-n", str(n), "main").split()
jobs = 0
alljobs = 0
nine = 0
differ = []
for c in commits:
    parents = git("rev-list", "--parents", "-n", "1", c).split()[1:]
    os.environ.pop("CHANGED_FILES", None)
    os.environ.update(EVENT="push", BEFORE=parents[0], HEAD=c)
    files = plan.changed_files()
    core, run, _ = plan.plan(files)
    k = int(core) + len(run)
    jobs += k
    alljobs += 9
    nine += k == 9
    nr = [p for p in git("diff", "--name-only", "--no-renames", parents[0],
                         c).splitlines() if p.strip()]
    core2, run2, _ = plan.plan(nr)
    if (core, [e["pkg"] for e in run]) != (core2, [e["pkg"] for e in run2]):
        differ.append((c[:8], core, [e["pkg"] for e in run], core2,
                       [e["pkg"] for e in run2]))
print("commits", len(commits), "jobs run", jobs, "of", alljobs,
      "; all nine:", nine, "; parents of each: 1 to",
      max(len(git("rev-list", "--parents", "-n", "1", c).split()) - 1
          for c in commits))
print("plans that --no-renames changes:", len(differ))
for d in differ:
    print("  ", d)

# constructed cases: file lists as `git diff --name-only` prints them for
# a pure rename (rename detection on, git's default), and as
# --no-renames prints them
cases = {
    "core R file moved to an extension": (
        ["extensions/frmtmb.eam/R/y.R"],
        ["R/y.R", "extensions/frmtmb.eam/R/y.R"]),
    "core test moved to dev/attic": (
        ["dev/attic/test-y.R"],
        ["dev/attic/test-y.R", "tests/testthat/test-y.R"]),
    "eam test moved to frmtmb.ode": (
        ["extensions/frmtmb.ode/tests/testthat/test-y.R"],
        ["extensions/frmtmb.eam/tests/testthat/test-y.R",
         "extensions/frmtmb.ode/tests/testthat/test-y.R"]),
    "testthat snapshot only": (
        ["tests/testthat/_snaps/print.md"], None),
}
for name, (renamed, plain) in cases.items():
    core, run, _ = plan.plan(renamed)
    line = "CASE %-36s rename-detected: core %s, ext %s" % (
        name, core, " ".join(e["pkg"] for e in run) or "none")
    if plain:
        core2, run2, _ = plan.plan(plain)
        line += " | --no-renames: core %s, ext %s" % (
            core2, " ".join(e["pkg"] for e in run2) or "none")
    print(line)

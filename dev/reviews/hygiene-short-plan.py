"""Short check of lane hygiene's fix round, F2: coverage-plan.py.

1. A scratch repository (created under <scratch>, never in frmtmb):
   a core file moved into an extension, a testthat snapshot alone, an
   extension snapshot alone, a README alone, and the rename case again
   with diff.renames=copies set in the repository. The lane's script is
   run verbatim as the workflow runs it (EVENT=push, BEFORE, HEAD).
2. The last N first-parent commits of main, read from a bare clone
   under <scratch>, through the script's own changed_files().

  uv run --no-project python -I dev/reviews/hygiene-short-plan.py \
      <lane-tree> <source-repo> <scratch> <N>
"""
import importlib.util
import os
import shutil
import subprocess
import sys

tree, src, scratch, n = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
script = os.path.join(tree, ".github", "scripts", "coverage-plan.py")
env0 = dict(os.environ, GIT_OPTIONAL_LOCKS="0")
for k in ("GIT_DIR", "GIT_WORK_TREE", "CHANGED_FILES"):
    env0.pop(k, None)


def git(cwd, *a, env=None):
    return subprocess.run(["git", "-c", "user.name=s", "-c",
                           "user.email=s@s", *a], cwd=cwd, env=env or env0,
                          capture_output=True, text=True,
                          check=True).stdout


def run_script(repo, before, head):
    e = dict(env0, EVENT="push", BEFORE=before, HEAD=head)
    r = subprocess.run([sys.executable, "-I", script], cwd=repo, env=e,
                       capture_output=True, text=True, check=True)
    return " | ".join(r.stderr.strip().splitlines()), r.stdout


def write(repo, path, text="x\n"):
    p = os.path.join(repo, *path.split("/"))
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        f.write(text)


## 1. scratch repository ----------------------------------------------------
stamp = str(os.getpid())
repo = os.path.join(scratch, "planrepo-" + stamp)
os.makedirs(repo)
git(repo, "init", "-q", "-b", "main")
body = "".join("line %d of a core R file\n" % i for i in range(40))
write(repo, "R/y.R", body)
write(repo, "R/z.R", "z\n")
write(repo, "extensions/frmtmb.eam/R/x.R")
write(repo, "README.md")
git(repo, "add", "-A")
git(repo, "commit", "-q", "-m", "base")


def step(name, change):
    before = git(repo, "rev-parse", "HEAD").strip()
    change()
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", name)
    head = git(repo, "rev-parse", "HEAD").strip()
    default = git(repo, "diff", "--name-only", before, head).split()
    err, out = run_script(repo, before, head)
    print("CASE %s" % name)
    print("  git diff --name-only (default):", " ".join(default))
    print("  script:", err)
    print("  stdout core:", out.splitlines()[0])


step("git mv R/y.R extensions/frmtmb.eam/R/y.R",
     lambda: git(repo, "mv", "R/y.R", "extensions/frmtmb.eam/R/y.R"))
step("add tests/testthat/_snaps/print.md only",
     lambda: write(repo, "tests/testthat/_snaps/print.md", "# snap\n"))
step("add extensions/frmtmb.ode/tests/testthat/_snaps/o.md only",
     lambda: write(repo, "extensions/frmtmb.ode/tests/testthat/_snaps/o.md"))
step("edit README.md only", lambda: write(repo, "README.md", "y\n"))
git(repo, "config", "diff.renames", "copies")
step("diff.renames=copies; git mv R/z.R dev/attic/z.R",
     lambda: (os.makedirs(os.path.join(repo, "dev", "attic")),
              git(repo, "mv", "R/z.R", "dev/attic/z.R")))

## 2. last N first-parent commits of main -----------------------------------
bare = os.path.join(scratch, "planbare-%s.git" % stamp)
subprocess.run(["git", "clone", "-q", "--bare", "--no-local", src, bare],
               env=env0, check=True)
spec = importlib.util.spec_from_file_location("plan", script)
plan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plan)
os.environ.update(GIT_DIR=bare, GIT_OPTIONAL_LOCKS="0")
os.environ.pop("CHANGED_FILES", None)
commits = git(bare, "rev-list", "--first-parent", "-n", str(n),
              "ec8178d3").split()
jobs = nine = 0
maxpar = 0
for c in commits:
    parents = git(bare, "rev-list", "--parents", "-n", "1", c).split()[1:]
    maxpar = max(maxpar, len(parents))
    os.environ.update(EVENT="push", BEFORE=parents[0], HEAD=c)
    files = plan.changed_files()
    if files is None:
        print("UNKNOWN files for", c[:8])
    core, run, _ = plan.plan(files)
    k = int(core) + len(run)
    jobs += k
    nine += k == 9
print("DRYRUN commits %d (tip %s) jobs run %d of %d; all nine: %d; "
      "parents up to %d" % (len(commits), commits[0][:8], jobs,
                            9 * len(commits), nine, maxpar))

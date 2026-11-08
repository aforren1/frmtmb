"""Apply the fix round's edits to dev/ci-findings.md (lane ci, review of
2026-11-06). Each old text must occur exactly once."""
import pathlib

p = pathlib.Path("dev/ci-findings.md")
t = p.read_text(encoding="utf-8")
crlf = "\r\n" in t
t = t.replace("\r\n", "\n")


def rep(old, new):
    global t
    assert t.count(old) == 1, old[:70]
    t = t.replace(old, new)


rep("""| 4 CRAN lag script | `dev/release/cran-lag.R`. No release checklist exists under `dev/release/` (or anywhere in the repository): needs the user. |""",
    """| 4 CRAN lag script | `dev/release/cran-lag.R`. Fixed in the fix round: it failed open when CRAN was unreachable (review, MAJOR 1); it now fails closed and reads a library stack. No release checklist exists under `dev/release/` (or anywhere in the repository): needs the user. |""")

rep("""| 6 other red or flaky | See section 6. |""",
    """| 6 other red or flaky | See section 6. brms-tiers can exceed 90 minutes on a cold cache and a slow mirror: limits raised in the fix round. |""")

rep("""Windows OpenBLAS build honors `OPENBLAS_CORETYPE` was not checked. So
the mechanism (which kernel, which operation) is not established; the
runner's CPU was not recorded. ubuntu-next now records it (below).""",
    """Windows OpenBLAS build honors `OPENBLAS_CORETYPE` was not checked. So
the mechanism (which kernel, which operation) is not established; the
runner's CPU was not recorded. ubuntu-next now records it (below).

Lead, from the review (minor 2), checked in the log: 26.04 installs
Ubuntu's x86-64-v3 build, not a generic one. `un-26-core.log`:

    Preparing to unpack .../3-libopenblas0-pthread_0.3.32+ds-5_amd64v3.deb ...
    Preparing to unpack .../4-libopenblas0_0.3.32+ds-5_amd64v3.deb ...

(334 lines of that log name `amd64v3`; the apt indexes are
`resolute-updates/main amd64v3 Packages` and so on.) Ubuntu compiles an
amd64v3 package for x86-64-v3, which includes AVX2 and FMA, so the
compiler may contract a multiply and an add into one FMA in OpenBLAS's
C code. `dev/ciharden-openblas.sh` uses upstream's generic Windows
build, `OpenBLAS-0.3.32-x64.zip`, which is a different binary. That can
explain why it did not reproduce. Upstream publishes no Windows build
for x86-64-v3; a local reproduction would mean building OpenBLAS from
source with `-march=x86-64-v3`, which is not cheap, so it was not
tried. 24.04 has no amd64v3 packages. The record step now also prints
`dpkg -l 'libopenblas*'` and `apt-cache policy libopenblas0-pthread`,
whose index line names the variant.""")

rep("""(26.04 job totals: 85 min on 10-19 and 10-26, 72 min on 11-02, against
a 90-minute limit.) `stan-next`'s `timeout-minutes` is now 150, about
1.6 times the worst job measured. The push job `brms-likelihood.yaml`
keeps 90: its cache hits, and its last five runs took 5 to 10 minutes.""",
    """(26.04 job totals: 85 min on 10-19 and 10-26, 72 min on 11-02, against
a 90-minute limit.) `stan-next`'s `timeout-minutes` is now 150. The
worst job was cut off at 90; from its identity tier of 88.3 minutes plus
about 7 to 10 minutes of setup and methods tier (the 26.04 jobs), its
full length would have been about 98, so 150 is about 1.5 times that.
The identity tier also has a step limit of 120 minutes now, as in
brms-likelihood.yaml (below).

Correction (review, minor 3). The first report said the push job
`brms-likelihood.yaml` keeps 90 minutes because "its cache hits, and
its last five runs took 5 to 10 minutes". That was wrong for the run on
3cd2a0c8: brms-tiers 42214681196 took 44 minutes, and its steps
(`gh api .../runs/42214681196/jobs`) were setup-r 25 minutes and
setup-r-dependencies 15 minutes (the slow mirror), then the identity and
methods tiers about 1 minute each from a warm cache. Its cache prefix
carries the rstan version, so a new rstan release gives a cold cache,
and ubuntu-next measured a cold identity tier on 24.04 at 49 to 88
minutes. A cold cache and a slow mirror together: about 41 + 88 + 3 =
132 minutes, which exceeds 90. So brms-likelihood.yaml now has a job
limit of 180 and a step limit of 120 on the identity tier (88 at the
slowest measured, a factor of 1.36), the test-coverage pattern: the
step limit catches a wedged run, the job limit leaves room for setup.
The job limit is 180 rather than 150 because 41 + 120 + the gp and
methods tiers exceeds 150. stan-next keeps 150 (it runs weekly and
reports only) with the same 120 on its identity tier.""")

rep("""- "Record the numerical libraries" prints the CPU model
  (`/proc/cpuinfo`) and runs with `OPENBLAS_VERBOSE=2`, which makes a
  DYNAMIC_ARCH OpenBLAS print its kernel ("Core: ...") at load. Not
  verified that Ubuntu's build prints it: CI only.""",
    """- "Record the numerical libraries" prints the CPU model
  (`/proc/cpuinfo`) and runs with `OPENBLAS_VERBOSE=2`, which makes a
  DYNAMIC_ARCH OpenBLAS print its kernel ("Core: ...") at load. Not
  verified that Ubuntu's build prints it: CI only. In the fix round it
  also prints `dpkg -l 'libopenblas*'` and `apt-cache policy
  libopenblas0-pthread` (each followed by `; true`, so a missing tool
  cannot fail the step).""")

p.write_text(t.replace("\n", "\r\n") if crlf else t, encoding="utf-8",
             newline="")
print("ok, crlf =", crlf)

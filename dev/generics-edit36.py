import io

P = "C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-findings.md"
s = io.open(P, encoding="utf-8", newline="").read()


def sub1(old, new):
    global s
    n = s.count(old)
    if n != 1:
        raise SystemExit("pattern found %d times: %r" % (n, old[:60]))
    s = s.replace(old, new)


sub1("""The block above is the replacement: both clocks timed on every
replicate, both ticks printed, 30 interleaved replicates one fresh
process each, two controls that must report about zero, and a
permutation test on the paired minima. The attribution survives and is
the part worth keeping: adding `loadNamespace("nlme")` and
`loadNamespace("generics")` to the BASE arm accounts for most of the
cost, and `ACTIVE - SWAP` is not distinguishable from zero on the same
run, which is the side-by-side the mechanism change needed.
""",
"""The block above is the replacement: both clocks timed on every
replicate, both ticks printed, 30 interleaved replicates one fresh
process each, all three builds as arms of ONE run, two controls, and a
permutation test on the paired minima. It was run with nothing else on
the machine, after both `R CMD check` runs had finished.

What it says, and what it does not:

* **The two new hard Imports cost +0.0168 s.** That is the BASE arm
  with `loadNamespace("nlme")` and `loadNamespace("generics")` added,
  minus BASE. The review's headline was +0.0167 s on the swap build and
  its own Imports arm +0.0183 s; an earlier run of this script, on a
  busier machine, gave +0.0158 s. Three runs, one figure: about 17 ms,
  and it is the Imports.
* **The whole change costs +0.0210 s on the active-binding build**, and
  +0.0241 s on the swap build in the same run. The residual after the
  Imports is +0.0042 s for ACTIVE and +0.0073 s for SWAP.
* **The two designs cannot be told apart at load.** `ACTIVE - SWAP` is
  -0.0031 s with p = 0.1251. Across this session's three runs the
  headline `ACTIVE - BASE` read +0.0280, +0.0229 and +0.0210 s, a
  spread of 7 ms that is larger than any difference between designs,
  so the page reports the Imports attribution as the finding and the
  headline as a range.
* **One control is not clean and it is reported rather than hidden.**
  `CONTROL nothing` gives p = 0.0293 at a magnitude of 0.0000 s. The
  arm is two `Sys.time()` reads, 0.0002 s, so its minima differ by a
  microsecond or two and a permutation test on MINIMA of values that
  close ranks by noise rather than by effect. The magnitude is the
  honest reading there; the p-value on that arm is not.
* `proc.time()` measured a minimum tick of 0.0100 s and a median of
  0.0200 s in this run, so the same difference on the old instrument
  reads as +0.0300 s, three ticks. That is the instrument that produced
  the withdrawn +0.040 s.
""")

sub1("""### 5. There is no new `R CMD check` NOTE

Round 1 introduced one, `unlockBinding`, and put a choice to the user
between keeping it and giving back a load order. **That choice was
false and is withdrawn**: the active binding removes the NOTE and
repairs more than the swap did. The remaining NOTE on core is the
pre-existing V8 math-rendering one.
""",
"""### 5. There is no new `R CMD check` NOTE

Round 1 introduced one, `unlockBinding`, and put a choice to the user
between keeping it and giving back a load order. **That choice was
false and is withdrawn**: the active binding removes the NOTE and
repairs more than the swap did. Confirmed on the final build:
`checking R code for possible problems ... OK`, and core's one NOTE is
the pre-existing V8 math-rendering one.
""")

sub1("""## What was run

<!--BLOCK:suite-->
""",
"""## What was run

The session driving this lane ended once, mid-round, with an
`R CMD check` of core half-built. Before anything else the resumed run
confirmed every R file this lane touched still parses
(`dev/generics-parsecheck.R`: 128 files, the one failure being the
review's own `dev/genrev-scale2.R`, which this lane never edited), that
the lane's library holds no hollow package directory, and that the
installed core is newer than its newest source edit. The check was
then re-run from scratch, alone on the machine.

The core suite below ran against the build with the `hypothesis()` fix
in; no source changed after it.

<!--BLOCK:suite-->
""")

sub1("""<!--BLOCK:check-->
""",
"""<!--BLOCK:check-->

The first core check of the active-binding build was
`Status: 1 ERROR, 1 NOTE`: the `hypothesis()` regression described
under BLOCKER A, found here and nowhere else. After the fix, the block
above.
""")

io.open(P, "w", encoding="utf-8", newline=chr(10)).write(s)
print("ok")

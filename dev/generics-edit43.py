# Punch round 2: the findings page. Markers are spliced by
# dev/generics-assemble.R afterwards, so every count stays generated.
P = "C:/Users/adf44/source/r/frmtmb-wt-generics/dev/generics-findings.md"
s = open(P, encoding="utf-8", newline="").read()


def sub1(old, new):
    global s
    n = s.count(old)
    if n != 1:
        raise SystemExit("pattern found %d times: %r" % (n, old[:70]))
    s = s.replace(old, new)


# ---- header -------------------------------------------------------------
sub1("""page did not survive and are corrected with their constructions.
""",
"""page did not survive and are corrected with their constructions.

Punch round 2, after the "Round 2 re-check" in the same review file.
One BLOCKER, E, which round 1's own `hypothesis()` fix caused in
`frmtmb.sample`, and five nits; what changed is marked **R3** below.
""")

# ---- nit 2: credit the refusals, not the binding -----------------------
sub1("""| every owner absent, frmtmb's own methods not frmtmb's | 10 of 26 | 10 of 26 | 0 of 26 |
""",
"""| every owner absent, frmtmb's own methods not frmtmb's **(R3: the refusals' win, not the binding's)** | 10 of 26 | 10 of 26 | 0 of 26 |
""")

sub1("""`dev/generics-out2/blocks/notice.txt` for the `hypothesis()` row, and
the collision test below for the rest.
""",
"""`dev/generics-out2/blocks/notice.txt` for the `hypothesis()` row, and
the collision test below for the rest.

**R3. Two rows in that table are not the binding's.** The
owner-absent row, 10 of 26 to 0 of 26, is BLOCKER B's ten refusal
methods: review registered those same ten on the SWAP build in-session
and got 0 of 26 there too (`dev/genrev-r2-attrib.R`). It holds on
either mechanism and says nothing about which one ships. The rows that
ARE the binding's are the stale owner, the non-generic owner, the
`frmtmb.sample` re-exports with brms only loaded, and the NOTE.
""")

# ---- nit 3: the 22 against 23 reason ------------------------------------
sub1("""**On 23 against the review's 22.** The review measured the SWAP
re-export residue at 22 of 26; this lane measures 23. The difference is
exactly one name and it is the BLOCKER C method difference again. Of
the 23, two are answered by a `.default` rather than by nothing,
`loo_compare` and `posterior_summary`. The review's probe counts a name
as lost only when NO method answers, so it does not count
`posterior_summary`; its name list also carries `expose_functions`
where this one carries `loo_compare`, which are both lost and cancel.
23 minus `posterior_summary` is 22. Both measurements agree on the
thing that matters: ACTIVE is 0 of 26 by either count.
""",
"""**On 23 against the review's 22. R3: the round-2 reason was wrong
and is withdrawn here.** Round 2 said the review's probe did not count
`posterior_summary` because a `.default` answered it. It did count it:
`posterior_summary` is on the review's own lost list
(`dev/genrev-out/p-sampleprop.txt`). The one-name difference is in the
NAME LISTS, not the method. The review's list omits `loo_compare` and
carries `expose_functions` in its place, and `frmtmb.sample` does not
re-export `expose_functions` (checked against its NAMESPACE), so with
only `frmtmb.sample` attached that name is not on the search path and
cannot be lost there. The review's 22 are exactly this lane's 23
without `loo_compare`. The count stood either way, and ACTIVE is 0 of
26 by both.
""")

# ---- BLOCKER E, its own section ------------------------------------------
sub1("""**What it costs.** One closure call per ACCESS instead of one per
session, so the cost model changed and was re-measured:
""",
"""### BLOCKER E. Round 1's `hypothesis()` fix removed the note from draws

**R3.** Round 1 moved the arming of the reserved-name shadowing note
out of core's `hypothesis()` generic and into core's two methods,
`hypothesis.frmtmb_fit()` and `hypothesis.frmtmb_multiple()`. That was
right for core. But there was a THIRD method that had relied on the
generic's arming: `frmtmb.sample`'s `hypothesis.frmtmb_draws()`, which
calls `hyp_env_vals()`, never armed the note, and could not, because
`hyp_shadow_arm()` was not on core's internal export list. So on draws
the note was gone in EVERY session, brms loaded or not, and no test in
either package covered it. Review measured it with
`dev/genrev-r2-drawsnote.R`: BASE 1 note on the fit and 1 on its draws,
the round-1 tree 1 and 0.

**The reason my guard missed it is precise, and it generalizes.** The
"no shared generic carries work in its own body" test stops work
coming BACK into a generic. It cannot see work that LEFT the generic
and did not reach every method. A guard that checks where work is must
also check where work went.

The fix, as review described it:

* `hyp_shadow_arm()` and `hyp_shadow_disarm()` join core's internal
  export list on `?frmtmb-sampling-api`, beside the `hyp_*` helpers
  `frmtmb.sample` already uses, with a paragraph saying every
  `hypothesis()` method must arm the note itself and why the generic
  cannot. `frmtmb.sample` uses `import(frmtmb)`, so it needs no
  NAMESPACE change of its own.
* `hypothesis.frmtmb_draws()` arms on entry and restores on exit, the
  same two lines as core's methods.
* Two tests in `frmtmb.sample`'s `test-draws-methods.R`, both seen
  failing first. One asserts the note appears exactly once on draws,
  with the same count on the fit as a control, so a construction that
  shadows nothing cannot pass it. It uses `fake_draws()`, because the
  note is emitted while the hypothesis is parsed against the fit and
  before any draw is read, so it needs no Stan build. The other is the
  WHERE-THE-WORK-WENT guard: every `hypothesis.*` method in both
  packages' namespaces must call the arm and disarm pair.

<!--BLOCK:blockere-->

Read the rows in order. On BASE the note test passes, because the
generic armed; the guard fails on all three methods, because on BASE
none of them armed, which is correct for a guard about where the work
lives. With core fixed and the draws method not, both tests fail and
the guard names exactly `hypothesis.frmtmb_draws`. With both fixed,
both pass, and the review's own probe on the real sampler gives 1 note
on the fit and 1 on draws, with and without brms.

`hypothesis` was the only borrowed generic with work in its body at
the base commit, 1 of 28 (`dev/genrev-r2-basebodies.R`), so this one
method is the whole blast radius.

**What it costs.** One closure call per ACCESS instead of one per
session, so the cost model changed and was re-measured:
""")

# ---- load cost: settled -----------------------------------------------------
sub1("""* **The two new hard Imports cost +0.0168 s.** That is the BASE arm
  with `loadNamespace("nlme")` and `loadNamespace("generics")` added,
  minus BASE. The review's headline was +0.0167 s on the swap build and
  its own Imports arm +0.0183 s; an earlier run of this script, on a
  busier machine, gave +0.0158 s. Three runs, one figure: about 17 ms,
  and it is the Imports.
""",
"""* **The two new hard Imports cost +0.0168 s.** That is the BASE arm
  with `loadNamespace("nlme")` and `loadNamespace("generics")` added,
  minus BASE. An earlier run of this script, on a busier machine, gave
  +0.0158 s.
* **R3, the disagreement is settled.** Review re-measured with 80
  interleaved rounds, arm order shuffled within each round, a
  high-resolution clock and paired bootstrap intervals
  (`dev/genrev-r2-loadcost.R`): IMPORTS minus BASE +0.0213 s
  [+0.0164, +0.0251], ACTIVE minus SWAP -0.0008 s [-0.0074, +0.0032],
  null control -0.0000 s. It withdrew its round-1 figure. That +0.0167 s
  was SWAP minus BASE on the minima of 30 runs, and it falls below the
  5th percentile of the same statistic resampled: **a low draw, not the
  cost of the whole change.** The finding on this page stands as it
  was: the two Imports are the cost, the two designs cannot be told
  apart at load, and the whole-change headline is a range.
""")

# ---- nit 1: truncation in the BLOCKER D bullets -------------------------------
sub1("""  "response")` reproduces `fitted()` exactly. Under truncation
  `fitted()` is the TRUNCATED mean and the formula is wrong, by a
  relative 0.1150 at the worst row on a `y | trunc(lb = 2000)` fit of
  this design. The review reports 0.4548 on its own construction; the
  magnitude depends on how much of the distribution is cut off, the
  finding is the same either way, and both are far above the machine
  precision the untruncated identity holds to.
""",
"""  "response")` reproduces `fitted()` exactly. Under truncation
  `fitted()` is the TRUNCATED mean and the formula answers a different
  question.

  **R3, and round 2's number was mislabeled.** Round 2 called 0.1150
  "the relative error at the worst row". It is `max |diff| /
  max(fitted)`; the worst row is 0.5163. And it is not a measurement at
  all but an IDENTITY: with `a = (log(lb) - mu) / sigma`, the per-row
  relative shortfall is `1 - pnorm(-a) / pnorm(sigma - a)`, because the
  truncated lognormal mean is `exp(mu + sigma^2 / 2) * pnorm(sigma - a)
  / pnorm(-a)`. It has no single value; it depends on where the bound
  sits in each row's distribution. The page now gives the formula and
  no number, and the test asserts the identity row by row to machine
  precision instead of a threshold on its size. That assertion is seen
  failing under mutants M3 and M6, and M1 fails the same block on its
  untruncated half. (Review's earlier 0.4548 was on its own
  construction and is superseded by this.)

<!--BLOCK:trunc-->
""")

# ---- nit 5: defect 3, narrowed --------------------------------------------------
sub1("""That is R's own `registerS3methods()` resolving
`S3method(posterior::as_draws, ...)`, which the base commit also
carries. In practice an older or partial `posterior` earlier on the
library path breaks frmtmb entirely rather than degrading. Not this
lane's, not fixed.
""",
"""That is R's own `registerS3methods()` resolving
`S3method(posterior::as_draws, ...)`, which the base commit also
carries. In practice an older or partial `posterior` earlier on the
library path breaks frmtmb entirely rather than degrading.

**R3, narrowed: "identical on BASE" is true of THAT construction and
not in general.** Review built a `posterior` that exports every generic
frmtmb registers on except `as_draws_rvars`, loaded before frmtmb. BASE
loads, SWAP loads, and this lane's build FAILS. Reproduced here:

<!--BLOCK:partialposterior-->

The cause is not the active binding. It is BLOCKER B's new
`S3method(posterior::as_draws_rvars, frmtmb_fit)` refusal, which is a
name no earlier NAMESPACE asked `posterior` for; `as_draws_list` is new
the same way. So this lane DID widen the exposure, by two names. The
realistic exposure is small: `posterior` 1.0.0, its first CRAN release,
already exports `as_draws_rvars` (the v1.0.0 tag's NAMESPACE, fetched
from the stan-dev repository), and the 1.7.0 installed here does. What
breaks is a partial or pre-CRAN install. `posterior (>= 1.0.0)` is now
declared in core's Suggests. That states the requirement and lets
`R CMD check` and installers see it; it does NOT stop a partial
`posterior` already on a library path from breaking `loadNamespace`,
which only removing the directives could do. The base-build failure
for a `posterior` lacking `as_draws` itself is unchanged and still not
this lane's.
""")

# ---- what was run -----------------------------------------------------------------
sub1("""The core suite below ran against the build with the `hypothesis()` fix
in; no source changed after it.
""",
"""**R3, punch round 2.** Round 2 changed core in comments, in two names
added to the internal export list with their `@aliases` and a
paragraph on `?frmtmb-sampling-api`, in one Suggests floor and in
prose, and changed `frmtmb.sample`'s code in one method and its tests.
As instructed, `R CMD check` was re-run on `frmtmb.sample`, whose code
changed, on a quiet machine with no other R process and without
`--no-manual`: `Status: OK`, tests OK. It was NOT re-run on core; the
core block below is the round-1 resumed run. Core's NAMESPACE did
change, by two exports, so the checks that change could break were run
directly on the installed package instead: `tools::undoc()`,
`tools::codoc()`, `tools::checkDocFiles()` and `tools::checkS3methods()`
all report nothing, on core and on `frmtmb.sample`. Both new Rd
changes were verified by rendering (`dev/generics-rdcheck.R`, and
`Rd2txt` on `?frmtmb-sampling-api`).

The core test files round 2 can reach, one process each:

<!--BLOCK:r2files-->

The full core suite below ran against the build with the round-1
`hypothesis()` fix in; round 2 changed no core code path it exercises
beyond the files just listed. The `frmtmb.sample` suite below is round
2's.
""")

# ---- version --------------------------------------------------------------------
sub1("""- A downstream package that registers a method on any of these generics
  must match the OWNER's formals. `frmtmb.sample` is edited here for
  that reason, so its `frmtmb (>= 0.55.1)` floor has to move to the
  bumped core version.
""",
"""- A downstream package that registers a method on any of these generics
  must match the OWNER's formals. `frmtmb.sample` is edited here for
  that reason, so its `frmtmb (>= 0.55.1)` floor has to move to the
  bumped core version.
- **R3.** The floor is now a hard requirement, not only a formals one:
  `frmtmb.sample`'s draws method calls `hyp_shadow_arm()` and
  `hyp_shadow_disarm()`, which no released core exports. Against the
  current release it would fail at the first `hypothesis()` on draws.
- **R3.** Core's Suggests now declares `posterior (>= 1.0.0)`; see
  defect 3. That is a dependency floor, not a version choice for either
  package.
""")

# ---- scripts table ------------------------------------------------------------------
sub1("""| `dev/generics-rdcheck.R` | the scale page verified by rendering it |
""",
"""| `dev/generics-rdcheck.R` | the scale page verified by rendering it |
| `dev/generics-trunc.R` | **R3** the truncation identity, swept over three bounds |
| `dev/generics-edit39.py` to `dev/generics-edit43.py` | **R3** BLOCKER E's test, its fix, the guard, the nits, this page |
""")

open(P, "w", encoding="utf-8", newline="\n").write(s)
print("ok")

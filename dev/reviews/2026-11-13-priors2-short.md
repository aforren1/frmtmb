# Short check of lane priors2's fix round (round of 2026-11-13, base ec8178d3)

Checks the fixes for F1 to F5 and L1 to L3 of
`2026-11-13-priors2.md` only. Lane library `wt-priors2-lib` (the build
of 15:19, after the `frm_stop()` change), base `rellib-r21`. Scripts and
outputs: `dev/reviews/priors2-short-*`. No git writes, no lane file
edited, at most 3 R processes (2 after the machine notice). Free space
on C: 145.4 GB at start, 145.0 GB at end. The scratch merge tree (item
5) was in the session scratchpad and is deleted.

The lane ran one `git diff --stat` in its fix round. That refreshes the
index's stat cache only; it does not write worktree files. The changed
file set here comes from `git hash-object` against `git ls-tree
ec8178d3` (content, not stat): 8 R files, 1 new R file
(`R/fit-open-ends.R`), 2 NEWS files, 2 sample R files, 5 sample test
files, 26 core test files (2 new).

## Verdict: FIX

F2, F4, F5, L1, L2, L3 and the merge hold. F1 is still not a proof
(case H), F3's "refused by name" is false on the user's path (`lb = 2`
on `zi` is dropped without a word), and the box-first design leaks a raw
optimizer warning on fits that now converge.

## 1. F2, the closed box first: PASS

- 33-fit survey, rerun on the lane library
  (`priors2-short-survey-lane.tsv`), compared in hex with
  `dev/priors2-survey-hex.R` against `priors2-out/survey-before.tsv`
  (rellib-r21): 29 of 33 bitwise identical
  (`priors2-short-survey-hex.txt`). The 4 others are as the lane says:
  `zi_unif_on` errored on 0.82.0, and the three `_push` fits.
- The review's 15 fits (`priors2-review-bitwise.R lane`): `cmp` with
  `priors2-review-bitwise-base1.tsv` is identical.
- 5 new fits with no bounds on fresh data, seed 23
  (`priors2-short-bitwise.R`: gaussian, glmm `(1 + x | g)`, 2-component
  mixture, cumulative, nonlinear with nlpar priors), one process per
  arm: `cmp` of the two tsv files is identical (code, objective, logLik
  and every outer parameter in hex). The glmm ends at code 1 on both
  builds.

## 2. F1, the infinite-density refusal: FIX

The refusal is positional, not a proof. `open_end_rescue()` sets
`opt$open_unbounded` when the kept point is within `open_end_at()`'s
slack of an infinite-density end, and `check_map_unbounded()` refuses
on that flag alone. The slack is `sqrt(.Machine$double.eps) *
max(1, |lb|, |ub|)`, about 1.5e-8 in absolute terms at `lb = 0`. The
objective is not read. A finite maximum within 1.5e-8 of the end is
refused with "has no maximum here".

Case H (`priors2-short-f1.R`, seed 5, n = 80): `poisson("identity")`,
`y ~ 0 + x`, `x = U(0.5, 2) * scale`, S counts of 1, the rest 0,
`gamma(0.5, 1)` with `lb = 0`. The log posterior is `(S - 0.5) log b -
b (sum x + 1) + c`. It goes to -Inf at b = 0, and the maximum is
`b* = (S - 0.5) / (sum x + 1)`. The negative log posterior by hand at
scale 1e7, S = 1: -4.672 at b* = 4.854e-10, -1.718 at b*/1e3, 329.5 at
1e-300. So the posterior vanishes at the end, and the maximum is
finite.

| scale, S | b* | lane, `lb = 0` | lane, `lb = b*/1e4` | 0.82.0, `lb = 0` |
| --- | --- | --- | --- | --- |
| 1e5, 1 | 4.854e-08 | code 0, 4.853823305e-08 | code 0, 4.853823206e-08 | optimizer failed, NA/NaN gradient |
| 1e5, 3 | 2.427e-07 | code 0, 2.426911653e-07 | code 0, 2.42691164e-07 | the same error |
| 1e7, 1 | 4.854e-10 | refused: "has no maximum here" | code 0, 4.853823775e-10 | the same error |
| 1e7, 3 | 2.427e-09 | refused | code 0, 2.426911886e-09 | the same error |
| 1e9, 1 | 4.854e-12 | refused | code 0, 4.853823615e-12 | the same error |
| 1e9, 3 | 2.427e-11 | refused | code 0, 2.426911888e-11 | the same error |

`beta(0.5, 2)` with `lb = 0, ub = 1` gives the same pattern: fitted at
scale 1e5, refused at 1e7 and 1e9. 0.82.0 errored on all, so it is not a
regression, but the refusal states a falsehood. Fix: refuse only on
evidence that the objective diverges. For example: the kept objective
is -Inf, or the data part of the objective is finite at the end, or the
objective still falls along the open coordinate (evaluate it at `q =
end + (q_stop - end) / 2^k`, k = 1..m, and refuse only if it keeps
falling). Add case H at scale 1e7 to `test-prior-bounds-ends.R`.

Second defect in the same message: under `beta(0.5, 2)` with `lb = 0,
ub = 1` the refusal names "x (beta(0.5, 2), ub = 1)", but the infinite
end is `lb = 0`. `box_end_density()` writes `end_labels[[nm]]` once per
side, so the upper (zero-density) end overwrites the lower label.

Reruns (`priors2-short-probe-lane.txt`, `dev/priors2-fix-probe.R`):
the output is identical to `fix-probe-after.txt`. A: code 0,
120.129961301298, b = 1.383301921. B: code 0, 88.9031549557894,
0.6478533563. C: code 0, 29.5619832612003, slope 0.3022615629 (the
0.82.0 local mode; see the note below). D `y`: code 0,
60.0663831131304. D `yneg`: refused (correct: the likelihood is finite
at 0). G: code 0, 36.5172433865701.

Note for the user, not a FIX: C's posterior is unbounded at 0 (gaussian
likelihood finite there), yet C returns the interior local mode with
code 0 and nothing said, as 0.82.0 did. The first build refused it.

### Warning on a correct fit (fix-round regression): FIX

The box run now comes first, and its raw nlminb warning "NA/NaN
function evaluation" escapes on fits that the rescue then brings to
the mode with code 0: A, B, the three pushed E fits, and every H row
(`priors2-short-probe-lane.txt`, `priors2-short-f1-lane.txt`). On
0.82.0, E stopped at code 1 with the gradient warning, which was
correct for that stop. The first lane build gave no warnings (review,
item 6). `test-prior-bounds-ends.R` allows the warning (lines 27, 86),
so the tests hide it. Fix: muffle the box run's NA/NaN warnings when
`open_end_rescue()` keeps the refit, and replay them when it does not.

## 3. F3, `bound_inside()`: FIX (user path), helper PASS

`priors2-short-f3c.R lane` (frmtmb's own links, `get_link()`):

| link | ends | outside | most linkinv calls |
| --- | --- | --- | --- |
| logit, `1 - 10^-seq(2, 5, .01)` | 602 | 0 | 23 |
| logit, `10^-seq(2, 12, .01)` | 2002 | 0 | 3 |
| logit, review ub, `1 - 1e-15`, `2^-1070`, `1e-300` | 8 | 0 | 98 |
| logit, `runif(5000)` | 10000 | 0 | 38 |
| probit, `runif(5000)` | 10000 | 0 | 6 |
| probit, `1 - 10^-seq(2, 15, .01)` | 2602 | 0 | 1 |
| cloglog, `runif(5000)` | 10000 | 0 | 2 |
| log, `exp(U(-700, 700))` | 10000 | 0 | 1 |
| log, `1 + k eps`, k in -2000:2000 | 8002 | 0 | 1 |
| sd scale (exp), `exp(U(-50, 50))` | 10000 | 0 | |

Termination: the doubling ends by `step` overflowing (at most about
2100 steps from 2^-1074) and the bisection by adjacent doubles. The
worst case seen is 98 calls. `bound_inside(0, 2, plogis, "lb")` refuses
by name at once.

Edges through `internal_bound()`: logit `ub = 1` gives Inf (no box);
logit `lb = 0` and log `lb = 0` or `-1` give -Inf (no box, brms's `lb =
0` on a log-linked dispersion, handled by name in the comment); logit
`lb = 1`, `ub = 0`, `ub = 2`, `ub = -1` and log `ub = 0` are refused
("outside the support ... empty box").

Defect: **logit `lb = 2` gives -Inf**. `qlogis(2)` is NaN, and
`internal_bound()` reads every NaN lower bound as "below the support"
before `bound_inside()` runs. On the user's path, `frm(cnt ~ z, family
= zero_inflated_poisson(), prior = set_prior("", class = "zi", lb =
2))` returns code 0, objective 148.316241056629, which is the fit with
no prior (`priors2-short-f3c-lane.txt`). It does the same on 0.82.0
(`priors2-short-f3-base.txt`), but the fix round and NEWS say "a bound
that no value maps inside is refused by name". The lane's test calls
`bound_inside()` directly, so it does not see this. Fix: in
`internal_bound()`, a NaN lower bound is "constrains nothing" only
where it is below the support (for example
`isTRUE(linkinv(-.Machine$double.xmax) >= v)`); above it, refuse as an
empty box. Test it through `frm()`.

Side note, no action: with `stats::make.link()` links, whose inverses
clamp at 2.2e-16, `bound_inside()` refuses `ub` below 2.2e-16
(`priors2-short-f3b-lane.txt`, 2363 of 10000 log ends). frmtmb's links
do not clamp, and `zi ub = 1e-20` and `sigma ub = 1e-20` through
`frm()` match 0.82.0, so no user path reaches this.

## 4. NEWS, compat, L1 to L3: PASS

- Core NEWS (development version): no "Through" and no version
  number. BREAKING marks on the `thetacov` rename, class `theta`,
  `cor` with `resp`, the centered-intercept bound, the
  infinite-density refusal and the step-start message. The refusal
  bullet says what remains ("where that second run goes to an end at
  which the prior's density is infinite"). The bullet "A bound that no
  value maps inside is refused by name" is false for `lb = 2` on `zi`
  (item 3).
- Sample NEWS: no version number; BREAKING on `thetacov`. Floor in
  words only; DESCRIPTION still `frmtmb (>= 0.82.0)` (consolidation
  must raise it, as the lane says).
- F5: `R/compat.R:2167` reads `class = \"thetacov\", coef =
  \"thetaac_1\"`.
- L1 (`priors2-short-l1.txt`): `set_prior()` with class `theta` and
  coef `theta_4`, `thetaac_1` or `thetar_2`, with a density, a bound or
  neither, is refused with "Write class = "thetacov", coef = ...". A row
  edited into a prior list carries the `thetacov` pointer through the
  dirichlet refusal. `set_prior("", class = "theta", coef = "a")`
  prints `theta_a ~ (flat)`.
- L2: the refusal ends "frm_sample(), which takes the bound on the
  class row, not on one coefficient" (`R/fit-open-ends.R:286`; D `yneg`
  output).
- L3: filed LOW under "### From lane priors2" in `dev/test-backlog.md`,
  and said in the `open_end_flags()` doc.

## 5. Merge with famend: PASS, one conflict

Scratch tree (session scratchpad, deleted): the lane's R/, with
`R/priors.R` and `R/fit.R` three-way merged with famend's current files
(`git merge-file`, base `ec8178d3`, CRLF to LF). `R/priors.R` 0
conflicts, `R/fit.R` 0 conflicts; both `prior_start_inside()` and
`open_end_rescue()` are present. `test-prior-table-brms.R`: 1 conflict,
the at0 block (priors2's lines vs famend's caller-start block). It was
resolved with the merge note's block (caller start, code 0, objective
within 1e-8 of m8).

Do not take famend's whole file, as the merge note says: priors2 also
changed the `cor`/`resp` test in that file (base line 180). famend's
copy keeps the base assertion that `cor` with `resp = "y2"` is
accepted, which priors2 now refuses. The three-way merge keeps both.

`pkgload::load_all()` of the scratch tree, then the merged file
(`priors2-short-merge.R`, `priors2-short-merge.txt`): pass 105, fail 0,
error 0, skip 0, warn 0.

## 6. Test files on the lane library: PASS

One file per process (`priors2-short-tests.sh`,
`priors2-short-tests-out/`). Every `lib:` line names
`wt-priors2-lib`.

| file | pass | fail | err | skip | warn |
| --- | --- | --- | --- | --- | --- |
| test-prior-bounds-ends.R | 66 | 0 | 0 | 0 | 0 |
| test-prior-mixture-theta.R | 55 | 0 | 0 | 0 | 0 |
| test-prior-table-brms.R | 101 | 0 | 0 | 0 | 0 |
| test-mixture-starts.R | 135 | 0 | 0 | 0 | 0 |
| test-conditions-census.R | 11 | 0 | 0 | 0 | 0 |
| sample test-default-priors-autocor.R | 101 | 0 | 0 | 0 | 0 |
| sample test-default-priors-dpar.R | 64 | 0 | 0 | 0 | 0 |
| sample test-flat-prior.R | 132 | 0 | 0 | 0 | 0 |
| sample test-me-prior-draws.R | 41 | 0 | 0 | 0 | 0 |
| sample test-prior-draws.R | 95 | 0 | 0 | 0 | 0 |

OpenBLAS 0.3.32 with LAPACK, `OPENBLAS_CORETYPE=Haswell`,
`OPENBLAS_NUM_THREADS=4`, run with the copy that already exists in
`dev/ciharden-out/Rob0.3.32-lapack` (not rebuilt):
`test-prior-bounds-ends.R` 66 and `test-prior-mixture-theta.R` 55, with
0 fail, error and warn. The copy has no corename query, so the kernel
it selected is not verified.

`dev/ciharden-out/` (0.3.26 and 0.3.32 builds, the review's and the
lane's) is 905.6 MB. It was left in place; delete it before
consolidation.

## What the lane must do

1. F1: make the refusal a proof (case H at scale 1e7 must fit at b* =
   4.854e-10). Fix the end label for a density with two open ends.
2. Muffle the box run's "NA/NaN function evaluation" when the rescue
   is kept. Drop the `allow_warnings()` for it in
   `test-prior-bounds-ends.R`.
3. F3: refuse a NaN lower bound above the support in
   `internal_bound()` (`zi lb = 2`). Test it through `frm()`.

## Re-check after fix round 2

Lane library build of 17:21. Scripts `dev/reviews/priors2-short-r2*.R`,
outputs `priors2-short-r2*`. At most 2 R processes. Free space on C:
146.6 GB at the end. `dev/ciharden-out/` is gone (the lane deleted it).

### Verdict: FIX (one item, or the user waives it)

Items 1, 3, 4, 5 and 6 pass. The probe's refusal is a proof in every
case built here (item 2). One case remains: a posterior shown to be
unbounded is returned at an interior local mode with code 0 and no
word.

### 1. Case H, D, the label: PASS

`dev/priors2-fix-probe.R` on the lane (`priors2-short-r2-probe-lane.txt`)
is identical to the lane's `fix2-probe-after.txt`.

- H, `gamma(0.5, 1)`: code 0 at b* on all six rows. For example,
  scale 1e7, S 1 gives 4.853823773e-10 (b* 4.853823771e-10), and
  scale 1e9, S 3 gives 2.426911934e-11 (b* 2.426911888e-11).
- H, `beta(0.5, 2)`: code 0 at b* on all six rows.
- D `yneg`: refused. D' (`beta(0.5, 2)`, lb 0, ub 1) is refused, and
  the message names "z (beta(0.5, 2), lb = 0)".
- No NA/NaN warning on A, B, E or H.

H at 1e7 and 1e9 still warns "Standard errors are not available ... a
bound holds it". That is the SE check's positional slack. The lane
filed it LOW.

### 2. Is the refusal a proof? PASS, with one silent case (FIX)

P1 (`priors2-short-r2.R`, seed 9, n = 5): exponential("identity"),
`y = rexp(5) * scale`, `gamma(0.5, 1)`, lb = 0. The objective `(n +
0.5) log b + S/b + b` falls for many decades and turns up only at b* =
S / 5.5.

| scale | b* | lane | 0.82.0 |
| --- | --- | --- | --- |
| 1e-4 | 1.026526439e-4 | code 0, 1.026507281e-4 | optimizer failed |
| 1e-9 | 1.026526439e-9 | code 0, 1.026526445e-9 | the same |
| 1e-12 | 1.026526439e-12 | code 0, 1.026526457e-12 | the same |
| 1e-30 | 1.026526439e-30 | code 0, 1.02652644e-30 | the same |
| 1e-100 | 1.026526439e-100 | code 0, 1.026526439e-100 | the same |
| 1e-300 | 1.026526439e-300 | code 0, 1.026526131e-300 (gradient warning) | the same |

At every scale the box run or the refit reached the turn, so the
probe has nothing to refuse. The grid steps two decades at a time down
to 2^-1074. A turn is missed only if it lies within the last two
decades above 2^-1074, which these data cannot produce.

P2: gaussian slope that the data put below 0 (likelihood finite at 0),
`gamma(a, 1)`, lb = 0. The posterior is unbounded.

- a = 0.5, 0.9 or 0.99: refused. The fall to 2^-1074 is (1 - a) * ~730,
  which is more than the margin of 1. This is correct.
- a = 0.999, 0.9999 and `beta(0.999, 2)`: the fall is less than 1. The
  fit is kept at z = 2^-1074, code 1, with "The fit stopped at the end
  of the support of a prior whose density is infinite there ... not a
  proven maximum (code 1)". This is acceptable: no fit of a diverging
  posterior comes back silent at code 0 along this path.

**Silent case (FIX).** P3, `y ~ x + z` (seed 1, n = 150), with
`gamma(0.5, 1)` and lb = 0 on x (gaussian, so the likelihood is finite
at 0 and the posterior is unbounded):

- 0.82.0: code 1, objective -Inf, x = 0. The box run reached the pole.
- Lane: code 0 at the interior local mode 0.3929088467, objective
  163.395461346473, with no warning.

The lane's own row D `y` (code 0, 60.0663831131304) is the same case.
C is the same posterior, but its box run never reached the pole. In
all three the MAP does not exist, and the fit says nothing. 0.82.0 at
least reported -Inf for P3 and D `y`.

Recommended fix: where the fit ends inside and the prior's density is
infinite at a box end, evaluate the objective at that end with the
other parameters held at the fit (the probe's last point). If it is
-Inf, or below the fit by the probe's margin, keep the fit and warn:
the posterior is unbounded toward that end, and the point is a local
mode. That claim is a proof, so the warning is correct. Otherwise, the
user must accept a silent local mode as the documented behavior.

### 3. The held NA/NaN warning: PASS

- Unit run (`priors2-short-r2c.R lane`): `open_end_hold()` with an open
  lower end holds the warning. `release(open_refit = FALSE)` raises it
  again ("NA/NaN function evaluation", class `frmtmb_warning`; 0.82.0
  raised `simpleWarning`). `release(open_refit = TRUE)` raises nothing.
  With no open end, the warning passes straight through.
- Real fits (`priors2-short-r2b.R`): poisson("identity") and
  exponential("identity") with lognormal, `gamma(2, 2)`, `inv_gamma(3,
  1)` or `beta(2, 2)` at lb = 0. On 0.82.0, 16 of 20 ended at code 1
  (objective Inf) or failed, most with the raw NA/NaN warning. On the
  lane, all 20 give code 0 with no warnings. The four with a start of
  1e-3 are bitwise identical on both builds, with no warning on
  either.
- I found no real fit whose box run raises NA/NaN beside an open end
  and then returns its own result. The 24 interior fits of
  `priors2-short-r2c.R` (identity links, intercept near 0, start (1, 1))
  raise no warning on either build. They are bitwise identical in hex
  between the builds.
- Interior fits beside an open end (P3: lognormal, `gamma(2, 2)`,
  `inv_gamma(3, 1)`, `beta(2, 2)` on x, review's C, and
  poisson("identity") with no prior and with `uniform(0, 5)`): hex
  and warnings are identical to 0.82.0. The exceptions are P3's
  `gamma(0.5, 1)` row (the silent case above) and poisson("identity")
  with `lognormal(0, 1)`, which goes from code 1 at Inf to code 0 at
  1.389530246.

### 4. Bounds above the range: PASS

- Through `frm()`, `zi lb = 2` is refused: "lb = 2 is outside the
  support ... empty box". `zi lb = -1` is bitwise the fit with no
  prior. `zi lb = 1` is refused, as before.
- `internal_bound()` with `get_link()`: logit, probit and cloglog
  refuse `lb = 2` and `1.5`, and give -Inf for `lb = -1` and `0`.

Side finding, LOW, and present on 0.82.0 too: frmtmb's cauchit
`linkfun` is `tan()`-based, so it returns finite values outside (0,
1). `lb = -0.5` maps to b = 1.2e-16, which is a box at probability 0.5.
`ub = 1.5` and `ub = 2.5` map to a box at 0.5 too. `zero_inflated_poisson(link_zi
= "cauchit")` does not exist, so I found no user path that reaches it.
Check `v` against the inverse link's range before `linkfun`, or file
it.

### 5. Bitwise: PASS

- The 33-fit survey (`priors2-short-r2-survey-hex.txt`): 29 of 33
  bitwise identical to rellib-r21. The other 4 are `zi_unif_on` and the
  three `_push` fits.
- The review's 15 fits: `cmp` is identical.
- The short review's 5 fits: `cmp` is identical.

### 6. Tests on the lane library: PASS

`priors2-short-r2-tests-out/`. Every `lib:` line names `wt-priors2-lib`.

- `test-prior-bounds-ends.R`: pass 80.
- `test-prior-mixture-theta.R`: pass 55.
- `test-prior-table-brms.R`: pass 101.

Fail, error, skip and warn are 0 on all three.

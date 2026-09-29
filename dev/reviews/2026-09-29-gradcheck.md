# Review of lane wt-gradcheck: the two-stage convergence criterion

Reviewer's brief was to FALSIFY. Read beside
`dev/gradcheck-findings.md`.

Everything below was measured on this box with R 4.6.1. "lane" is
`C:/Users/adf44/source/r/wt-gradcheck-lib`; "base" is the read-only
0.64.0 reference at `C:/Users/adf44/source/r/rellib-r3`. Every reviewer
script is `dev/gradcheck-rev-*.R` and every reviewer log is under
`dev/gradcheck-rev-log/`. Each relation below is marked IDENTITY (true
by construction, checked with `identical()` or by reading the code) or
MEASUREMENT.

## 0. Is the worker's build the worktree?

`dev/gradcheck-rev-00-env.R`. The deparsed body of every function the
diff touches (`fit_outer_box`, `grad_bound_active`, `grad_headroom`,
`grad_verdict`, `grad_warning_msg`, `check_convergence`, `diagnose`,
`fit_assembled`, `autoscale_prefit`) is `identical()` between the lane
library and a fresh `sys.source()` of the worktree `R/`: 9 of 9.
IDENTITY. No rebuild was needed.

One mismatch that does not affect R code: the library was built at
23:10:58 and `man/diagnose.Rd` and `man/frmtmb_control.Rd` were
regenerated at 23:36:14, so the INSTALLED help is older than the
worktree's. Every Rd figure below is rendered from the worktree source,
which is newer than every R source and is what a release would ship.

## 1. The two-stage criterion

### 1.0 The structural fact that bounds the whole risk

IDENTITY, from the code. `check_convergence()` reaches
`grad_verdict()` only inside

    } else if (is.finite(g) && g > control$grad_tol) {

where `g` is `max(abs(gvec * (fit$par_units %||% 1)))`: character for
character the old test, on the same quantity, against the same
threshold. Inside the verdict, `proj <= gmax` always, because the
projection only zeroes components, and the warning needs `proj > tol`
too. Therefore

    warn_lane(fit)  implies  warn_base(fit)

on every fit, with no exception. The lane build CANNOT raise a gradient
warning base does not.

This settles task 1(c) as posed. "An unusable Hessian counts as a
confirmation" cannot manufacture a new false alarm, because the
unchanged gradient gate stands in front of it. The only failure mode it
has is failing to REMOVE a base false alarm, which is harmless.

MEASUREMENT consistent with that identity: lane-only gradient warnings
observed over 180 fits of 9 new designs, 15 rows of the worker's own
true-positive set, 17 rows of my miss probe and 17 model shapes:
**0**.

**Verdict: holds.**

### 1(a) The KKT sign test

`dev/gradcheck-rev-02-kkt.R` (log `kkt-lane.txt`) and
`dev/gradcheck-rev-03-probe.R` (log `probe-lane.txt`). All MEASUREMENT.

| construction | seed | numbers | verdict |
|---|---|---|---|
| A2 `ub = 0.1` on `b` binds, `sigma` interior | 101 | `max_grad` 121.25, `grad_proj` 4.860e-07, `held = [x]`, silent | holds |
| A3 `ub` on `x1` AND `lb` on `x2`, both bind | 9103 | `max_grad` 147.09, `grad_proj` 2.583e-07, `held = [x1,x2]`, silent | holds |
| A4 `lb = 2` on class `sd`, stored as a log | 9104 | box lower on `theta_1` is 0.69314718 = log 2; `par[theta_1]` is 0.69314718; gradient +39.12 at a LOWER bound; `held = [theta_1]`; `grad_proj` 2.110e-04; silent | holds |
| A5 autoscaled column, spread 1e-4, no bound | 9105 | `par_units` = (1, 9947.68, 1); silent | see 1(a)-units |
| A6 autoscaled column AND `ub = 100` on `b` | 9105 | `max_grad` 0.06836, `grad_proj` 1.898e-04, `held = [xs]`, silent | holds |
| B1 a parameter ON a bound whose gradient points INTO the box | 9201 | below | holds |

A4 is the piece worth naming: the class-`sd` bound is spelled on the sd
scale and `resolve_bounds()` plus the class transform put it into the
box as `log(2)`, in the same space as `fit$opt$par`. The KKT test
therefore compares like with like on the log scale. It is not tricked.

B1 is the case the brief asks for and the one the worker's test file
only unit-tests through `grad_bound_active()`. No ordinary fit reaches
it: with `ub = 0.1` and a true slope of -2 the optimizer stops at
`x = -2.0210112`, well inside. So the fit object was poked onto the
bound (`opt$par[x] <- 0.1`, cache rebuilt from `fit_outer_box()`) and
the verdict recomputed. There the gradient is **+932.662**, which at an
UPPER bound points down into the feasible set.
`grad_bound_active()` returns `FALSE FALSE FALSE`, `grad_proj` is
1978.19, the headroom is `NA` and `warn` is `TRUE`. The true shortfall
at that point is 989.09307 log-likelihood units. The sign test does not
swallow it.

A1 was the reachable attempt at the same thing (`lb = 0.1`, start AT
the bound, `iter.max = 1`): nlminb steps off the bound to
`x = 1.0107075`, so it does not exercise the sign, but it does warn,
251.66332 units short, with "the curvature there is unusable".

**Verdict: holds on every piece the brief names.**

### 1(a)-units: which scale the test runs on, and a falsified findings claim

Sound part. `grad_bound_active()` receives the RAW gradient and the RAW
parameter vector. `par_units` is nlminb's `scale` argument and not a
reparameterization (`R/fit.R:1836`), autoscale fits a standardized COPY
and maps back (`R/autoscale.R:1-7`), and the class transforms put every
bound in the internal space. `par_units > 0`, so the sign test is
invariant to it. Correct.

Unsound part, and it is a documentation defect.
`diagnose()$max_grad` is the RAW largest component while
`diagnose()$grad_proj` is the largest component in `par_units` units.
On an autoscaled fit those are different scales, and the printed block
puts them on adjacent lines
(`dev/gradcheck-rev-03-probe.R` row B2, seed 9202,
`par_units` = (1, 1.01016e+05, 1.02608, 1)):

    Max |gradient|: 1001 at z
      1 parameter is held by a bound (z); the largest gradient no bound
      holds is 1.685

`dev/gradcheck-findings.md` section 7 says "the new `grad_proj` is in
natural units and the difference is documented on the manual page". It
is not. The rendered `man/diagnose.Rd` (lines 50-60,
`dev/gradcheck-rev-log/rd.txt`) says only that `grad_proj` "is the
largest component over the rest", and the words `par_units`,
`autoscale` and "units" appear nowhere near it. **The findings claim is
falsified.**

The same mismatch makes `expect_identical(d$grad_proj, d$max_grad)` in
the test block "with the bound gone, nothing is excluded" pass only
because that fit is not autoscaled and its argmax has
`par_units == 1`. It is not the identity its comment asserts.

### 1(b) The Newton decrement on a Laplace model

`dev/gradcheck-rev-12-hess.R`, log `hess-lane.txt`, seed 8501, poisson
`y ~ x + (1 | g)`, 50 groups, 1000 rows, 3 outer parameters, 50 inner.

- `fit$obj$he(p)` on a Laplace object **ERRORS**: "Hessian not yet
  implemented for models with random effects." So every random-effects
  fit takes the fallback. IDENTITY (the error is deterministic).
- The fallback is `stats::optimHess()` differencing `fit$obj$gr`, the
  MARGINAL gradient, so the Hessian is the marginal one. Checked
  against a hand-written central difference of the same marginal
  gradient at step 1e-5: max relative gap **7.683e-08**. MEASUREMENT.
- The headroom by the two routes: `grad_headroom()` 4.844167884e-13,
  the central-difference Hessian 4.844167468e-13. Agreement to 8
  significant figures. MEASUREMENT, two routes.

Direct cost of the stage. Blocks grown past 1.2 s, minimum of 3 rounds,
`proc.time()` resolution measured at 10.0 ms in the same process.

| fit | np | one marginal gradient | `grad_headroom()` | whole fit | stage / fit |
|---|---|---|---|---|---|
| poisson `(1 \| g)`, 5000 groups, n = 20000 | 3 | 0.00719 s | 0.18750 s | 1.5500 s | 0.121 |
| poisson `(x \| g)`, 5000 groups, n = 20000 | 5 | 0.01289 s | 0.45000 s | 2.5500 s | 0.177 |
| gaussian `s(x) + z`, n = 20000 | 5 | 0.02547 s | 0.44250 s | 5.1000 s | 0.087 |
| gaussian `gp(t, k = 40)`, n = 700 | 4 | 0.00922 s | 0.13187 s | 1.1600 s | 0.114 |
| gaussian `(1 + x1..x5 \| g)` us block, 200 groups, n = 4000 | 29 | 0.00412 s | 0.39500 s | 2.0400 s | 0.194 |

Both 5000-random-effect fits TRIP the trip-wire (`gmax` 7.429e-02 and
9.956e-03) and are correct (headroom 1.269e-07 and 2.186e-08), so they
do pay the stage, and what they pay is 12.1 and 17.7 percent of one
fit. That is the "a large-n GLMM that trips pays it" case the brief
asks to see. Across five large-latent designs the stage is 8.7 to 19.4
percent of a fit. The worker's own worst case (0.726 at np = 63)
remains the worst recorded; nothing here exceeds it.

`gp(t)` with no `k=` is refused by the package ("gp() without k= builds
a dense 600-point covariance; use k= for the Hilbert-space
approximation"), so a dense-covariance gp cannot be built to test.

A DISSOLVED measurement, recorded rather than deleted. A second
instrument (`dev/gradcheck-rev-05-cost.R`, log `cost-lane.txt`) timed
whole fits at `grad_tol = 1e-3` against `grad_tol` above the fit's own
gradient, interleaved, with a control arm. Its control read 0.992 and
1.048 in one run and 0.927 and 1.032 in another, while the effect read
0.880 to 1.167. At 3 rounds on this loaded box that instrument cannot
resolve the stage. Use the direct table above.

**Verdict: holds.**

### 1(c) Can an unusable Hessian create a new false alarm?

No, and not by luck: section 1.0.

Empirically, the two benign-indefinite constructions the brief names do
not even reach stage 2. `B4a`, a variance component truly at zero (60
groups, 2400 rows): `max_grad` 5.408e-06, under `grad_tol`, `pdHess`
TRUE, silent. `B4b`, `student()` on gaussian data at n = 20000 with
`nu` running off its link: `max_grad` 5.189e-04, under `grad_tol`,
`pdHess` TRUE, silent. Over the 180-fit extended grid of section 3, no
fit warned with `headroom = NA`.

The one construction where an unusable Hessian decided the verdict is
`M3d mixture rel.tol = 1e-2` (`dev/gradcheck-rev-09-miss.R`, seed 505):
`headroom = NA`, the fit warns, and it really is 0.281479
log-likelihood units short. A correct call.

**Verdict: holds; the hazard is unreachable by construction.**

## 2. Bitwise identity of fits

`dev/gradcheck-rev-04-ident.R` and `-07-ident2.R`, compared by
`-06-cmp.R` with `identical()` on `coef`, `opt$par`, `logLik`,
`obj$fn(par)`, `opt$objective`, `opt$convergence`, `max|grad|` and the
FULL gradient vector. IDENTITY claims, tested with `identical()` and
not with a printed comparison.

13 model shapes ran in both arms and are bitwise identical on all 8
fields: a gaussian mixture, a collapsed mixture, two-response
`set_rescor(TRUE)`, `cumulative` with `thres(gr = grp)`, an `mi()`
two-part model, `gp(t)`, `s(x)`, `REML = TRUE`,
`zero_inflated_poisson`, complete separation, `cumulative` with a rare
top category, `student()`, and a 2000-group poisson GLMM.

Four more shapes in `-07-ident2.R`: a nonlinear `nl = TRUE` decay, the
same fit stopped short, an importance-corrected GLMM
(`importance = 1000L`) and a 300-group importance fit. Objectives agree
to all 15 printed digits on both builds.

Base raised 3 gradient warnings the lane build does not, all on correct
fits: rescor at `gmax` 1.069e-03, ordinal `thres(gr=)` at 1.009e-03,
the 2000-group GLMM at 6.805e-03. Lane-only warnings: 0.

The worker's own 15-row true-positive set was rerun through their own
script on both arms (`dev/gradcheck-06-truepos2.R`, logs
`repro-truepos2-lane.txt`, `repro-truepos2-base.txt`). Every objective
column is byte-identical between arms; base warned on 12 rows, lane on
7.

**Verdict: holds.**

## 3. False alarms and misses

### 3.1 An extended design set

`dev/gradcheck-rev-08-extend.R`, logs `extend-lane.txt`,
`extend-base.txt`. Nine designs the worker did not have, 2 sample sizes,
10 replicates, one seed per cell and none shared: 180 rows per build,
**180 distinct seeds, 0 fit errors** on either build. Counts generated
by the script, pasted verbatim:

```
       design     n  reps  base  lane  conv!=0        max gmax
           gp  2000    10     0     0        0    2.358316e-04
           gp 20000    10     0     0        0    2.018227e-04
 mix_collapse  2000    10     2     0        0    1.889428e-03
 mix_collapse 20000    10     8     0        0    1.916065e-02
     near_sep  2000    10     0     0        0    5.696882e-05
     near_sep 20000    10     0     0        0    3.639846e-05
     rare_top  2000    10     0     0        0    2.295882e-04
     rare_top 20000    10     4     0        0    1.694572e-03
         reml  2000    10     1     0        0    2.851078e-03
         reml 20000    10     5     0        0    2.458426e-03
   separation  2000    10     0     0       10   2.131861e-243
   separation 20000    10     0     0       10   5.911833e-228
       smooth  2000    10     0     0        0    9.847308e-04
       smooth 20000    10     1     0        0    1.429894e-03
   student_nu  2000    10     0     0        0    7.567016e-04
   student_nu 20000    10     2     0        0    1.834089e-03
    zero_infl  2000    10     0     0        0    9.430324e-04
    zero_infl 20000    10     5     0        0    8.941436e-03
TOTAL                 180    28     0       20
```

base 28 / 180 = 0.156; lane 0 / 180 = 0. The 20 `convergence != 0` rows
are the two separation cells, where the MLE is at infinity; both builds
speak there through the optimizer's own status, not the gradient.

`frmtmb.sample` draws objects do NOT route through the new code.
`frm_sample()` on a formula takes `fit_assembled(objective_only = TRUE)`,
which returns at `R/fit.R:987`, before `bounds` is resolved (991) and
before `check_convergence()` (1114). `frm_sample()` on a fit reuses a
fit whose verdict is already cached. Read from the code; confirmed by
`test-draws-methods.R` and `test-stan-control.R` being identical
between builds (section 4).

**Verdict: holds, extended.**

### 3.2 The misses

`dev/gradcheck-rev-09-miss.R`, logs `miss-lane.txt`, `miss-base.txt`.
17 rows, each with a reference fit of the same data optimized hard
(`rel.tol = x.tol = 1e-14`, 8 restarts). Base warned on 14, lane on 9.
Every warning base raises and lane does not is on a fit within
**8.24e-05** log-likelihood units of its reference:

| row | shortfall | base | lane |
|---|---|---|---|
| poisson `rel.tol = 1e-5` | 8.23811e-05 | warns | silent |
| poisson `rel.tol = 1e-6` | 1.73227e-05 | warns | silent |
| collinear `eps = 1e-8` | 2.01089e-09 | warns | silent |
| near-separation `rel.tol = 1e-6` | 3.90061e-07 | warns | silent |
| nonlinear decay from a far start | 1.25155e-10 | warns | silent |

The headroom's calibration is a MEASUREMENT with two independent sides.
shortfall / headroom, over every row where both are finite and the
reference is trustworthy: 1.001, 1.000, 0.9997, 1.000, 1.000 (the
poisson `rel.tol` sweep), 0.9984 (collinear 1e-5), 1.054, 0.9997,
0.9999 (near separation), 0.9956 (nonlinear), 0.6503 (mixture). Worst
overstatement 1.54x, on the mixture.

The worker's 2.256e-04 miss reproduces EXACTLY at their seed 402:
base warns, lane is silent, shortfall 0.000225669
(`repro-truepos2-*.txt`).

**Could not construct** a short fit with a LARGER shortfall that the
new criterion misses through a near-singular Hessian, after nine
attempts over four mechanisms: a collinear ridge at four
conditionings, near separation at three tolerances, a nonlinear plateau
from two far starts, and a collapsing mixture at two tolerances. Every
one either warns or is inside 8.3e-05 units. Absence of a construction
is not proof; the attempts and their scripts are recorded so the next
reader does not redo them.

A RETRACTION, made on the page as the rules require. My first pass
called the worker's section 5 collinear rows wrong: at `eps = 1e-6` my
hard reference recovers only 1.68536e-08 units while the headroom reads
0.529293, a factor of 3e7. **I was wrong; the worker is right.**
`dev/gradcheck-rev-13-step.R` settles it by taking the step, which needs
no optimizer:

```
eps=1e-06  rcond(H)=2.881e-13  |step|=42585  predicted=0.529293
           best drop=0.528027 at t=1  ratio=0.9976
   objective       : 1363.65765091862  -> 1363.12962372381
   independent nll : 1363.65765091862  -> 1363.12962372381  drop 0.528027
```

The "independent nll" line is a hand-written gaussian log likelihood
over `cbind(1, x1, x2)` computed outside RTMB, and it agrees with the
tape to 15 digits at BOTH points, so the 0.528 is not a cancellation
artifact of the tape. The lower point exists, so the fit is genuinely
short and the reference optimization failed to find it. The
"shortfall" column on my M3a rows at `eps` 1e-6 and 1e-7 is the
optimizer's failure, not the truth.

### 3.3 The trip-wire's gap, and the "never less diagnostic" rule

Reproduced in kind at my own seed. `dev/gradcheck-rev-09-miss.R` row
M1, seed 501: `y ~ xs` with `sd(xs) = 1e-7` and `autoscale = FALSE`
stops with `convergence = 0`, `max_grad` 1.601e-05, and
**297.175** log-likelihood units below the autoscaled fit of the same
data, with NO warning of any kind on EITHER build.

So the brief's question has a clean answer: the old build did not warn
there either, the new check is not less diagnostic at that point, and
the rule is not broken by the code. The worker's 1609.6 figure is a
property of their design; the mechanism reproduces. Reaching it needs
`autoscale = FALSE`, a non-default.

The rule IS broken by the manual page. Section 7.1.

**Verdict: holds.**

## 4. The extension suites, which the worker did not run

Every extension test file mentioning `diagnose`, `max_grad`, `grad`,
`convergence` or `check_convergence` (grep over
`C:/Users/adf44/source/r/frmtmb/extensions/*/tests`) was run against
the lane core with the extensions from `rellib-r3`, and against the
reference core. One file per R process, `NOT_CRAN=true` and
`FRMTMB_SCALE_TESTS=true`, runner `dev/gradcheck-rev-runfile.R`.

The test environment's parent is the extension's own namespace, as
`test_check()` does. Without that, four files report artificial errors
on their own internals; that is what the first-pass `ext-*` logs show
and why they were rerun as `e2-*`. The `e2-*` figures are the ones to
read.

| file | package | lane | base |
|---|---|---|---|
| test-cross-pairs.R | coupling | pass=52 fail=0 error=0 skip=0 | pass=52 fail=0 error=0 skip=0 |
| test-contaminant.R | eam | pass=32 fail=0 error=0 skip=0 | pass=32 fail=0 error=0 skip=0 |
| test-ndt-bound.R | eam | pass=79 fail=0 error=0 skip=0 | pass=79 fail=0 error=0 skip=0 |
| test-units.R | eam | pass=32 fail=0 error=0 skip=0 | pass=32 fail=0 error=0 skip=0 |
| test-gddm-reference.R | eam | pass=166 fail=0 error=0 skip=0 | pass=166 fail=0 error=0 skip=0 |
| test-hmm-starts.R | latent | pass=84 fail=0 error=0 skip=0 | pass=84 fail=0 error=0 skip=0 |
| test-hmm.R | latent | pass=118 fail=0 error=0 skip=0 | pass=118 fail=0 error=0 skip=0 |
| test-lca.R | latent | pass=105 fail=0 error=0 skip=0 | pass=105 fail=0 error=0 skip=0 |
| test-rlddm-ndt.R | learn | pass=53 fail=0 error=0 skip=0 | pass=53 fail=0 error=0 skip=0 |
| test-ode.R | ode | pass=70 fail=0 error=1 skip=0 | pass=70 fail=0 error=1 skip=0 |
| test-draws-methods.R | sample | pass=148 fail=0 error=0 skip=0 | pass=148 fail=0 error=0 skip=0 |
| test-stan-control.R | sample | pass=47 fail=0 error=0 skip=0 | pass=47 fail=0 error=0 skip=0 |
| test-difference.R | spline | pass=59 fail=0 error=0 skip=0 | pass=59 fail=0 error=0 skip=0 |

The single `test-ode.R` error is identical in both arms and is an
artifact of running outside `pkgload` ("No packages loaded with
pkgload" from `local_mocked_bindings`), not a finding about either
build.

Every one of those 13 logs is BYTE-IDENTICAL between the two arms after
the header lines, so the messages agree and not only the counts.

`frmtmb.eam`'s bound-hitting fits are exactly what `test-ndt-bound.R`
(the `ndt` bound) and `test-contaminant.R` (the collapsing contaminant
`lambda`) cover; both are unchanged. `test-contaminant.R` captures the
whole `diagnose()` printout and matches with `expect_match`, so the two
new printed lines do not break it.

SCALE TIER. The seven `test-scale.R` files call `diagnose()$max_grad`
on every row through `helper-scale.R`, which is the read the brief is
worried about. Run with `FRMTMB_SCALE_TESTS=true`, two at a time:

| file | package | lane | base |
|---|---|---|---|
| test-scale.R | coupling | pass=6 fail=0 error=0 skip=0 | pass=6 fail=0 error=0 skip=0 |
| test-scale.R | eam | pass=18 fail=0 error=0 skip=0 | pass=18 fail=0 error=0 skip=0 |
| test-scale.R | latent | pass=6 fail=0 error=0 skip=0 | pass=6 fail=0 error=0 skip=0 |
| test-scale.R | learn | pass=14 fail=0 error=0 skip=0 | pass=14 fail=0 error=0 skip=0 |
| test-scale.R | sample | pass=3 fail=0 error=0 skip=0 | pass=3 fail=0 error=0 skip=0 |
| test-scale.R | spline | pass=2 fail=0 error=0 skip=0 | pass=2 fail=0 error=0 skip=0 |
| test-scale.R | ode | pass=2 fail=0 error=0 skip=0 | pass=2 fail=0 error=0 skip=0 |

Seven of seven agree exactly, and all seven logs are byte-identical
between the arms after the header lines. The `frmtmb.ode` row took
about 52 minutes per arm (one `test_that` block, an `RTMBode` fit),
which is why an earlier draft of this section recorded it as
unfinished. It landed before the review was filed.

Two rows are worth quoting because they are fits the change affects.
`frmtmb.learn`:
`diag=conv=0,maxgrad=0.00257,pdHess=TRUE,nbadse=0,nflat=0` at 20000
rows, which base warns on and lane does not. `frmtmb.ode`:
`diag=conv=1,maxgrad=0.00235,pdHess=TRUE,nbadse=0,nflat=0`, where base
raises TWO warnings and lane raises only the optimizer's own. Both
assertions pass on both builds, because they read the string rather
than the warning.

The FUZZ TIER, which the brief did not name and which is the sharpest
risk in the change. `helper-fuzz.R:1761` MUTES every metamorphic
invariant on a fit whose warnings match
`"Large maximum absolute gradient"`. Removing 40 percent of those
warnings therefore puts MORE fits under the invariants than before, and
any latent invariant failure among them would go red on the lane build
while base was green. Run with `FRMTMB_FUZZ=true`
(`dev/gradcheck-rev-17-fuzz.R`):

| plan size | lane | base |
|---|---|---|
| 120 | pass=2 fail=0 error=0 skip=0 | pass=2 fail=0 error=0 skip=0 |
| 300 (the file's own default) | pass=2 fail=0 error=0 skip=0 | pass=2 fail=0 error=0 skip=0 |

No invariant fires from the unmuting. The warning's first clause is
byte-preserved, so the pattern still matches at all.

**Verdict: holds. 20 extension test files, both builds, every log
byte-identical between the arms.**

## 5. `R/confint.R` and `R/autoscale.R`

`R/autoscale.R`: two comment edits, no code. Verified by `git diff`.

`R/confint.R`: changes are confined to `diagnose()`. Three new return
fields, two new printed lines, one roxygen paragraph, and the
clean-line condition moving from a hardcoded `out$max_grad < 1e-3` to
the fit's own verdict. `confint()` itself is untouched.

`dev/gradcheck-rev-11-profile.R`, compared by `-15-cmp2.R` (log
`profile-cmp.txt`). On a gaussian GLM, a 40-group poisson GLMM and a
bounded gaussian fit, `identical()` is TRUE for all 13 saved objects:
`confint()` Wald, `confint(parm =, method = "profile")`, `vcov()`, the
sdreport `pdHess` flag, and every pre-existing `diagnose()` field
(`convergence`, `message`, `max_grad`, `worst_grad`, `pdHess`,
`bad_se`, `flat`, `singular`, `separation`, `unbounded_dpar`,
`predictor_scale`, `nonfinite_trials`, `scale`). IDENTITY claims,
tested with `identical()`.

The clean-line change is a real behaviour change and it is in the
direction the worker says. Base's condition read the RAW `max_grad`
against a hardcoded `1e-3`, ignoring `frmtmb_control(grad_tol =)`
entirely, so base could print "No convergence problems detected" on a
fit whose own warning had just fired. `dev/gradcheck-rev-03-probe.R`
row B3 swept `grad_tol` over 1e-8, 1e-3 and 1e-1 on a correct gaussian
fit at n = 20000 and the lane build prints the clean line in all three,
consistently with its own silence.

**Verdict: holds.**

## 6. Tests

`dev/gradcheck-rev-runfile.R`, one file per process, counts taken from
the testthat results object so an aborted file cannot read clean.

- `test-grad-verdict.R` on the reference build: **pass=23 fail=10
  error=7**. Reproduces the worker's figure exactly.
- On the lane build: **pass=46 fail=0 error=0 skip=0**. Reproduces.

Of the 10 base failures, 7 are behavioural and 3 are a new field
reading `NULL`; the 7 errors are `NULL` arithmetic inside `expect_lt` /
`expect_gt`. The behavioural ones, seen to fail on the unfixed code:

- `expect_false(r$grad)` returns TRUE on the bounded fit, on the
  gaussian GLM at n = 50000, and on the ordinal fit at n = 6000.
- `diagnose()` prints neither "held by a bound" nor "the largest
  gradient no bound holds".
- the warning carries neither "curvature there is unusable" nor
  "in log-likelihood".

Absolute numeric tolerances: ONE.
`expect_lt(abs(d$grad_headroom / realized - 1), 0.05)` in "the headroom
predicts what a Newton step actually gains". 0.05 is a bare constant,
not a ratio to anything the run measures, and the file's own header
claims "No absolute tolerance appears below". Every other numeric
factor (1e3, 100, 1, 1e-3) is a margin against `grad_tol`, a standard
error, or a quantity measured twice by two routes. NIT.

Also a NIT: "the bounded fit really is the constrained optimum" hard
codes `0.1` twice rather than reading `cap`, so changing the fixture's
default silently breaks the block.

### The block pinned to FLIP

`dev/gradcheck-rev-10-flip.R`, log `flip-lane.txt`. The block is
`expect_false(diagnose(r$fit, quiet = TRUE)$pdHess)` with a comment
saying the covariance machinery is not bound-aware and to flip the
block when it is. Constructing the case where that guarded condition is
ABSENT: the same data and seed, a different cap.

```
ub=0.1   x=0.1   pdHess=FALSE held=[x]  ev 693.40  73.10 -29.63  -> pass
ub=0.5   x=0.5   pdHess=FALSE held=[x]  ev 715.95 102.86 -26.23  -> pass
ub=1     x=1     pdHess=TRUE  held=[x]  ev 738.69 164.41   4.70  -> FAIL
ub=1.5   x=1.5   pdHess=TRUE  held=[x]  ev 705.0  256.8  118.8   -> FAIL
ub=1.8   x=1.8   pdHess=TRUE  held=[x]  ev 631.5  305.3  233.6   -> FAIL
ub=1.9   x=1.9   pdHess=TRUE  held=[x]  ev 609.2  314.2  262.7   -> FAIL
ub=1.95  x=1.95  pdHess=TRUE  held=[x]  ev 602.5  316.6  271.2   -> FAIL
ub=1.99  x=1.99  pdHess=TRUE  held=[x]  ev 600.1  317.4  274.1   -> FAIL
```

In 6 of 8 caps the bound still holds `x`, the covariance machinery is
still not bound-aware, and the assertion FAILS. So the block does NOT
pin bound-awareness; it pins that this one design at `cap = 0.1` has an
indefinite UNCONSTRAINED Hessian. The assertion is fine and the defect
record is fine; the comment overstates what it guards, which matters
because a future reader will flip it for the wrong reason. NIT, with
the construction.

That table is also the clearest single piece of evidence FOR the
change: at `ub >= 1` a bounded fit is now fully clean (`pdHess` TRUE,
`bad_se` empty, no gradient warning) where base warned.

**Verdict: the counts reproduce; the FLIP block's comment is
falsified.**

## 7. Docs

`dev/gradcheck-rev-01-rd.R` renders both pages with `tools::Rd2txt`
(log `rd.txt`), as the rules require, rather than reading the source.

### 7.1 `man/frmtmb_control.Rd` asserts a guarantee this lane's own
### measurement contradicts

The rendered `grad_tol` entry ends:

    ... the gradient alone warned on 289 of them and the two readings
    together on none, while every fit stopped short of its optimum by
    more than 'grad_tol' still warned.

That clause is false as written. Counterexample, section 3.3: a fit
**297.175** log-likelihood units short of its optimum at
`grad_tol = 1e-3` with no warning of any kind
(`dev/gradcheck-rev-09-miss.R` row M1, seed 501). The worker's own
findings section 7 records the same shape at 1609.6 units, so this is
not a new discovery, only a sentence that contradicts the record it
cites. Read the other way, the clause is scoped to "720 CORRECT fits",
where it is vacuous because that grid holds no short fits. Either
reading misinforms.

The fix is one sentence: name the trip-wire's gap beside the clause, or
drop the clause. The `autoscale` default is what closes the gap in
practice and that belongs in the same sentence.

### 7.2 `man/diagnose.Rd` omits the units caveat the findings say it carries

Section 1(a)-units. `grad_proj` and `max_grad` are on different scales
under `autoscale`, the printed block puts them on adjacent lines, and
the page says nothing about it while the findings assert that it does.

### 7.3 `vignettes/diagnostics.Rmd` cites a measurement that was not made

The new paragraph reads:

    the same gaussian design warned on 0 of 20 replicates at 400 rows
    and 20 of 20 at 20,000.

The design grid is 200 / 1000 / 5000 / 20000 rows
(`dev/gradcheck-03-design.R:100`). There is no 400-row cell anywhere in
`dev/gradcheck-log/`; the only `400` in the lane's logs is the
`n = 400` bootstrap control row in `12-base.txt` and `12-lane.txt`. The
gaussian figure at 200 rows IS 0 of 20, so the sentence is right with
"200" and unsourced with "400". One character.

### What is right

The "Convergence problems" section is true of the new behaviour. Step 1
now lists the three readings; step 7 sends the reader to the
log-likelihood number rather than to a gradient band, which is the
right advice and is what the criterion makes possible. US English, no
em dashes, no emoji, no spaced hyphen standing in for one.

`NEWS.md` sits under `# frmtmb (development version)`, carries no
version number, says plainly which warning stops appearing, names the
three new `diagnose()` fields and the hardcoded `1e-3` it replaces, and
makes the bitwise-identity claim that section 2 confirms with
`identical()`.

Style on the whole diff: 0 added lines over 80 columns in `R/fit.R`,
`R/confint.R`, `R/autoscale.R`, `tests/testthat/test-grad-verdict.R`,
`NEWS.md` and `vignettes/diagnostics.Rmd`.

**Verdict: 7.1, 7.2 and 7.3 falsified. The rest holds.**

## 8. Cost of a refit loop

`dev/gradcheck-rev-14-boot.R`, log `boot-lane.txt`. The worker timed the
two BUILDS in two processes, which the instrument rules do not allow.
This times the two READINGS inside ONE process on the lane build:
`restarts = 0` in both arms so `grad_tol` governs only the warning
path, arm "off" sets `grad_tol` to ten times the fit's own gradient so
no Hessian is built, arms interleaved per round, minimum of 5 rounds,
and a second "off" arm as the control.

| design | with stage 2 | without | ratio | control |
|---|---|---|---|---|
| gaussian n=20000, `frm_bootstrap(nsim = 20)` | 2.860 s | 1.560 s | **1.833** | 1.064 |
| gaussian n=400, `frm_bootstrap(nsim = 20)` | 0.060 s | 0.070 s | 0.857 | 0.857 |

The n = 20000 row is a valid measurement: blocks of 1.5 to 3.1 s, well
past 1.2 s, control 1.064. The penalty is **1.833x**, not the 1.65x the
worker recorded. Direction and mechanism are theirs; the magnitude is
worse than recorded.

The n = 400 row is NOT a control. At 0.06 to 0.12 s it is 6 to 12 ticks
of a 10.0 ms clock, and it reads 0.857 in the effect arm AND in the
control arm, which is the instrument and not the code. The worker's
"the n = 400 row is the control ... it must report 1.00, and it does"
rests on 0.13 s against 0.13 s, 13 ticks, below resolution. The claim
about the code is right; the evidence offered for it is not evidence.
NIT, with the measurement.

**Verdict: holds in direction and mechanism; the recorded ratio
understates the cost and the offered control is under-resolved.**

## 9. Findings the worker did not record

### 9.1 The warning names one quantity and prints another, and its
### remedy pointer names the wrong parameter

`dev/gradcheck-rev-16-msg.R`, seed 9301: a poisson fit with
`set_prior("", class = "b", coef = "x", ub = 0.1)` that binds and a
loosened optimizer, so the free components really are short.

    Large maximum absolute gradient at the optimum (52.5) (1 component
    held by a bound excluded): one Newton step from there would still
    gain 0.131 in log-likelihood, more than frmtmb_control(grad_tol =
    0.001). ... diagnose() names the offending parameter

    Max |gradient|: 5971 at x
      1 parameter is held by a bound (x); the largest gradient no bound
      holds is 52.49

The maximum absolute gradient is 5971.25, not 52.5.
`grad_warning_msg()` formats `v$proj` inside a clause the worker kept
verbatim for grep-compatibility, so the sentence now names one quantity
and prints another. Worse, the warning promises that `diagnose()` names
the offending parameter, and `diagnose()` names `x`, which is the ONE
parameter the verdict deliberately excluded. The parameter the headroom
actually comes from is never named anywhere. Reachable at default
settings on any bounded fit that warns.

### 9.2 `diagnose()` runs the verdict on an importance-corrected fit,
### which the fit's own check refuses to judge

`dev/gradcheck-rev-09-miss.R` rows M4, seed 506, a 60-group poisson
GLMM with `importance = 500L`:

| rel.tol | max_grad | grad_proj | grad_headroom | clean line | fit warnings |
|---|---|---|---|---|---|
| 1e-3 | 0.8236 | 0.8236 | 0.002337 | FALSE | 0 |
| 1e-2 | 3.843 | 3.843 | 0.06561 | FALSE | 0 |
| 1e-1 | 6.231 | 6.231 | 0.2312 | FALSE | 0 |

`check_convergence()` takes the `!is.null(fit$importance)` branch and
says nothing about the gradient, for the reason the worker's section 8
gives: a Monte Carlo objective's gradient carries an O(N^-1/2) error, so
no gradient criterion of any shape applies. `diagnose()` calls
`grad_verdict()` unconditionally, so it now publishes a
`grad_headroom` computed from exactly that gradient and a Hessian
differenced from it, and its clean-line condition consults `gv$warn` on
such a fit.

This is NOT a regression: base withheld the clean line too, through
`max_grad < 1e-3`, and nothing warns on either build. It is a new
number the package's own reasoning calls inapplicable, plus a Hessian
built on every `diagnose()` call on such a fit. The consistent fix is
one guard: skip the verdict, or return `NA`, when
`!is.null(fit$importance)`, matching `check_convergence()`.

### 9.3 `diagnose()` builds a Hessian it did not build before

On any fit whose projected gradient exceeds `grad_tol`, `diagnose()`
now costs one Hessian. On a fit that already warned this is free (the
verdict is cached on `fit$cache`), but `fit_set_outer()` in
`R/brms-shapes.R` replaces the cache, so `diagnose()` on a perturbed
fit pays it, and so does every `diagnose()` inside a refit loop. Cost
per call is the `grad_headroom()` column of section 1(b): 0.13 to 0.45
s on the large-latent designs. NIT; worth one sentence in the findings
beside the fit-path cost, which is the only one the worker measured.

### 9.4 `fit$lower` and `fit$upper` are unreachable, so
### `fit_outer_box()`'s fallback always returns an open box

`frm()` has no `lower` or `upper` argument (`R/fit.R:537`), and every
internal `fit_assembled()` call passes `NULL` or `fit$lower`, which is
that same `NULL` (`R/fit.R:616`, `R/allfit.R:94`, `R/confint.R:1581`,
`R/autoscale.R:259`). So the `tryCatch(resolve_bounds(fit, fit$lower,
fit$upper))` branch can only ever yield the open box the line above it
already returns. Not a bug: it errs toward warning, which is the safe
direction, and it is what makes the perturbed-fit case the worker
recorded harmless. But the function's comment says "`fit$lower` and
`fit$upper` carry only the bounds a CALLER passed" and no caller can
pass any. NIT.

### 9.5 No cache-key hazard

Checked, not a finding, recorded so it is not rechecked. `fit$cache` is
an environment, so `$` does not partial-match, and the two new slot
names collide with nothing: the existing slots across `R/` are `sdr`,
`flat_pars`, `Vjoint`, `Qjoint`, `features`, `tbl`, `acc`,
`warned_nonfinite_cov` and `warned_singular_precision`. Nothing in `R/`
enumerates cache names. `test-bracket-access.R` passes identically on
both builds.

## 10. Core test files rerun

Every core file that can reach a convergence warning or that captures
`diagnose()` output, one per process, both builds. All IDENTICAL,
counts and messages.

| file | lane | base |
|---|---|---|
| test-grad-verdict.R | pass=46 fail=0 error=0 skip=0 | pass=23 fail=10 error=7 skip=0 |
| test-diagnostics-ux.R | pass=127 fail=0 error=0 skip=0 | pass=127 fail=0 error=0 skip=0 |
| test-predfix.R | pass=92 fail=0 error=0 skip=0 | pass=92 fail=0 error=0 skip=0 |
| test-review-v28.R | pass=74 fail=0 error=0 skip=0 | pass=74 fail=0 error=0 skip=0 |
| test-autoscale.R | pass=46 fail=0 error=0 skip=0 | pass=46 fail=0 error=0 skip=0 |
| test-conditions.R | pass=150 fail=0 error=0 skip=0 | pass=150 fail=0 error=0 skip=0 |
| test-verbose.R | pass=40 fail=0 error=0 skip=0 | pass=40 fail=0 error=0 skip=0 |
| test-thres.R | pass=67 fail=0 error=0 skip=0 | pass=67 fail=0 error=0 skip=0 |
| test-confint-anova.R | pass=35 fail=0 error=0 skip=0 | pass=35 fail=0 error=0 skip=0 |
| test-importance.R | pass=226 fail=0 error=0 skip=0 | pass=226 fail=0 error=0 skip=0 |
| test-bracket-access.R | pass=33 fail=0 error=0 skip=0 | pass=33 fail=0 error=0 skip=0 |
| test-fuzz.R (gated off) | pass=0 fail=0 error=0 skip=1 | pass=0 fail=0 error=0 skip=1 |
| test-fuzz.R (`FRMTMB_FUZZ=true`, size 300) | pass=2 fail=0 error=0 skip=0 | pass=2 fail=0 error=0 skip=0 |

`test-grad-verdict.R`'s base column is the seen-to-fail run, not a
regression.

## 11. What I could not settle

- The worker's whole-core-suite figure (PASS 11924, FAIL 0, ERROR 0,
  SKIP 152) was not reproduced. A 180-file run did not fit beside the
  falsification work; the 12 files a convergence warning can reach were
  rerun instead (section 10) and each matches the reference build.
- The worker's `R CMD check --as-cran` run (1 NOTE) was not repeated;
  the lane rules put the authoritative one on the consolidating
  session.
- No construction was found in which the new criterion misses a short
  fit by more than 8.24e-05 log-likelihood units through a
  near-singular Hessian, after nine attempts over four mechanisms
  (section 3.2). Absence of a construction is not proof.

## Verdict

NOT MERGEABLE. The criterion itself survives every attempt to break it,
and nothing in the R logic needs to change. What blocks is that three
pieces of user-facing text state things this lane's own measurements
contradict, and a fourth misdirects the user the warning exists to
help.

Blocking, in order:

1. **7.1** `man/frmtmb_control.Rd` tells the user that "every fit
   stopped short of its optimum by more than `grad_tol` still warned".
   A fit 297.175 log-likelihood units short does not warn
   (`dev/gradcheck-rev-09-miss.R` M1, seed 501), and the worker's own
   findings section 7 records 1609.6 units for the same mechanism. This
   is the sentence a user would rely on to decide whether silence means
   convergence. One sentence to fix.
2. **9.1** The warning prints the PROJECTED gradient inside a clause
   that says "maximum absolute gradient" (52.5 against an actual
   5971.25), and then tells the user "diagnose() names the offending
   parameter" while `diagnose()` names the parameter the verdict
   excluded. The parameter the headroom comes from is never named.
   `dev/gradcheck-rev-16-msg.R`, seed 9301, default settings.
3. **7.3** `vignettes/diagnostics.Rmd` cites "0 of 20 replicates at 400
   rows". The grid has no 400-row cell
   (`dev/gradcheck-03-design.R:100`); the figure belongs to 200 rows.
   One character.
4. **7.2** `man/diagnose.Rd` does not say that `grad_proj` and
   `max_grad` are on different scales under `autoscale`, while
   `dev/gradcheck-findings.md` section 7 asserts that it does. The
   printed block puts the two numbers on adjacent lines
   (1001 and 1.685 on one fit).

Nits, in order:

5. **9.2** `diagnose()` publishes a `grad_headroom` on an
   importance-corrected fit, which `check_convergence()` deliberately
   refuses to judge. One guard would make the two agree.
6. **6** The FLIP block's comment claims to pin bound-aware standard
   errors; six of eight caps on the same data falsify that reading
   (`dev/gradcheck-rev-10-flip.R`). Reword the comment.
7. **8** The recorded refit-loop ratio is 1.65x; a within-process
   instrument with a resolved control reads **1.833x**. The `n = 400`
   row offered as the control is 13 clock ticks and cannot serve as
   one.
8. **6** One absolute tolerance in `test-grad-verdict.R`
   (`abs(headroom / realized - 1) < 0.05`) against a header that says
   there are none, and `0.1` hard coded twice where `cap` is available.
9. **9.3** The new `diagnose()` Hessian (0.13 to 0.45 s on large-latent
   fits) is not recorded beside the fit-path cost.
10. **9.4** `fit_outer_box()`'s comment describes a caller-passed bound
    that no caller can pass.

---

# Re-check, 2026-09-29, after punch round 1

Rebuilt library at `C:/Users/adf44/source/r/wt-gradcheck-lib`, built
01:08:01, every R source 00:58 to 00:59 and both edited Rd files
01:07:53, so all of them predate the build. `dev/gradcheck-rev-00-env.R`
re-run: the deparsed body of all 9 touched functions is `identical()`
between the library and a fresh source of the worktree. Only what
changed is covered below.

## (a) Every number and name in the new warning and block

`dev/gradcheck-rev-18-recheck.R`, log `recheck-a-lane.txt`. The script
recomputes the five quantities from the fit's own tape, the box and base
R, without calling `grad_verdict()`, then compares with `identical()`.

Seed 9301, `ub = 0.1` on `x` binds and the free set is short:

| quantity | independent | `diagnose()` | `identical()` |
|---|---|---|---|
| raw max and its parameter | 5971.250447 at `x` | 5971.250447 at `x` | TRUE |
| bound-held set | `[x]` | `[x]` | TRUE |
| projected max and its parameter | 52.48866159 at `z` | 52.48866159 at `z` | TRUE |
| headroom | 0.1308330663 | 0.1308330663 | TRUE |

The warning now reads, verbatim:

    Large maximum absolute gradient at the optimum (5971 at x), of which
    1 component is held by a bound (x), leaving 52.5 at z as the largest
    no bound holds: one Newton step over the parameters no bound holds
    would still gain 0.131 in log-likelihood, more than
    frmtmb_control(grad_tol = 0.001). The fit may not have converged.
    diagnose() reports all three numbers; ...

Blocking finding 9.1 is **fixed**. The headline number is the raw
maximum with its parameter, the projection is named with ITS parameter,
and the closing clause no longer sends the reader to the excluded
parameter.

Seed 9202, autoscaled and `ub` on `z`: raw max 1000.889458 at `z`,
natural-unit max 1026.993245, `grad_proj` 1.685440903 at `xs`, all four
`identical()` TRUE. The printed block now carries the line

    This fit was standardized internally, so 'Max |gradient|' is on the
    raw scale and every number below it is in the per-parameter units
    the fit was judged in

Blocking finding 7.2 is **fixed in the output as well as on the page**.

Seed 9103, two bounds bind: `grad_proj_par` is `sigma_(Intercept)`, the
correct free argmax, and not `NA`. `diagnose()$grad_headroom` is `NA`
while my unconditional recomputation gives 5.303e-17; that is NOT a
disagreement, because `grad_proj` is 2.583e-07, under `grad_tol`, so the
package does not measure it and my script does regardless.

**(a) holds.**

## (b) The fuzz grep, on real warnings

`dev/gradcheck-rev-19-recheck2.R`. The worker's
`dev/gradcheck-15-fuzzgrep.R` feeds hand-built lists to
`grad_warning_msg()`, which tests the formatter rather than the package.
I built five REAL fits, captured the warning the package actually
raised, and matched it against `FUZZ_NONCONVERGENCE` read out of
`tests/testthat/helper-fuzz.R`:

| construction | shape reached | matched |
|---|---|---|
| plain, short, finite headroom | held=0 headroom=finite | TRUE |
| one bound held, short free set | held=1 headroom=finite | TRUE |
| far start, curvature unusable | held=0 headroom=NA | TRUE |
| two bounds held, short free set | held=2 headroom=finite | TRUE |
| bound plus far start | held=0 headroom=NA | TRUE |

4 of the 4 reachable message shapes (held 0 or more, headroom finite or
NA), every one matched. The ABSENT case: the same pattern against
"Model fitted successfully with no problems" returns FALSE, so the grep
is not passing because it matches everything.

The worker's fifth case, `no_names` (`gmax_par` NA), is a defensive
formatter branch that no fit reaches, so testing it through the
formatter is the only way to reach it, which is what they do. Their
script also re-runs green here.

**(b) holds, and is stronger than the worker's own assertion.**

## (c) The importance withholding, and one REGRESSION it introduced

The withholding itself works. On my M4 seed 506 fit at three optimizer
tolerances, `grad_proj`, `grad_proj_par` and `grad_headroom` are all
`NA`, `grad_bound_held` is empty, neither the "Newton step" nor the
"held by a bound" line is printed, `fit$cache$grad_verdict` is NULL (so
no Hessian is built, which also removes the cost), and `nwarn` is 0 on
both builds, so `check_convergence()` is unchanged there. The same
design with no `importance =` still prints the Newton line and warns, so
the withholding is about the branch and not about the design.

BUT the guard failed CLOSED, and this is a new false alarm.
`R/confint.R:1251` sets `gv <- NULL` on any importance fit, and line
1482 reads

    (degenerate || (!is.null(gv) && !isTRUE(gv$warn)))

With `gv` NULL and `degenerate` FALSE that conjunct is FALSE for EVERY
importance fit, so "No convergence problems detected" is withheld
unconditionally rather than only when the verdict would have
complained. Base gated the same line on `out$max_grad < 1e-3`, which a
well-converged importance fit passes.

Constructed, `dev/gradcheck-rev-20-impclean.R`, seed 9503, a 20-group
poisson GLMM with `importance = 2000L` and `se = TRUE`:

```
                              conv max_grad     pdHess nbadse under_tol CLEAN
lane  importance nd=2000 9503   0  0.00094754   TRUE   0      TRUE      FALSE
base  importance nd=2000 9503   0  0.00094754   TRUE   0      TRUE      TRUE
```

`singular`, `separation`, `unbounded_dpar` and `predictor_scale` are all
clean too, and the same data with no `importance =` prints the clean line
on both builds (`max_grad` 5.2669e-05). So a correct importance fit that
satisfies all seven remaining conjuncts can no longer say so, on any
importance fit at all. Two further rows (seeds 9501, 9502) sit above
`grad_tol` and withhold on both builds, which is why the fixture had to
be constructed rather than taken from the earlier probe.

The seven other conjuncts stay meaningful on an importance fit. The fix
is one token: `(degenerate || is.null(gv) || !isTRUE(gv$warn))`.

**(c) the withholding holds; the clean-line conjunct is a regression.**

## (d) The quadraticity yardstick

It is a MEASURED quantity, and I checked the law it rests on rather than
taking it. On the block's own fit (seed 402):

```
realized                 = 0.3824926653
half-step drop/realized  = 0.7500187394   (0.75 exactly for a quadratic)
quarter-step/realized    = 0.43751082     (0.4375 exactly for a quadratic)
nonquad (from the fit)   = 2.49858e-05
yard                     = 2.49858e-05    (the sqrt(eps) floor is NOT active)
headroom error           = 2.39689e-05
error / yard             = 0.959298       assertion (< 10) passes
```

The quarter-step ratio is a third measurement and it is not in the test:
it confirms the objective really follows `D(2t - t^2)` along the step,
so `0.75` is the exact analytic value for a quadratic and not a fudge
factor. `error / yard` is 0.959, which means the headroom's error IS the
objective's third-order term, exactly as the comment claims. Against the
worker's recorded range of 0.96 to 1.29 over six constructions, the
bound of 10 leaves 7.8x of headroom over the worst observed case. It is
a ratio to something the run measures, which is what the house rule
asks. Nit 8 is addressed.

Worth stating: the assertion is now very tight in absolute terms. It
admits a headroom wrong by a factor of only 1.00025 on this fit. That is
a property of the fit being nearly quadratic there rather than of the
construction, and it is discriminating rather than fragile: a headroom
wrong by 2x would read `error/yard` of about 4e4.

It does still fail on the reference build, but as an ERROR (`expect_lt`
on `NULL/realized`), which the lane rules call the weak form. The
behavioural pins in the file are the `expect_false(r$grad)` blocks and
the new units and importance blocks, which fail properly.

**(d) holds; it is a measured yardstick, not a disguised constant.**

## (e) Docs

Rendered with `tools::Rd2txt` (`dev/gradcheck-rev-01-rd.R`, log
`rd2.txt`).

- **7.1 fixed.** The guarantee clause is gone. `?frmtmb_control` now
  carries a "WHAT IT DOES NOT COVER" paragraph naming the mechanism (the
  first reading gates the second) and my construction verbatim:
  `autoscale = FALSE`, column spread 1e-7, convergence reported, maximum
  absolute gradient 1.6e-5, 297 units below the standardized fit, no
  warning of any kind, with the `autoscale` default as the remedy.
- **7.2 fixed.** `?diagnose` now says `max_grad` is raw while
  `grad_proj`, `grad_proj_par` and the warning are in the `par_units`
  the fit was judged in, "not comparable: 1001 against 1.685 on one
  measured example", that the printed block says so when it applies, and
  that `grad_headroom` is invariant to the scaling. The gradient
  paragraph is now "reported four ways" and states that `worst_grad` and
  `grad_proj_par` are DIFFERENT parameters on a bounded fit.
- **7.3 fixed.** The vignette reads "0 of 20 replicates at 200 rows",
  and gains a blind-spot paragraph carrying the same 297-unit
  construction.
- **Nit 9 fixed.** `?diagnose` documents the Hessian cost, including
  that a fit whose cache was replaced pays 0.13 to 0.45 s.
- **Nit 5 documented.** A "NOT ON AN IMPORTANCE-CORRECTED FIT"
  paragraph, pointing the reader at `fit$importance$grad`.
- NEWS states the message change ("no longer prints the projected
  gradient under a phrase that says maximum", "no longer sends the
  reader to the one parameter the verdict excluded"), the four
  `diagnose()` fields, the importance withholding and the not-covered
  gap. No version number, still under `# frmtmb (development version)`.

Style: 0 added lines over 80 columns in `R/fit.R`, `R/confint.R`,
`NEWS.md`, `vignettes/diagnostics.Rmd` or
`tests/testthat/test-grad-verdict.R`; no em dash, no emoji. Three
roxygen-GENERATED lines in `man/diagnose.Rd` are 81 to 82 columns, which
is roxygen's wrapping and not hand-written text.

**(e) holds. All three blocking documentation findings are fixed.**

## (f) Test files

| file | build | result |
|---|---|---|
| test-grad-verdict.R | lane | pass=63 fail=0 error=0 skip=0 |
| test-grad-verdict.R | base | pass=29 fail=18 error=8 skip=0 |
| test-diagnostics-ux.R | lane | pass=127 fail=0 error=0 skip=0 |
| test-autoscale.R | lane | pass=46 fail=0 error=0 skip=0 |
| test-importance.R | lane | pass=226 fail=0 error=0 skip=0 |
| test-verbose.R | lane | pass=40 fail=0 error=0 skip=0 |
| test-confint-anova.R | lane | pass=35 fail=0 error=0 skip=0 |
| test-fuzz.R, `FRMTMB_FUZZ=true`, `FRMTMB_FUZZ_N=300` | lane | pass=2 fail=0 error=0 skip=0 |

Both `test-grad-verdict.R` counts match what the coordinator reported.
The five other lane files are unchanged from the first pass. The gated
fuzz tier is still green at the default plan size, so the unmuting still
costs nothing.

Of the 18 base failures, the ones that matter are behavioural and new
this round: 4 in "an importance-corrected fit gets no verdict at all",
1 in "grad_proj is in the units the verdict uses, max_grad is raw", 1 in
"with the bound gone" on `grad_proj_par`, and 1 in the rewritten
covariance block on `grad_bound_held`. That block now MEASURES its
premise and `skip_if_not`s when the full Hessian is positive definite,
then asserts the free-set Hessian is not, so nit 6 is addressed: the
comment names my cap sweep and cites `dev/gradcheck-rev-10-flip.R`.

## (g) The refit cost, for the record

The worker remeasured interleaved with a control: 1.914 and 1.463 with
controls 1.022 and 0.989. My own within-process figure was 1.833 with
control 1.064. The two agree, and both are above the 1.65 first
recorded. Nit 7 is addressed.

## Re-check verdict

NOT MERGEABLE, on one finding, and it is new this round.

**Blocking. The importance clean-line conjunct fails closed.**
`R/confint.R:1482` reads
`(degenerate || (!is.null(gv) && !isTRUE(gv$warn)))`, and `gv` is NULL on
every importance fit, so `diagnose()` can never print "No convergence
problems detected" on one however well it converged. Constructed at seed
9503: `conv = 0`, `max_grad` 0.00094754 under `grad_tol`, `pdHess` TRUE,
no bad standard errors, all four remaining conjuncts clean, and base
prints the clean line where lane does not
(`dev/gradcheck-rev-20-impclean.R`). This is a regression against base
introduced by punch item 5, and it is the shape the standing rule warns
about: a check that fires on a correct model is worse than no check. One
token fixes it: `(degenerate || is.null(gv) || !isTRUE(gv$warn))`.

Everything else from the first pass is resolved. Blocking findings 7.1,
9.1, 7.3 and 7.2 are fixed and verified by construction; nits 5, 6, 7,
8, 9 and 10 are addressed. No new defect was found in the criterion, the
warning text, the fuzz grep, the yardstick or the documentation.

---

# Re-check 2, 2026-09-29, after punch round 2

Library rebuilt 01:50:25, `R/confint.R` 01:49:59 and
`tests/testthat/test-grad-verdict.R` 01:41:51, both predating the build.
All 9 touched function bodies `identical()` between the library and a
fresh source of the worktree (`dev/gradcheck-rev-00-env.R`). One finding
was open; it is closed.

## The spelling taken, and why it is the better one

The worker did not take my `is.null(gv) ||` suggestion, and was right not
to. That spelling would have called every importance fit clean on the
gradient, including one whose gradient is large. What landed instead is

    (degenerate ||
       if (is.null(gv)) {
         isTRUE(out$max_grad < (fit$control$grad_tol %||% 1e-3))
       } else {
         !isTRUE(gv$warn)
       })

the 0.64.0 gate with the setting honoured instead of hardcoded. Absence
of a verdict now decides neither way: it hands the question back to the
check that existed before this round.

## The regression is closed

`dev/gradcheck-rev-20-impclean.R` re-run on the rebuilt library, log
`rc2-impclean-lane.txt`, same three seeds:

| fit | max_grad | under grad_tol | lane CLEAN | base CLEAN |
|---|---|---|---|---|
| importance nd=2000 seed 9503 | 0.00094754 | TRUE | TRUE | TRUE |
| importance nd=2000 seed 9501 | 0.0034419 | FALSE | FALSE | FALSE |
| importance nd=4000 seed 9502 | 0.013354 | FALSE | FALSE | FALSE |

Seed 9503, the construction that carried the regression, now prints the
line on both builds. 9501 and 9502 withhold on both. `conv = 0`,
`pdHess` TRUE, `nbadse` 0 on all three, and the three matched
no-importance fits still print the line. Lane and base now agree on all
six rows, where before they disagreed on 9503.

The seen-to-fail evidence for the first of the two new test blocks is in
this review's own round-1 log, `impclean-lane.txt`: lane CLEAN = FALSE
against base CLEAN = TRUE at seed 9503, which is the behavioural failure
that block asserts against. The pre-fix library is gone, so that log is
the record.

## The absent case, constructed across the setting

`dev/gradcheck-rev-21-impgate.R`, log `rc2-impgate-lane.txt`. The
coordinator asked for one cell. A gate that ignored `grad_tol` would
give the same answer in every column, so the whole sweep is what
separates the two, and three fits with three different max absolute
gradients give three different flip points. Every cell also asserts the
verdict is still withheld (`grad_proj`, `grad_headroom` NA and
`grad_bound_held` empty) and that no other conjunct went unclean.

```
seed   max|grad|    gt=1e-04   gt=1e-03   gt=1e-02   gt=1e-01
9501   0.00344187   FALSE      FALSE      TRUE       TRUE
                    expected: FALSE      FALSE      TRUE       TRUE
9502   0.0133543    FALSE      FALSE      FALSE      TRUE
                    expected: FALSE      FALSE      FALSE      TRUE
9503   0.000947541  FALSE      TRUE       TRUE       TRUE
                    expected: FALSE      TRUE       TRUE       TRUE
```

12 of 12 cells follow `max_grad < grad_tol`, no cell was flagged for an
unclean conjunct, and the flip point moves with the fit: between 1e-3
and 1e-2 for 9501, between 1e-2 and 1e-1 for 9502, between 1e-4 and 1e-3
for 9503. The fallback honours the setting, in both directions, and the
verdict stayed withheld at all 12 tolerances. This is the claim's
complement as well as the claim: at 1e-4 even seed 9503 loses the line.

## Tests

| file | build | result |
|---|---|---|
| test-grad-verdict.R | lane | pass=72 fail=0 error=0 skip=0 |
| test-grad-verdict.R | base | pass=38 fail=18 error=8 skip=0 |

Both match what the coordinator reported. The 18 base failures still
include the 4 behavioural ones in "an importance-corrected fit gets no
verdict at all" and the units block.

The two new blocks are built both ways round, which is what the standing
rule asks. "an importance fit under grad_tol still prints the clean
line" asserts each of the other seven conjuncts separately and
`skip_if_not`s on the premise, so it cannot pass because something else
was wrong. "an importance fit over grad_tol does not print the clean
line" is the absent case and carries a `hit` flag, so it skips rather
than passing vacuously if no draw lands above the tolerance.

Worth recording plainly: both new blocks PASS on the reference build,
because base's gate is `max_grad < 1e-3` and `grad_tol` is 1e-3 by
default, so base cannot distinguish the two spellings. That is correct
rather than a weakness. The regression was lane-only, so the pre-fix
lane library is the only build that could show it, and it did
(pass=71 fail=1, and independently this review's round-1 measurement).

## Re-check 2 verdict

MERGEABLE.

The one open finding is closed by a spelling better than the one I
proposed, verified on the construction that found it and on a
twelve-cell sweep that moves `grad_tol` across three fits in both
directions. Nothing new was found. Every blocking finding and every nit
from the first pass is resolved: the criterion, the bitwise identity of
fits, the extension and core suites, the fuzz tier, `confint()`, the
warning text, the two-scales reporting, the importance withholding, the
measured yardstick, both manual pages, the vignette and NEWS.

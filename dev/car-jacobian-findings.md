# Lane wt-car-jacobian: findings

Worktree `C:/Users/adf44/source/r/frmtmb-wt-car-jacobian`, branch
`wt-car-jacobian`, branched at 68d6782 (0.52.0). Main checkout not
touched. Private library
`<scratch>/cj-lib`, the worktree core installed into it.

Inputs read end to end: `dev/reviews/2026-09-05-esicar.md` (the
reviewer's two rounds, residual R2) and `dev/esicar-findings.md` (the
wt-esicar lane, section 8 and punch P2).

## 1. The defect, restated from the code

`expand_b()` (`R/covstruct.R:1781`) sends an esicar block through
`car_center()` (`:1307`), so the coefficient vector the linear
predictor multiplies is `c = P b` with

    P = I - sum_j (1 / n_j) s_j s_j'

the per-component centering projection. Every delta-method path pairs
the `Z` columns with `b` instead, through `dc/db = I`:

- `lp_delta_A()` `R/predict.R:1485`, the `add_b_cols` else branch
  `:1497` (`A <- cbind2(A, Zc)`), reached whenever
  `frame$has_rr` is FALSE, which is every esicar fit that has no `rr`
  block;
- `rr_jacobians()` `R/predict.R:637` fills IDENTITY entries for every
  non-`rr` block (`:661`), so an `esicar` + `rr` fit leaks the same way
  through the `has_rr` branch;
- `ranef(condVar = TRUE)` `R/methods-fit.R:604` reads
  `sdr$diag.cov.random` (`:611`), the conditional variance of `b`, and
  displays it against `cvec = P b`.

`b` splits orthogonally into the field `f = P b` and the component
means `m_j`, and the esicar precision `tau L + P0` is block diagonal on
that split with `var(m_j) = con_sd^2` exactly, independent of `n_j` and
of `tau`. So the identity Jacobian adds exactly `con_sd^2` to the
variance of every prediction and every condSD.

(sections below are filled as they are measured)

## 2. What landed

| file:line | change |
| --- | --- |
| `R/covstruct.R:1332` | `car_center_jacobian()`, the derivative of `car_center()`: the per-component projection `P`, as triplets in within-block positions. A singleton emits NO entries. |
| `R/covstruct.R:1364` | `car_center_condsd()`, `sqrt(diag(V) - con_sd^2)`, which IS `sqrt(diag(P V P'))` in closed form. |
| `R/predict.R:632` | `rr_jacobians()` roxygen: the function is no longer rr-only. |
| `R/predict.R:672` | `rr_jacobians()` esicar branch, splicing `P` through `c_idx`/`b_idx`. |
| `R/predict.R:1496` | `lp_delta_A()` roxygen. |
| `R/predict.R:1516` | `lp_delta_A()` derives the need for a Jacobian from `frame_needs_expand()`, with the caller's `has_rr` as a fast path, and pairs the Z columns through `d cvec/d b` whenever b-space is not coefficient space. |
| `R/methods-fit.R:611` | `ranef()` reads `diag.cov.random` as VARIANCES; the esicar branch projects, every other block takes the sqrt as before. |

NO call site of `lp_delta_A()` was touched: it fills in the Jacobian
when the caller passed none, so `predict()`, `frm_lp_basis()`,
`lp_basis_nl()` and the conditional-effects path all inherit the fix
without an edit in a sibling lane's region.

### The form chosen, and the one not chosen

The exact delta method, NOT pinning `kappa`. `P` annihilates the inert
coordinate whatever its scale, so the standard error is independent of
`con_sd` for ANY `kappa`, and pinning would add a second behavior
change (a silently ignored argument, and an esicar log-likelihood that
moves at 1e-12 away from the default `con_sd`) without moving a
reported number. `con_sd` therefore still sets the inert coordinate's
prior scale, which keeps `car_cov()` consistent with the density that
`draw_b()` samples, and `icar` and `bym2` are untouched by
construction.

## 3. The projection, verified three ways on three graphs

`<scratch>/cj-verify2.R`. `Jb` is `rr_jacobians()`'s
`d cvec/d b`; `P` is built from `W` outside the package.

| check | 4x4 lattice | disconnected, c = 2 | singleton, nj = (4,4,1) |
| --- | --- | --- | --- |
| `Jb` vs analytic `P`, max abs | 0 | 0 | 0 |
| `Jb` vs central difference of `expand_b()`, max abs | 1.578e-10 | 1.056e-10 | 1.172e-10 |
| `predict()` se vs `[X, Z P] V [X, Z P]'` from `frm_joint_cov()` | BITWISE | BITWISE | BITWISE |
| the identity-Jacobian rebuild, off by, in the VARIANCE | 1.000000e-06 | 1.000000e-06 | 1.000000e-06 |
| `condSD` vs `sqrt(diag(P V P'))`, max abs | 5.024e-15 | 7.258e-15 | 1.121e-14 |
| unprojected `sqrt(diag(V))`, off by | 2.687e-06 | 4.263e-06 | 1.000e-03 |

The fourth row is the whole finding in one number: the old design's
variance excess is `con_sd^2 = 1e-6` on every graph, exactly, and the
new one reproduces an independent from-scratch rebuild bit for bit.

The singleton is the cleanest signature. Its field is exactly 0, and
its condSD was exactly `con_sd` before (0.01, 0.001, 1e-4 at those
three settings) and is exactly 0 now.

### `diag.cov.random` is the MARGINAL b block, not the conditional one

Worth recording because it cost me a wrong reference. `<scratch>/cj-dbg.R`:

    max |diag.cov.random - diag(solve(Q))[b]|   1.77e-15   (marginal)
    max |diag.cov.random - diag(solve(Q[b,b]))| 3.33e-03   (conditional)

so `ranef(condVar = TRUE)` reports the joint-posterior variance of `b`.
Pre-existing and unchanged by this lane; only the projection is new.
The closed form holds against BOTH covariances, which is the point of
using it:

    max |diag(P V P') - (diag(V) - con_sd^2)|, marginal     2.08e-17
    max |diag(P V P') - (diag(V) - con_sd^2)|, conditional  2.78e-17

## 4. The invariance table

`<scratch>/cj-baseline.R` fits the four car types and an esicar
`con_sd` sweep on `car_lattice_data(42)`, before and after, and
`<scratch>/cj-compare.R` joins the two snapshots.

### esicar prediction SEs and ranef condSDs against con_sd

Relative to the `con_sd = 1e-3` fit, maximum over all 96 rows / 16
levels:

| con_sd | se rel BEFORE | se rel AFTER | condSD rel BEFORE | condSD rel AFTER |
| --- | --- | --- | --- | --- |
| 1e-1 | 1.267482e-01 | 1.798561e-14 | 1.351923e-01 | 1.809664e-14 |
| 1e-2 | 1.333574e-03 | 2.731149e-14 | 1.427998e-03 | 1.687539e-14 |
| 1e-3 | 0 | 0 | 0 | 0 |
| 1e-4 | 1.334472e-05 | 3.843592e-12 | 1.429027e-05 | 3.647971e-12 |

In the VARIANCE, which is where the leak was exactly `con_sd^2` and
where the test now asserts the equality:

| pair | max abs var diff, se | max abs var diff, condSD | the old leak |
| --- | --- | --- | --- |
| 1e-2 vs 1e-3 | 2.4217e-15 | 1.3184e-15 | 9.9000e-05 |
| 1e-4 vs 1e-3 | 3.7789e-13 | 2.8656e-13 | 9.9000e-07 |

**Over the required 1e-2 to 1e-4 range: 3.78e-13 (se) and 2.87e-13
(condSD), both under 1e-12.**

### The residual is the optimizer, not con_sd

3.8e-12 relative rather than 0 because the esicar OBJECTIVE is itself
con_sd-invariant only to about 1e-11 (the `log|K|` against `log|H|`
cancellation is exact on paper, not in floating point), so the optimum
itself moves: `max |theta(1e-4) - theta(1e-3)| = 1.28e-12`, and the
standard errors follow it with a sensitivity of about 3. Tightening
the optimizer does not help, because there is nothing wrong with the
optimizer:

| optimizer setting | theta max diff | se var diff |
| --- | --- | --- |
| default | 1.282e-12 | 3.779e-13 |
| `restarts = 3` | 1.282e-12 | 3.779e-13 |
| `rel.tol = 1e-14, x.tol = 1e-12, restarts = 3` | 1.282e-12 | 3.779e-13 |

Pinning the OUTER parameter vector with a no-op optimizer, so that
`con_sd` is the only difference between the fits
(`<scratch>/cj-fixedpar.R`):

| graph | con_sd 1e-1 | 1e-2 | 1e-4 |
| --- | --- | --- | --- |
| 4x4 lattice, se max rel | 1.916e-13 | 1.863e-13 | 4.985e-11 |
| disconnected, se max rel | 2.859e-12 | 2.890e-12 | 9.229e-11 |
| singleton, se max rel | 1.438e-12 | 1.450e-12 | 1.882e-10 |

Note the DIRECTION. The residual SHRINKS as `con_sd` grows, because it
is the conditioning of `tau L + P0` (which degrades as `con_sd^-2`)
feeding the sparse solve. The leak did the opposite: it grew a
hundredfold per decade of increasing `con_sd`. Nothing that behaves
like `con_sd^2` is left.

### escar, icar and bym2

`identical()` on the fitted `logLik`, `theta`, `beta`, `b`, the fitted
values, `predict(se.fit = TRUE)$se.fit`, `ranef()` and
`ranef(condVar = TRUE)`'s condSD, before against after:

| type | logLik | theta | beta | b | fitted | se.fit | ranef | condSD |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| escar | IDENT | IDENT | IDENT | IDENT | IDENT | IDENT | IDENT | IDENT |
| icar | IDENT | IDENT | IDENT | IDENT | IDENT | IDENT | IDENT | IDENT |
| bym2 | IDENT | IDENT | IDENT | IDENT | IDENT | IDENT | IDENT | IDENT |
| esicar | IDENT | IDENT | IDENT | IDENT | IDENT | moved | IDENT | moved |

Exactly two things moved and they are the two that had to. `icar`
keeps its documented `con_sd` sensitivity untouched.

### What moved, on the esicar fit at the default

    before se range   0.192596926 .. 0.253203174
    after  se range   0.192594330 .. 0.253201199
    max |var difference|            1.000000e-06  = con_sd^2, exactly
    max relative drop in the se     1.347970e-05
    condSD max |var difference|     1.000000e-06

1.347970e-05 is the reviewer's independently measured
1.347969552e-05.

## 5. Scope: what was left out, and why

### Pinning `kappa` was NOT taken

The reviewer's fallback (`dev/reviews/2026-09-05-esicar.md`, P2 and
R2) was to pin esicar's `kappa_j` to the default constant so that
"con_sd is ignored" became literally true. It is unnecessary once the
Jacobian is exact, and it costs more than it buys:

- `P` annihilates the inert coordinate whatever its scale, so pinning
  cannot change any reported number by more than the conditioning
  residual measured in section 4;
- it would make `con_sd` a SILENTLY ignored argument on one type,
  which is a different documentation defect from the one being fixed;
- it would move esicar's log-likelihood at the 1e-12 level away from
  the default `con_sd`, so esicar would stop being bit-identical to
  0.52.0 for no gain;
- `con_sd` still buys something real: it is the prior scale of the
  coordinate `car_cov()` hands `draw_b()`, so pinning would put the
  density and the draw covariance on different footings.

The one genuine argument FOR pinning is conditioning: a user who sets
`con_sd = 1e-8` on an esicar term gets a needlessly ill-conditioned
`tau L + P0` (section 4's residual grows as `con_sd^-2`). That is a
separate, pre-existing hazard shared with `icar`, and it is a
different change from this one.

### `R/compat.R:1231` is now FALSE and is a sibling's file

`frm_compat_rules()` still says

    con_sd changes no esicar estimate, but it does scale that type's
    prediction standard errors and ranef() conditional SDs as con_sd^2
    in the variance, so leave it at the default there.

`R/compat.R` belongs to wt-protocol this round, so this lane did not
touch it. Proposed replacement for that clause, to be applied by
whoever holds the file:

    con_sd changes nothing an esicar term reports: the constraint is
    exact and the delta method uses the centering projection, so the
    estimates, the prediction standard errors and the ranef()
    conditional SDs are all invariant to it.

The rest of that rule (per-component against brms's global sum) stays
as written.

## 6. Combinations

`<scratch>/cj-combo.R`. The fix is in one seam, so every path that
reaches `lp_delta_A()` inherits it without a call-site edit.

| check | result |
| --- | --- |
| `predict(newdata =)` se against the in-sample se, same rows | 0.000e+00 max rel |
| `frm_lp_basis()`'s `sqrt(diag(A V A'))` against `predict()` se | 2.220e-16 max rel |
| `conditional_effects()` on an esicar fit | returns its frame |
| esicar + `(1 \| g)`: the `g` block's Jacobian entry | exactly `I` |
| esicar + `(1 \| g)`: the car block's Jacobian entry | exactly `P` |
| esicar + `rr(it + 0 \| sub, d = 2)`: car block | exactly `P` |
| esicar + rr: the rr block keeps its loadings and 7 theta columns | yes |
| esicar + rr: `predict(se.fit = TRUE)` finite | yes |

The last three close the review's "not a defect, recorded" note that
`esicar` + `rr` inherited the same approximation because
`rr_jacobians()` filled identity entries for every non-rr block. It
no longer does.

## 7. Cost

`<scratch>/cj-cost.R`, lattices with 4 rows per location.

| graph | type | fit | `predict(se.fit)` | `rr_jacobians()` | `Jb` nnz |
| --- | --- | --- | --- | --- | --- |
| 6x6, n = 36 | esicar | 1.9 s | 0.20 s | 0.001 s | 1296 |
| 6x6, n = 36 | icar | 0.1 s | 0.04 s | - | - |
| 10x10, n = 100 | esicar | 0.4 s | 0.10 s | 0.002 s | 10000 |
| 10x10, n = 100 | icar | 0.2 s | 0.11 s | - | - |

Building the projection costs 1 to 2 ms and `predict(se.fit = TRUE)`
is not measurably slower than `icar`'s. `P` is dense within a
component, so `Jb` carries `sum_j n_j^2` entries; at 5000 locations in
one component that would be 2.5e7, but `get_joint_cov()` already
materializes a dense `p x p` covariance at that size, so the Jacobian
is not what binds first.

## 8. Verification

One `test_file()` per process, `NOT_CRAN=true`,
`test_file(package = "frmtmb")`, private library `<scratch>/cj-lib`
with the worktree core installed into it. Baselines measured, not
recalled: `git archive HEAD` (68d6782) into `<scratch>/cj-base-src`,
installed into its own `<scratch>/cj-libbase`, and the files run
against it.

### The named files

| file | at 68d6782 | this lane |
| --- | --- | --- |
| test-car-spde.R | 22 tests, 124 assertions | 23 tests, 135 assertions, 0 fail |
| test-predict-lp-basis.R | 5 tests, 30 assertions | 5 tests, 30 assertions, 0 fail |
| test-predict-newdata.R | 6 tests, 12 assertions | 6 tests, 12 assertions, 0 fail |
| test-importance.R | - | 31 tests, 185 assertions, 0 fail |
| test-message-uniqueness.R | - | 1 test, 6 assertions, 0 fail |
| test-bracket-access.R | - | 3 tests, 8 assertions, 0 fail |

`test-car-spde.R` gains 1 test and 11 assertions: the new
"the esicar delta method is Z P V P' Z', not Z V Z'" block, plus the
singleton's condSD, the disconnected graph's per-component projection
rebuild, and the third `con_sd` fit in the invariance block. Nothing
was removed; the `con_sd` block was rewritten in place, keeping its
name's first clause and flipping the second.

The two `test-predict*` files are unmoved at their baselines, which is
the check that matters for `lp_delta_A()`: every other block type
still goes through the identity.

### The gated brms log-density tier

`FRMTMB_BRMS_FIT_TESTS=true`, `NOT_CRAN=true`, `FRMTMB_STAN_CACHE` and
`R_MAKEVARS_USER` pointing at my own copy of
`C:/Users/adf44/source/r/frmtmb/dev/stan-cache` at
`<scratch>/cj-stan-cache`.

    test-brms-likelihood.R   33 tests, 373 assertions,
                             0 fail, 0 error, 0 skip

373, unmoved, as the brief says it should be.

One caveat, with its control. Main's cache holds 53 `.rds` programs,
and the first run compiled ONE
(`c17f0e928eb7ed5b9fcd26784c8221d1.rds`), leaving 54. That is main's
cache being one program short of the tier, not this lane: the second,
warm run went 54 to 54 with nothing recompiled and the same 33/373.
Nothing in this lane touches a density or a generated program.

### Full core suite, one file per process

    117 files, 1161 tests, 6647 assertions,
    0 fail, 0 error, 92 skip

By-name audit against
`ls tests/testthat/test-*.R`: 117 on disk, 117 in the log, nothing on
disk missing from the log, nothing in the log that is not on disk, no
file logged twice. Every file reported `fail = 0 error = 0`.

Skips: test-brms-methods 45, test-brms-likelihood 30, test-brms-priors
11, test-brms-agreement 2, test-rl-example 2, test-fuzz 1, test-ps 1.

NOTE on the file count. The brief says 116 files at 68d6782;
`git ls-tree --name-only HEAD tests/testthat/` counts 117 tracked
`test-*.R` files there, and the working tree has the same 117 with one
of them modified by this lane. The audit is by name against what is on
disk, so the count is stated rather than assumed.

### Against the base's own full suite, file by file

`git archive HEAD` into `<scratch>/cj-base-src`, installed into
`<scratch>/cj-libbase`, and the whole suite run against it the same
way, one file per process.

    BASE  117 files, 1160 tests, 6637 assertions, 91 skip
    LANE  117 files, 1161 tests, 6648 assertions, 91 skip

A by-name join reports three files differing. Two of them are the
shared machine, not this lane, and both were re-measured until they
agreed:

| file | base | lane | verdict |
| --- | --- | --- | --- |
| test-car-spde.R | 22 / 124 | 23 / 135 | THIS LANE, intended |
| test-perf.R | 2 / 3, **1 fail** | 2 / 3, 0 fail | a timing flake in the BASE run: it ran while `R CMD check` had the machine, and passes twice on the base and once on the lane when re-run alone |
| test-ps.R | 9 / 51, 0 skip | 9 / 50, 1 skip | `skip_if_not_installed("brokenstick")`. brokenstick was absent from the SHARED user library when the lane suite ran and present when the base suite ran, because a sibling lane installed it in between. Re-run now, both give 9 / 51 / 0 skip |

So exactly ONE file moved, and it is the one this lane edits. The
+1 test and +11 assertions are its whole delta. The raw lane figure of
6647 assertions and 92 skips above is the same run with test-ps.R
caught on the wrong side of that library race; 6648 / 91 is what it
reads once brokenstick is there for both.

### R CMD check --as-cran

`--as-cran --no-manual`, `_R_CHECK_CRAN_INCOMING_=false`, pandoc 3.8.3
on PATH from `RSTUDIO_PANDOC`, checked against `<scratch>/cj-lib`,
tarball built WITH vignettes.

    Status: OK          no WARNINGs, no NOTEs
    examples ............................... [64s] OK
    examples with --run-donttest ........... [66s] OK
    tests .................................. [30m] OK
      [ FAIL 0 | WARN 1 | SKIP 97 | PASS 6538 ]
    re-building of vignette outputs ........ [391s] OK

The vignette rebuild is the step that knits the edited
`vignettes/frmtmb.Rmd`, and it passed.

`_R_CHECK_FORCE_SUGGESTS_=false` was needed and is recorded rather than
hidden. A first run without it stopped at
"Packages suggested but not available: 'brokenstick', 'frmtmb.spline'",
which is this machine's library and not the package: `frmtmb.spline`
lives in `extensions/` and was installed into `<scratch>/cj-lib` from
the worktree (read only; that extension belongs to wt-spline-span), and
`brokenstick` is a CRAN Suggests that was absent at that moment.

### roxygen

`roxygen2::roxygenise()` twice on the worktree leaves `man/` and
`NAMESPACE` untouched: `git status` after both passes lists only the
six files this lane edited. Every roxygen block added here is `@noRd`;
nothing exported, no `@param` or `@rdname` touched.

## 9. Files touched

`R/covstruct.R`, `R/predict.R`, `R/methods-fit.R`,
`tests/testthat/test-car-spde.R`, `vignettes/frmtmb.Rmd`, `NEWS.md`,
plus this file. 232 insertions, 43 deletions. `man/` and `NAMESPACE`
unchanged. Nothing committed; the worktree HEAD is still 68d6782.

No sibling lane's file was edited. In `R/predict.R` the three hunks sit
at 629, 669 and 1494, well clear of wt-spline-span's `frm_lp_basis()`
region at about 3068, and no call site of `lp_delta_A()` was touched -
which is why `lp_basis_nl()` at 3229, inside that region, picks the fix
up without an edit. `R/compat.R:1231` needs the correction quoted in
section 5 and belongs to wt-protocol.

## 10. Summary

The identity Jacobian is gone from every reported quantity on an
esicar term.

- `predict(se.fit = TRUE)` is bitwise equal to an independent
  `[X, Z P] V [X, Z P]'` rebuilt outside the package on three graphs;
- `ranef(condVar = TRUE)` matches `sqrt(diag(P V P'))` to 1e-14, and a
  singleton reports exactly 0 where it used to report exactly
  `con_sd`;
- both are invariant to `con_sd` over 1e-2 to 1e-4 to 3.8e-13 in the
  variance, against a leak that was exactly `con_sd^2`;
- the residual shrinks as `con_sd` grows, where the leak grew; it is
  the esicar objective's own con_sd-invariance limit, which moves
  `theta` by 1.3e-12 between those fits;
- `escar`, `icar` and `bym2` are bit-identical on every reported
  quantity, and so is esicar's `logLik`, `theta`, `beta`, `b`,
  `ranef()` and `VarCorr()`.

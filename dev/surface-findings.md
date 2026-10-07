# Lane surface: brms post-processing calls on an ML fit

## Punch round 1 (2026-10-07)

The review is `dev/reviews/2026-10-07-surface.md` (NOT MERGEABLE on
B1 and B2). What changed:

- **B1, `add_criterion(overwrite = TRUE)` on draws did not recompute.**
  `add_criterion.frmtmb_draws()` now runs `x$criteria[new] <- NULL`
  before computing, as brms's does, so `loo()`'s `use_stored` path
  cannot hand back the object being replaced. `test-add-criterion.R`
  now stores a loo from `ndraws = 100`, overwrites it, and asserts the
  new one's `dims` (all draws) and its estimates against `loo(ds)`.
  Seen to fail on the build before the fix
  (`dev/surface-p1-tests/addcrit-before-fix.log`):

  ```text
  FAIL [ add_criterion() stores loo and waic, and loo() reads them ]: Expected `attr(dsb$criteria$loo, "dims")[1L]` to equal `ndraws(ds)`.
    `actual`: 100
  `expected`: 200
  RESULT frmtmb.sample test-add-criterion.R pass=39 fail=2 error=0 skip=0 warn=0
  ```

- **B2, `pp_mixture()`'s `Est.Error`.** It is now the standard
  deviation of the law the logit Wald interval is built from,
  `plogis(Z)` with `Z ~ N(logit p, lse^2)`, integrated numerically
  (`logitnormal_sd()`): 100-point Gauss-Hermite where the logit SD is
  at most 4 and `stats::integrate()` row by row where it is wider,
  because Gauss-Hermite converges slowly there (`dev/surface-p1-gh.R`,
  `dev/surface-out/p1-gh.txt`: against `integrate()` over a grid of
  logit means -20 to 0 and SDs 0.1 to 10, 40 nodes are off by up to
  4.2%, 100 by 0.65%, 200 by 0.058%; at SD at most 4, 100 nodes by
  4.1e-6). The logit and its standard error are read off the smaller
  of `p` and its complement (the sum of the other components), so a
  probability near 1 does not lose the digits of `1 - p`: before that,
  the two components of a two-component fit disagreed by 0.1% on a row
  near the edge. Unfiltered, all 300 rows, against brms 2.23.0's
  `pp_mixture()` on the same data (the reviewer's brms fit,
  `dev/surface-rev-out/brms-ppmix.rds`; `dev/surface-p1-ppbrms.R`,
  `dev/surface-out/p1-ppbrms-lane.txt`):

  ```text
                  n to_brms  to_mc brms_sd
  [-Inf,0.0001)  26   4.750 18.200 0.00715
  [0.0001,0.001) 35   2.370  5.960 0.01470
  [0.001,0.01)   79   1.270  2.590 0.04220
  [0.01,0.05)    60   1.030  1.530 0.07780
  [0.05,0.2)     57   0.843  1.090 0.14300
  [0.2,0.5)      43   0.914  0.953 0.21400
  all 300 rows: median Est.Error / brms SD 1.11 | / MC SD 1.81
  rows with p in [0.001, 0.999]: 239 | Est.Error / brms SD range of the bin medians 0.843 to 1.27
  interval width / brms's: median 1.04
  estimates inside brms's 95% interval: 300 of 300
  p 0.503333 lse 0.987: Gauss-Hermite 0.206309, integrate() 0.206309, ratio 1.00000000
  p 6.53861e-08 lse 5.936: Gauss-Hermite 0.0450179, integrate() 0.0450179, ratio 1.00000000
  p 0.999683 lse 3.63: Gauss-Hermite 0.0979515, integrate() 0.0979515, ratio 1.00000006
  p 0.00953235 lse 1.679: Gauss-Hermite 0.0625583, integrate() 0.0625583, ratio 1.00000000
  p 0.999006 lse 2.691: Gauss-Hermite 0.0656417, integrate() 0.0656417, ratio 1.00000000
  ```

  These reproduce the reviewer's figures. Within 0.001 of an edge the
  law overstates brms's SD (2.4 and 4.8, bin medians); the help page
  and NEWS say so, with the unfiltered numbers and the brms comparison,
  and the help page no longer cites a dev script (m8). The delta
  method's figure (0.36 of brms's) is given for contrast. A new test
  checks `Est.Error` against `integrate()` of the same law at four rows
  and in the log case; seen to fail on the build before the change
  (`dev/surface-p1-tests/postfit-before-fix.log`: `pass=60 fail=6`).

- **m1.** The findings' reason for not merging `update(prior =)` is
  withdrawn (item 3 below); NEWS says "replaces, where brms merges";
  filed in `dev/test-backlog.md`.
- **m2.** `dev/lane-rules.md` now says 57 active bindings (33 + 24);
  `R/generic-owners.R`'s header counts the owner table (33) and the
  four static imports.
- **m3.** `plot_fit_brms_args`: `plot` is "whether the pages are drawn
  or only returned", `newpage` "whether the first page starts a new
  graphics page", as `plot.brmsfit()` uses them.
- **m4.** The `+Inf`/`NaN` scope gap (`us` with two or more
  coefficients, `cs`, `ar1`) is filed in `dev/test-backlog.md`, and the
  "sampled non-centered" claim is corrected below.
- **m5.** NEWS now says the mice agreement is the estimates (4.1e-6),
  that the SEs are 0.93 of mice's (ML within-imputation variances), and
  that the complete-data df is `df.residual()` of the fit, which counts
  `sigma` (20 against lm's 21), as `anova()` of a `frm_multiple()`
  already documents. `dfcom` was left as it is: one convention for
  `fixef()`, `summary()` and `anova()` of a pooled fit.
- **m6.** `pp_mixture()` of draws of a non-mixture model stops in
  brms's words, as the fit method does. This makes the ported-suite row
  `brmsfit-methods:719` (recorded "pass in frmtmb's own words") stale:
  frmtmb.sample's `test-brms-suite-methods.R` reports "STALE
  OWN-WORDS: brms's assertion holds as written", the expected flip for
  the ledger regeneration (verdicts not edited).
- **m7.** `test-covstruct-floor.R` builds each of the eight block types
  the floor touches (`gp`, `ou`, `homcs`, `homtoep`, `exp`, `gau`,
  `mat`, one-coefficient `gr(cov =)`) from a `dry_run = "frame"` and
  asserts no `+Inf` at log sd -1137.64; on rellib-r6 all eight fail
  that assertion (`dev/surface-p1-tests/floor-base.log`, `pass=13
  fail=23 error=1`).
- **Versions.** Core's DESCRIPTION is `0.68.1.9000` (the development
  marker) and frmtmb.sample's floor is `frmtmb (>= 0.68.1.9000)`, so
  the package cannot be installed against 0.68.1, which lacks the
  generics it imports. The release numbers are left to consolidation.

Reruns on the final build (one file per process, every gate on):

```text
core:  test-brms-postfit 66 pass | test-covstruct-floor 39 | test-multiple-methods 39 |
       test-message-uniqueness 6 | test-arg-refusal 114 | test-generic-collision 61 |
       test-update-prior-drop 9 | test-portability 96; 0 fail, 0 error, 0 warn each
frmtmb.sample full suite (dev/surface-suite-sum.R p1sample):
       49 files, 2655 pass, 1 fail, 0 error, 1 skip, 0 warn
       the fail: test-brms-suite-methods.R, brmsfit-methods:719, STALE OWN-WORDS (m6)
       the skip: test-scale.R (the scale tier's gate)
```

R CMD check `--as-cran` of frmtmb.sample (`dev/surface-check/
frmtmb.sample.Rcheck/00check.log`): `Status: OK`.

### Re-check notes (MERGEABLE; n1 fixed, n2 and n3 noted)

- **n1, fixed.** `logitnormal_sd()`'s adaptive path stopped with "the
  integral is probably divergent" at 3 of the reviewer's 168 grid
  points (p 1e-3 and 1e-5 at logit SD 3000, p 1e-30 at SD 1e4), which
  would abort `pp_mixture()` for the whole fit. The error is now caught
  and the integral redone in pieces (`logitnormal_pieces()`: breaks at
  the normal's centre and 1, 3, 6 and 14 SDs either side, and at the
  logistic's transition `z = 0` and 3, 10 and 40 either side), with
  `stop.on.error = FALSE` per piece so no input stops the call. A test
  at the three points (finite, at most 0.5, above 0.95 of it) failed
  on the build before (`dev/surface-p1-tests/postfit-n1-before.log`:
  "the integral is probably divergent") and passes after:
  `test-brms-postfit.R` pass=76 fail=0 error=0 warn=0. The reviewer's
  scripts rerun on the lane library (`dev/surface-n1-lnsd*.R`, copies
  with this lane's library paths; `dev/surface-out/n1-lnsd*.txt`):
  "errors: 0 of 168"; the NaN/Inf/0 inputs, the log case, the
  three-component mixture and the large-SD limits read as in the
  reviewer's logs.
- **n2, noted, not fixed.** At p = 1e-30 with logit SD 10 the single
  `integrate()` over `mu +- 14 s` (which does not fail, so the pieced
  fallback is not reached) misses the logistic's transition 7 SDs out:
  2.21e-6 against a careful 1.41e-6, relative 0.568
  (`n1-lnsd2.txt`). The absolute size is negligible.
- **n3, noted.** At logit SD 0 the result is rounding noise, 3.05e-16
  at p 0.2, not 0.

Files touched in this round: `R/brms-postfit.R` (`pp_mixture.frmtmb_fit`,
`mixture_prob_band`, new `gh_std_normal`, `logitnormal_sd`, roxygen),
`R/conditional-effects.R` (`plot_fit_brms_args`),
`R/generic-owners.R` (header), `DESCRIPTION` (Version), `NEWS.md`,
`man/pp_mixture.Rd`, `tests/testthat/test-brms-postfit.R`,
`tests/testthat/test-covstruct-floor.R`; frmtmb.sample `R/loo.R`
(`add_criterion.frmtmb_draws`), `R/methods-draws.R`
(`pp_mixture.frmtmb_draws`), `DESCRIPTION` (floor), `NEWS.md`,
`tests/testthat/test-add-criterion.R`; `dev/lane-rules.md`,
`dev/test-backlog.md`, this file, `dev/surface-p1-ppbrms.R`,
`dev/surface-p1-gh.R`; `dev/surface-check.sh` (`CHECK_PKGS`).

---

Date 2026-10-06. Worktree `frmtmb-wt-surface`, branch `wt-surface`,
base 4f5ea39f (frmtmb 0.68.1). Base build `rellib-r6` (frmtmb 0.68.1,
frmtmb.sample 0.16.0); lane build `wt-surface-lib` (core and
frmtmb.sample from this worktree; the other six extensions are the
rellib-r6 builds behind it, loading the lane's core). brms 2.23.0,
R 4.6.1, Windows.

Every number below is pasted from a log under `dev/surface-out/`,
`dev/surface-tests/`, `dev/surface-suite-*/`, `dev/surface-port-out/`
or `dev/surface-bv-out/` (all gitignored), with the script that wrote
it. "base" is rellib-r6, "lane" is the lane library first.

## Summary

| item | outcome |
|---|---|
| 1 `stancode()`, `standata()`, `pp_mixture()` on a fit | FIXED: refusals by name for the first two; `pp_mixture()` computed at the estimates in brms's array, ordinal mixtures included. Generics moved to core, as `parnames` did |
| 2 `plot()` of a fit and `N`, `variable`, `regex` | FIXED: each of brms's 11 `plot.brmsfit()` arguments refused by name with the reason; and `plot()` of DRAWS now draws brms's display, so the route the refusal names works |
| 3 `update()` and a stale `lkj` prior | FIXED: dropped as brms drops it, with a message where brms is silent |
| 4 `frm_multiple()` pooled post-processing | FIXED: `fixef()`, `summary()`, `conditional_effects()` pool by Rubin's rules; pooled table uses brms's names. `plot()`, `as_draws_array()`, `nchains()` stay refused (no draws) |
| 5 `frm_sample()` repeats a frame-build warning | NOT A DEFECT: the record's "6 times" was all six warnings of the call; the frame-build warning appears once. Guard test added |
| 6 `(cs(1) | g)` with brms attached; categorical message | FIXED, both |
| 7 sum-to-zero component lists no Intercept rows | DOCUMENTED DIVERGENCE, with brms's Stan code as the reason; no change |
| 8 `fitted()` with no fixed `disc` column | FIXED |
| 9 `frm_sample(fit)` on an exact `gp()` does not move | FIXED at the cause: the field's log density was `+Inf` at an underflowed sd |
| 10 `summary(waic =)`, `add_criterion()` | REFUSED BY NAME on a fit (no honest ML meaning); `add_criterion()` on draws implemented with brms's semantics |

## Item 1: stancode(), standata(), pp_mixture()

Repro: `dev/surface-repros.R` (one line per defect), logs
`dev/surface-out/repros-base.txt` and `repros-lane.txt`.

Base (rellib-r6):

```text
1 stancode(fit), brms not loaded                   ERROR: no applicable method for 'stancode' applied to an object of class "frmtmb_fit"
1 pp_mixture(gaussian mixture fit)                 ERROR: no applicable method for 'pp_mixture' applied to an object of class "frmtmb_fit"
1 stancode(fit), brms loaded                       ERROR: Data must be specified using the 'data' argument.
1 standata(fit), brms loaded                       ERROR: Data must be specified using the 'data' argument.
```

Lane:

```text
1 stancode(fit), brms not loaded                   ERROR: stancode() has no meaning for a frmtmb fit: there is no Stan program. ...
1 pp_mixture(non-mixture fit)                      ERROR: Method 'pp_mixture' can only be applied to mixture models.
1 pp_mixture(gaussian mixture fit)                 ok
   dim 120 4 2 | dimnames[[2]] Estimate Est.Error Q2.5 Q97.5 | [[3]] P(K = 1 | Y) P(K = 2 | Y)
   max |Estimate - mixture_probs()| 0
1 pp_mixture(ordinal mixture fit)                  ok
   dim 300 4 2
1 stancode(fit), brms loaded                       ERROR: stancode() has no meaning for a frmtmb fit: ...
1 brms::stancode(fit)                              ERROR: stancode() has no meaning for a frmtmb fit: ...
1 brms::pp_mixture(fit)                            ok
```

Decisions.

- **Where the generics live.** frmtmb.sample defined `stancode`,
  `standata` and `pp_mixture` (owner brms). A fit method there would
  leave a core-only user with brms loaded at brms's default. They moved
  to core's `frm_generic_owners` with owner brms, and frmtmb.sample
  re-exports them, exactly as `parnames` moved. `add_criterion` was
  added to the same table. frmtmb.sample keeps its draws methods,
  registered in both tables.
- **`standata(fit)` refuses rather than returning the frame.** The
  frame is not brms's data list (no `Z_1_1`, `J_1`, ...).
  `brms_standata_view()` in the ported-suite helper maps part of it and
  says it leaves the rest out "rather than numbered by a guess"; a
  partial list under brms's names would be read as brms's.
- **`pp_mixture()` on a fit.** `Estimate` is `mixture_probs()`;
  `Est.Error` the delta method (`fit_fd_se()`, group effects
  included); quantile columns a Wald interval on the logit. brms's
  names: `P(K = k | Y)`, rows numbered (brms's `pp_mixture.brmsfit()`
  source, `dev/surface-out/brms-src.txt`). `summary = FALSE`,
  `robust`, `ndraws`, `draw_ids`, `newdata` and a non-NULL
  `re_formula` refused by name. frmtmb.sample's draws method now uses
  the same component names (it used `class1`, `class2`).
- **How good is `Est.Error`.** `dev/surface-ppmix-mc.R`
  (`dev/surface-out/ppmix-mc.txt`): a two-gaussian mixture, 300 rows
  (data seed 4), against the SD of the probability over 4000 draws
  (seed 1) from N(estimate, `vcov(full = TRUE)`):

```text
delta Est.Error / Monte Carlo SD (R = 4000, seed 1), rows with MC SD > 0.01:
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.03177 0.47225 0.69395 0.68823 0.91083 1.12780
the same, rows with MC SD > 0.05:
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 0.5963  0.7665  0.8967  0.8937  1.0439  1.1278
rows whose logit-Wald interval overlaps the MC 95% range: 1
interval width, Wald / MC, rows with MC width > 0.02:
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 0.9186  1.0456  1.4812  1.6555  2.1249  3.2033
```

  WITHDRAWN at punch round 1: these figures filter out the 56 rows
  where the ratio is worst, and the delta SE they describe is no longer
  `Est.Error`. See "Punch round 1", B2, for the unfiltered numbers and
  the replacement.

  The first-order SE understates the spread near 0 and 1, where the
  probability is far from linear; the logit interval is at least 0.92
  of the Monte Carlo range everywhere and usually wider. The help page
  says to read the interval near the edges. Against the posterior SD of
  `frm_sample()` draws (2 chains, 1500 iterations, seed 11, default
  priors) the ratio has median 0.53 (`dev/surface-check.R`,
  `dev/surface-out/check.txt`), which adds prior and posterior
  nonlinearity to the same effect.

## Item 2: plot() of a fit

brms 2.23.0's `plot.brmsfit(x, pars, combo, nvariables, N, variable,
regex, fixed, bins, theme, plot, ask, newpage, ...)`
(`dev/surface-out/brms-src.txt`): `N` is an alias of `nvariables`, and
all of them choose and lay out per-parameter histogram and trace panels
of draws. A fit has none, so there is no honest ML meaning short of a
profile per parameter, which `plot(profile(fit, parm))` already is.
Each of the 11 is refused by name (`plot_fit_brms_args`, passed to
`frm_check_dots(.unsupported =)`), and the message names
`plot(frmtmb.sample::frm_sample(fit), variable = , regex = , N = )`.

That route was itself a refusal at 0.16.0 (`plot.frmtmb_draws()`
named `mcmc_plot()`). It now is brms's display: `bayesplot::mcmc_combo()`
pages with brms's arguments and defaults, `N` honored with brms's
deprecation warning. Two hand-translation rows ("plot() of draws") were
this refusal.

Base and lane, `dev/surface-repros.R`:

```text
base: 2 plot(fit, N = 2, ask = FALSE)  ERROR: plot() has no argument `N`. Did you mean `x`?. It takes: x, which, ask
lane: 2 plot(fit, N = 2, ask = FALSE)  ERROR: plot() cannot honor `N`: brms sets how many parameters' posterior panels go on one page with it. ...
lane: 2 plot(fit) control              ok
```

## Item 3: update() and a stale prior

brms 2.23.0's `update.brmsfit()` (`dev/surface-brms-src2.R`,
`dev/surface-out/brms-src2.txt`) sets
`attr(dots$prior, "allow_invalid_prior") <- TRUE` on the stored prior
(and on a new one merged with the old user priors), so every
specification that matches no parameter of the new model is dropped,
silently. The vigport reviewer's run confirmed it samples.

frmtmb now does the same for the stored prior: `update()` passes
`object$prior` by value marked `allow_invalid_prior`, and `frm()` drops
(`prior_drop_unmatched()`) each specification that does not resolve
alone on the new model, the test `resp_missing_refusal()` uses.

Message or silence: the tiebreaker is brms unless brms is clearly
wrong. The drop is brms's and is kept. A message is added: the
behavior it replaces was a refusal, and the standing rule is that a
default must never be less diagnostic than what it replaces; a silent
drop of a user's prior would be. The message is suppressible and
changes no result.

`dev/surface-check.R` (`check.txt`):

```text
update() message: update() drops the prior specification of the original fit that matches no parameter of the updated model, as brms drops it: lkj(2) on cor
update() vs direct fit with the two kept priors: logLik diff 0
update(prior = <with cor>): No random-effect correlations match class=cor. ...
```

Not done: brms MERGES a prior passed to `update()` with the stored
user priors (new rows win per slot) and lets both be invalid; frmtmb's
`update(prior =)` replaces the stored prior and checks it as `frm()`
does. Left as it is, as it was on 0.68.1; filed in
`dev/test-backlog.md` ("Filed by lane surface"). (The reason given
here first, that merging would mean dropping a specification the user
types, was wrong, as the review says: frmtmb could merge the stored
rows and still check the new rows strictly.)

## Item 4: frm_multiple()

The brms_missings vignette's calls on a `brm_multiple()` fit
(`dev/surface-port-out/r6/results-spell/brms_missings.log`):
`summary()`, `plot(variable = "^b", regex = TRUE)`,
`as_draws_array()`, `nchains()`, two cascades on `subset_draws()`, and
`conditional_effects(fit, "age:chl")`. On 0.68.1 `summary()` "ran" as
`summary.default`, printing the list layout (the hand translation
labels it BEHAVIOR), and `fixef()` had no method.

Added: `fixef.frmtmb_multiple()` (Rubin over `fixef()` of each fit:
brms's rows, names and order, ordinal thresholds included, t interval
on the Barnard-Rubin df); `summary.frmtmb_multiple()` (brms's blocks
with `df` and `fmi` where brms prints R-hat and ESS; distributional
parameters without a formula pooled on their link and transformed
back; the pooled variance components); `conditional_effects()` pooled
point by point on the first imputation's grid (brms's
`brm_multiple()` keeps the first data set too), on the scale its Wald
band is symmetric on (`ce_pool_scale()`, an attribute the fit method
now sets). `$pooled` rows carry brms's names (`Intercept`,
`sigma_Intercept`).

Kept refused: `plot()`, `as_draws_array()`, `nchains()` read draws or
chains; the per-imputation convergence check they serve has no ML
counterpart beyond the m stored fits. `plot()`'s message now names
`summary()`, `fixef()` and `conditional_effects()`.

Checks (`dev/surface-check.R`, `check.txt`; nhanes, m = 5, mice seed 1):

```text
pooled CE vs hand pooling: max |estimate diff| 0  max |se diff| 0  rows 300
pooled CE se / first imputation's se: median 1.139016
fixef(fm) Estimate vs mice::pool(lm): max |diff| 4.114721e-06
```

The hand pooling is `frm_linpred()` of each fit on the grid and
`rubin_pool()`, so the zero is an identity on the identity link (same
numbers through two routes), not an independent measurement; the mice
comparison is independent.

## Item 5: frm_sample() and a frame-build warning

`dev/test-backlog.md` filed "6 times" from `dev/rel068-cs-probe.R`.
Rerun on rellib-r6 (`dev/surface-cs-probe-copy.R`, the probe with this
worktree's Stan cache, `dev/surface-out/cs-probe-r6.txt`):

```text
frm_sample                         ok [6 warning(s): Category specific effects for this family should be considered experim]
```

The probe's `step()` prints the COUNT of all warnings and the FIRST
one's text. Counting by message (`dev/surface-sample-repros.R`,
`dev/surface-out/sample-base.txt`), same data and seed:

```text
5 frm_sample(formula route)              ok; 6 warning(s)
     1 x Bulk Effective Samples Size (ESS) is too low, indicating pos
     1 x Category specific effects for this family should be consider
     1 x Examine the pairs() plot to diagnose sampling problems
     1 x Tail Effective Samples Size (ESS) is too low, indicating pos
     1 x The largest R-hat is 1.08, indicating chains have not mixed.
     1 x There were 1 divergent transitions after warmup. See
5 frm_sample(fit route)                  ok; 5 warning(s)
5 brms::make_standata (brms's own count) ok; 1 warning(s)
```

So the frame-build warning appears once, as in `frm()` and brms. No
code change; a guard test in frmtmb.sample's `test-add-criterion.R`
asserts the count of 1.

## Item 6: (cs(1) | g) with brms attached; the categorical message

Cause: `parse_linpred()`'s `sub_specials()` rewrites a call to a
special name as `.frm_cs()` when a FUNCTION of that name is visible.
Without brms there is no `cs` function, the call stays `cs(1)`, and
`check_cs_in_bar()` refuses it. With brms attached, `brms::cs` is
found, the call becomes `.frm_cs(1)`, the check missed it, and
`model.frame()` died. The check now recognizes both spellings and
reports the term as written. brms 2.23.0 fits the model (stancode ok);
frmtmb has no group-level category-specific term and refuses it, as
before.

```text
base: 6 (cs(1) | g), brms attached after frmtmb   ERROR: variable lengths differ (found for '.frm_cs(1)')
lane: 6 (cs(1) | g), brms attached after frmtmb   ERROR: A category-specific effect inside a group-level term, (cs(1) | g), is not supported: ...
base: 6 categorical CE predict  ERROR: method = "predict" has no meaning on an ordinal family: ...
lane: 6 categorical CE predict  ERROR: method = "predict" has no meaning on family 'categorical', whose response is a category: ...
```

## Item 7: sum-to-zero component, class Intercept rows

Not changed; the reason, run: `dev/surface-brms-stz.R`
(`dev/surface-out/brms-stz.txt`) shows brms declaring
`ordered[nthres] Intercept_mu1`, computing
`Intercept_mu1_stz = Intercept_mu1 - mean(Intercept_mu1)`, putting the
prior on the raw `Intercept_mu1`, and passing only the centered vector
to the likelihood. The prior is on a parameter whose common location
the likelihood never sees; frmtmb estimates the `nthres - 1` free
directions, so there is no parameter for those rows to address. This
is the decision 0.68.0 made for one sum-to-zero family, and the
mixture component follows it. `dev/surface-stz.R`
(`dev/surface-out/stz-lane.txt`):

```text
one family, set_prior(class = 'Intercept')         ERROR: Prior target not found (class=Intercept): with threshold = 'sum_to_zero' there is no threshold vector for class "Intercept" to address. brms puts that prior on ...
mixture, set_prior(class = 'Intercept', dpar = 'mu1') ERROR: Prior target not found (class=Intercept, dpar=mu1): with threshold = 'sum_to_zero' there is no threshold vector ...
mixture, set_prior(class = 'Intercept', dpar = 'mu2') (control) ok
```

brms lists 5 Intercept rows for the stz component (one class row and
one per threshold; the brief said 4), frmtmb 0, for one family and as
a component alike (`repros-lane.txt`, item 7).

## Item 8: fitted() with no fixed disc column

`lp_eta_design()` formed `X %*% est$betad[idx]` with `betad` absent;
`lp_fixed_part()` returns the zero vector of the design's rows for a
zero-column design, as `build_objective()` does on the tape, and is
used at the three call sites that formed it (two in
`lp_eta_design()`, one in `simulate-newdata.R`).

## Item 9: frm_sample(fit) on an exact gp() fit

Diagnosis, all on rellib-r6, the gpby review's construction (60
points, data seed 5):

1. Not the tape and not the prior. `dev/surface-gp-diag.R`: the joint
   density and all 42 gradient entries are finite at the mode init.
2. The mode init is the trigger. `dev/surface-gp-diag2.R`:

```text
mode init, prior obj               first stepsizes 2.44e-04 1.16e+01 8.84e+00 1.03e+00 | post-warmup accept 0.000
random init, prior obj             first stepsizes 1.22e-04 1.44e+01 2.43e+00 2.40e-01 | post-warmup accept 0.895
mode + N(0, 1e-3) init             first stepsizes 0.000488 7.814780 1.102219 0.098049 | post-warmup accept 0.886
```

   `dev/surface-gp-diag3.R`: retaping elsewhere does not help; a 1e-12
   offset does not either (accept 0.000, 150 of 150 divergent each).
3. The cause. `dev/surface-gp-diag4.R`: the second warmup iteration
   (step size 11.6, adapted up from 2.4e-4) moved the chain to a point
   with `energy__ = -Inf`, a log density of `+Inf`, and it never left:

```text
      accept_stat__ stepsize__ treedepth__ n_leapfrog__ divergent__  energy__
 [1,]        0.8813     0.0002          10         1023           0 -141.7842
 [2,]        0.6667    11.5922           1            3           1      -Inf
 [3,]        0.0000     8.8450           0            1           1      -Inf
 [4,]           NaN     1.0326           2            3           0      -Inf
 [5,]        0.0000        NaN           0            1           1      -Inf
distinct draws over all 300 iterations: 2
the point the chain moved to:
      beta      betad      theta      theta
  110.6830   350.6820 -1137.6400    96.6115
gp field log density there: Inf
```

   At log sd -1137.64, `exp(2 * log sd)` is 0, the field's covariance
   is the zero matrix, and RTMB's `dmvnorm()` of a nonzero field
   against it returns `+Inf`.

Fix: `sd2_floored(log_sd) = exp(2 * log_sd) + 1e-300` for the variance
of every one-scale dense block (`gp`, `ou`, `homcs`, `homtoep`, the
spatial entries, `gr(cov =)` with one coefficient; the same expression
in six places). The floor is below half an ulp of any variance above
9e-285 (log sd above about -326), so it is the variance bit for bit
there: `test-covstruct-floor.R` asserts `identical()` over 200001
points of [-326, 50], and the formula route's chain on this data is
the same on both builds (accept 0.966, step size 0.00166, both logs).
The ou block had the same `+Inf` on base (its test fails there).

`dev/surface-sample-repros.R`, fit route, sampler seed 4:

```text
base: 43 columns, NA with sd 0; accept 0.000, stepsize NaN, divergent 300
lane: 43 columns, 0 with sd 0; accept 0.886, stepsize 0.00363, divergent 0
```

The chain moves but mixes poorly (R-hat 1.89 at 600 iterations; the
formula route's is 1.56 on both builds): the centered field is a
funnel, and frmtmb.sample has no non-centered form for `gp()`
(`ncp_reasons`). That is a sampler property, not this defect.

Sampler seeds 1 to 8 on both builds (`dev/surface-gp-sweep.R`,
`dev/surface-out/gp-sweep-base.txt`, `gp-sweep-lane.txt`):

```text
arm base
seed 1: accept 0.770, stepsize 0.00349, divergent 0, distinct draws of b_Intercept 300
seed 2: accept 0.925, stepsize 0.00208, divergent 0, distinct draws of b_Intercept 300
seed 3: accept 0.908, stepsize 0.00144, divergent 0, distinct draws of b_Intercept 300
seed 4: accept 0.000, stepsize NaN, divergent 300, distinct draws of b_Intercept 1
seed 5: accept 0.856, stepsize 0.0042, divergent 0, distinct draws of b_Intercept 300
seed 6: accept 0.959, stepsize 0.00119, divergent 0, distinct draws of b_Intercept 300
seed 7: accept 0.570, stepsize 0.00397, divergent 0, distinct draws of b_Intercept 297
seed 8: accept 0.922, stepsize 0.00309, divergent 0, distinct draws of b_Intercept 300
arm lane
seed 1: accept 0.770, stepsize 0.00349, divergent 0, distinct draws of b_Intercept 300
seed 2: accept 0.925, stepsize 0.00208, divergent 0, distinct draws of b_Intercept 300
seed 3: accept 0.908, stepsize 0.00144, divergent 0, distinct draws of b_Intercept 300
seed 4: accept 0.886, stepsize 0.00363, divergent 0, distinct draws of b_Intercept 300
seed 5: accept 0.856, stepsize 0.0042, divergent 0, distinct draws of b_Intercept 300
seed 6: accept 0.959, stepsize 0.00119, divergent 0, distinct draws of b_Intercept 300
seed 7: accept 0.570, stepsize 0.00397, divergent 0, distinct draws of b_Intercept 297
seed 8: accept 0.922, stepsize 0.00309, divergent 0, distinct draws of b_Intercept 300
```

So the defect was 1 seed in 8 here (the review's seed 4 is the one),
and the floor changes no other chain: the seven that moved on base
are the same on the lane to every printed digit, as the bit-identity
argument says they must be. The OpenBLAS emulator was not run: the
change is one scalar addition whose result is the operand itself, so
it cannot depend on the BLAS.

## Item 10: summary(waic =) and add_criterion()

brms 2.23.0's `summary.brmsfit()` has no `waic` (its dots swallow it)
and `add_criterion.brmsfit()` accepts `loo`, `waic`, `kfold`,
`loo_subsample`, `bayes_R2`, `loo_R2`, `marglik`
(`dev/surface-out/brms-src.txt`). Every one is a posterior quantity;
`vignettes/brms-migration.Rmd` has `AIC()`/`BIC()` replace `loo()` on
a fit, and frmtmb refuses `loo()`, `waic()` and `bayes_R2()` on a fit.
So on a fit: `summary(waic =)` and each criterion of `add_criterion()`
are refused by name with the reason, naming `AIC()`/`BIC()` and the
draws route. An "aic" criterion was not invented: brms has no such
name, and `AIC()` needs nothing stored.

On draws, `add_criterion.frmtmb_draws()` follows brms: compute,
store in `x$criteria`, keep unless `overwrite`, save to `file`; and
`loo()`/`waic()` with no further argument return the stored object
(brms's `use_stored`). `check.txt`:

```text
add_criterion stored: loo waic | loo(ds2) is the stored object: TRUE
```

## Tests

New files, each seen to fail on rellib-r6 first. Base, the final
versions of the files (`dev/surface-tests/final-*-base.log`):

```text
RESULT frmtmb.sample test-add-criterion.R pass=5 fail=9 error=5 skip=0 warn=0
RESULT frmtmb test-brms-postfit.R pass=1 fail=2 error=7 skip=0 warn=0
RESULT frmtmb test-covstruct-floor.R pass=4 fail=7 error=1 skip=0 warn=0
RESULT frmtmb test-cs-bar-brms.R pass=4 fail=5 error=0 skip=0 warn=0
RESULT frmtmb test-disc-no-fixed.R pass=2 fail=0 error=3 skip=0 warn=0
RESULT frmtmb test-multiple-methods.R pass=0 fail=1 error=6 skip=0 warn=0
RESULT frmtmb test-update-prior-drop.R pass=4 fail=0 error=1 skip=0 warn=0
```

The behavioral failures, not only missing symbols: `stancode()` with
brms loaded reached brms's "Data must be specified"; `update()` stopped
with "No random-effect correlations match class=cor"; `fitted()` with
"requires numeric/complex matrix/vector arguments"; the child process
with "variable lengths differ (found for '.frm_cs(1)')"; the gp and ou
log densities were `Inf`; `summary()` of a `frm_multiple()` was not a
`summary.frmtmb_multiple`; `plot()` of draws refused. The passes on
base are the controls each file carries (the case where the guarded
thing is absent: a prior that still applies, a disc predictor with a
fixed column, the refusal without brms, the floor's bit identity, the
item-5 guard).

Lane, in the full run below; the passes of the absent-case controls are
the same tests.

Existing tests edited: `test-portability.R` (FN-11 expected the
`conditional_effects()` refusal; it now expects the pooled object and
the refusal of `band = "boot"`), core and frmtmb.sample
`test-generic-collision.R` (the moved generics), frmtmb.sample
`test-loo.R` (expected `plot()` of draws to refuse).

## The suites

All eight, one file per process, every gate on (`NOT_CRAN`,
`FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_DRMTMB_FIT_TESTS`, `FRMTMB_FUZZ`), 24
processes at a time (`dev/surface-suite.sh`, runner
`dev/surface-runtest.R`), then the files the last two edits reach run
again on the final build (`dev/surface-suite-rerun1/`). Generated by
`dev/surface-suite-sum.R full1 rerun1` (`dev/surface-out/suite-sum.txt`):

```text
Files run again in dev/surface-suite-rerun1 :
- frmtmb--test-arg-refusal: first pass 113 fail 1 err 0 warn 0, rerun pass 114 fail 0 err 0 warn 0
- frmtmb--test-brms-postfit: first pass 59 fail 0 err 0 warn 0, rerun pass 59 fail 0 err 0 warn 0
- frmtmb--test-generic-collision: first pass 61 fail 0 err 0 warn 0, rerun pass 61 fail 0 err 0 warn 0
- frmtmb--test-message-uniqueness: first pass 5 fail 2 err 0 warn 0, rerun pass 6 fail 0 err 0 warn 0
- frmtmb--test-multiple-methods: first pass 39 fail 0 err 0 warn 0, rerun pass 39 fail 0 err 0 warn 0
- frmtmb--test-perf: first pass 1 fail 2 err 0 warn 0, rerun pass 3 fail 0 err 0 warn 0
- frmtmb--test-portability: first pass 94 fail 1 err 0 warn 0, rerun pass 96 fail 0 err 0 warn 0
- frmtmb.sample--test-add-criterion: first pass 38 fail 0 err 0 warn 0, rerun pass 38 fail 0 err 0 warn 0
- frmtmb.sample--test-generic-collision: first pass 61 fail 0 err 0 warn 0, rerun pass 61 fail 0 err 0 warn 0

files: 361 | with a RESULT line: 361 | whose frmtmb is the lane's: 361

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 221 | 17390 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 49 | 2653 | 0 | 0 | 1 | 0 |
| frmtmb.coupling | 11 | 542 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 29 | 1743 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 10 | 360 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 15 | 501 | 0 | 0 | 2 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.spline | 15 | 590 | 0 | 0 | 1 | 0 |
| **total** | **361** | **24328** | **0** | **0** | **15** | **0** |

files with a failure, an error, a warning or no RESULT line: 0
files with a skip: 7 (the seven test-scale.R files, the scale tier's gate)
```

What the first-pass failures were: `test-arg-refusal.R` caught
`add_criterion.frmtmb_fit()` with a `...` it never read (it now takes
brms's `model_name`, `overwrite`, `file`, `force_save` and checks its
dots); `test-message-uniqueness.R` caught the pooled summary's `mc_se`
refusal duplicating the fit's text; `test-portability.R` FN-11 was the
deliberate change above; `test-perf.R` is the wall-clock bound filed
at 0.68.0 (20.02 s against 20 s and 3.1 against 2.0 under 24 parallel
processes), 3 of 3 alone. The ported brms-suite files are all green,
so no ledger row flips.

## The ports

The vigport record's counts are from 0.67.0; 0.68.0 did not rerun the
ports. So the before arm is rellib-r6 (0.68.1), measured here with the
same scripts (`dev/surface-port.sh r6` and `lane`: `run-all.R 120 raw`
and `spell`, `summarize.R`, `_run-all.sh`), and it reproduces the
record's 0.67.0 numbers exactly (post CLEAN 63, post spell OK 71 of
102, hand ok 448 of 539). Generated by `dev/surface-portcmp.R`
(`dev/surface-out/portcmp.txt`):

```text
## mechanical port, spell pass
model spell OK: r6 39 -> lane 39 of 42; class CLEAN: r6 36 -> lane 36
post  spell OK: r6 71 -> lane 72 of 102; class CLEAN: r6 63 -> lane 64
other spell OK: r6 97 -> lane 97 of 106; class CLEAN: r6 97 -> lane 97
  post  CLEAN     63 ->  64
  post  FAIL      25 ->  24
rows whose spell status or message changed: 11
regressions (spell OK on r6, not on the lane): 0

## hand translation
rows: 539 | ok: r6 448 -> lane 452
rows whose outcome or message changed: 16
regressions (ok on r6, not on the lane): 0
```

The rows that run now: mechanical `brms_missings.7.1`
(`conditional_effects(fit_imp1, "age:chl")`); hand
`pp_mixture(fit_mix2)`, `conditional_effects(fit_imp1, 'age:chl')`,
SAMPLE `plot(fit1, N = 2, ask = FALSE)` and SAMPLE
`add_criterion(fit1, 'loo')`. The other changed rows are refusals by
name where the message was generic or wrong: the 5 `plot()` calls with
`N`, `variable`, `regex` (gap rank 3), the 3 `add_criterion()` calls
(rank 5), `summary(waic = TRUE)` (rank 10), `stancode()` and
`standata()` ("Data must be specified" before). The count of failing
calls these refusals leave is deliberate: each has no maximum-likelihood
meaning.

One hand row moved from one refusal to another: SAMPLE
`plot(fit1, variable = 'simo', regex = TRUE)` now runs `plot()` of
draws and stops with brms's own "No valid variables selected.",
because frmtmb.sample names the `mo()` simplex `zeta1_1`, `zeta1_2`
(the internal softmax coordinates) where brms has
`simo_moincome1[1..3]` (`dev/surface-simo.R`). That naming gap is not
this lane's and is filed below.

## R CMD check

`dev/surface-check.sh` (built with vignettes, `--as-cran`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH,
`R_LIBS` = lane library, rellib-r6, user library), logs
`dev/surface-check/<pkg>.Rcheck/00check.log`:

- frmtmb: `Status: 2 NOTEs`. The V8 math-rendering NOTE on the HTML
  manual (expected), and the examples-timing NOTE for
  `frm_hazard_reads` (user 1.64 s, elapsed 13.66 s), `VarCorr` (1.49,
  9.37) and `dharma_residuals` (1.10, 21.39): elapsed many times user
  is load, other lanes running beside it, and none of the three is an
  example this lane touched.
- frmtmb.sample: `Status: OK`.

## Files touched

Core, shared files and the functions edited in them:

- `R/conditional-effects.R`: `plot.frmtmb_fit()` (and its roxygen),
  new `plot_fit_brms_args`; `ce_display_kind()` (the categorical
  message); `conditional_effects.frmtmb_fit()` (sets the
  `pool_scale` attribute), new `ce_pool_scale()`.
- `R/methods-fit.R`: `brms_summary_args` (the `waic` entry).
- `R/multiple.R`: `frm_multiple()` (row names, roxygen),
  `plot.frmtmb_multiple()` (message), `conditional_effects.frmtmb_multiple()`
  (refusal replaced by pooling); new `multiple_row_names()`,
  `fixef.frmtmb_multiple()`, `multiple_pool_fixef()`,
  `multiple_block()`, `summary.frmtmb_multiple()`,
  `print.summary.frmtmb_multiple()`.
- `R/predict.R`: `lp_eta_design()`; new `lp_fixed_part()`.
- `R/covstruct.R`: the `nll` of `ou`, `homcs`, `homtoep`,
  `spatial_entry()`, `gr_cov` (one coefficient) and `gp`; new
  `sd2_floored()`.
- `R/confint.R`: `update.frmtmb_fit()` (and its roxygen).
- `R/fit.R`: `frm()` (one call before `fit_assembled()`).
- `R/priors.R`: new `prior_drop_unmatched()`.
- `R/frame.R`: `check_cs_in_bar()`.
- `R/simulate-newdata.R`: one line of `xb`, through `lp_fixed_part()`
  (in the function that builds the per-dpar designs there).
- `R/generic-owners.R`: `frm_generic_owners` (four names).
- New `R/brms-postfit.R`: `stancode()`, `standata()`, `pp_mixture()`,
  `add_criterion()` and their fit methods, `mixture_prob_band()`,
  `add_criterion_reasons`, `add_criterion_names()`.
- `NAMESPACE`, `man/` (new `stancode.Rd`, `pp_mixture.Rd`,
  `add_criterion.Rd`, `frm_multiple-methods.Rd`; changed
  `frm_multiple.Rd`, `plot.frmtmb_fit.Rd`, `update.frmtmb_fit.Rd`),
  `NEWS.md`, `vignettes/brms-migration.Rmd`.
- Tests: new `test-brms-postfit.R`, `test-multiple-methods.R`,
  `test-update-prior-drop.R`, `test-disc-no-fixed.R`,
  `test-cs-bar-brms.R`, `test-covstruct-floor.R`; edited
  `test-generic-collision.R`, `test-portability.R`.

frmtmb.sample:

- `R/methods-draws.R`: the `pp_mixture` and `stancode`/`standata`
  generics removed; `pp_mixture.frmtmb_draws()` (brms's names);
  `plot.frmtmb_draws()` rewritten (brms's display).
- `R/loo.R`: `loo.frmtmb_draws()`, `waic.frmtmb_draws()` (stored
  criterion); new `draws_stored_criterion()`, `add_criterion_options`,
  `add_criterion.frmtmb_draws()`.
- `R/generic-owners.R`: `sample_generic_owners` (three names out),
  header count.
- `R/reexports.R`: four re-exports.
- `NAMESPACE`, `man/` (new `add_criterion.frmtmb_draws.Rd`,
  `plot.frmtmb_draws.Rd`, `pp_mixture.frmtmb_draws.Rd`; removed
  `pp_mixture.Rd`; changed `frmtmb-draws-refusals.Rd`, `reexports.Rd`),
  `NEWS.md`, `vignettes/brms-posterior.Rmd`.
- Tests: new `test-add-criterion.R`; edited `test-generic-collision.R`,
  `test-loo.R`.

dev/: `surface-*` scripts and this file.

## What I did not do

- `update(prior =)` does not merge with the stored prior as brms does
  (item 3).
- The pooled `conditional_effects()` does not pool `method = "predict"`
  or the bootstrap and profile bands, nor a reported dpar whose band
  scale is chosen per grid (a mixing weight); all refused by name.
- `pp_mixture()` on a multivariate fit or with `newdata` is refused.
- No version numbers were changed. frmtmb.sample's DESCRIPTION still
  says `frmtmb (>= 0.68.0)`; it needs the next frmtmb (it imports the
  moved generics), so the floor must be raised at consolidation.
  Bumps: frmtmb minor (new exported generics and methods, breaking
  name changes), frmtmb.sample minor (generics moved, `plot()` and
  `add_criterion()` new), floor the new frmtmb.
- `x$rhats` of a `brm_multiple()` fit (the vignette rounds it) has no
  counterpart; the m stored fits carry their own convergence codes.
- The floor was not extended to blocks with several scales (`us`,
  `diag`, `cs`, `toep`, `hetar1`). CORRECTED at punch round 1: the
  first version said frmtmb.sample samples them non-centered by
  default; on the fit route every prior is flat and `ncp_plan()` keeps
  them centered. The review measured `us` with two or more
  coefficients and `cs` at `+Inf` and `ar1` at `NaN` at log sd
  -1137.64 on both builds; latent on 16 sampler seeds. Filed in
  `dev/test-backlog.md`.
- The gap ranking of `dev/vigport-gaps.R` was not regenerated; the
  port comparison above lists every row that moved.

## Defects found, not fixed

- **frmtmb.sample names the `mo()` simplex `zeta1_1`, `zeta1_2`**, the
  internal softmax coordinates, where brms's draws have
  `simo_moincome1[1]` to `[3]` on the simplex
  (`dev/surface-simo.R`). `plot(ds, variable = "simo", regex = TRUE)`
  and any brms script that selects `simo_` find nothing. The hand
  translation's 0.42.0 label already said so for `summary()`.
- **frmtmb.sample's exact `gp()` sampling mixes poorly from either
  route** (R-hat 1.56 to 1.89 at 600 iterations, 197 to 215 transitions
  over the tree depth, `dev/surface-out/sample-*.txt`): the field is
  sampled centered, which is a funnel; brms samples the standardized
  `zgp`. Known (gpby r4); a non-centered form for `gp()` would need the
  kernel's Cholesky factor per draw.

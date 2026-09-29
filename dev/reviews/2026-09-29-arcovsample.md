# Review of lane wt-arcovsample, 2026-09-29

Adversarial review. Worktree
`C:\Users\adf44\source\r\frmtmb-wt-arcovsample`, base `9b4bb650`. Lane
build `C:/Users/adf44/source/r/wt-arcovsample-lib`, read only; shared
base build `C:/Users/adf44/source/r/rellib-r3`, read only. R 4.6.1,
brms 2.23.0, rstan 2.32.7, StanHeaders 2.39.1, tmbstan 1.2.1, loo
2.10.1, bayesplot 1.16.0, mvtnorm 1.4.2. No `pinlib` on any path;
`R_MAKEVARS_USER` set to `C:/Users/adf44/Documents/.R/Makevars.win`,
whose `-std=gnu++17` line was asserted before every sampled or fitted
result. No Stan cache, so every brms fit below compiled fresh.

Every number carries the script that made it and its seed. The
reviewer's scripts are `dev/arcovsample-rev-*` and their output is in
`dev/arcovsample-rev-log/`. No package file was edited.

## Verdict

**MERGEABLE**, with one blocking documentation fix of one sentence and
one blocking item for the consolidating session, not for the lane.

The substance holds. I could not break the row density: on my own
design and seeds, 15 constructions agreed with the taped objective or
its R composition to 7.1e-14 absolute at worst (8.5e-16 relative), and
6 of them agreed with brms 2.23.0 cell by cell to 5.3e-15 against cells
whose own spread is 0.5 to 1.0. The Student-t `rescor` fix is real, and
larger than the lane reported: the reference build is wrong at 8 of 15 `nu`
values, by up to 34504 log units, and `NaN` at `nu = 1e306`. Every test
count the lane reported reproduced exactly, on both builds, in both
tiers.

Blocking:

1. `?sample-log_lik`'s new section says each group's first `max(p, q)`
   rows take the unshifted mean. Only the FIRST row of each group does,
   at every order measured, including `max(p, q) = 3`. New text, and it
   is the sentence that says what a row conditions on.
2. For the consolidating session: `frmtmb.sample/DESCRIPTION` still
   reads `Imports: frmtmb (>= 0.64.0)`, and that combination installs,
   loads and then dies at the first call. It is not the lane's job to
   pick the number, but the artifact as it stands permits a broken
   install and the failure reaches `posterior_predict()` as well as
   `log_lik()`.

Nits are at the end. Two PRE-EXISTING defects the lane did not file are
recorded in section 9.

## 0. The instrument

`dev/arcovsample-rev-00-integrity.R`. Every function the diff touches
was deparsed from the installed namespace and compared with a
`deparse()` of the same definition parsed out of the worktree source:
`arma_cond_resp`, `arma_cond_dpars`, `rescor_row_loglik`,
`draws_chain_id`, `draws_loglik_factors`, `draws_row_loglik`,
`posterior_predict.frmtmb_draws`. **7 same, 0 differing.** The lane
build is the worktree source, so I did not build my own.

Both builds carry `frmtmb 0.64.0` and `frmtmb.sample 0.12.0`, so a
version string cannot tell them apart and a wrong `.libPaths()` order
would be invisible. Every reviewer script therefore prints whether
`arma_cond_resp` is in `getNamespaceExports("frmtmb")` and
`stopifnot()`s it against the arm it was asked for. In a clean process
the reference build has neither new name, exported or internal.

## 1. Row density identity. HOLDS.

`dev/arcovsample-rev-03-rows.R`, seed 5701, log
`dev/arcovsample-rev-log/03-rows.txt`. N = 42 in 6 groups of UNEQUAL
length (4, 7, 11, 5, 9, 6) with two interior time points removed from
each, rows shuffled. 200 NUTS draws per model, one chain, EVERY draw
checked.

Three instruments, and they are not the same relation:

- **(1) IDENTITY, tape form.** For a model with no random effect
  `fit$obj$fn(p)` IS the taped negative log likelihood, so the row sum
  must equal `-obj$fn(p)`. The internal parameter vector is rebuilt from
  each draw by inverting `parList()`, and the inversion is checked by a
  round trip on `last.par.best` (`identical()`) before it is used. This
  is the SAMPLING API on one side and the compiled tape on the other.
- **(2) IDENTITY, R form.** `build_objective(frame)(estimates)` is the R
  composition RTMB tapes, and for a random-effect model it is the JOINT
  nll, so the row sum must equal it minus each block's own prior. Used
  where `obj$fn` is a Laplace marginal instead.
- **(3) MEASUREMENT, and the only one that is not frmtmb's own code.**
  brms 2.23.0's generated model block for `arma()`, transliterated in
  this review from `make_stancode()` output read here, driven by the
  draws matrix's own brms-scale columns. `J_lag` and the `order(gr,
  time)` sort were confirmed against `make_standata()` rather than
  assumed. The aterm cases add their own closed forms.

```
construction                        (1) tape   (2) R form  (3) brms   cell sd
gaussian arma(2,2), p=q=2, no RE   2.84e-14   0            4.44e-15   0.628
gaussian ar(2), sigma ~ z          2.84e-14   0            0          0.663
gaussian ma(2), weights(w)         2.84e-14   0            3.55e-15   1.010
gaussian arma(2,2), cens(cc)       2.13e-14   0            6.22e-15   0.705
gaussian ar(2), trunc(lb = 0)      2.13e-14   0            4.44e-16   0.657
student arma(2,2) + (1|g), nu fit  (Laplace)  7.11e-15     3.38e-14   1.034
rescor + ar(1) both, gaussian      7.11e-14   -            2.22e-15*  0.882
rescor + ar(1) both, student       5.68e-14   -            see below  0.854
```

Relative residual against the tape is 5.0e-16 to 8.5e-16 throughout.
No non-finite cell anywhere. `loo()` ran on all of them with finite
pointwise `elpd_loo`. The row sums move across draws (sd 1.60 to 3.71),
so the identity is not satisfied by a constant.

`*` the rescor cell reference is `mvtnorm::dmvnorm` / `mvtnorm::dmvt`
at the same `mu`, `sigma` and `C`, with `mu` from the review's own
transliteration. For the STUDENT rescor row the reference disagreed by
up to 5.23, and the reference is what is wrong: the sampled `nu` reached
1.012e304 and `mvtnorm::dmvt` forms `lgamma()` as written, which
section 2 measures as the failing side. The tape identity holds there at
5.68e-14.

### Edge constructions. HOLDS.

`dev/arcovsample-rev-10-edges.R`, seed 501,
`dev/arcovsample-rev-log/10-edges.txt`. Residual against instrument (2)
at every draw:

- a group with ONE row, a group with two, a group with interior gaps, a
  group of nine, at `(p, q)` = (1,0), (2,0), (0,2), (2,2): **0** in all
  four.
- two responses, `ar(2)` on ONE of them, no rescor: **0**. The other
  response is left alone.
- `ar(1)` on BOTH responses, no rescor: 3.55e-15.
- one response `cov = FALSE` and the other `cov = TRUE`: refused, and
  the message names `y2` and `ar(t, g, cov = TRUE)`, the response that
  is the problem. Correct.

### Against brms 2.23.0, cell by cell. HOLDS.

`dev/arcovsample-rev-09-brms.R`, seed 9301, N = 39 in 6 unequal groups
with interior gaps, rows shuffled. brms is fitted (1 chain, 300
iterations, 100 kept) and then every one of its 100 draws is
transplanted into a frmtmb draws object by parameter name, so both
packages report a row density at the SAME parameter vector and nothing
rests on two samplers agreeing. The transplant STOPS if any frmtmb
column is left unfilled, and `posterior_epred()` is compared first
because it is a function of every transplanted parameter. brms's
`log_lik()` undoes its own sort, so both matrices are in the user's row
order; that is checked rather than assumed, by epred.

These are the constructions the lane did NOT run.

```
model                        epred max|d|   log_lik max|d|  rel to cell sd
gaussian arma(2,2)           1.11e-15       5.33e-15        7.3e-15
student arma(1,1) + (1|g)    8.88e-16       3.11e-15        3.6e-15
gaussian ma(2), weights(w)   4.44e-16       3.11e-15        3.7e-15
gaussian arma(2,2), cens(cc) 8.88e-16       1.78e-15        2.5e-15
gaussian ar(2), trunc(lb=0)  0              0               0
gaussian ar(2), sigma ~ z    0              0               0
```

All six brms fits compiled Stan fresh in this run.

`elpd_loo` from `loo::loo()` on the two matrices agrees to nine printed
digits in every row, which is an IDENTITY given the cell agreement and
not a second result.

The groups' FIRST rows, where brms's `e_s = 0` convention could have
differed silently: max 1.11e-15 over the 6 columns in the student case
and **0** in the other four.

`identical()` on the two matrices is FALSE even where max|d| is 0,
because brms's matrix carries dimnames and frmtmb's does not. The VALUES
are bitwise equal there.

## 2. The Student-t `rescor` fix. HOLDS, and is larger than reported.

### The arithmetic, against a third party

`dev/arcovsample-rev-01-tnu.R`, seed 99, K = 3, n = 6,
`dev/arcovsample-rev-log/01-tnu.txt`. Two references, because neither
covers the whole range: `mvtnorm::dmvt(log = TRUE)` is the same density
but is itself built on a naive `lgamma()`, so it is a reference only up
to about `nu = 1e6`; `mvtnorm::dmvnorm` is the `nu -> Inf` limit, which
the t approaches as O(1/nu). MEASUREMENTS, both outside frmtmb.

```
     nu   new vs dmvt   new vs dmvnorm   new vs mvt_std   OLD vs the
                                            _loglik       better ref
  1.5e0   8.88e-16      0.610            0                8.88e-16
  3.0e0   1.78e-15      0.367            0                1.78e-15
  1.0e1   3.55e-15      0.140            0                1.78e-15
  1.0e2   3.38e-14      0.0259           0                2.04e-14
  1.0e4   7.34e-12      2.81e-4          0                3.73e-14
  1.0e6   9.31e-11      2.81e-6          0                2.55e-10
  1.0e8   8.21e-08      2.81e-8          0                4.38e-08
  1.0e12  1.55e-03      2.81e-12         0                5.61e-04
  1.0e20  3.33          1.15e-14         0                68.04
  1.0e306 NaN           1.11e-13         0                NaN
```

The new branch equals the objective's own `mvt_std_loglik()` to **0**
at every `nu`, which is an IDENTITY (same two helpers, same order). It
matches `mvtnorm::dmvt` where `dmvt` is trustworthy, and thereafter it
is `dmvt` that leaves the truth: the residual against the gaussian
limit is 2.81/nu to four digits at every `nu` from 1e4 to 1e10, which is
the O(1/nu) convergence the density has.

**The small-nu end did not regress.** At `nu` in 1.05, 1.5, 2, 2.5, 3,
4, 7, 10, 30, 50 the new form against `mvtnorm::dmvt` has relative
residual 0 to 9.36e-16, that is 0 to 4 ulp.

### Through the exported function, on both builds

`dev/arcovsample-rev-13-rrl.R`, seed 777, K = 2 with a distributional
`sigma` on one response. ONE parameter vector per `nu`, through the
`logm1` link, given to the objective and to `rescor_row_loglik()`.

```
     nu     reference build resid      lane build resid
  3.0e0     2.84e-14                   2.84e-14
  1.0e2     3.13e-13                   2.84e-14
  1.0e3     1.83e-11                   2.84e-14
  1.0e5     1.53e-10                   0
  1.0e7     6.65e-07                   0
  1.0e9     1.43e-05                   0
  1.0e11    1.05e-02   BAD             2.84e-14
  1.0e15    107.7      BAD             2.84e-14
  1.0e20    2267.9     BAD             0
  1.0e50    5721.8     BAD             0
  1.0e100   11478      BAD             0
  1.0e200   22991      BAD             0
  1.0e300   34504      BAD             2.84e-14
  1.0e306   NaN        BAD             0
```

Reference build: **8 of 15** `nu` values disagree with the objective by
more than 1e-6 relative or are not finite. Lane build: **0 of 15**,
residual at most 2.84e-14 on a value of about 149.6, that is 1.9e-16
relative. The lane's own figures (5.6e-6 at 1e10, 227 log units at
1e20) are its `tnu.txt` design, K = 2 and n = 5; the magnitude is
design dependent and mine is larger. Not a contradiction, but the NEWS
figures are specific to a construction the NEWS does not name.

`logLik()` and the fit are unaffected: `-obj$fn` in the table above is
the same on both builds.

### The pin was seen to fail

`dev/arcovsample-rev-run.R`, one file per process, `NOT_CRAN=true`:

- `tests/testthat/test-mv-gaps.R` reference build: **pass=53 fail=4**.
  The four are behavioral, not a missing symbol: 1.34e-4 at `nu = 1e10`,
  5443 at 1e20, 82810 at 1e300, plus the convergence assertion at
  77367. Quoted in `dev/arcovsample-rev-log/mv-gaps-ref.txt`.
- Lane build: **pass=57 fail=0 err=0 skip=0**.

Exactly the lane's numbers.

### The downstream consequence

`dev/arcovsample-rev-12-rescornan.R`, seed 777, a Student-t `rescor`
model with NO autocorrelation and `nu` walked over 15 values to 1e306.
On the lane build `log_lik()` is finite at all 15, 0 non-finite cells of
750, and the row sum converges to the gaussian limit and stays there.
I could not run the same file on the reference build, for a reason that
is itself the lane's other fix: see section 3.

## 3. No behavior change elsewhere. HOLDS.

`dev/arcovsample-rev-04-pp.R` runs the same code on the same
DETERMINISTIC draws objects (`stanfit = NULL`, a fixed parameter matrix
built from a fixed seed, every column filled or the script stops) on
each build and saves the result;
`dev/arcovsample-rev-04b-cmp.R` compares them. Each
`posterior_predict()` call is wrapped in its own `set.seed(707)`.
Seed 808, N = 30 in 5 unequal groups, 6 draws.

```
case                          posterior_predict   posterior_epred
cov = FALSE arma(1,1)         identical()         identical()
cov = TRUE ar()               identical()         identical()
no autocorrelation            identical()         identical()
rescor, two responses         identical()         identical()
a laplace draws object        identical()         identical()
(1 | g), full draws           identical()         identical()
```

`identical()`, not a tolerance. The laplace case is `identical()` on
both builds and all NaN on both builds; see section 9, that is
pre-existing and not this lane's.

`dev/arcovsample-rev-log/draws-methods-{lane,ref}.txt`:
`test-draws-methods.R` is **pass=148 fail=0** on BOTH builds, which is
the same file that exercises `posterior_predict()` on the shapes the
diff touched.

### The refusals

Same script. Lane messages, checked against the form they name:

- `ar(t, g, cov = TRUE)`: quotes the label, says "an R-side residual
  correlation MATRIX", and adds "brms's default cov = FALSE form of
  ar() does factor ... refit without cov = TRUE". True: section 1
  measures that form against brms.
- `cosy(t, g)` and `unstr(t, g)`: same message with the label swapped
  and NO `cov = FALSE` sentence. Asserted by grep on the message text,
  both ways.
- `mi()` with an `ar()` term: the in-model-imputation message, not the
  autocorrelation one. On the reference build the same model gave the
  autocorrelation message, because `cov = FALSE` was refused one step
  earlier. Both messages were true when they fired; the lane's is the
  more specific.
- `loo()` reaches all of these through `log_lik()` and reports the same
  text. Checked for all five models.
- a `cov = FALSE` model: NO error on the lane build, the old
  autocorrelation refusal on the reference build.

`test-message-uniqueness.R` lane: **pass=6 fail=0 err=0 skip=0**.
`test-conditions-census.R` lane: **pass=11 fail=0**. Core
`test-conditions.R` lane: **pass=150 fail=0**.

### `draws_chain_id()`

The `NULL`-`stanfit` guard is not cosmetic and I saw it fail. On the
reference build, `dev/arcovsample-rev-12-rescornan.R` dies before it
reaches any density:

```
Error in x$stanfit@sim :
  no applicable method for `@` applied to an object of class "NULL"
Calls: log_lik -> log_lik.frmtmb_draws -> draws_chain_id -> %||%
```

A behavioral failure, on a construction written for another purpose.

## 4. `pp_check()`'s `loo_*` types. HOLDS. A shipped defect.

`dev/arcovsample-rev-07-ppcheck.R`, seed 31, THREE models so that
"every model" is a measurement: a plain `y ~ x` gaussian fit, one with
`(1 | g)`, and one with `ar(t, g)`. Draws from a fixed matrix, 20 draws,
`ndraws = 10`.

On the REFERENCE build and on the lane build, character for character
identical, on all three models:

```
[dens_overlay]    OK ggplot2::ggplot
[stat]            OK ggplot2::ggplot
[loo_pit_overlay] rlang_error: One of 'lw' and 'psis_object' must be
                  specified.
[loo_pit]         getvarError: argument "lw" is missing, with no default
[loo_intervals]   getvarError: argument "psis_object" is missing, with
                  no default
[loo_ribbon]      getvarError: argument "psis_object" is missing, with
                  no default
```

4 of 6 types tried, 3 of 3 models, 2 of 2 builds. The cause is
structural and I confirmed it by reading
`pp_check.frmtmb_draws`: it resolves the bayesplot function and calls
it with the response and the predictions, and nothing anywhere in it
computes `lw` or a `psis_object`, which `bayesplot::ppc_loo_*` require
and which `brms:::pp_check.brmsfit` computes from `log_lik()`. Filed,
not fixed, correctly: PRE-EXISTING, four types, every model, outside
this lane's gap. The lane saw the real thing.

## 5. Tests. HOLDS. Every count reproduced.

`dev/arcovsample-rev-run.R`, ONE test file per R process,
`NOT_CRAN=true`, gated runs also `FRMTMB_BRMS_FIT_TESTS=true`. Every
run prints which library it loaded and asserts the new exports are
present or absent as the arm requires. Skips were read, not counted.

```
arm    tier    package        file                       pass fail err skip
ref    plain   frmtmb.sample  test-loo.R                   77    0   4    2
lane   plain   frmtmb.sample  test-loo.R                  103    0   0    2
ref    gated   frmtmb.sample  test-loo.R                   83    0   5    0
lane   gated   frmtmb.sample  test-loo.R                  111    0   0    0
ref    plain   frmtmb         test-autocor-cond.R          46    0   1    0
lane   plain   frmtmb         test-autocor-cond.R          56    0   0    0
ref    plain   frmtmb         test-mv-gaps.R               53    4   0    0
lane   plain   frmtmb         test-mv-gaps.R               57    0   0    0
ref    plain   frmtmb.sample  test-draws-methods.R        148    0   0    0
lane   plain   frmtmb.sample  test-draws-methods.R        148    0   0    0
lane   plain   frmtmb.sample  test-message-uniqueness.R     6    0   0    0
lane   plain   frmtmb.sample  test-conditions-census.R     11    0   0    0
lane   plain   frmtmb         test-autocor.R              192    0   0    0
lane   plain   frmtmb         test-conditions.R           150    0   0    0
floor  plain   frmtmb.sample  test-loo.R                   32    0  16    2
```

Identical to the lane's figures in every row it reported. The two skips
on the ungated `test-loo.R` are the two `FRMTMB_BRMS_FIT_TESTS` blocks,
read from the skip reason, and the gated run executes both with
`skip=0`. The gated runs compiled brms Stan models fresh, so the
`Makevars` assertion is load bearing and it fired.

The four reference-build errors in `test-loo.R` are behavioral, on the
old refusal, not a missing symbol: the blocks build their draws from
`frmtmb::brms_par_labels()`, which exists on 0.64.0. The gated fifth is
the brms row-by-row block.

The one WEAK pin is `test-autocor-cond.R`, which errors on the
reference build because `arma_cond_resp` does not exist there. The lane
says so itself and carries the behavioral content INSIDE the block:
the shifted mu must differ from the unshifted one by more than
0.05 sd(y), and on a two-response model only the response with the term
may move. I confirmed both assertions are there and that they are not
vacuous: `dev/arcovsample-rev-06-firstrows.R` measures the shift as
nonzero on every row but the first of each group.

### Tolerances

No new absolute numeric tolerance. Every new assertion is a ratio:
`1e3 * .Machine$double.eps` and `1e4 * .Machine$double.eps` inside
`expect_equal(tolerance =)` and `expect_lt()` (relative in waldo, and
`.Machine$double.eps` is read at run time), `1e-10 * stats::sd(llb)`,
`1e-10 * stats::sd(d$y)`, `0.05 * abs(sum(plain))`,
`0.05 * stats::sd(d$y)`, `64 * .Machine$double.eps * abs(vals[4])`,
`0.01 * abs(vals[4])`. The bare `tolerance = 1e-10` and `1e-12` lines in
`test-loo.R` are all pre-existing; `git show HEAD:` confirms.

### Gating

The new brms block calls `skip_unless_brms_fit()` then `skip_sampler()`,
which is the same pair, in the same order, as the existing brms block
20 lines below it. Consistent. `skip_sampler()` is redundant there,
because the block builds its draws without a sampler; harmless, and
listed as a nit.

I did NOT rerun the whole frmtmb.sample suite or `R CMD check --as-cran`.
The consolidating session runs the authoritative check, and the lane's
section 7 records its own.

## 6. Docs. ONE BLOCKING ERROR, otherwise holds.

`dev/arcovsample-rev-05-rd.R` renders the four topics with
`tools::Rd2txt` and greps the RENDERED text;
`dev/arcovsample-rev-05b-rd.R` does the same for the diff's ADDED lines
alone, because whole-file counts carry pre-existing text.

Over 583 added lines: **0** em dashes, **0** en dashes, **0** curly
quotes, **0** British spellings against a 17-word list. The 21 hits for
` - ` are all arithmetic or Stan code (`Y[n] - mu[n]`,
`lgamma(nu / 2)`, `log(nu - 1)`), not a hyphen standing in for an em
dash. Rendered output for all four topics: 0 of each bad character, 0
British spellings.

### BLOCKING: the sentence is not true of the code

`?sample-log_lik`, "Autocorrelation, and what a row conditions on",
rendered:

> A row's `mu` is the one-step conditional mean: the linear predictor
> plus a regression on the OBSERVED earlier residuals of that row's
> group, with each group's first `max(p, q)` rows taking the unshifted
> mean.

`dev/arcovsample-rev-06-firstrows.R`, seed 6161, 4 groups of 7 rows.
For each `(p, q)` it measures which within-group positions have
`arma_cond_dpars()`'s mu equal to `eval_dpars()`'s:

```
term                       max(p, q)  positions with NO shift
ar(t, g, p = 1)            1          1
ar(t, g, p = 2)            2          1
ma(t, g, q = 2)            2          1
arma(t, g, p = 2, q = 2)   2          1
arma(t, g, p = 3, q = 1)   3          1
```

Position 2 of every group is shifted in 4 of 4 groups at every order.
brms initializes `Err` to zero and then fills `Err[n + 1, i]` from
`J_lag[n]`, so row 2 of a group already carries row 1's residual in
column 1; rows 2 to `max(p, q)` are PARTIALLY shifted, not unshifted.
Only the first row of each group takes the unshifted mean. The source is
`extensions/frmtmb.sample/R/loo.R:237` and the rendered
`man/sample-log_lik.Rd:91`. One sentence.

This matters because it is exactly the statement the section exists to
make, and because `p = 2` and `q = 2` are the orders the change makes
newly reachable through `log_lik()`.

### The rest of the docs

Rendered and read in full:

- `?sample-log_lik`'s "Likelihoods with no per-observation column" now
  names "a residual correlation MATRIX (`cov = TRUE`, `cosy()`,
  `unstr()`)". True of the code: `autocor_is_cond()` is
  `isFALSE(ac[["cov"]])`, `parse_autocor_call()` sets `cov = TRUE`
  unconditionally for `cosy` and `unstr`, and `arma_cond_resp()`
  returned `character(0)` on both in section 1's edge script.
- `?sample-loo`'s "A time series, and what is left out": every clause
  checks out. "this is brms's elpd" is section 1's brms table, nine
  digits. "not the row's influence on its neighbors, because every
  retained column still reads `y_t` through its own lagged residuals"
  is a true description of the construction.
- `?frmtmb-sampling-api`'s new paragraph names the call order
  (`arma_cond_dpars()` before `row_lpdf()` and before
  `rescor_row_loglik()`) and says it is safe to call unconditionally.
  Both verified: `arma_cond_dpars()` on a `cov = TRUE` fit and on a fit
  with no block returned `dpv` `identical()`ly unchanged, and the rescor
  path in section 1 gets the shift, which it can only do if the order
  is as documented.
- `?frmtmb-autocor`'s added bullet is true and correctly scoped. It does
  create a forward reference: core now asserts what `frmtmb.sample`
  does, so an installation with core bumped and the extension not bumped
  has core documentation that is false. See blocking item 2. A nit, not
  a blocker on its own.

### NEWS

Both entries sit under a `(development version)` heading, which is the
FIRST heading in each file. **0** lines matching
`[0-9]+\.[0-9]+\.[0-9]+` inside either block, **0** em dashes, **0**
lines over 80 columns. Both say what a user sees change: `log_lik()`,
`loo()`, `waic()` and `psis()` accept a class of model they refused; one
refusal message changed wording; a Student-t `rescor` `log_lik()` was
`NaN` and is not; a draws object with no `stanfit` no longer errors. The
frmtmb.sample block states the dependency in words in its first line.
All true as measured.

## 7. Floors. HOLDS, and is sharper than reported.

`dev/arcovsample-rev-08-floor.R`. The worker's `frmtmb.sample` was
installed into a scratch library
`C:/Users/adf44/source/r/wt-arcovsample-rev-floorlib` with only
`rellib-r3`'s frmtmb 0.64.0 behind it.

- `R CMD INSTALL` **succeeds**, exit 0, including "testing if installed
  package can be loaded".
- `library(frmtmb.sample)` **succeeds**. `NAMESPACE` has
  `import(frmtmb)`, a whole-namespace import, and only 30 names in an
  explicit `importFrom()`, so nothing names the two new functions and
  nothing fails at load.
- The first call fails LOUDLY:
  `simpleError: could not find function "arma_cond_resp"`, from
  `log_lik()` on a `cov = FALSE` `ar()` fit.
- `posterior_predict()` on the same fit fails the SAME way. That is
  wider than the lane's account: the floor bites `posterior_predict()`
  too, which worked on 0.64.0 for this model class, so a mixed install
  breaks a method that was not part of the gap.
- `test-loo.R`: **pass=32 fail=0 err=16 skip=2**. 16 of the file's
  blocks error, including every pre-existing `log_lik()` and `loo()`
  block, not only the new ones.

So the failure is loud, immediate and total rather than silent. The
consolidating session must raise `Imports: frmtmb (>= ...)`; nothing in
the repository stops the broken pairing today.

## 8. What the lane reported that I did not re-derive

- `R CMD check --as-cran` for either package (cost; the consolidating
  session runs the authoritative one).
- the whole 36-file frmtmb.sample suite, and the gated brms-suite ports.
- the lane's own seeds 4021, 4022, 4023 and its `dev/arcovsample-*.R`
  scripts. I wrote my own designs instead, which is the stronger test of
  the same claims, and every shared count matched.

## 9. Two pre-existing defects the lane did not file

Both found while probing this lane, both identical on the reference
build, so neither is caused by the change. Recorded so they are not
re-found.

1. **`posterior_predict()` and `posterior_epred()` on
   `frm_sample(laplace = TRUE)` draws return NaN silently.**
   `dev/arcovsample-rev-11-laplace.R`, seed 1212, N = 30 in 6 groups.
   On BOTH builds: `y ~ x + ar(t, g) + (1 | g)` gives 4500 non-finite
   cells of 4500, `y ~ x + (1 | g)` gives 3000 of 4500, and the only
   signal is a repeated base warning "NAs produced". `log_lik()` on the
   same object refuses with a clear frmtmb message, which is the right
   shape; the two predict methods should refuse the same way rather than
   returning a matrix of NaN.
2. **`frm_sample(laplace = TRUE)` on a model with no random effect dies
   with an internal error.** Same script, same on both builds:
   `Error in -obj$env$random : invalid argument to unary operator`. A
   bare R error from a public argument combination, not an frmtmb
   message.

## Nits

1. `tests/testthat/test-mv-gaps.R`: the new `test_that()` title line is
   86 columns. The house rule is 80 in R sources.
2. `tests/testthat/test-mv-gaps.R`: no blank line between the new
   block's closing `})` and the `test_that()` that follows it. Every
   other pair in the file has one.
3. The new brms block in `test-loo.R` calls `skip_sampler()`, which it
   does not need: it builds its draws with `stanfit = NULL`. Harmless,
   and consistent with its neighbour, but it makes the block skip on a
   machine whose tmbstan is broken even though tmbstan has nothing to
   do with it.
4. The core NEWS entry quotes 5.6e-6 at `nu = 1e10` and 227 log units at
   `nu = 1e20` without naming the construction they came from. Those
   figures are design dependent: on my seed 777 design the same
   quantities are 1.4e-5 at 1e9 and 2268 at 1e20. The `@noRd` comment on
   `rescor_row_loglik()` does cite `dev/arcovsample-log/tnu.txt`; the
   NEWS does not.
5. `?frmtmb-autocor` in CORE now asserts what `frmtmb.sample` does. See
   blocking item 2: with core bumped and the extension not, that
   sentence is false.
6. `loo()` on a `cov = FALSE` ARMA fit returns a number with no message
   about what the leave-one-out quantity is, and the lane's own test
   asserts that silence (`expect_no_message`). It is documented in
   `?sample-loo` and brms is silent too, so parity argues for it; noted
   only because the covariance form gets a paragraph of explanation and
   this one gets none at the console.

## Final line

**MERGEABLE** once the `max(p, q)` sentence in
`extensions/frmtmb.sample/R/loo.R` and the rendered
`man/sample-log_lik.Rd` is corrected to name the first row of each
group, and on the understanding that the consolidating session raises
`frmtmb.sample`'s floor on `frmtmb` before either package ships.

---

# Re-check after punch round 1, 2026-09-29

Scope as set by the coordinator. Same harness, same read-only rules; new
scripts `dev/arcovsample-rev-14-recheck.R`, `-15-eqn.R`, `-16-eqn2.R`,
`-17-nuthresh.R`, logs `dev/arcovsample-rev-log/14-recheck.txt`,
`15-eqn.txt`, `16-eqn2.txt`, `17-nuthresh.txt`, `r2-*.txt`. No package
file was edited.

## Instrument, rebuilt library

`dev/arcovsample-rev-14-recheck.R` part (a): 5 of 5 touched functions
deparse identically between the rebuilt lane namespace and the worktree
source (`arma_cond_resp`, `arma_cond_dpars`, `rescor_row_loglik`,
`draws_row_loglik`, `draws_loglik_factors`). `R/compat.R` carries the
sentence as DATA, so it was compared as a string: the new rule is
present. Every run below printed `newexports=TRUE` from
`getNamespaceExports("frmtmb")`.

## (1) The corrected wording, rendered. RESOLVED.

`tools::Rd2txt` on `extensions/frmtmb.sample/man/sample-log_lik.Rd` and
`man/frmtmb-autocor.Rd`. The retracted claim is gone in all four
spellings I searched for: **0 hits in each topic**. The replacement, as
rendered:

```
A group starts with no residuals behind it, so lag i first reaches
the row at within-group position i + 1: the FIRST row of a group gets
no lagged term, a row at position k gets only the lags up to k - 1,
and from position max(p, q) + 1 onward a row gets all of them.
```

That matches my seed-6161 measurement of the previous round exactly
(positions with no shift: position 1 only, at every order tried) and it
matches the worker's rule.

The two-argument `\eqn` is in `man/frmtmb-autocor.Rd`;
`sample-log_lik.Rd` has no `\eqn` at all, it uses `\code`. All three
back-ends, `dev/arcovsample-rev-log/16-eqn2.txt`:

- `Rd2txt`: "from position max(p, q) + 1 onward", the ASCII second
  argument. Clean.
- `Rd2HTML`: the LaTeX first argument inside R's math span
  (`class="reqn"`), which is what every other `\eqn` in the same file
  produces.
- `Rd2latex`: the first argument kept for math mode, second dropped.
  Parses, 287 lines.
- `tools::checkRd()`: returns no condition.

The brace-or-backslash counts my sweep reported (2 in
`sample-log_lik.txt`, 9 in `frmtmb-autocor.txt`) are all PRE-EXISTING:
`\deqn` math at rendered lines 15, 16, 43, 44, 117, 122, 123 and example
code at 163, 174, 278, 280. None is in the new text. 0 em dashes, 0
British spellings in either topic.

## (2) Unit-length and short groups through `log_lik()`. HOLDS.

`dev/arcovsample-rev-14-recheck.R` part (d), seed 6162, N = 26 in groups
of 9, 7, 3, 1, 6, rows shuffled, 3 parameter vectors per model. Row sum
against `build_objective(frame)(estimates)`:

```
term                  max |rowsum - objective|   rel   finite cells
ar(p = 3)             0                          0     TRUE
arma(p = 3, q = 2)    0                          0     TRUE
ma(q = 3)             0                          0     TRUE
ar(p = 1)             0                          0     TRUE
```

Exact zero, not a rounded print. The 1-row group and the 3-row group at
`p = 3` are the cases where a group is shorter than the order; they
carry no special handling and need none.

## (3) The stronger form of the rule, tested as stated. HOLDS.

The new wording claims more than "only position 1 is unshifted": it
claims position `k` carries only the lags up to `k - 1`. I tested that
by falsification rather than by re-measuring whether a shift is present.
On `ar(t, g, p = 3)`, seed 6161, 4 groups of 7, each coefficient was
changed ALONE by 0.11 and the rows that moved were recorded:

```
lag 1 moves exactly the rows at within-group positions >= 2   TRUE
lag 2 moves exactly the rows at within-group positions >= 3   TRUE
lag 3 moves exactly the rows at within-group positions >= 4   TRUE
```

`identical()` against the predicate `pos >= i + 1` over all 28 rows, for
each of the three lags. The stronger sentence is true of the code.

## (4) Test counts. All reproduced.

`dev/arcovsample-rev-run.R`, one file per process, `NOT_CRAN=true`.

```
package        file                       pass fail err skip   worker said
frmtmb.sample  test-loo.R                  103    0   0    2   103
frmtmb.sample  test-loo.R (gated)          111    0   0    0   111
frmtmb         test-mv-gaps.R               57    0   0    0    57
frmtmb         test-autocor-cond.R          56    0   0    0    56
frmtmb         test-compat.R               639    0   0    0   639
frmtmb.sample  test-compat-preflight.R      81    0   0    0    55, 5 skips
frmtmb.sample  test-message-uniqueness.R     6    0   0    0     6
```

The gated run compiled brms Stan fresh and reported skip=0, so nothing
hid behind a skip.

`test-compat-preflight.R` is the one row that differs, and it is a
library-path effect, not a disagreement. That file lives in
`extensions/frmtmb.sample/tests/testthat/`, and the review path puts
`rellib-r3` behind the lane library, so the suggested siblings resolve.
With the lane's own two-entry path (`LANE`, then the user library) I
reproduced **pass=55 fail=0 err=0 skip=5** exactly and read all five
reasons: four `{frmtmb.ode} is not installed`, one `{frmtmb.eam} is not
installed`. Both counts are honest; the review path runs 5 more blocks
and they pass.

## (5) The backlog and the NEWS rewrite. Accurate, three record nits.

`dev/test-backlog.md`'s new section states every figure correctly
against my own logs: 4500 non-finite of 4500 and 3000 of 4500 for the
laplace predictive methods, identical on both builds; the bare
`-obj$env$random` error; 4 of 6 `pp_check()` types on 3 of 3 models and
2 of 2 builds, with the error text as I measured it; and the
`nchains.frmtmb_draws()` item with the `draws_chain_id()` failure I saw
cited by a log path that exists.

The rewritten core NEWS bullet on the opening rows is accurate and
correctly scoped ("The likelihood was always this; only the sentences
were wrong"), and the `rescor_row_loglik()` bullet's decision to quote
no single drift figure and to name two constructions instead is the
right call, since my own design gave a different magnitude.

Record nits, all in prose rather than in code:

1. **`dev/test-backlog.md` says "all three are outside that lane's
   gap"** and the section lists FOUR items (1 high, 2 medium, 1 low). It
   also says "Items 1 and 2 are wrong ANSWERS, not refusals", but item 2
   is the bare `-obj$env$random` error, which is neither. This is the
   typed-count failure `dev/lane-rules.md` names under "Generate counts
   into the document; do not type them".
2. **`NEWS.md` names `nu = 3.6e305` as where the as-written head
   "overflows to `Inf - Inf`". Measured 5.0654746e305**, bisected in log
   space (`dev/arcovsample-rev-log/17-nuthresh.txt`; the same figure at
   K = 1, 2, 3 and 5, because `lgamma()` itself overflows above argument
   2.5327373e305). At the 3.6e305 the NEWS names, the as-written head
   returns **0** where the truth is about 702.88, so that point is in
   the SILENT regime, not the `NaN` one. The silent regime starts far
   lower still: my own `13-rrl-{ref,lane}.txt` has the reference build
   107.7 log units out at `nu = 1e15`.
3. **`NEWS.md` cites `dev/arcovsample-rev-log/13-rrl.txt`, which does
   not exist.** The reviewer's files are `13-rrl-ref.txt` and
   `13-rrl-lane.txt`. A shipped NEWS file pointing at a missing path.

## The five nits from the first pass

Fixed, checked one at a time. `tests/testthat/test-mv-gaps.R` has no
line over 80 columns and has a blank line before the `test_that()` that
follows the new block. The new brms block no longer calls
`skip_sampler()` and says why in a comment. `?frmtmb-autocor` now ends
its new sentence with "Whether the installed `frmtmb.sample` does is
that package's own question", which removes the assertion that a mixed
install would make false. The core NEWS no longer quotes a single drift
figure without its construction. Over the whole punch-1 diff the only
added lines past 80 columns are the one `R/compat.R` data string, whose
file convention is one long string per row and whose line was already
351 characters before the edit, and four generated `.Rd` lines from
roxygen. No hand-written R source line is over 80.

## What is still open, unchanged from the first pass

`extensions/frmtmb.sample/DESCRIPTION` still reads
`Imports: frmtmb (>= 0.64.0)`, which permits the pairing that installs,
loads and then dies at the first call
(`dev/arcovsample-rev-log/08-floor.txt`: `could not find function
"arma_cond_resp"`, `test-loo.R` err=16, and `posterior_predict()` broken
as well). That is the consolidating session's number to pick, not the
lane's, and the lane states the dependency in words in both NEWS files.

## Re-check verdict

**MERGEABLE.** The blocking sentence is fixed, the corrected rule is
true of the code in its stronger form as well, the `\eqn` renders in all
three back-ends and `checkRd()` is silent, unit-length and short groups
give an exact-zero residual, and all seven test counts reproduce. The
three remaining items are record nits in `dev/test-backlog.md` and
`NEWS.md`, none of which changes a code path, plus the standing
DESCRIPTION floor for the consolidating session.

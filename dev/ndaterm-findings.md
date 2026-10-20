# Lane ndaterm: addition terms on newdata

Base 47cf659d (frmtmb 0.71.0). Base build `rellib-r9`, lane build
`wt-ndaterm-lib`. Every number below has its script and log; logs are
under `dev/ndaterm-out/` and `dev/ndaterm-suite-*/` (gitignored).

## The defect, verified

The reviewer of lane bin2port was right about the cause. Core
`predict_simulate()` (`R/predict-brms.R`) and frmtmb.sample's
`posterior_predict.frmtmb_draws()` built the newdata addition terms
only under `trunc()`, `rate()` or `thres(gr = )` and passed `list()`
otherwise. The binomial-type simulators read a missing `trials()` as
`%||% 1`. `arma_cond_fill_dpars()` (`R/autocor.R`) read `trunc()` only.

`dev/ndaterm-probe.R r9` against `mine` (seed 7; logs
`probe-r9.txt`, `probe-mine.txt`). Newdata `x = 0`, `N = 1000, 2`; for
`se()`, `s = 0.05, 5`; 400 draws unless stated.

Each case reads r9 (0.71.0), then the lane, then the reference.

- binomial `predict(summary = FALSE)` column means: r9 0.555, 0.545
  (max 1, 1); lane 563.30, 1.08 (max 614, 2). `simulate()` 562.987,
  1.095; `fitted()` 563.1103, 1.1262; brms 562.72, 1.15 (max 617, 2).
- binomial `predict()` Estimate: r9 0.5825, 0.5750; lane 563.1450,
  1.2125.
- zero_inflated_binomial column means: r9 0.375, 0.420; lane 417.260,
  0.765; `fitted()` 410.86, 0.822.
- beta_binomial column means: r9 0.440, 0.415; lane 401.785, 0.815;
  `fitted()` 393.22, 0.786.
- gaussian `se(s)` column sd: r9 0.029457, 0.029457; lane 0.058689,
  4.767247; brms 0.06104, 4.94885.
- custom family, `vint()` as the trial count: r9 ERROR "invalid
  arguments"; lane means 563.30, 1.08; `simulate()` 562.987, 1.095.
- custom family, `vreal()` as the known sd: r9 ERROR "invalid
  arguments"; lane sd 0.058689, 4.767247, as `se()`.
- `multinomial(K = 3)`: r9 ERROR "needs trials()"; lane row sums 1000
  and 2, newdata's N.
- wiener `dec(upper)`, bias 0.3, column means: r9 ERROR "needs the
  boundary"; lane 0.86737, 0.61116, 0.81447; `simulate()` 0.84857,
  0.60550, 0.81396.
- binomial newdata without `N`: r9 Estimate 0.55 and no condition;
  lane refused, naming `trials(N)`; brms refuses ("can neither be
  found in 'data'").

brms 2.23.0 numbers: `dev/ndaterm-brms.R` (seed 7, one chain of 300
draws), log `brms.txt`. brms also refuses newdata without `s` for the
`se(s)` model.

frmtmb.sample, `dev/bin2port-rev-newdata2.R`'s case reproduced in the
new test (seed 7, one chain, 200 draws; log
`t-sample-newdata-r9.txt`): `posterior_predict(newdata)` means 0.595
and 0.520 where `posterior_epred(newdata)` gives 576.28 and 1.153; the
`se(s)` draws have sd 0.0329 in both columns where 0.0598 and 5.0001
are due.

### Not reproducible: se() with cov = FALSE ARMA

The brief asks for an `se()` newdata test on a `cov = FALSE` ARMA
model. That model cannot be fitted: `check_autocor_cond()` refuses it
at the fit, with brms's own words ("Please set cov = TRUE in ARMA
structures when including known standard errors"; probe log, all four
builds). `cov = FALSE` is accepted for gaussian and student only, and
of their terms only `trunc()` reaches the simulator, which the old
code read. So `arma_cond_fill_dpars()` had no live defect. It now
takes every term all the same, so that the next term a family reads
there is not lost. A test pins the fit-time refusal instead
(`test-newdata-aterms.R`, "se() and cov = FALSE ARMA are refused").
The `trunc()` fill is unchanged: minimum draw -2.992 over bound -3 on
both builds.

## How old the defect is

`dev/ndaterm-probe.R <arm>` on rellib-r6 (0.68.1), r7 (0.69.1), r8
(0.70.0) and r9 (0.71.0), one process each, logs
`probe-<arm>.txt`. All four give the same numbers on every case of
the r9 column above (binomial col means 0.555, 0.545, max 1, 1;
`se(s)` sd 0.029457, 0.029457; vint and vreal "invalid arguments";
wiener "needs the boundary").

From the history (`git log -S`):

- frmtmb.sample's `posterior_predict(draws, newdata =)` passed `list()`
  from the method's first version, 65fed080 (2026-08-31, frmtmb
  0.15.0, before the extraction). d20e105f (0.22.0) added `trunc()`,
  38701fdc `thres(gr = )`, 4594a82e `rate()`.
- Core `predict(newdata =)` has drawn since 862032e3 (2026-09-22),
  when `predict()` became brms's predictive summary; that commit wrote
  the `has_trunc()` gate. Corrected in fix round 1 (review finding 5):
  862032e3 is not an ancestor of the 0.60.0 release commit df8af533
  (`git merge-base --is-ancestor` exits 1) and is one of the 0.61.0
  release commit bac1e5f6 (exits 0). The DESCRIPTION of 862032e3 still
  read 0.60.0, which is what misled the first version of this record.

So it is in every release from 0.61.0 (core) and from 0.15.0 (draws).

## What changed

Core:

- `predict_simulate()` (`R/predict-brms.R`) and
  `arma_cond_fill_dpars()` (`R/autocor.R`) call `aterms_for_newdata()`
  on every newdata call. `aterms_for_newdata()` was already the shared
  construction (`simulate()`, `fitted()`, `predict_mean_se()`,
  `ord_prob_aterms()`, emmeans' `emm_epred_part()`); the defect was the
  gates in front of it, so no new construction was needed.
- `newdata_unread_aterms` (`R/predict.R`): the terms no newdata path
  evaluates, `cens()`, `weights()`, `mi()`, the `se_sigma` flag, the
  `thres()` count, `subset()` and `index()`, as brms's
  `validate_newdata()` fills or ignores them. Shared by
  `aterms_for_newdata()` and the guard below.
- `aterm_coerce()` (`R/predict.R`): a term's newdata value coded as
  the frame coded it (registered coercion, `thres(gr = )` codes). Used
  by `aterms_for_newdata()` and by `ce_aterms()`, which had no
  coercion.
- `aterm_trials()` (`R/families.R`): binomial, zero_inflated_binomial,
  beta_binomial and zero_inflated_beta_binomial refuse a missing
  `trials()` by name in `lpdf`, `valid_y`, `init_dpars`, `mean_fn`,
  `var_fn`, `dev_fn` and `sim`, instead of `%||% 1`.
  `check_trials_given()` already refuses such a model at the fit, so
  only a post-fit caller that drops the term reaches it.
  `mixture_mu_start()` keeps its `%||% 1`: it serves beta components,
  which have no trials. `sep_certificate()` (`R/se-check.R`) keeps its
  default for the same reason (bernoulli).
- `check_family_aterms()` (`R/families.R`), called first in
  `sim_response()` and `response_mean()`: a family's declared
  `required_aterms` must be in the list it is handed, else a refusal
  naming the term. Terms in `newdata_unread_aterms` are not asked for.
- `ce_aterms()` (`R/conditional-effects.R`): `aterm_coerce()`, and
  `thres(gr = )` is strict (a grid row without its group is refused,
  not drawn without it).

frmtmb.sample: `posterior_predict.frmtmb_draws()` calls
`aterms_for_newdata()` on every newdata call. It needs nothing new
from core (`aterms_for_newdata()` is exported at 0.71.0).

frmtmb.eam: `ddm_indicator_binary()` in `ddm_indicator_mean()` and
`ddm_indicator_sim()` refuses a finite boundary value other than 0 or
1. `NA` passes.

## conditional_effects() and dec(): two more silent answers

`dev/ndaterm-ce.R <arm>` (seed 7, bias 0.3, 300 trials; logs
`ce-r9.txt`, `ce-mine.txt`). The CE estimate is the median of 2001
draws at the grid point x = -0.013.

- A `dec()` factor with levels `up, down`: the fit codes `up` as 0, the
  lower boundary (`ddm_coerce_dec()`). r9 CE at `side = up` 0.7013,
  `side = down` 0.4316, while `simulate()` medians at x = 0 are 0.4282
  and 0.6852: swapped, because `ce_aterms()` took `as.numeric()` of the
  factor. Lane: `up` 0.4316, `down` 0.7013, consistent with
  `simulate()`.
- A numeric `dec(upper)`: the grid holds `upper` at its mean, 0.5733,
  and every wiener consumer reads anything but 1 as the lower boundary.
  r9 default grid 0.4859, the same as `upper = 0` (0.4859), against
  0.7454 at `upper = 1`. Lane: refused, "the boundary indicator must be
  0 (lower) or 1 (upper) on every row, and 100 row(s) hold 0.5733",
  with the remedy. Not fixed by holding the variable elsewhere: brms
  holds it at the mean too, but brms's wiener `posterior_predict()`
  does not condition on `dec`, so it has no such defect.

## Tests, each seen to fail on r9

`dev/ndaterm-runfile.R <pkg> <file> <arm>`, one process per file.

| file | r9 | lane |
|---|---|---|
| `tests/testthat/test-newdata-aterms.R` | pass=17 fail=14 error=4 | pass=39 |
| `frmtmb.sample/.../test-newdata-aterms-draws.R` | pass=4 fail=5 | pass=9 |
| `frmtmb.eam/.../test-newdata-dec.R` | pass=3 fail=4 error=2 | pass=10 |

Logs `t-newdata-<arm>.txt`, `t-sample-newdata-<arm>.txt`,
`t-eam-dec-<arm>.txt`. Each `test_that()` block that pins a defect
fails on r9; the blocks that construct the absent case ("nothing new
fires where no term is due", the in-sample and `NA` cases) pass on
both, as they must. Tolerances are multiples of Monte Carlo standard
errors measured in the run (mean of draws: 6 se of the difference; sd
of draws: 6 x 1/sqrt(2(n - 1)) relative; a median: 6 binomial se of
the share below it).

Absent cases constructed: poisson, bernoulli and gaussian fits with no
addition term through `predict(newdata)` and `predict()` with
`expect_no_condition()`; binomial in-sample draws bounded by the
frame's `N`; `posterior_predict()` of draws in sample; wiener with
`upper = 1` given, an `NA` row, and in-sample `predict()`.

## Suites

All eight packages, one file per process, lane library first, 8 at a
time (`bash dev/ndaterm-par.sh all 8 <all eight>`, summarized by
`dev/ndaterm-summ.sh all`; per-file logs and `RESULTS.txt` in
`dev/ndaterm-suite-all/`). Every `lib:` line names the lane's core.

```
ran 385 of 385
frmtmb          files=237 pass=17479 fail=0 error=0 skip=185 warn=0
frmtmb.sample   files=54  pass=2706  fail=0 error=0 skip=4   warn=0
frmtmb.eam      files=31  pass=1818  fail=0 error=0 skip=3   warn=0
frmtmb.coupling files=11  pass=542   fail=0 error=0 skip=5   warn=0
frmtmb.latent   files=10  pass=360   fail=0 error=0 skip=2   warn=0
frmtmb.learn    files=15  pass=430   fail=0 error=0 skip=13  warn=0
frmtmb.ode      files=11  pass=549   fail=0 error=0 skip=1   warn=0
frmtmb.spline   files=16  pass=594   fail=0 error=0 skip=1   warn=0
```

Gated tier, `FRMTMB_BRMS_FIT_TESTS=true`, on the 46 files that skipped
above or are brms-suite files (`dev/ndaterm-gated-jobs.lst`, logs
`dev/ndaterm-suite-gated/`): ran 46 of 46, pass=3732 fail=0 error=0
skip=14. The 14 skips are `test-drmtmb-agreement.R` (13, gated on
`FRMTMB_DRMTMB_FIT_TESTS`) and `test-fuzz.R` (1, `FRMTMB_FUZZ`, the
grammar fuzz tier). `test-drmtmb-agreement.R` with both variables set:
pass=131 fail=0 (`dev/ndaterm-out/t-drmtmb.txt`). The brms-suite files
pass with no stale verdict, so this lane flips no row of
`dev/brmsport-verdicts.tsv` (core `test-brms-suite-methods.R` 162
passes, frmtmb.sample's 61).

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`;
a 2000 x 2000 `%*%` takes 0.09 s there against 9.39 s on the reference
BLAS), core suite one file per process (`dev/ndaterm-suite-ob026/`):
ran 237 of 237, pass=17480 fail=0 error=0 skip=185. The one extra pass
is `test-v11.R` (23 against 22), whose `if (is.finite(ll_ref))` branch
rounding decides; the base build gives 23 under the same OpenBLAS too.

One session per changed package (`dev/ndaterm-onesession.R`, logs
`dev/ndaterm-out/onesession-<pkg>.txt`):

```
ONESESSION frmtmb.eam pass=1818 fail=0 error=0 skip=3 warn=0
ONESESSION frmtmb.sample pass=2706 fail=0 error=0 skip=4 warn=0
ONESESSION frmtmb pass=17474 fail=9 error=0 skip=185 warn=0
```

The nine core failures are in `test-compat-register.R` (2) and
`test-compat.R` (3 blocks): "expects = accepts a name another package
already supplies", "a declared but unresolved rule side is reported,
not dropped", "every rule wins at least one pair", "rule patterns name
real kinds, groups, and features", "declared post-fit methods are
callable". These are the files lane-env.md names as failing in one
session on 0.71.0 (lane testiso fixes them); this lane does not touch
the compatibility registry, and both files pass one per process.

R CMD check --as-cran (`dev/ndaterm-check.sh <pkg>`, built and checked
in `dev/ndaterm-check/<pkg>/`, `R_LIBS` lane, rellib-r9, user):

- frmtmb.sample: Status: OK
- frmtmb.eam: Status: 1 NOTE (the HTML manual's V8 math rendering, the
  expected one)
- frmtmb: Status: 1 NOTE (the same V8 NOTE; tests OK in 22 min)

## Fix round 1 (review `dev/reviews/2026-10-14-ndaterm.md`)

The sections above are the first round's record. Where this section
says otherwise, it supersedes them. The suite and check counts above
are the first round's; the counts for the current tree are under
"Round 1 suites" below.

### 1. Class-wise simulators read other rows' terms (HIGH): fixed

The mixture simulator, the grouped mixture's and frmtmb.latent's hmm
`state_sim` drew each class (state) on its rows `idx` but handed the
component the whole term list, so the j-th row of a class read the
`trials()` of the j-th row of the data. Each now passes the class's
rows: `subset_obs(aterms, idx, n)` at `R/families.R` (mixture `sim`,
`mixture_sim_groups()`), and `hmm_rows_of()`, a local copy of
`subset_obs()` (not exported), in frmtmb.latent's `R/hmm.R`.

`dev/ndaterm-r1-probe.R <arm>` (the reviewer's constructions, seed 7;
logs `r1-probe-r9.txt`, `r1-probe-mine.txt`):

- `mixture(binomial(), binomial())`, N in {2, 50}, `simulate()`: share
  of draws above N 0.25033 on r9, 0 on the lane; mean draw on the
  N = 2 rows 11.767 on r9, 0.98606 on the lane, against `fitted()`
  0.98305.
- the same model, `mixture(groups = ~id)`: share 0.21263 on r9, 0; mean
  on N = 2 rows 9.2056 on r9, 0.83504, against `fitted()` 0.8388.
- `hmm(K = 2, binomial())`: share 0.25973 on r9, 0; mean on N = 2 rows
  12.204 on r9, 0.86749, against `fitted()` 0.9069.

The sweep for the same pattern (every per-class or per-state call of a
component's `sim`, `lpdf`, `mean_fn`, `var_fn`, `comp_lpdf`,
`state_lpdf`, `state_mean`, in core and the seven extensions; grep
listed in the round's transcript, the sites by name here):

- the three sites above were the only ones that subset the rows of the
  dpars and not of the terms;
- `sim_response()`'s trunc() rejection loop, `sim_autocor_cond()` and
  `predict_sim_rows()` already subset both;
- the mixture `lpdf`, `comp_lpdf`, `mean_fn`, the hmm `state_lpdf`,
  `state_mean`, `hmm_var_response()` and `hmm_mean_response()`
  evaluate every row at once, with full-length dpars and terms, so
  they are aligned (the hmm block's `y` is in frame row order,
  `hmm.R` "y = yv");
- `mvn_sim_rows()` (mixture_mvn) and `lca()` read no addition term; the ordinal
  per-category loop (`ord_cumulative_logpmf_cs()`) and the grouped
  thresholds (`thres.R`) read none either; frmtmb.coupling has no
  class-wise draw.

The same probe found a second defect in the path the lane had opened.
`predict_dpar_values()` (core) and frmtmb.sample's `posterior_predict()`
read each dpar on newdata, or under `re_formula`, with
`frm_linpred(type = "response", dpar = )`. That is the REPORTING
scale: for a mixture's `theta` it is already the softmax probability,
and the simulator applied the softmax again. Both now take the link
inverse of the link-scale predictor, the scale `dpars_natural()` gives
and `simulate()` uses. `dev/ndaterm-r1-mixgap.R` (log
`r1-mixgap.txt` before, `r1-mixgap-fix.txt` after, 4000 draws):

```
fitted 483.9209 0.9678417  Est.Error 21.03141 0.04206282
before: propagate_error FALSE means 443.1785 0.88225 (se 3.60, 0.0122)
after:  propagate_error FALSE means 486.6025 0.9695  (se 3.70, 0.0123)
simulate(newdata) means 480.249 0.9565
```

`dev/ndaterm-r1-theta.R` (gaussian mixture, weight 0.19 on a component
6 higher; logs `r1-theta-<arm>.txt`): `predict(newdata)` means 0.778,
1.818, 2.876 on r9 against `fitted()` 0.150, 1.165, 2.180; the lane
0.104, 1.194, 2.185. On draws (`dev/ndaterm-r1-theta-draws.R`):
`posterior_predict(newdata)` means 1.253, 1.753, 2.792 on r9 against
`posterior_epred()` 0.150, 1.165, 2.181.

### 2. Default conditional_effects() on numeric dec() refused: fixed

`conditional_effects()` holds the numeric variables of `thres()`,
`thres(gr = )`, `vint()` and every registered term at
`round(median())` unless `conditions` sets them, with a message
("holding the addition-term variable(s) upper at 1, ..."), beside the
trials() holding (`ce_int_vars()`, `R/conditional-effects.R`). brms
2.23.0's `get_int_vars()` holds trials, thres and vint variables so
(read and run by the reviewer, `ndaterm-rev-brms.R`), whether or not a
predictor reads the variable too, and frmtmb follows that rather than
excluding shared variables. A registered term is added because
frmtmb's families read it as data. The default calls
`conditional_effects(fit)` and `conditional_effects(fit, "x")` run on
a numeric-`dec()` wiener fit, epred and predict, and the default
equals `conditions = data.frame(upper = 1)`. The refusal stays for an
explicit newdata value (`fitted(newdata = data.frame(upper = 0.5))`).
This also closes the first round's "not done" item: a custom family
reading `vint(N)` as a count now runs in CE predict (N held at 10),
where r9 gave 400 "NAs produced" warnings and an error
(`test-newdata-aterms.R`, 53 escaped warnings on r9).

### 3. dec() on newdata coded by newdata's own levels: fixed

At frame assembly, a term coded by its levels keeps the table of what
each value became (`resp[["aterm_codes"]]`, carried to the fit's spec
by `carry_finalized_responses()`): a factor for any term
(`as.numeric()` gives level codes), and a factor, character or logical
value of a registered term, whose coercion may build levels from its
input. `aterm_coerce()` reads newdata (a factor, character or logical
value) through that table and refuses a value the fit did not see by
name; a refit (a subset) is coded through it too. When every variable
of a term is in newdata, `aterms_for_newdata()` now reports the
evaluation's own error ("Addition term dec(side) on newdata: ...")
instead of "could not be evaluated". The mechanism is generic: the core
test uses `vint()` of a factor, with no registered term.

Covered in eam's `test-newdata-dec.R` (each block seen to fail on r9):
relevel(), character, a one-level factor, one character value, a value
the fit never saw (refused), and a one-row logical `dec(lg)`.

### 4. The refusal's remedy: fixed

It ends "; supply the variable". All five assertions match the prefix
"could not be evaluated on newdata" and pass.

### 5. Age: corrected

0.61.0 (bac1e5f6), not 0.60.0; see "How old the defect is".

### 6. NEWS marks

The strict `thres(gr = )` of the CE grid is under Breaking changes in
NEWS.md. frmtmb.eam's DESCRIPTION floor (`frmtmb (>= 0.70.0)`) must
rise at consolidation to the core release that carries this lane: its
`test-newdata-dec.R` fails on frmtmb 0.71.0 (below). eam's NEWS says
so.

### 8. zero_inflated_binomial at the zi boundary: not the same class

`dev/ndaterm-r1-probe.R`, last block (seed 7, data with no
inflation): the zi intercept is at its boundary (link-scale standard
error 4926.7). In-sample `predict()` averages 3.043 against
`fitted()` 5.735 on both r9 and the lane. `predict(propagate_error =
FALSE)` gives 5.7354 and `simulate()` 5.7358, so the simulator reads
the right rows and parameters; the gap is the parameter-draw law,
whose draws of a link-scale zi with sd 4926 put zi near 0 or 1 half
the time each. Filed, not fixed (below).

### Tests of round 1, each seen to fail

| file | r9 | lane |
|---|---|---|
| `test-newdata-aterms.R` (core) | pass=18 fail=25 error=5 warn=53 | pass=54 |
| `frmtmb.sample/.../test-newdata-aterms-draws.R` | pass=4 fail=8 | pass=12 |
| `frmtmb.eam/.../test-newdata-dec.R` | pass=6 fail=11 error=4 | pass=23 |
| `frmtmb.latent/.../test-hmm-aterms.R` (new) | pass=0 fail=2 | pass=2 |

Logs `t-newdata-r1-<arm>.txt`, `t-sample-r1-<arm>.txt`,
`t-eam-dec-r1-<arm>.txt`, `t-latent-hmm-<arm>.txt`. Every new block
fails on r9 for its behavioral reason. The sample's theta block passes
on r9 under the first spelling (a mean comparison, 6 se = 1.06 against
a gap of 1.1); it was rewritten as the share of draws in the second
component against `1 - theta1` (0.125 off, tolerance 0.097 on r9).
Its fix is the sample's own: with the lane's core alone over r9's
frmtmb.sample (arm `coreonly` of `dev/ndaterm-runfile.R`, log
`t-sample-r1-coreonly.txt`), that block still fails (0.125 >= 0.097).
The same arm shows the guard of the first round at work: r9's
frmtmb.sample on the lane's core refuses the binomial newdata
(`aterm_trials()`) where it used to draw 0/1 counts.

### Round 1 suites

All eight packages, one file per process, lane library first, 8 at a
time (`dev/ndaterm-par.sh all2`, `dev/ndaterm-summ.sh all2`; logs in
`dev/ndaterm-suite-all2/`). Every `lib:` line names the lane's core.

```
ran 386 of 386
frmtmb          files=237 pass=17494 fail=0 error=0 skip=185 warn=0
frmtmb.sample   files=54  pass=2709  fail=0 error=0 skip=4   warn=0
frmtmb.eam      files=31  pass=1827  fail=0 error=1 skip=3   warn=0
frmtmb.coupling files=11  pass=542   fail=0 error=0 skip=5   warn=0
frmtmb.latent   files=11  pass=362   fail=0 error=0 skip=2   warn=0
frmtmb.learn    files=15  pass=430   fail=0 error=0 skip=13  warn=0
frmtmb.ode      files=11  pass=549   fail=0 error=0 skip=1   warn=0
frmtmb.spline   files=16  pass=594   fail=0 error=0 skip=1   warn=0
```

The one error was frmtmb.eam's `test-ndt-bound.R:376`, "a group the fit
never saw cannot fail open": `ndt_group()` is a registered term, so an
unseen label `s9` on newdata is now refused by core's table lookup,
"holds value(s) 's9' the fit did not see", before eam's own lookup
says "was not fitted to". The refusal is the right one (it is finding
3's rule, for every registered term); the assertion was changed to
core's message, and core's message was changed to list the values the
fit saw instead of their codes (ndt_group codes are label hashes such
as 253570877197). eam rerun one file per process after that change
(`dev/ndaterm-suite-eam3/`): ran 31 of 31, pass=1831 fail=0 error=0
skip=3. `ndt_time()` keeps eam's own message.

Gated tier (`FRMTMB_BRMS_FIT_TESTS=true`, the same 46 files,
`dev/ndaterm-suite-gated2/`) on the final build: ran 46 of 46,
pass=3732 fail=0 error=0 skip=14 (13 drmTMB-gated, 1 fuzz-gated, as in
round 0). No brms-suite verdict goes stale.

OpenBLAS 0.3.26 with LAPACK, core one file per process, final build
(`dev/ndaterm-suite-ob026b/`): ran 237 of 237, pass=17495 fail=0
error=0 skip=185; the one extra pass against the reference BLAS is
`test-v11.R` again (23 against 22).

One session per changed package, final build
(`dev/ndaterm-out/onesession2-<pkg>.txt`):

```
ONESESSION frmtmb pass=17489 fail=9 error=0 skip=185 warn=0
ONESESSION frmtmb.sample pass=2709 fail=0 error=0 skip=4 warn=0
ONESESSION frmtmb.eam pass=1831 fail=0 error=0 skip=3 warn=0
ONESESSION frmtmb.latent pass=362 fail=0 error=0 skip=2 warn=0
```

The nine core failures are the same five `test-compat-register.R` and
`test-compat.R` blocks as in round 0, which lane-env.md names as failing
in one session on 0.71.0 (lane testiso's).

R CMD check --as-cran on the final build (`dev/ndaterm-check.sh`,
`dev/ndaterm-out/check2-<pkg>.txt`): frmtmb 1 NOTE, frmtmb.eam 1 NOTE,
frmtmb.latent 1 NOTE (each the HTML manual's V8 math rendering, the
expected one); frmtmb.sample OK.

## Not done, and defects found but not fixed

- **To file: `predict()` at a boundary parameter.** A minimal repro,
  `dev/ndaterm-r1-probe.R` last block:
  `frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial())` on
  `rbinom(n, N, plogis(0.3 + 0.5 x))`, n = 200, seed 7. `predict()`
  averages 3.04 over the rows where `fitted()` says 5.74, with no
  condition raised. Cause: `propagate_error = TRUE` draws the zi
  intercept from N(estimate, 4926.7^2) on the logit scale. The draw
  law at a boundary, older than this lane; the boundary message and
  `frmtmb_boundary_fit` are the seam a fix would read.
- `predict()` of `se(s, sigma = TRUE)` where sigma is estimated near 0
  overflows the parameter draws (sd 1.6e141; r9 Inf with the "draw was
  not finite" warning). The same draw-law class as the item above.
- frmtmb.sample's `posterior_predict()` on newdata does not add a
  shared `nu` (`with_shared_nu()`, a rescor student fit) as core's
  `predict_dpar_values()` does; not reached by this lane's cases.
- First-round items now closed: the CE `vint()` count (finding 2), the
  coercion by newdata's own levels (finding 3), the refusal's remedy
  (finding 4).

- **To file (R1-c): `mixture(cumulative(), cumulative())` draws NA.**
  `predict(propagate_error = TRUE)` gives NA in 42.75% of draws, in
  sample and on newdata, with no condition raised; the same on r9.
  `propagate_error = FALSE` and a plain `cumulative()` give none.
  Repro: the reviewer's `dev/ndaterm-rev-r1ord.R`. Older than this
  lane.
- **To file (R1-d): `mixture_mvn()` `fitted(newdata =)` labels.**
  `mixture_mvn(K = 2, D = 2)` labels its two response dimensions
  "P(Y = 1)" and "P(Y = 2)", and reports quantiles of 1 and 0 beside
  Estimates of 1.36 and -1.04; the same on r9. Repro: the reviewer's
  `dev/ndaterm-rev-r1dpar.R`. Older than this lane.

## Fix round 2

The re-check of round 1 (MERGEABLE, two text fixes) asked for:

- R1-a. frmtmb.eam's `ddm_indicator_binary()` message and roxygen no
  longer say that `conditional_effects()` holds `dec()` at its mean
  (it holds it at the rounded median since round 1). The message is
  now "wiener: the boundary indicator on newdata must be 0 (lower) or
  1 (upper) on every row, and 1 row(s) hold 0.5." The unseen-value
  refusal names the term once: "Addition term dec(side) holds
  value(s) 'left' the fit did not see, so they have no code; the fit
  saw 'down', 'up'" (`aterms_for_newdata()` adds its "on newdata:"
  lead only to a message that does not already name the term).
  Texts from `dev/ndaterm-r2-msg.R`, log `r2-msg.txt`.
- R1-b. The reporting-scale NEWS bullets carry their ages, verified
  read-only. Core: the per-dpar `type = "response"` read of
  `predict_dpar_values()` was written by 862032e3, an ancestor of the
  0.61.0 release commit bac1e5f6 and not of 0.60.0's df8af533.
  frmtmb.sample: its per-dpar `type = "response"` read predates
  afd00a43 (the parent already has it); afd00a43 (2026-09-05) gave a
  mixture's weights their softmax reporting scale ("A mixture's
  weights have a softmax response scale"), which made that read
  wrong. afd00a43 has frmtmb.sample 0.1.0 in its DESCRIPTION and is an
  ancestor of 21f947bd, the release of frmtmb 0.52.0 and frmtmb.sample
  0.2.0 (`git merge-base --is-ancestor` exits 0). So every
  frmtmb.sample release from 0.2.0.
- R1-c and R1-d filed above.

Reruns on the final build, one file per process (`dev/ndaterm-runfile.R`,
logs `dev/ndaterm-out/t-r2-*.txt`): core `test-newdata-aterms.R`
pass=54; eam `test-newdata-dec.R` pass=23, `test-ndt-bound.R` pass=79,
`test-message-uniqueness.R` pass=4, `test-conditions-census.R` pass=11;
all fail=0 error=0 warn=0.

## Bumps and floors

- frmtmb: minor. BREAKING: newdata without `trials()` or `se()` is
  refused; binomial-type families refuse a missing `trials()`;
  `required_aterms` checked after the fit; a level-coded term's unseen
  value refused; CE holds whole-valued term variables at
  `round(median())`; the strict CE `thres(gr = )`.
- frmtmb.sample: minor (BREAKING refusal of newdata without
  `trials()`/`se()`). Floor unchanged, `frmtmb (>= 0.71.0)`: it calls
  only `aterms_for_newdata()`, exported at 0.71.0, and its own link
  inverse. Measured in round 0 (`dev/ndaterm-floor.R`, pass=9); the
  round-1 additions need nothing more from core.
- frmtmb.eam: minor (BREAKING 0-or-1 refusal). Floor must move from
  `frmtmb (>= 0.70.0)` to the core release carrying this lane, at
  consolidation.
- frmtmb.latent: patch (a simulator fix, nothing refused). Floor
  unchanged, `frmtmb (>= 0.63.0)`: `hmm_rows_of()` is local. Measured
  on frmtmb 0.71.0 (the oldest build kept): the lane's frmtmb.latent
  alone over rellib-r9 passes `test-hmm-aterms.R`, pass=2 fail=0
  (`dev/ndaterm-floor-latent.R`, log `floor-latent.txt`).

## Files

Changed: `R/predict-brms.R`, `R/predict.R`, `R/autocor.R`,
`R/families.R`, `R/conditional-effects.R`, `R/frame.R`, `R/parse.R`
(roxygen of `frmtmb_register_aterm(coerce =)`),
`man/frmtmb_register_aterm.Rd`, `NEWS.md`,
`extensions/frmtmb.sample/R/methods-draws.R`,
`extensions/frmtmb.sample/NEWS.md`,
`extensions/frmtmb.eam/R/wiener-post.R`,
`extensions/frmtmb.eam/NEWS.md`,
`extensions/frmtmb.latent/R/hmm.R`, `extensions/frmtmb.latent/NEWS.md`,
`extensions/frmtmb.eam/tests/testthat/test-ndt-bound.R` (one assertion,
round 1 suites).

Tests: `tests/testthat/test-newdata-aterms.R`,
`extensions/frmtmb.sample/tests/testthat/test-newdata-aterms-draws.R`,
`extensions/frmtmb.eam/tests/testthat/test-newdata-dec.R`,
`extensions/frmtmb.latent/tests/testthat/test-hmm-aterms.R`.

Lane scripts: `dev/ndaterm-probe.R`, `dev/ndaterm-ce.R`,
`dev/ndaterm-brms.R`, `dev/ndaterm-runfile.R`, `dev/ndaterm-par.sh`,
`dev/ndaterm-summ.sh`, `dev/ndaterm-onesession.R`,
`dev/ndaterm-install.R`, `dev/ndaterm-check.sh`, `dev/ndaterm-floor.R`,
`dev/ndaterm-r1-probe.R`, `dev/ndaterm-r1-mixgap.R`,
`dev/ndaterm-r1-theta.R`, `dev/ndaterm-r1-theta-draws.R`,
`dev/ndaterm-floor-latent.R`, `dev/ndaterm-r2-msg.R`. The
reviewer's `dev/ndaterm-rev-*` were read, not edited.

Shared with other lanes: `R/predict.R`, `R/families.R` and
`R/predict-brms.R` (lane catmulti), `R/frame.R`, `R/parse.R`,
`extensions/frmtmb.sample/R/methods-draws.R` (lane sampleprior),
`NEWS.md` and the three extension NEWS files.

Process notes: `dev/ndaterm-par.sh` was edited (one line) while the
first round's run was live, so that run's own RESULTS.txt was garbled;
its summary was rebuilt from the per-file logs by
`dev/ndaterm-summ.sh`. Scratch libraries live under the session
scratchpad and are deleted after use; reinstall to rerun
`dev/ndaterm-floor.R` or the `coreonly` arm.

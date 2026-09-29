# Review of lane wt-thresrefit

Round of 2026-09-28, reviewed 2026-09-29. Worktree under review
`C:\Users\adf44\source\r\frmtmb-wt-thresrefit`, branch `wt-thresrefit`.
Arms:

- `lane`: the worker's build, `C:/Users/adf44/source/r/wt-thresrefit-lib`.
- `base`: the shared reference build of 0.64.0,
  `C:/Users/adf44/source/r/rellib-r3`.

Every reviewer script is `dev/thresrefit-rev-*.R` with its log beside it.
Nothing in the worktree outside those files was changed, and no git
operation was run.

## The build under review is the worktree

`dev/thresrefit-rev-00-env.R`. The installed `thres_pin_of_fit()` and
`thres_pin_apply()` deparse to the worktree's `R/thres.R` bodies,
`assemble_frame()` carries the `thres_pin` formal,
`influence.frmtmb_fit()` passes it, and both `draw_prior_entry()` and
`prior_entry_label()` carry their new branches. No reinstall was needed.
`base` has none of these symbols, confirmed in the same script.

## Claim 1: no in-package refit recounts thresholds. HOLDS.

`dev/thresrefit-rev-01-boot.R`, log `-after-log.txt`. Data seed 202,
n = 40, four categories, `table(y) = 9/16/11/4`. `simulate(nsim = 60,
re_formula = NA)` at seed 11 gives exactly one replicate (31) whose
response stops at category 3, on all four families. For each family, the
`refit()` the bootstrap performs against a `thres(3)`-pinned `frm()` on
the replicate's own data:

| family | refit n_tau | pinned n_tau | logLik relative difference | coef names identical |
|---|---|---|---|---|
| cumulative | 3 | 3 | 1.835e-11 | TRUE |
| sratio | 3 | 3 | 2.351e-11 | TRUE |
| cratio | 3 | 3 | 2.351e-11 | TRUE |
| acat | 3 | 3 | 8.980e-12 | TRUE |

Grouped, data seed 503, n = 96, three levels, `nthres = 3/2/2`,
`nraw = 7`; at simulate seed 13, 10 to 11 of 40 replicates lose a level's
top category. First such replicate, `refit()` against
`thres(k, gr = g)`:

| family | refit n_tau | pinned n_tau | logLik relative difference |
|---|---|---|---|
| cumulative | 7 | 7 | 6.086e-11 |
| sratio | 7 | 7 | 3.514e-11 |
| cratio | 7 | 7 | 3.514e-11 |
| acat | 7 | 7 | 2.642e-11 |

`frm_bootstrap(FUN = c(fixef, tau_raw))` returns 60 x 4 (ungrouped) and
40 x 8 (grouped) with 0 NA cells, on `re_formula = NA` and on
`re_formula = NULL` alike, all four families.

The path a USER writes, `simulate()` then `frm()` by hand, DOES fit K-1
thresholds: 2 instead of 3 ungrouped and 6 instead of 7 grouped, on every
family, while reaching the same log-likelihood to 7.4e-11 to 1.3e-10
relative. That is not this lane's defect, and the worker says so; it is
recorded here because a reader of `?frm_bootstrap` will meet it.

The remaining paths the findings table names, measured rather than read
off the code (`dev/thresrefit-rev-15-otherrefits.R`, log beside it, data
seed 202, `table(y) = 12/10/16/2`, fitted `n_tau = 3`):

| path | result |
|---|---|
| `frm_allfit()` | n_tau 3 on all four arms (nlminb, optim, bobyqa, nloptr_lbfgs) |
| `anova(refit = TRUE)` | runs, Df 4 and 5, Chisq 0.0834 |
| `confint(method = "profile")` on `tau_raw_3` | (0.4778, 1.5159), est 1.0104 |
| `conditional_effects(band = "boot")` | 400 rows, 0 NA in `lower__`/`upper__` |
| autoscale pre-fit | n_tau 3 |

The structural reason is checkable in one grep: `assemble_frame()` is
called from `frm()`, `influence()`, `par-template.R`, `priors.R` and
`frm_simulate()` and from nowhere else, and `refit()` replaces
`frame[["y"]][[1L]]` in place. Verdict: HOLDS.

## Claim 2: an interior loss is "left unpinned". FALSIFIED as stated.

`dev/thresrefit-rev-02-interior.R`, logs `-after-log.txt` and
`-before-log.txt`. Seed 901, n = 60, `table(y) = 30/1/18/11`, so exactly
one row (row 2) sits in category 2 and deleting it empties an interior
category.

The two logs are IDENTICAL line for line except the library banner, so
nothing here is a regression. What they show is that "left unpinned" does
not mean "left alone": for an ordered-factor response the unpinned refit
recounts and the row of the table is MISALIGNED, which is the same silent
wrong answer the lane set out to remove.

| response coding | NA in the deleted row | the row |
|---|---|---|
| integer | 0 | 0.282734278  0.0886563053  -21.8788664  0.384413086 |
| ordered factor | 1 | 0.282731951  0.0886511169  0.384408733  NA |
| unordered factor | frm() refuses (brms parity) | |
| character | 0 | 0.282734278  0.0886563053  -21.8788664  0.384413086 |

The `thres(3)` reference on the subset is
`0.282734278  0.0886563053  -21.4148732  0.384413086`. The ordered
factor's `tau_raw_2` column therefore holds `tau_raw_3`'s value
(0.384408733 against 0.384413086) and `tau_raw_3` is NA.
`cooks.distance()` at that row is NA. Grouped (seed 902, `nthres = 3/2/2`,
the single interior row in level b) is worse in the same way: 3 NA cells
and every coefficient from position 3 onwards one slot early,
`0.0991004049 1.09910277 -0.168122769 1.2476021 1.47989557 NA NA NA`
against the reference
`0.0991082092 1.09908382 -23.1346723 -0.168138177 1.24759251 -21.9741049 1.47990302 -20.8608018`.

Integer and character codings are unaffected, because they carry no
`y_levels`, so the pin applies and `max(y)` is unchanged.

What brms does, measured rather than assumed
(`dev/thresrefit-rev-04-brms.R`, brms 2.23.0, `make_standata()` so no
Stan is compiled):

| construction | brms |
|---|---|
| integer 1..3, cumulative | `nthres = 2` |
| ordered factor, 4 declared levels, 3 used, cumulative / sratio / acat | `nthres = 2` |
| integer 1,3,4 (interior absent), cumulative / sratio | `nthres = 3` |
| ordered factor, 4 levels, interior level 2 unused, cumulative | `nthres = 2` |
| `thres(3)` pinned, either coding | `nthres = 3` |
| `thres(gr = g)` with one level topping out lower, integer | `nthres = 2/2` |
| `thres(gr = g)` with an ordered factor of 4 declared levels | ERROR, "Number of thresholds is smaller than required by the response" |

So the premise in the assignment is wrong: brms does NOT error when the
response leaves a category unused. It recounts silently, exactly as
frmtmb's base build does, in the interior case as well as the top one.

Verdict: leaving the MODEL unpinned there is defensible and matches brms.
Leaving the influence ROW misaligned is a hole, and the new documentation
denies it (see claim 7). The cheap correct behavior is to refuse the unit
when the remaining labels are not an initial segment, so the row is all
NA instead of shifted. The worker's own test asserts the unpinned branch
at `thres_pin_apply()` but never measures what `influence()` then
returns, which is why the misalignment survived.

## Claim 3: identical rows across codings. HOLDS, with one exception.

`dev/thresrefit-rev-03-coding.R`, logs `-after-log.txt` and
`-before-log.txt`. Worker's seeds 501 (n = 50, `table(y) = 20/11/18/1`)
and 502 (n = 96, per-level max 4/3/3, one row in level a's category 4).

| model | integer | ordered factor | unordered factor | character |
|---|---|---|---|---|
| `bf(y ~ x)` | 0 NA, row filled | `identical()` TRUE to integer | frm() refuses | `identical()` TRUE to integer |
| `bf(y \| thres(gr = g) ~ x)` | 0 NA | `identical()` TRUE | frm() refuses | `identical()` TRUE |
| `bf(y \| thres(3) ~ x)` | 0 NA, row filled | 4 of 4 cells NA in the deleted row | frm() refuses | `identical()` TRUE |

`identical()` TRUE is bitwise, on the whole table, not a tolerance. The
unordered-factor refusal is pre-existing and is brms parity (brms 2.23.0
refuses it too).

The third row is the finding. `thres_pin_apply()` returns at
`if (is.null(pin) || !is.null(av[["thres"]])) return(out)`, so a count the
user wrote skips the LABEL restoration as well as the count injection.
The result is that the one case where the user has already pinned the
count by hand is the one case `influence()` cannot refit: the ordered
factor's subset has 3 levels, `thres(3)` asks for 4 categories, and
`frame.R`'s ordered-factor refusal fires inside the per-unit `tryCatch`.
Base and lane are identical here (4 NA both), so it is not a regression,
and moving the label work above that early return would close it.

## Claim 4: agreement to 2.4e-5 relative at worst. HOLDS, as noise.

`dev/thresrefit-rev-05-agree.R`, log beside it. The worker's seed 501
reproduces to the digit: `x` 0.000627432768754419 against
0.000627417852989218, relative 2.38e-05, and 3.12e-07 of one column
standard deviation. Four fresh seeds (1501..1504) on cumulative and
sratio give a worst identified-coefficient relative difference of
9.4e-12 to 3.6e-07, so 2.4e-5 is the worst of six cases and is inflated
by the small denominator, as the findings say.

The question the assignment asked, settled on log-likelihoods. The
leave-one-out refit `influence()` performs was reproduced through the same
internal path (`assemble_frame(thres_pin = )` then `fit_assembled()`), its
coefficients are BITWISE identical to the influence row in all ten cases,
and its log-likelihood against the `thres(3)` reference:

| seed | family | logLik loo-refit | logLik reference | relative |
|---|---|---|---|---|
| 501 | cumulative | -52.38100757488011 | -52.38100757368454 | 2.28e-11 |
| 501 | sratio | -52.34085133678953 | -52.34085133372233 | 5.86e-11 |
| 1501 | cumulative | -53.07071216707108 | -53.07071216701262 | 1.10e-12 |
| 1501 | sratio | -53.23155609707752 | -53.23155609546244 | 3.03e-11 |
| 1502 | cumulative | -52.22040333172659 | -52.22040333387736 | 4.12e-11 |
| 1502 | sratio | -52.08626202779279 | -52.08626202981490 | 3.88e-11 |
| 1503 | cumulative | -51.13064696036680 | -51.13064695834981 | 3.95e-11 |
| 1503 | sratio | -51.46531501426501 | -51.46531501518086 | 1.78e-11 |
| 1504 | cumulative | -49.28730542568737 | -49.28730542631965 | 1.28e-11 |
| 1504 | sratio | -49.10842818985080 | -49.10842818852837 | 2.69e-11 |

Same parameter count, same coefficient names, same maximized
log-likelihood to 1e-11 or better, and the only coefficient that moves is
`tau_raw_3`, whose deviation is 0.00012 to 0.026 relative. That is a flat
ridge, not a different model. Verdict: HOLDS.

An instrument of my own failed and is recorded rather than deleted: I
also evaluated one fit's objective at the other's parameter vector using
`fit$obj$par`, which in RTMB is the START vector and not the optimum, so
the "gap" column that script prints is meaningless. The log-likelihood
comparison above is the evidence.

## Claim 5: the recycling is unreachable. HOLDS, census incomplete.

`dev/thresrefit-rev-06-priorb.R`, logs `-before-log.txt` and
`-after-log.txt`. The public path on the base build, for all four ordinal
families with `set_prior("normal(0, 2)", class = "Intercept")` and K = 4,
seed 702:

    base : simpleError  values must be length 1,
           but FUN(X[[1]]) result is length 3
    lane : frmtmb_error A prior on tau_raw[1:3] is a density on 3
           parameters at once, ... (cumulative keeps its older refusal,
           "would not produce an ordered one")

The complement, which the guard must not close on: K = 2 (one threshold)
on cratio, acat and sratio draws on both arms and the `pars` report is
bit-identical, `tau_raw_1 = 3.2534387, 0.9516278` and
`b_x = 0.9604788, -0.8118185`. A multivariate model, the `rescor` route,
is refused earlier by name ("frm_simulate() supports univariate models").

The fit path and the prior reports are untouched. With the same
`class = "Intercept"` plus `class = "b"` priors, base and lane print the
same 17-digit log-likelihood and coefficients on cumulative
(-61.26023248446688) and sratio (-61.38101822691157), and
`prior_summary()` and `default_prior()` print the same rows.

Where the census is wrong. `dev/thresrefit-rev-07-multiidx.R` reads the
entries the resolver builds and finds TWO more multi-index routes, both
reachable through `frm_simulate()`:

| model | entry | base | lane |
|---|---|---|---|
| `(1 + x + z \| g)`, `set_prior("lkj(2)", class = "cor")` | `theta` idx 4,5,6 | `simpleError: 'length = 3' in coercion to 'logical(1)'` | named refusal, "A prior on theta[4:6] is a density on 3 parameters at once" |
| `ar(time = t, gr = gg, p = 2)`, `set_prior("normal(0, 0.3)", class = "ar")` | `thetaac` idx 1,2 | `simpleError: values must be length 1, but FUN(X[[4]]) result is length 2` | named refusal, "A prior on thetaac[1:2] is a density on 2 parameters at once" |

The fit path on both of those models is bit-identical between arms
(logLik -175.98786419096427 and -123.63837585947768, theta and thetaac
printed at 12 digits in the log). So there is NO regression, and the
change improves two more paths. But the findings say `rescor` is "the
other multi-index entry the resolver can build", which is wrong, and the
NEWS bullet names only the ordinal case, so a user who writes
`set_prior("lkj(2)", class = "cor")` inside `frm_simulate()` meets a new
message that NEWS does not mention.

frmtmb.sample, run with the WORKER's core underneath and frmtmb.sample
from `rellib-r3` (`dev/thresrefit-rev-13-sample.R`; arm separation
verified in `dev/thresrefit-rev-14-armcheck.R`, which prints
`wt-thresrefit-lib/frmtmb` against `rellib-r3/frmtmb`):

| file | lane | base |
|---|---|---|
| test-prior-route.R | pass=9 fail=0 err=0 skip=0 | pass=9 fail=0 err=0 skip=0 |
| test-prior-update.R | pass=7 fail=0 err=0 skip=0 | pass=7 fail=0 err=0 skip=0 |
| test-default-priors-brms.R | pass=25 fail=0 err=0 skip=0 | pass=25 fail=0 err=0 skip=0 |
| test-reparam.R | pass=292 fail=0 err=0 skip=0 | pass=292 fail=0 err=0 skip=0 |
| test-conditions.R | pass=6 fail=0 err=0 skip=0 | pass=6 fail=0 err=0 skip=0 |

No install was needed, and `grep -rn "prior_entry_label|draw_prior"` over
every extension's `R/` and `tests/` matches nothing.

## Claim 6: the recorded base counts, and the guard's absent case. HOLDS.

`dev/thresrefit-rev-runtest.R` counts from `as.data.frame(test_file())`,
not from the reporter's capped summary, and prints the block count so an
aborted file cannot pass as clean.

`test-thres-refit.R` on `base`: `pass=43 fail=15 err=2 skip=0`, exactly
the recorded figure, with the same 15 failures and 2 errors
(`dev/thresrefit-rev-testlog-base/test-thres-refit.R.log`). One error is
the weak form ("object 'thres_pin_of_fit' not found"); the worker flags
it as such. The grouped positional shift is caught directly, as recorded:
`abs(got[[unident]])/scale_id` is 0.6 against a required 3.

Thirteen other files run on both arms give identical counts:

| file | lane | base |
|---|---|---|
| test-frame.R | 22 | 22 |
| test-influence-plot.R | 21 | 21 |
| test-thres.R | 67 | 67 |
| test-ordinal.R | 109 | 109 |
| test-ordinal-fitted.R | 129 | 129 |
| test-sratio-thresholds.R | 29 | 29 |
| test-simulate-ergonomics.R | 50 | 50 |
| test-boot.R | 49 | 49 |
| test-setprior.R | 27 | 27 |
| test-unpinned-seams.R | 30 | 30 |
| test-multiple-pooling.R | 33 | 33 |
| test-autoscale.R | 46 | 46 |
| test-ce-bands.R | 167 | 167 |

all with fail=0 err=0 skip=0.

The ABSENT case, by construction (`dev/thresrefit-rev-09-absent.R`, seeds
601 to 603). `thres_pin_of_fit()` returns NULL on gaussian, gaussian with
`(1 | g)`, bernoulli, poisson and categorical, and returns
`grouped=FALSE nthres=3`, `grouped=FALSE nthres=3` and
`grouped=TRUE nthres=2/2` on the three ordinal models that have nothing
to lose. Diffing the whole log between arms leaves ONLY the added
`thres_pin_of_fit ->` lines and the library banner: every influence
table's dimension, NA count, digest of the rounded matrix, 17-digit sum
and first four Cook's distances are unchanged. So the pin is inert where
it should be.

No absolute numeric tolerance is written in `test-thres-refit.R`.
Coefficient agreement is `|difference| / sd(column of the influence
table)`, the log-likelihood one is `|difference| / |log-likelihood|`
against `sqrt(.Machine$double.eps)`, and the unidentified threshold is
compared to the largest identified threshold of the same fit. The three
bare numbers that remain (1e-3 and 3) are cutoffs on dimensionless
ratios whose denominators the run measures, which is what the rule asks
for.

The `accepts_aterms` half of the gate in `thres_pin_of_fit()` cannot be
false for any family in the package: all four ordinal families declare
`c("weights", "thres")`. It is defensive for extension families, which
is right, but the test exercises only the `type != "ordinal"` branch.

## Claim 7: the docs. RENDER CLEAN, TWO STATEMENTS FALSE.

`dev/thresrefit-rev-08-rd.R` renders `man/frm.Rd`,
`man/frm_bootstrap.Rd` and `man/influence.frmtmb_fit.Rd` with
`tools::Rd2txt`. All three render, the new "Ordinal thresholds" section
and the new `\details{}` appear, and a byte scan of every line the diff
adds finds no em dash, no en dash, no emoji, no British spelling and no
unescaped `%` in any Rd. NEWS additions are all within 80 columns.

Two rendered statements are false in cases I measured above.

- `?influence.frmtmb_fit`: "Each deletion refits the same model, so an
  ordinal response keeps the threshold count of the full-data fit and the
  per-level counts of `thres(gr = )`." Not so when the deletion empties
  an INTERIOR category of a factor response (claim 2), when the user
  wrote `thres(K)` on a factor response (claim 3), or when
  `influence(groups = )` deletes a whole `thres(gr = )` level (which the
  findings record and the docs do not).
- `?frm`: the same sentence, plus "`influence()` and `cooks.distance()`
  rebuild the design from a subset of the data and carry the count over",
  with the same exceptions.

One nit: `?frm_bootstrap` states the replicate's log-likelihood matches a
`thres(K)`-pinned fit "to about 2e-11 relative". Measured over four
families, grouped and ungrouped, the range is 8.98e-12 to 6.09e-11. A
number in user documentation that a different optimizer or box will move.

The NEWS entry is true as far as it goes. It does not mention the
`class = "cor"` and `class = "ar"` routes the same change alters (claim
5), and it does not mention the `influence(data = )` change of failure
mode (finding 1 below).

## Claim 8: the 40 reachable files, and the whole suite. HOLDS.

`dev/thresrefit-rev-suite.sh` with `dev/thresrefit-rev-files.txt`, one
file per R process, 6 at a time, `NOT_CRAN=true`, arm `lane`, log
`dev/thresrefit-rev-suite-lane-log.txt`, per-file logs in
`dev/thresrefit-rev-testlog-lane/`.

<!-- generated by dev/thresrefit-rev-suite.sh; pasted verbatim -->
    launched 40 files  arm=lane
    files=40 pass=5322 fail=0 err=0 skip=0
    files with no REVRESULT line (aborted): none

The per-file pass counts were diffed against the worker's generated block
in `dev/thresrefit-findings.md` line by line: IDENTICAL for all 40, and
the totals agree at 5322.

The whole suite, all 180 `test-*.R`, one per process, 6 at a time, arm
`lane`, log `dev/thresrefit-rev-suite-all-log.txt`, per-file logs in
`dev/thresrefit-rev-testlog-all/`:

<!-- generated by dev/thresrefit-rev-suite.sh; pasted verbatim -->
    launched 180 files  arm=lane
    files=180 pass=11949 fail=0 err=0 skip=152
    files with no REVRESULT line (aborted): none

Both figures are the worker's exactly. The 152 skips are the gated tiers,
which this run does not enable, and the per-file breakdown matches the
one in the findings: test-brms-methods.R 46, test-brms-likelihood.R 37,
test-drmtmb-agreement.R 13, test-brms-priors.R 12, the thirteen
test-bcm-*.R files 39 between them, test-brms-agreement.R 2,
test-rl-example.R 2, test-fuzz.R 1.

## Findings, in severity order

### 1. BLOCKING: the new documentation is false in three cases

One of those cases returns a misaligned row.
`?influence.frmtmb_fit` now reads "Each deletion refits the same model,
so an ordinal response keeps the threshold count of the full-data fit and
the per-level counts of `thres(gr = )`", and `?frm` says the same. Three
measured cases contradict it, and one of them is a silent wrong number
rather than an NA.

The wrong number, from claim 2 (`dev/thresrefit-rev-02-interior.R`, seed
901, `table(y) = 30/1/18/11`, ordered-factor response). The deleted row
is

    0.282731951  0.0886511169  0.384408733  NA

and the `thres(3)` reference for the same subset is

    0.282734278  0.0886563053  -21.4148732  0.384413086

so the `tau_raw_2` column of that row reports `tau_raw_3`. Grouped (seed
902) the shift runs over five coefficients. The only visible sign is a
trailing NA, which is exactly the signature the lane set out to remove,
and the documentation now says it cannot happen.

The other two cases the sentence covers and gets wrong: a hand-written
`thres(K)` on a factor response (finding 3, whole row NA), and
`influence(groups = )` deleting a whole `thres(gr = )` level (7 of 24
cells NA in `dev/thresrefit-rev-11-edges.R` section 3, which the findings
record and the documentation does not).

Neither the misalignment nor the NA is a regression: base and lane are
identical line for line, and brms 2.23.0 recounts the same way
(`dev/thresrefit-rev-04-brms.R`), so leaving the MODEL unpinned is
defensible. The documentation is not. Fix by refusing the unit when the
remaining labels are not an initial segment, so the row is honestly all
NA, or by naming the three exceptions in both help pages. One sentence
each would do.

### 2. BLOCKING: `influence(data = )` gives an all-NA table in silence

0.64.0 raised an error on the same two inputs. Two constructions, both
measured on both arms.

`dev/thresrefit-rev-11-edges.R` section 1, seed 1101, n = 50. The fit has
three categories and `n_tau = 2`. Two rows of the `data` argument are
moved into category 4, which the fit never saw.

    base : subscriptOutOfBoundsError, "subscript out of bounds"
    lane : returns a 50 x 3 frmtmb_influence, 150 of 150 cells NA,
           50 of 50 rows entirely NA, no warning and no message

`dev/thresrefit-rev-12-newlevel.R`, seed 1103, k = 96,
`thres(gr = g)` over levels a/b/c, `nthres = 3/2/2`. The `data` argument
carries a fourth level "d".

    base : subscriptOutOfBoundsError, "subscript out of bounds"
    lane : returns a 96 x 8 frmtmb_influence, 768 of 768 cells NA,
           no warning and no message

The mechanism: with the pin, `thres_finalizer()`'s "Number of thresholds
is smaller than required by the response" and `thres_pin_apply()`'s own
"cannot add a level" both fire INSIDE the per-unit `tryCatch` in
`influence.frmtmb_fit()`, which returns NULL and leaves the row NA.
Without the pin, the recounted refit produced coefficient names that are
not columns of the table and `fe[i, names(fe_i)] <- fe_i` threw. So a
loud error became a silent empty answer, on a documented argument, in a
lane whose purpose is to remove silent wrong answers from ordinal
refits.

`influence()` has no failure counter and warns about nothing, so a
partially failed table is already silent; this change makes a totally
failed table silent too. `dev/thresrefit-findings.md` records the first
construction as "that unit is skipped", which understates it by a factor
of the sample size, and does not record the second at all.

In mitigation: `?influence.frmtmb_fit` describes `data` as "the original
model data", so data whose response reaches a category the fit never saw
is a misuse, and the result is an empty answer rather than a wrong one.
It is still a loud failure turned quiet, in a lane about quiet failures.

Remedy, small: count the units whose refit returned NULL and
`frm_warning()` with that count, or let the pin's refusal propagate for
`data = ` rather than swallowing it. Either one also improves the
pre-existing partial-NA case.

### 3. A hand-written `thres(K)` plus an ordered factor is still skipped

Numbers under claim 3: 4 of 4 cells NA in the deleted row, on both arms,
while the integer and character codings give the full row. Caused by
`thres_pin_apply()` returning before the label restoration when
`av[["thres"]]` is already set. One line of reordering closes it. Claim 3
as written ("identical rows across codings") is true only for a formula
with no hand-written `thres()`.

### 4. The reachability census for defect B misses two routes.

Numbers under claim 5. `class = "cor"` on a block with more than one
correlation and `class = "ar"` (or `ma`, `cortime`) above order 1 build
multi-index entries that `frm_simulate()` reaches. No regression, an
improvement in both, but the findings call `rescor` "the other" such
entry and NEWS names only the ordinal case.

### Nits

- `?frm_bootstrap` pins "about 2e-11 relative"; measured 8.98e-12 to
  6.09e-11 over four families, grouped and ungrouped.
- `prior_entry_label()`'s new branch fires on `length(e$idx) != 1L`,
  which includes `integer(0)`, and the label then reads `comp[NA:NA]`. No
  resolver I found builds one, but the condition is wider than the
  comment beside it describes.
- `test-thres-refit.R` opens every block with `skip_on_cran()`, so the
  pin does not run on CRAN. Consistent with the rest of the suite; noted
  only.
- The `accepts_aterms` half of the `thres_pin_of_fit()` gate has no
  construction that makes it false, because all four ordinal families
  declare `thres`.

## Test files run, with counts

Arm `lane`, `NOT_CRAN=true`, one file per R process:
`dev/thresrefit-rev-testlog-lane/` for the 40, and
`dev/thresrefit-rev-testlog-all/` for all 180. Arm `base`:
`dev/thresrefit-rev-testlog-base/`. frmtmb.sample:
`dev/thresrefit-rev-samplelog/`. The tables above carry the per-file
numbers.

## Claim-by-claim verdicts

| claim | verdict |
|---|---|
| 1. No in-package refit recounts | HOLDS |
| 2. Interior loss is "left unpinned" | FALSIFIED as stated: the row is misaligned, not merely unpinned. Not a regression, and brms recounts the same way, so the model choice is defensible; the documentation is not |
| 3. Identical influence rows across codings | HOLDS for integer, ordered factor and character with no hand-written `thres()`; FALSIFIED for `bf(y \| thres(K) ~ x)` with an ordered factor, where the row is all NA |
| 4. Agreement is optimizer noise | HOLDS, log-likelihoods agree to 1e-11 relative or better on ten constructions |
| 5. Defect B unreachable from outside | HOLDS; the census of multi-index routes is incomplete by two |
| 6. Base counts and the guard | HOLDS, 43/15/2 reproduced exactly, guard inert on eight constructions |
| 7. Docs | Render clean and pass the style scan; two new statements are false |
| 8. Suite counts | HOLDS, 40 files 5322 pass and 180 files 11949 pass, both identical to the record |

## Verdict

NOT MERGEABLE.

**Superseded on 2026-09-29 by the "Re-check after punch round 1" section
at the end of this file, which is MERGEABLE. Both blocking items below
were fixed and the fixes were verified by construction. The re-check also
corrects one statement made above: the `groups = ` deletion of a
`thres(gr = )` level was recorded as "7 of 24 cells NA, a limitation", and
its rows were in fact filled and misaligned, so it was a third instance of
this lane's own defect.**

Blocking, in order:

1. The new `?influence.frmtmb_fit` and `?frm` paragraphs say
   unconditionally that a deletion keeps the fitted threshold count.
   Three measured cases contradict them, and in one the deleted row
   reports `tau_raw_3` in the `tau_raw_2` column with no sign but a
   trailing NA. Refuse the unit, or name the exceptions in both pages.
2. `influence(data = )` returns a table of nothing but NA, silently,
   where 0.64.0 raised an error: 50 of 50 rows on a response reaching a
   higher category, 96 of 96 on a new `thres(gr = )` level, no warning in
   either. Introduced by this change. Fix with a failure count and a
   warning, or by letting the pin's refusal propagate.

Nits and record corrections:

3. `bf(y | thres(K) ~ x)` with an ordered-factor response is the one
   coding the pin skips, because a user-written count short-circuits the
   label restoration as well. One line of reordering.
4. The defect-B census misses the `class = "cor"` and `class = "ar"`
   routes, which the same change also turns from an internal error into a
   named refusal; NEWS should say so.

Nits: the "about 2e-11" figure in `?frm_bootstrap`;
`prior_entry_label()`'s `!= 1L` admitting `integer(0)`; the
`accepts_aterms` gate with no false construction.

Everything else measured in this review holds, including the whole 180
file suite, the 40 reachable files per file, the base-arm failure counts,
the bit-identical influence tables on five non-ordinal families, and the
frmtmb.sample prior suites on both cores.

# Re-check after punch round 1, 2026-09-29

Same harness, same two arms. `lane` was reinstalled by the worker;
`dev/thresrefit-rev-00-env.R` confirms the installed `thres_pin_recode()`,
`thres_pin_apply()` and `thres_pin_of_fit()` deparse to the worktree's
`R/thres.R`, and that `influence()` carries the failure counter. New
scripts `dev/thresrefit-rev-20-recode.R` through `-rev-24-docclaims.R`,
each with a `-before-log.txt` and an `-after-log.txt`. Only what changed
is below.

## Item 1: the reviewer's own interior and grouped constructions

`dev/thresrefit-rev-20-recode.R`. My rev-02 generators, unchanged.

| construction | arm | integer | ordered factor | character | `identical()` to integer |
|---|---|---|---|---|---|
| seed 901, ungrouped, `table(y) = 30/1/18/11`, one interior row | base | 0 NA | 1 NA | 0 NA | ordered FALSE, character TRUE |
| | lane | 0 NA | 0 NA | 0 NA | both TRUE |
| seed 902, grouped, `table(y) = 74/1/17/4`, `nthres = 3/2/2` | base | 0 NA | 3 NA | 0 NA | ordered FALSE, character TRUE |
| | lane | 0 NA | 0 NA | 0 NA | both TRUE |
| seed 904, grouped, `table(y) = 61/1/30/4` | base | 0 NA | 3 NA | 0 NA | ordered FALSE, character TRUE |
| | lane | 0 NA | 0 NA | 0 NA | both TRUE |

`identical()` is over the WHOLE table, bitwise, not a tolerance. On `lane`
every construction reports 0 warnings and `cooks.distance()` with no NA.
Against the hand-pinned reference on the same subset, the deleted row's
IDENTIFIED coefficients agree to 6.87e-13 to 1.44e-08 of one column
standard deviation, and the columns that disagree are exactly the
log-increments of the emptied category, which bound a category nobody
chose: seed 901 `tau_raw_2` -21.88 against -21.41, seed 904 `tau_raw_2`
-21.59 against -28.77, `tau_raw_5` -20.41 against -19.77, `tau_raw_7`
-21.23 against -28.22. Those are minus-infinity parameters stopped at
different places on a flat ridge, and they sit in their OWN columns on
both sides, which is the whole point. HOLDS.

My seed 902 does reproduce the misalignment on `base` (3 NA), unlike the
worker's own seed 902; my generator empties category 2 globally, theirs
does not. Not a contradiction, and both of mine are fixed.

## Item 2: label sets that are not "1".."K"

`dev/thresrefit-rev-21-labels.R`, seed 901 with one row in the interior
category 2 AND one in the top category 4, so both emptyings happen in one
table. Every coding gives the same FULL-fit log-likelihood,
-48.83081807764106, and the same four coefficients at 14 digits.

| coding | `y_levels` in the fit | lane `identical()` to integer | base |
|---|---|---|---|
| ordered factor "1".."4" | 1,2,3,4 | TRUE | FALSE |
| lo < mid < hi < top (alphabetical order would be hi, lo, mid, top) | lo,mid,hi,top | TRUE | FALSE |
| d < c < b < a, declared in descending alphabetical order | d,c,b,a | TRUE | FALSE |
| `levels = 1:5` with level 5 taken by no row even in the FULL fit | 1,2,3,4 | TRUE | FALSE |
| character "1".."4" | NULL | TRUE | TRUE |
| grouped `thres(gr = g)` with lo/mid/hi/top | lo,mid,hi,top | TRUE | FALSE (3 NA) |
| `bf(y \| thres(3) ~ x)` with lo/mid/hi/top | lo,mid,hi,top | TRUE, 0 NA | FALSE, 8 NA |

On `lane` all of these are 0 NA and 0 warnings. The declared-but-unused
fifth level is dropped in the FULL fit too, so the pin's label set is the
fit's OBSERVED labels and is self-consistent; that matches brms, which
recounts from the data (measured last round). The recode uses the fitted
labels rather than any sort order, which the descending case pins. HOLDS.

A character vector of NON-numeric labels (lo, mid, ...) is refused by
`frm()` on BOTH arms with "missing value where TRUE/FALSE needed". That is
pre-existing and nothing in this lane touches it, but see the nits.

## Item 3: the recode is inert where it should be

`dev/thresrefit-rev-09-absent.R` rerun on the new build and diffed against
the stored `base` log. The diff is again ONLY the library banner and the
added `thres_pin_of_fit ->` lines: every influence table's dimension, NA
count, digest of the rounded matrix, 17-digit sum and first four Cook's
distances are unchanged for gaussian, gaussian with `(1 | g)`, bernoulli,
poisson and categorical (pin `NULL` on all five) and for the three ordinal
integer models with every category present (pin present and inert).

Added, because rev-09 had no such case
(`dev/thresrefit-rev-23-preflight.R`): an ordered-factor ordinal response
at seed 2301 with `table(y) = 26/15/20/19`, so no deletion empties
anything. 0 NA, 0 warnings, and `identical()` TRUE across integer, ordered
and character, sum -35.366433494278262 on all three. So the recode's
`identical(levels, lv_fit)` early return is exercised and inert. HOLDS.

## Item 4: the four failure cases

`dev/thresrefit-rev-22-failcount.R`, both arms.

| case | base | lane |
|---|---|---|
| `data = ` with TWO rows in an unseen category (seed 1101, n = 50) | `subscriptOutOfBoundsError` | ERROR "influence(): all 50 deletion refits failed ... The first was observation 1: Number of thresholds is smaller than required by the response ... reaches category 4" |
| `data = ` with ONE row in an unseen category, row 7 | `subscriptOutOfBoundsError` | WARNING "49 of 50 deletion refits failed and their rows are NA. The first was observation 1: ...", table 50x3 returned, row 7 filled, 49 rows all NA |
| the same as an ordered factor, so the RECODE refuses rather than the count guard | `subscriptOutOfBoundsError` | WARNING with "the response holds category '4' that the fitted model never saw, whose categories are '1', '2', '3'" and row 7 filled |
| `data = ` with a fourth `thres(gr = )` level "d" (seed 1103, k = 96) | `subscriptOutOfBoundsError` | ERROR "all 96 deletion refits failed ... 'd' ... cannot add a level" |
| `groups = 'g'` where g IS the threshold factor | 3 x 8 table, 7 of 24 cells NA, SILENT, and every row misaligned (row a: 0.1466 -1.0725 0.8713 -0.7922 0.2610 NA NA NA) | ERROR "all 3 deletion refits failed ... The first was 'g' level a: ... have no rows left, so the 3 threshold(s) the fitted model gives them cannot be estimated ... cannot drop a level" |

Every message names the unit (observation 1, 'g' level a) and the reason.
The 49-of-50 count is right by construction: only the deletion of row 7
removes the offending row. HOLDS.

This also upgrades my earlier characterization of the `groups = ` case. I
had recorded it as "7 of 24 cells NA, a limitation"; the base rows are
FILLED and WRONG, not merely short, so that case was a third instance of
the lane's own defect and is now refused.

## Item 5: the control

Same script. `influence(groups = "id")` on
`bf(y | thres(gr = g) ~ x + (1 | id))`, seed 1201, n = 120, where `id` is
NOT the threshold grouping factor: 8 x 8, 0 NA cells, 0 warnings,
`cooks.distance()` with no NA, IDENTICAL on both arms. The smallest number
of rows any `g` level keeps across the eight `id` deletions is 35.
Observation-wise on the same model: 120 x 8, 0 NA, 0 warnings, both arms.
Seed 1101 with `data` untouched: 50 x 3, 0 NA, 0 warnings, both arms.
`influence(groups = "sub")` on a `thres(gr = g)` model and a bivariate
model with two ordinal responses
(`dev/thresrefit-rev-11b-after-log.txt`): 12 x 8 and 60 x 7, 0 NA. The new
refusals do not fail closed. HOLDS.

`dev/thresrefit-rev-01-boot.R` rerun on the new build diffs EMPTY against
its punch-round-0 output, so nothing in the bootstrap and `refit()` arm
moved.

## Item 6: test-data2.R, and whether the error should come earlier

| file | base | lane |
|---|---|---|
| test-data2.R | pass=29 fail=2 err=0 skip=0 | pass=31 fail=0 err=0 skip=0 |

Both are the worker's figures exactly, and the 2 base failures are the two
new `expect_error()` calls, seen to fail. The block that ASSERTED the
silent all-NA table now asserts the error and its cause.

Is the new error right for a non-ordinal fit whose `data2` is lost? Yes,
and a pre-flight would not be better. Measured
(`dev/thresrefit-rev-23-preflight.R`): one `assemble_frame()` on the FULL
data refuses for three of the four failure cases and ASSEMBLES for the
fourth, `groups = 'g'`, where only the per-level subsets fail. So a
pre-flight cannot replace the counter. It would also destroy a useful
result: in the one-offending-row case the single successful unit is the
deletion of that row, which is the informative row, and the current code
hands it back with a warning where a pre-flight would refuse the call. The
saving would be small: 50 failing units cost 0.040 s against 0.230 s for
the same 50 units succeeding, because a failing unit dies at frame
assembly before any optimization. `proc.time()` ticks at 10 ms here, so
0.040 s is four ticks and the figure is coarse; the ORDER is what matters
and it is not in doubt. The message names `influence()`, the count, the
unit and the underlying condition, and says nothing ordinal-specific,
which is right because the data2 cause is not. HOLDS as delivered.

## Item 7: the reachable files

`dev/thresrefit-rev-suite.sh` with the 44 names in
`dev/thresrefit-suite.sh`, one file per R process, `NOT_CRAN=true`, 6 at a
time, arm `lane`, logs in `dev/thresrefit-rev-testlog-lane2/`.

<!-- generated by dev/thresrefit-rev-suite.sh; pasted verbatim -->
    launched 44 files  arm=lane
    files=44 pass=5618 fail=0 err=0 skip=2
    files with no REVRESULT line (aborted): none

The 2 skips are the gated brms tier inside test-brms-agreement.R, which an
ungated run does not enable. Per-file pass counts were diffed against the
worker's generated block: IDENTICAL on all 40 files that block lists,
including test-thres-refit.R at 118 and test-data2.R at 31.

`test-thres-refit.R` on `base`: **pass=68 fail=23 err=5 skip=0**, the
worker's figure exactly. Of the 5 errors, 2 are the weak form (object
'thres_pin_of_fit' not found) and 3 are behavioral
(`subscriptOutOfBoundsError` in `[<-`, and the `vapply` death). No
expectation in the file uses an absolute numeric tolerance: every
comparison is a ratio to `spread`, `scale_id` or `|logLik|`, and the
coding-invariance pins use `expect_identical()` on the whole table.

## Item 8: the three Rd topics

`dev/thresrefit-rev-08-rd.R`, log `dev/thresrefit-rev-08b-rd-log.txt`. All
three render; no em dash, no en dash, no British spelling and no unescaped
percent sign in any line the diff adds; NEWS additions within 80 columns.
The earlier "about 2e-11" figure is gone, replaced by "to optimizer
tolerance", and the hand-written `simulate()` then `frm()` behavior I
recorded is now documented in `?frm_bootstrap`.

The statement that blocked last round is now true as measured. "Each
deletion refits the same model, so an ordinal response keeps the threshold
count of the full-data fit, and the per-level counts of `thres(gr = )`,
whether the response is coded as integers, as a character vector or as an
ordered factor" is pinned by items 1, 2 and 3 above. "Deleting the last
observation in a category, top or interior, leaves that category's
threshold in the model and unidentified, which shows as a large
displacement in its column, never as a shorter row or a coefficient under
the next column's name" is pinned by the same runs.

The two new "cannot be refit" statements are true of the code, with two
imprecisions, both measured and both in the direction of promising less
than the code does. They are in the nits.

## New findings this round, all nits

1. **The `data = ` refusal is described more widely than it fires.** Both
   pages say `data = ` "holding a response category ... that the fit never
   saw describes a different model". Falsified in one construction
   (`dev/thresrefit-rev-24-docclaims.R`, seed 2401): a fit that pins
   `thres(4)` by hand on data reaching only category 3 declares five
   categories, and `influence(data = )` with one row in category 4 refits
   every unit, 60 x 5, 0 NA, 0 warnings, on both arms. The same fit with a
   row in category 6 does warn, 59 of 60 failed. So the condition is
   "outside the fitted model's threshold layout", not "a category the fit
   never saw". The sentence costs a reader a capability rather than giving
   them a wrong number, which is why it is a nit and not a blocker. The
   ordered-factor spelling of that model cannot be fitted at all
   ("thres(x = 4) asks for 5 categories, and the response is an ordered
   factor with 3 levels in the data"), which is pre-existing.
2. **"top or interior" is an incomplete enumeration.** The BOTTOM category
   behaves the same way and is also fixed, on a construction the worker's
   evidence does not have (seed 2402, `table(y) = 1/20/13/26`, one row in
   category 1). Lane: 0 NA, 0 warnings, `identical()` TRUE across integer,
   ordered and character, and the deleted row agrees with the `thres(3)`
   reference to 2.21e-05 and 5.09e-05 of a column standard deviation on
   the identified coefficients. Base: the ordered factor gives 1 NA and
   the misaligned row 0.941044305 -0.85430991 0.0631241488 NA. Say "any
   category" or name the bottom.
3. **The docs now invite a character coding that `frm()` refuses for
   non-numeric labels**, with "missing value where TRUE/FALSE needed",
   identically on both arms. Pre-existing and outside this lane, but the
   new sentence points readers at it, so it is worth a filed defect.
4. **Record.** The findings say "118 expectations in 13 blocks";
   `grep -c "^test_that("` gives 14, and my runner reports `blocks=14`.
   Separately, the generated count block still lists 40 files while
   `dev/thresrefit-suite.sh` now names 44: test-brms-agreement.R (168, 2
   skips), test-mv-gaps.R (46), test-tabular-inputs.R (9) and test-v07.R
   (25) are run but not in the block.
5. **Cosmetic.** `thres_pin_recode()` reads `p[["levels_y"]]` before it
   tests `is.null(p)`. It works because a `[[` read on `NULL` is `NULL`,
   but it reads as a bug at the point where a reader is checking the
   guard.
6. **My own stale script, recorded so nobody mis-reads its log.**
   `dev/thresrefit-rev-11-edges.R` calls `thres_pin_apply()` directly with
   a synthetic pin carrying the OLD field name `groups`, now `gr_levels`,
   so its log shows "whose levels are ''". That is my script, not the
   package.

## Re-check verdict

MERGEABLE.

Both blocking items are fixed, and the fix is the right one: the codes are
put back on the fitted scale rather than the labels being papered over, so
the refit is the integer-coded refit bit for bit. Verified on my own
constructions and on five the worker did not have: the BOTTOM category, a
non-alphabetical label set, a descending label set, a declared level that
no row takes even in the full fit, and a hand-pinned count that exceeds
the observed categories. Nothing regressed: the bootstrap arm diffs empty
against punch round 0, the inertness log diffs only by the added pin
lines, the controls are identical on both arms, and 44 files give 5618
pass, 0 fail, 0 err.

Remaining, none blocking: two documentation sentences that promise less
than the code delivers (the `data = ` refusal condition, and "top or
interior" where the bottom works too); a pre-existing `frm()` refusal for
non-numeric character responses that the new text now points at; a block
count of 13 that is 14; and a generated count block covering 40 of the 44
files the runner names.

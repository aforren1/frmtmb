# Review of item 2.5f, lane `brmsmatch`

Reviewer's worktree: `C:\Users\adf44\source\r\frmtmb-wt-brmsmatch`, read
only except `dev/bmrev-*` and this file. No git operation was run. Base
`64237e2`. Private library `C:/Users/adf44/source/r/bmrev-lib`, built by
`dev/bmrev-install.R`; core 0.57.0 and the BEFORE arm come from the
shared reference `C:/Users/adf44/source/r/rellib-r3`, read only.
StanHeaders 2.32.10 from `pinlib`, rstan 2.32.7, brms 2.23.0,
posterior 1.7.0, bayesplot 1.16.0.

## Verdict

**Mergeable with named fixes.** Two BLOCKERs, six nits. The three things
item 2.5f asked for are done and I could not falsify any of the three
headline measurements. Both BLOCKERs are about brms-matching claims that
the lane made and did not carry far enough, not about the work it did.

## Scripts this review added

| what | script |
| --- | --- |
| private install | `dev/bmrev-install.R` |
| brms's own bodies and `extract_pars` | `dev/bmrev-bodies.R` |
| the 68 methods, recounted and re-compared | `dev/bmrev-generics.R` |
| independent draws, brms's body on our array | `dev/bmrev-diag.R` |
| behavioural probes | `dev/bmrev-behave.R` |
| `summary.brmsfit`'s own diagnostics | `dev/bmrev-summary-brms.R` |
| last probes | `dev/bmrev-probe2.R` |

Draws saved at `dev/stan-cache/bmrev-draws.rds`: `y ~ x + (1 | g)`,
gaussian, n = 120, **data seed 17**, `frm_sample(chains = 4, iter =
1000, seed = 31337)`, 2000 x 11. Deliberately not the lane's seeds.

## 1. The headline: `identical()` TRUE against brms. NOT FALSIFIED

I tried three ways to break it.

**The lane is re-implementing brms's body, not calling brms.** I printed
the bodies myself (`dev/bmrev-bodies.R`). They are what the lane says:

    brms:::rhat.brmsfit        posterior::summarise_draws(
                                 as_draws_array(x, variable = pars, ...),
                                 rhat = posterior::rhat)
    brms:::neff_ratio.brmsfit  matrixStats::rowMins(
                                 cbind(ess_bulk, ess_tail)) /
                                 ndraws(draws)

`rhat.frmtmb_draws` and `neff_ratio.frmtmb_draws` evaluate the same
composition on `as_draws_array(x)`, with `pmin` for `rowMins` over two
columns. So the claim is "same quantity, computed the same way", not
"same function called". The findings say exactly this, under *An
identity, not a measurement*, and the wording is accurate. The three
places where a re-implementation could still be wrong are the array fed
in, the `rowMins`/`pmin` substitution, and whether `ndraws` is taken
before or after variable selection. All three are right.

**I ran brms's OWN compiled bodies on my own draws.** `dev/bmrev-diag.R`
registers an `as_draws_array` method for a shim class in posterior's
namespace, stubs `brms:::contains_draws` (the only other thing those two
bodies touch `x` with), and calls `brms:::rhat.brmsfit(shim)` and
`brms:::neff_ratio.brmsfit(shim)` so that brms's own bytecode computes
on our array. Against `rhat(ds)` and `neff_ratio(ds)`:

```
rhat        names equal TRUE   identical() TRUE   max |diff| 0
neff_ratio  names equal TRUE   identical() TRUE   max |diff| 0
```

Names on both sides: `Intercept x sigma_Intercept b[1] ... theta_1 lp__`.

**The control says the test is not empty.** On the same fit, rstan's
numbers do differ: sorted `Rhat` by 0.0040022254 relative, sorted
`n_eff/N` by 0.33629868 relative. `max|rhat - 1|` is 0.011884667 (ours)
against 0.012549052 (rstan's), so the two disagree by about a third of
the whole excess over 1 on this fit, the same order the lane reports on
its own.

**The cache-provenance argument.** The caller asked whether reproducing
one recorded figure to eight digits establishes that two draws objects
are the same. It does not, in principle: 0.00563108 is one scalar
functional of a 2000 x 11 array and a match is evidence, not proof. Here
it is moot, and I settled it by a stronger construction: I refit the
lane's construction from scratch (data seed 9, `frm_sample(chains = 4,
iter = 1000, seed = 20260915)`, 3.7 s off the Stan cache) and compared
to `dev/stan-cache/brmsmatch-draws.rds`:

```
identical(old$draws, new$draws): TRUE
max |old - new| over the matrix: 0
```

So the construction is bitwise reproducible and the rebuilt cache is
exactly what the documented recipe produces. What remains unprovable is
that the ORIGINAL `sgrev-draws.rds` used that recipe, and for that the
8-digit match is the only bridge. The findings sentence "so the two are
the same draws" claims more than that bridge carries. **Nit 1: reword to
"so the rebuild reproduces the recipe bitwise, and the recorded figure
agrees to eight digits", and cite `identical()` rather than the figure.**

## 2. BLOCKER: `summary(ds)` should have moved too

Verified true, verified stated in NEWS, and the reasoning behind leaving
it is wrong.

The disagreement is real. On my independent fit:

```
             var  summary Rhat    rhat(ds)           diff
       Intercept   1.01254905   1.0118847   0.00066438446
               x   1.00181275   1.0023347  -0.00052198749
 sigma_Intercept   0.99990396   1.0037623  -0.00385836438
         theta_1   1.00094355   1.0028598  -0.00191621770

             var  summary n_eff   neff_ratio(ds) * N        diff
       Intercept       292.63443            338.67715   -46.04272
               x      1600.37053           1028.59770   571.77283
 sigma_Intercept      1358.72017           1389.70257   -30.98241
         theta_1       491.43531            501.19789    -9.76257
```

`identical()` FALSE. A user reading `summary(ds)` sees `Rhat` below 1 for
`sigma_Intercept` (0.99990) while `rhat(ds)` says 1.00376, and sees 1600
effective draws for `x` where `neff_ratio(ds)` says 1029, a **55.6%
overstatement of the exact number the lane's own NEWS says a user reads
to decide whether to sample longer**.

NEWS does say so, plainly, under *The trade, stated plainly*, and the
vignette says it at more length. The disclosure is honest.

**But brms does not make this trade.** `dev/bmrev-summary-brms.R` prints
`brms:::summary.brmsfit`, which contains:

    c(measures) <- list(Rhat = posterior::rhat,
                        Bulk_ESS = posterior::ess_bulk,
                        Tail_ESS = posterior::ess_tail)

In brms, `summary(fit)`'s `Rhat` column and `rhat(fit)` are the same
number by construction, and brms's ESS columns are `Bulk_ESS` and
`Tail_ESS`, never `n_eff`. So "following brms moves these two away from
agreeing with `ds$stanfit`" is true of `rhat()` but it is not a
consequence of following brms that the package must now disagree with
itself: brms followed all the way leaves no disagreement. "because they
come from `rstan::summary()` by construction" names the current
construction, not a constraint.

Under the user's rule that frmtmb.sample matches brms down to output,
this package now prints a column headed `Rhat`, which is a brms header,
carrying a definition brms does not put behind that header, next to no
label, in the one place a user reads both.

**BLOCKER 1. Move `summary.frmtmb_draws()`'s convergence columns to
posterior's, as brms's `summary.brmsfit` does, and rename them
`Rhat`, `Bulk_ESS`, `Tail_ESS`.** If that is judged out of 2.5f's scope,
the minimum acceptable alternative is to rename the current column so it
does not carry brms's header with a different definition (for example
`Rhat_rstan`), and to say in `?draws-diagnostics` that brms's own
`summary()` reports the posterior quantity. Leaving a `Rhat` column that
means something else is the failure mode the item exists to remove.

## 3. BLOCKER: `pars` on `rhat()`/`neff_ratio()` is not brms's `pars`

This is a divergence the lane INTRODUCED while fixing others, and it is
invisible to the lane's own criterion because the argument NAME matches
at position 2.

brms's `rhat.brmsfit(x, pars = NULL, ...)` does not reach
`extract_pars` at all. It passes `variable = pars` to
`as_draws_array.brmsfit(x, variable, regex = FALSE, ...)`, so in brms
that slot is an EXACT-NAME selector defaulting to NULL. This package
routes it through `draws_select_variables()`, i.e. through brms's
`extract_pars` rule, which is the rule for `as.mcmc`, `mcmc_plot` and
`posterior_interval` but not for these two. Measured
(`dev/bmrev-behave.R`):

| call | brms | here |
| --- | --- | --- |
| `rhat(x)` | all 11 | all 11 |
| `rhat(x, "x")` | 1 variable | 1 variable |
| `rhat(x, "^x$")` | ERROR, missing variable | 1 variable |
| `rhat(x, NULL)` | **all 11** | **ERROR** |
| `rhat(x, NA)` | ERROR | all 11 |

`posterior::subset_draws(arr, variable = "^x$")` errors with "The
following variables are missing in the draws object: {'^x$'}", which is
the brms path. `neff_ratio()` behaves identically to `rhat()` here.

Three of five probed values diverge, and the worst is `rhat(ds, NULL)`:
NULL is brms's own default for that argument, so a ported script that
spells the default out gets "Argument 'pars' must be NA or a character
vector."

The Rd makes this worse by asserting the opposite. `?draws-diagnostics`,
a page whose subject is brms-matching, says of `pars`: "brms refuses a
`pars` that is neither `NA` nor character, and so does this." That is
true of `mcmc_plot()` and false of `rhat()` and `neff_ratio()`, which
share the page.

**BLOCKER 2. Either accept `NULL` on `rhat()`/`neff_ratio()` as "every
variable", brms's own default, or route those two through
`variable = pars` as brms does. Either way, correct the `@param pars`
sentence so it does not claim a brms rule these two methods do not
follow.** The regex-by-default behaviour is more permissive than brms
and I would keep it; it needs a sentence saying it is deliberately more
permissive, not a sentence saying it is brms's.

## 4. The two that answered silently. NOT FALSIFIED

Reproduced end to end.

`brms:::extract_pars(TRUE, ...)` and `extract_pars(0.9, ...)` both stop
with `Argument 'pars' must be NA or a character vector.`; so do `NULL`
and `1L`. Ours:

```
as.mcmc(ds, TRUE)            Argument 'pars' must be NA or a character vector.
posterior_interval(ds, 0.9)  Argument 'pars' must be NA or a character vector.
identical() to brms's string: TRUE
```

`draws_extract_pars()` is a faithful transcription of brms's body, minus
`na_value` and the deprecated `exact_match`.

**No over-refusal.** The calls brms answers are answered:

```
as.mcmc(ds, "^x$")             mcmc.list, chain varnames "x"
posterior_interval(ds, "^x$")  1x2, rownames "x", colnames 2.5% 97.5%
as.mcmc(ds, NA)                mcmc.list, all variables
as.mcmc(ds, c("x","sigma_Intercept"), TRUE)   mcmc.list, fixed = TRUE
posterior_interval(ds, NA, "x", 0.9)          1x2
mcmc_plot(ds, "^x$")           ggplot
```

## 5. The 68 and the six. NOT FALSIFIED, and the criterion is the
deeper scope error

`dev/bmrev-generics.R` counts from the INSTALLED `NAMESPACE`, not from
the lane's script: 126 `S3method()` lines mentioning `frmtmb_draws`,
**68 distinct generic names**, 28 in this package's own methods table.
Re-running the criterion resolved from wherever each method is
registered gives the lane's table character for character:

```
bayes_R2               4    robust             probs
conditional_effects    3    conditions         resp
hypothesis             3    class              alpha
pairs                  2    pars               variable
posterior_summary      2    pars               probs
pp_check               4    prefix             re_formula

diverge: 6   agree as far as both go: 62   no brmsfit method: 0
```

**But "agree as far as both go" is scored by truncating each side at its
own `...`, so any generic where OUR side reaches `...` first scores as
agreeing no matter what brms puts there.** That is the real scope error,
and it is one level below the one the lane found. Measured
(`dev/bmrev-behave.R`); every row below is a positional brms call that
this package ANSWERS with a different thing, silently:

```
call                             brms gives      here
as.matrix(ds, "b_x")             1 column        all 11, dim == as.matrix(ds)
as.array(ds, "x")                1 variable      500x4x11
as_draws_array(ds, "x")          1 variable      500x4x11
as_draws_df(ds, "x")             1 variable      2000x14
fixef(ds, FALSE)                 raw draws       summary, identical() TRUE
ranef(ds, FALSE)                 raw draws       summary, identical() TRUE
coef(ds, FALSE)                  raw draws       summary, identical() TRUE
VarCorr(ds, NULL, FALSE)         raw draws       summary, identical() TRUE
summary(ds, NULL, 0.5)           50% CI cols     identical() to summary(ds)
bayes_R2(ds, NULL, TRUE, TRUE)   robust med/MAD  Estimate Est.Error Q100
```

"identical() TRUE" above means identical to the same call with no
second argument at all, so nothing tells the caller.

The `summary = FALSE` family is the answer to the caller's question
about a DIFFERENT RETURN SHAPE under the same name: `fixef`, `ranef`,
`coef` and `VarCorr` each return a summary where brms returns the raw
draws, and the call is `identical()` to the no-argument call, so nothing
tells the caller. `bayes_R2` is worse than the lane's table suggests:
it is one of the six, and it ANSWERS, with a `Q100` column, rather than
erroring. The item singled out "answers silently" as the severe
subclass; at least one of the six is in it.

The rest of the surface refuses loudly and well: `posterior_samples`,
`nsamples`, `plot`, `loo(ds, ds)` and `waic(ds, ds)` all stop with a
message naming the replacement.

**Nit 2. Record the above in `dev/brmsmatch-findings.md`.** The section
*A NEW finding: the ten were not the whole of the defect* currently
reads as if six is the remaining count. Six is the count the criterion
can see. Extend `dev/brmsmatch-beyond.R` to report, beside each
"agree", whether our side truncated shorter than brms's, and list the
ten rows above as the higher-severity residue. This is out of 2.5f's
scope and correctly not fixed here; it is the filing that is incomplete.

**Nit 3. `pairs(ds, <name not present>)` stops with base R's "subscript
out of bounds".** `pairs.frmtmb_draws` indexes `a[, , keep]` with no
check, where `mcmc_plot()` on the same page gives "variable = names
nosuchvariable, which the draws do not contain. variables() lists what
is there". Pre-existing, unchanged by this lane, and `pairs` is already
filed as one of the six; worth fixing when that item is taken.

## 6. `nuts_params()` and `log_posterior()`. NOT FALSIFIED

    brms:::nuts_params.brmsfit    bayesplot::nuts_params(object$fit, pars, ...)
    brms:::log_posterior.brmsfit  bayesplot::log_posterior(object$fit, ...)

Both printed from the installed brms by `dev/bmrev-bodies.R`. Ours
delegate to the same bayesplot functions on `object$stanfit`, and on my
draws:

```
identical(nuts_params(ds),   bayesplot::nuts_params(ds$stanfit))   TRUE
identical(log_posterior(ds), bayesplot::log_posterior(ds$stanfit)) TRUE
nuts_params(ds) 12000x4, Parameter levels accept_stat__ stepsize__
  treedepth__ n_leapfrog__ divergent__ energy__
```

Leaving them is matching brms, not an omission. As a bonus the
positional `nuts_params(ds, "stepsize__")` works, because our `...`
forwards it into bayesplot's own `pars` slot; it returns the one level.

## 7. The tests, seen failing. NOT FALSIFIED

All four counts reproduce exactly, via the lane's own runner (which does
sum `error`), one file per process, against `rellib-r3` for the base arm
and `bmrev-lib` for mine.

```
file                       blocks  base rellib-r3         my build
test-draws-methods.R       21      PASS 99 FAIL 11 ERR 2  PASS 117 0 0
test-draws-spellings.R     10      PASS 50 FAIL  6 ERR 3  PASS  75 0 0
```

The base failures are BEHAVIOURAL, not missing symbols. Samples:

- `rhat(cs$ds)` returned names `beta[1] beta[2] betad ...`
- `is.na(rhat(cs$ds)["x"])` was TRUE
- `all.equal(rhat(cs$ds), bayesplot::rhat(cs$ds$stanfit))` was TRUE
  where the new test asserts FALSE
- `dimnames(pr)[[2L]]` was `Estimate Est.Error Q2.5 Q97.5` where the
  test wants `Q10 Q90`
- the signature check named all ten with position and both names
- the two base ERRORs are `check_pars(allpars, pars): no parameter ^x$`
  and `if (!summary) return(out): the condition has length > 1`, both a
  positional argument landing on the wrong slot, not a missing symbol

**Both guards carry an inverse case and both fire.**

- The signature guard asserts `length(brms_ten) == 10L` and then asserts
  that `first_divergence()` returns `2L` for `c("object", "ndraws",
  "resp")`, the signature `log_lik()` actually shipped with. That is a
  real inverse case: it proves the comparison machinery can name a wrong
  signature. Its scope is the hardcoded ten, so it cannot catch a
  regression on the other 58; that is the same limit as nit 2.
- The Rd guard prints a control file with one unescaped `%` and scores
  it 1. Run: every one of the seven Rd files scores `source = 2,
  rendered = 0`, the two being roxygen's header. **Nit 4:** the summary
  line "files with an UNESCAPED percent sign in the source: 7" is
  unconditional, since every roxygen Rd has 2. The per-file count is the
  informative number and would move to 3 on an author `%`, so the
  instrument works, but the findings describe it as a guard that "cannot
  fail open" and it is an instrument a human reads. Either subtract the
  header or drop the summary line.

## 8. `R CMD check --as-cran`. CONFIRMED, and the artifact is fresh

I did not re-run it. I verified the existing artifact corresponds to the
tree under review: `md5sum` of every file inside
`brmsmatch-check/frmtmb.sample_0.5.0.tar.gz` against the worktree is
IDENTICAL for `R/methods-draws.R`, `R/loo.R`, `NEWS.md`, both changed
test files and all 17 `man/*.Rd`. `DESCRIPTION` differs only by
`R CMD build` normalization (field wrapping, `Packaged:`, `Author:`,
`Maintainer:`).

```
* checking S3 generic/method consistency ... OK
* checking Rd line widths ... OK
* checking examples with --run-donttest ... [20s] OK
* checking tests ... [208s] OK
* checking re-building of vignette outputs ... OK
* checking PDF version of manual ... OK
* checking HTML version of manual ... OK
Status: OK
```

`grep -c "NOTE\|WARNING" 00check.log` is **0**. Run without
`--no-manual`, so the PDF and HTML manual sections did run. The S3 line
is the one that matters here and it passes: twelve changed signatures
still carry their owner's generic's arguments.

## 9. `is_arg_unset()`. Acceptable, no export wanted

`frmtmb:::is_arg_unset` is, in `frmtmb/R/utils.R:63`, exactly
`inherits(x, "frmtmb_arg_unset")`. The workaround in
`draws_refuse_newdata()` is therefore not an approximation of the
predicate, it IS the predicate, inlined in one place with a comment
saying why. Exporting a one-line `inherits()` wrapper to remove one
`inherits()` call is not worth a public name. Leave it. If core ever
gives the class more structure, the comment is where the next reader
looks.

## 10. `posterior_interval()`'s default variable set. A divergence

Measured: ours returns 4 of 11 rows, dropping `b[1] ... b[6]` and
`lp__`. brms's `posterior_interval.brmsfit` is
`as.matrix(object, pars = NA, variable = NULL)` into
`rstantools::posterior_interval()`, so brms returns all of them.

Under the stricter standard this is a divergence, not "a default rather
than a divergence". The distinction the findings draw, that a default is
not a positional slot, is the right distinction for item 2.5f's scope
and the wrong one for "matches brms down to output": a default IS
output. A ported script doing
`posterior_interval(fit)["r_g[1,Intercept]", ]` gets a subscript error
here.

The lane's body text already concedes this ("It is a real remaining
difference from brms and a candidate for a later item"), which is the
honest version. **Nit 5: change the heading sentence from "That is a
default, not a positional slot" to call it a divergence that is out of
this item's scope.** Correctly deferred either way; the argument that a
large fit has thousands of group-level rows is real, and brms's own
`summary()` drops them too.

## 11. House style

No new line over 80 columns in any package file the lane changed
(checked over the diff's added lines only). No em dash, no emoji, no
spaced hyphen standing in for one in the added text. **Nit 6:** three
lines over 80 in `dev/brmsmatch-findings.md`, lines 67 to 69, at 99, 91
and 99 columns. They are inside the generated block, pasted verbatim, so
the fix belongs in `dev/brmsmatch-summary.R`, which should wrap the name
lists.

## What I could not falsify, and what I tried

- The two diagnostics equal brms's. Tried: printing brms's bodies
  independently; running brms's own bytecode on an independently sampled
  array through a shim; a control showing rstan's numbers do differ on
  the same fit. `identical()` TRUE both, max |diff| exactly 0.
- The frmtmb names. `variables(ds)` and `names(rhat(ds))` are equal,
  11 of 11, on my own draws.
- The refusal, character for character, with no over-refusal.
- 68 methods, six divergences, reproduced from the NAMESPACE by an
  independently written script.
- `nuts_params()` / `log_posterior()` delegate, `identical()` to
  bayesplot on the stanfit.
- All four test counts, and the base failures are behavioural.
- `R CMD check` Status: OK, zero NOTEs, S3 consistency OK, artifact
  verified fresh by md5.
- The cache. Stronger than the lane's argument: bitwise identical on
  refit.

## Fixes, ranked

1. **BLOCKER.** `summary(ds)`'s `Rhat` and `n_eff` must move to
   posterior's, as `brms:::summary.brmsfit` does, or be renamed so they
   do not carry brms's header with another definition.
2. **BLOCKER.** `rhat()` and `neff_ratio()` must accept `pars = NULL`,
   brms's own default for that slot, and the `@param pars` sentence
   must stop claiming a brms rule these two do not follow.
3. Nit. Reword the cache-provenance sentence: `identical()` on refit,
   not "the two are the same draws" from one 8-digit match.
4. Nit. Record the ten silent positional divergences the
   truncate-at-dots criterion cannot see, and that `bayes_R2` is one of
   the six that ANSWERS rather than erroring.
5. Nit. Call `posterior_interval()`'s default variable set a
   divergence, not a default.
6. Nit. `pairs(ds, <absent name>)` gives "subscript out of bounds";
   give it the message its siblings on the same page give.
7. Nit. The Rd `%` summary line is unconditional; subtract roxygen's
   two or drop the line.
8. Nit. Three lines over 80 in the generated block of
   `dev/brmsmatch-findings.md`; wrap them in
   `dev/brmsmatch-summary.R`.

---

# Re-check, 2026-09-16: the two BLOCKERs only

Narrow re-check at the coordinator's request. The six nits are not
reopened and nothing below revisits them. Package reinstalled into
`C:/Users/adf44/source/r/bmrev-lib` from the current worktree. Every
number is on the REVIEWER's draws, `dev/stan-cache/bmrev-draws.rds`
(data seed 17, `frm_sample(chains = 4, iter = 1000, seed = 31337)`),
never the lane's. Scripts: `dev/bmrev-recheck.R`, `dev/bmrev-identity.R`.

**Both BLOCKERs are CLOSED.** Two nits, no new BLOCKER.

## BLOCKER 1, `summary(ds)`: CLOSED

```
colnames now:  mean sd 2.5% 97.5% Rhat Bulk_ESS Tail_ESS
colnames base: mean sd 2.5% 97.5% n_eff Rhat
n_eff column still present: FALSE
```

**1. The three columns are brms's measures.** Against
`posterior::summarise_draws(a, Rhat = posterior::rhat, Bulk_ESS =
posterior::ess_bulk, Tail_ESS = posterior::ess_tail)`, the call shape
`brms:::summary.brmsfit` literally contains, on the same subset array:

```
Rhat      identical() TRUE   max|diff| 0
Bulk_ESS  identical() TRUE   max|diff| 0
Tail_ESS  identical() TRUE   max|diff| 0
```

For `Rhat` I went further and ran **brms's own `rhat.brmsfit` body** on
my draws through the shim from the first pass:
`identical(summary(ds)[, "Rhat"], brms:::rhat.brmsfit(shim)[rownames])`
is TRUE, max |diff| 0. So `Rhat` is brms's own computation, not
posterior's called some other way.

The limit, stated: brms exposes `Bulk_ESS` and `Tail_ESS` through no
method of its own. `neff_ratio.brmsfit` returns
`pmin(bulk, tail) / ndraws` and `summary.brmsfit` needs a real
`brmsfit` to build its table. So for those two the strongest available
check is the `summarise_draws` call shape lifted from
`summary.brmsfit`'s body, which is what I used. `Rhat` carries the
end-to-end check for all three.

**The disagreement my first pass measured is gone.**

```
summary Rhat vs rhat(ds):  identical() TRUE, max|diff| 0
variables on opposite sides of 1:      0   (was 1: sigma_Intercept,
                                            0.99990396 against 1.0037623)
rows addressable by rhat(ds)[rowname]: 4 of 4
```

The row COUNT is unchanged at 4 of 11: `lp__` and the six `b[]` modes
are still dropped, so this did not quietly move the deferred item in
section 10.

**2. The `pmin` identity: an IDENTITY, and the direction matters.**

Both sides call `posterior::ess_bulk` and `posterior::ess_tail` on the
same per-variable draws with the same `ndraws`, so this is arithmetic,
not agreement between two estimates. What `identical()` adds is that the
construction has not drifted: same array, same draw count, no subsetting
artifact. That is worth having and it is not an independent measurement.

The direction is not cosmetic:

```
neff_ratio(ds)[keep] == pmin(B, T) / ndraws(ds)
  identical() TRUE    max|diff| 0

pmin(B, T) == neff_ratio(ds) * ndraws(ds)
  identical() FALSE   max|diff| 5.6843418860808015e-14   0.51 ulp
  (one row, theta_1: 501.19788626479237 against 501.19788626479232)
```

`NEWS.md` states the DIVISION direction and the new test asserts the
division direction with `expect_equal`, so both are correct and the test
is not flaky. The multiplication phrasing is the one that is not
bitwise. **Nit A: keep the division wording; do not let "`pmin` is
exactly `neff_ratio * ndraws`" reach a test or a NEWS line, because
`x/N*N` is a round trip and loses up to half an ulp.**

**3. Nothing else moved.** `summary(ds)` captured from the BASE build in
a child process against `rellib-r3` and compared column for column:

```
mean    identical() base vs now TRUE   max|diff| 0
sd      identical() base vs now TRUE   max|diff| 0
2.5%    identical() base vs now TRUE   max|diff| 0
97.5%   identical() base vs now TRUE   max|diff| 0
Rhat    identical() base vs now FALSE  max|diff| 0.0038583643830625292
rownames identical: TRUE
only in base: n_eff      only in now: Bulk_ESS Tail_ESS
```

All four distribution columns are also `identical()` to the statistic
recomputed straight off `ds$draws`. Exactly one column changed
definition, `Rhat`, which is the one that was supposed to. The
positional consequence is worth knowing and is not a defect: the table
goes from 6 columns to 7 and `Rhat` moves from position 6 to position 5,
so `summary(ds)[, 5]` changes from `n_eff` to `Rhat`. NEWS marks the
change BREAKING and says the `n_eff` column is gone.

The new test block was SEEN FAILING. Against `rellib-r3` with the
current test file, "summary() reports brms's three diagnostics, agreeing
with rhat()" gives 3 failures and 1 error, all behavioural: `n_eff`
present, `colnames` wrong, `s[, "Rhat"]` unequal to `rhat(ds)`, and
`s[, "Bulk_ESS"]` a subscript-out-of-bounds because the column does not
exist.

## BLOCKER 2, `pars`: CLOSED

```
here  rhat(x, pars = NULL, regex = FALSE, ...)
here  neff_ratio(object, pars = NULL, regex = FALSE, ...)
brms  rhat(x, pars = NULL, ...)
brms  neff_ratio(object, pars = NULL, ...)
```

Nine values probed side by side, ours against `brms:::rhat.brmsfit` run
on the same array. Both the outcome and the message agree on all nine.

```
pars          brms                               here
missing       OK n=11                            OK n=11
"x"           OK n=1  [x]                        OK n=1  [x]
"^x$"         ERROR: variables are missing ...   same
NULL          OK n=11                            OK n=11
NA            ERROR: Assertion ... Contains ...  same
TRUE          ERROR: Assertion ... not logical   same
0.9           ERROR: Assertion ... not double    same
"^b" regex=T  OK n=6  [b[1],b[2],b[3]]           OK n=6
c("x","lp__") OK n=2  [x,lp__]                   OK n=2
                                      matching 9 of 9
```

`neff_ratio()` on the first five: matching 5 of 5. `rhat(x, NULL)`
returns all 11 and `rhat(x, NA)` errors, both as asked.

**The refusals are `posterior::subset_draws()`, not a look-alike.**

```
rhat(ds, TRUE)  msg:  Assertion on 'variables' failed: Must be of type
                      'character', not 'logical'.
                call: check_existing_variables(variable, x, regex = regex,
                        exclude = exclude, scalar = scalar)
```

`conditionCall()` names posterior's own `check_existing_variables()`.
Calling `posterior::subset_draws(arr, variable = TRUE)` directly gives
the identical string, as do the `0.9` and the missing-name cases. And
the two method bodies plus `draws_diag_array()` no longer mention
`draws_extract_pars`, `draws_select_variables`, or the string "must be
NA or a character vector"; they mention `subset_draws`. So these two
have genuinely left the `extract_pars` rule rather than reproducing it.

**Dropping `variable` and `fixed` breaks no call brms answers.** The
lane's reasoning was plausible; it is now measured:

```
brms  rhat(x, variable = "x")     ERROR: formal argument "variable"
                                  matched by multiple actual arguments
here  rhat(ds, variable = "x")    OK n=11
brms  rhat(x, "x", fixed = TRUE)  OK n=1 [x]
here  rhat(ds, "x", fixed = TRUE) OK n=1 [x]
brms  rhat(x, "x", inc_warmup=F)  OK n=1 [x]
here  rhat(ds, "x", inc_warmup=F) OK n=1 [x]
```

brms answers NO call that spells `variable=` on these two, so removing
it costs nothing. `fixed=` and `inc_warmup=` are not over-refused: both
sides drop them into `...` and return the same thing.

**Nit B, in the opposite direction from over-refusal.**
`rhat(ds, variable = "x")` returns all 11 SILENTLY where brms errors,
because `...` absorbs it and `draws_diag_array()` takes no `...`. The
shared Rd still carries `\item{...}{Passed to the bayesplot function.}`,
which is true of `mcmc_plot()`, `pairs()`, `nuts_params()` and
`log_posterior()` on that page and false of these two, whose `...` goes
nowhere. Either refuse an unused named argument on these two or split
the `...` documentation so it stops promising a forwarding they do not
do. This is a nit and not a regression: 0.5.0 ships `rhat(x, ...)` and
`neff_ratio(object, ...)`, so `variable=` was never a released argument
on either.

The new `pars` test block was also seen failing: against `rellib-r3`,
"rhat() and neff_ratio() take brms's OTHER `pars` rule" errors with
`check_pars(allpars, pars): no parameter x`, the base build handing
`pars` to rstan.

## Tests after the change

One file per R process, never more than two at once.

```
file                     libs        blocks  result
test-draws-methods.R     bmrev-lib   22      PASS 127 FAIL 0 ERROR 0
test-draws-methods.R     rellib-r3   22      PASS  99 FAIL 14 ERROR 3
test-draws-spellings.R   bmrev-lib   10      PASS  75 FAIL 0 ERROR 0
test-sample-direct.R     bmrev-lib   17      PASS 136 FAIL 0 ERROR 0
```

`test-sample-direct.R` is included because it reads
`summary(...)[, "Bulk_ESS"]`, a column that did not exist before this
change; it passes.

## One factual note for the merge

The `R CMD check --as-cran` artifact in
`C:/Users/adf44/source/r/brmsmatch-check` is now **STALE**: `md5sum`
differs from the worktree for `R/methods-draws.R`, `NEWS.md` and
`tests/testthat/test-draws-methods.R`. Section 8's "Status: OK, zero
NOTEs" was verified against the earlier source and does not carry to
this one. The check must be re-run before merge. Not a review finding,
just a fact about the artifact.

## Re-check verdict

Both BLOCKERs closed. Two nits, A and B above, neither blocking. With
the eight nits from the first pass this is now **mergeable**.

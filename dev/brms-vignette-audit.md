# brms vignette audit: the hand-translated port

Audit date 2026-09-03, frmtmb v0.42.0 (worktree `wt-brmsvig`), brms
2.23.0, R 4.6.1 on Windows.

This audit replaces the mechanical measurement in
`dev/brms-vignette-port.md`. That one applied an AST transform to every
vignette chunk and counted what survived. This one is a HAND
TRANSLATION: a person read each vignette, wrote the frmtmb spelling of
every call, ran it, and labeled every divergence. The two answer
different questions. The mechanical audit answers "what does
find-and-replace get you". This one answers "what does a competent
porter meet, and what does the output look like when it runs".

Nothing in `R/` or `tests/` was changed. Edges are recorded, not fixed.
One known-broken area, the prior interface, is deliberately exercised at
its CURRENT spelling so the findings stay honest while a sibling lane
fixes it.

**Revision at v0.43.0 (consolidation, same date).** The prior lane
landed: `priors =` is renamed `prior =`, `prior()`/`prior_()`/
`prior_string()` are exported, `set_prior()` takes `nlpar`/`resp`, and
a `brmsprior` is translated rather than refused. The five scripts that
passed priors, plus the coexistence instrument, were swept to the new
spelling and re-run green at v0.43.0; their labels otherwise still
describe the v0.42.0 measurement. Also applied at consolidation, so no
longer open below: the `conditional_effects` guard fix (rank 5's
`method = "predict"` refusal), and the refusal pass over the draws and
fit surfaces (the silent `rescor_matrix()` `NULL`, two-model `loo()`,
draws `update()` and `plot()`, bare `loo`/`waic`/`bayes_R2`/
`expose_functions` on a fit). A full re-measurement against the new
API is the next audit round, not this file.

**Revision after v0.43.0 (the plot lane).** Break 2 below is fixed:
`plot.frmtmb_conditional_effects` now draws the condition sets as small
multiples on one page (tinyplot when installed, a `par(mfrow)` grid of
base panels otherwise), honors the `ncol` that used to be swallowed by
`...`, and shows in each panel only the observations belonging to that
condition, following brms's `make_point_frame()`. The rows below still
describe the v0.42.0 measurement and are not relabeled.

## Method

One runnable script per vignette under `dev/brms-vignettes/`. Each
script covers BOTH inference paths, in two banner-marked sections:

- **PATH 1, ML**: `frm()` and the frequentist post-processing surface.
- **PATH 2, SAMPLE**: `frm_sample()` and the `frmtmb_draws` surface.

The second path matters because the brms vignettes show POSTERIOR
output. A `frmtmb_draws` object is the closer analogue, and several
calls the ML path has to give up come back there.

Every translated call goes through `bv()` in
`dev/brms-vignettes/_harness.R`, which runs it, catches the error, and
writes one CSV row. The tables below come out of
`dev/brms-vignettes/_scoreboard.R` reading those rows, and the edge
inventory out of `_summarize.R`. A failing call does not stop the
script: the point is to record, not to block.

Each row carries one of five labels, applied to the TRANSLATION and not
to the outcome:

| label | meaning |
|---|---|
| CLEAN | the brms line runs unchanged apart from `brm` to `frm` and the removal of MCMC-only arguments |
| SPELLING | works, under a different name or argument |
| BEHAVIOR | runs, but the output differs from what brms documents |
| MISSING | no frmtmb path |
| REFUSAL | frmtmb refuses on purpose |

The path is carried in the label prefix, `ML:` or `SAMPLE:`.

frmtmb is loaded with `pkgload::load_all()`. brms is NOT attached in any
script but `brms-coexistence.R`: a post-processing name has to resolve
to frmtmb or fail, and that is the measurement. brms data sets are
reached with the `brms::` prefix.

Sizes are modest by instruction. Sampling is one chain of at most 500
iterations against brms's four chains of 2000, and bootstrap bands take
25 refits against the 200 default. Every place a short run costs the
vignette's illustrated output is recorded as a BEHAVIOR row rather than
papered over.

One sizing caveat worth stating, because it is not obvious and it cost
this audit real time. On a model with CROSSED CORRELATED random-effect
blocks, sampling cost tracks the geometry of those blocks and not the
row count. On `brms_multivariate`'s bivariate model a random 30-level
subset of `dam` (about 220 rows) samples in about 8 seconds, while the
FIRST 30 levels of `dam` (201 rows, fewer rows) took 30 to 60 seconds.
Subset a grouping factor at random, not by taking a prefix.

Reproduce:

```sh
export BV_OUT=/some/dir
for f in dev/brms-vignettes/brms_*.R dev/brms-vignettes/brms-*.R; do
  Rscript "$f" > "$BV_OUT/bv-$(basename "$f" .R).log" 2>&1
done
Rscript dev/brms-vignettes/_summarize.R    # every edge, one per line
Rscript dev/brms-vignettes/_scoreboard.R   # the tables in this file
```

## Headline

**567 translated calls across twelve scripts, 508 seconds of compute.**
Ten scripts are the shipped brms vignettes.
`dev/brms-vignettes/brms-coexistence.R` is this audit's own instrument
and is scored apart. `brms_threading` has no model to translate.

Every count below is generated by `dev/brms-vignettes/_scoreboard.R`
from the recorded CSVs, not typed by hand.

The claim under test, from `README.md` and `docs/index.md`: *a measured
audit of the brms vignettes puts about 7 of 10 of their model calls
through that transform unchanged.*

**The claim holds for model calls and is the wrong number to lead
with.**

On the ML path, across the ten vignettes:

| | calls | port unchanged | share |
|---|---|---|---|
| model calls | 56 | 38 | 68% |
| post-processing calls | 228 | 80 | 35% |
| **all calls** | **284** | **118** | **42%** |

"About 7 of 10" is close to right for MODEL calls, and it no longer
needs the earlier audit's concession about the family argument: `frm()
DOES default to gaussian` now, on the plain, `bf()`, `nl = TRUE` and
multivariate paths alike, so a vignette line with no `family` argument
ports verbatim. FN-1 is fixed and the claim stands on its own.

It is still not defensible as a statement about porting a VIGNETTE,
because a vignette is mostly post-processing and only about a third of
that ports unchanged. **The honest headline is 4 of 10 calls, not 7 of
10.** The 7-of-10 figure should be qualified as a model-call statistic
wherever it appears, or replaced.

The earlier figure also came from that audit's PROJECTION pass, which
stood in for two fixes not yet shipped. Both shipped. This number is
measured on v0.42.0.

On the sampling path, the same ten vignettes:

| | calls | port unchanged | share |
|---|---|---|---|
| model calls | 46 | 11 | 24% |
| post-processing calls | 140 | 55 | 39% |
| **all calls** | **186** | **66** | **35%** |

The model-call share falls because ONE structural change hits every
model: `brm(formula, data, chains =, iter =)` becomes two calls,
`frm()` then `frm_sample(fit, chains =, iter =)`, or one
`frm_sample(formula, data, ...)`. That is a single well-understood
edit, not 35 different ones. Post-processing does BETTER here than on
the ML path, and that is the finding that matters: the draws surface is
the closer analogue of what the vignettes actually show.

Of 56 distinct model calls, 29 port unchanged on at least one path.

## Per vignette

Model calls, then post-processing calls, as
`total | clean spelling behavior missing refusal`.

### ML path

| vignette | model | post |
|---|---|---|
| brms_customfamilies | 2 \| 1 1 0 0 0 | 16 \| 4 5 2 4 1 |
| brms_distreg | 7 \| 5 1 0 0 1 | 17 \| 7 1 7 2 0 |
| brms_families | 1 \| 1 0 0 0 0 | 43 \| 28 4 3 8 0 |
| brms_missings | 4 \| 4 0 0 0 0 | 13 \| 2 3 4 0 4 |
| brms_monotonic | 10 \| 5 1 0 3 1 | 15 \| 5 1 5 4 0 |
| brms_multilevel | 10 \| 8 1 1 0 0 | 22 \| 4 5 7 3 3 |
| brms_multivariate | 3 \| 3 0 0 0 0 | 21 \| 5 2 5 5 4 |
| brms_nonlinear | 6 \| 2 4 0 0 0 | 29 \| 13 6 5 2 3 |
| brms_overview | 7 \| 3 1 1 2 0 | 18 \| 4 1 7 6 0 |
| brms_phylogenetics | 6 \| 6 0 0 0 0 | 26 \| 8 2 13 2 1 |
| brms_threading | 0 | 8 \| 0 0 1 7 0 |
| *brms-coexistence* | 4 \| 2 1 0 0 1 | 29 \| 15 9 2 3 0 |

### SAMPLE path

| vignette | model | post |
|---|---|---|
| brms_customfamilies | 2 \| 2 0 0 0 0 | 9 \| 5 1 1 1 1 |
| brms_distreg | 4 \| 0 4 0 0 0 | 23 \| 12 1 8 1 1 |
| brms_families | 0 | 0 |
| brms_missings | 3 \| 1 1 0 1 0 | 12 \| 8 0 1 0 3 |
| brms_monotonic | 7 \| 0 4 1 1 1 | 11 \| 5 2 3 1 0 |
| brms_multilevel | 8 \| 4 2 2 0 0 | 18 \| 7 1 6 1 3 |
| brms_multivariate | 4 \| 0 3 1 0 0 | 17 \| 2 3 7 3 2 |
| brms_nonlinear | 6 \| 2 3 1 0 0 | 14 \| 5 5 3 0 1 |
| brms_overview | 5 \| 2 1 2 0 0 | 15 \| 5 2 6 0 2 |
| brms_phylogenetics | 7 \| 0 5 1 1 0 | 21 \| 6 4 10 1 0 |
| brms_threading | 0 | 0 |
| *brms-coexistence* | 2 \| 1 0 1 0 0 | 34 \| 29 0 0 0 5 |

Notes:

- **brms_phylogenetics** is the cleanest: six of six ML model calls port
  unchanged, including `gr(phylo, cov = A)` with `data2`, the `se()`
  meta-analysis, the observation-level Poisson model and `update()`.
  It loses ground only in post-processing, and almost every loss is one
  of two things: `sigma` printed on the log link scale, and
  `hypothesis()` without an evidence ratio.
- **brms_missings** ports all four ML model calls, `mi()` models
  included, and `frm_multiple` post-processing now REFUSES with
  messages that name the workaround. The v0.34 audit's FN-11 complaint,
  that `plot()` failed with a base-graphics message, is fixed.
- **brms_multilevel** ports eight of ten, including `mm()` with weights
  (FN-2, fixed), the `ult + omega + theta ~ 1 + (1|ID1|AY)` shorthand
  (FN-8, fixed) and `update()` with a new formula (FN-9, fixed).
- **brms_nonlinear** scores worst on model calls (two of six) for a
  single reason: every nonlinear model needs `start`, because frmtmb
  begins at zero where the gradient is `NaN` and brms's priors did that
  job. Nothing else about those four calls changes.
- **brms_multivariate** ports all three ML model calls unchanged,
  `mvbind()`, `set_rescor()`, `bf() + bf()`, per-response families,
  `lf()` and a smooth in one response included. It loses ground only
  after the fit: `pp_check()`, `fitted()` and `bayes_R2()` all refuse a
  multivariate object, and `add_criterion()` has no analog on either
  path.
- **brms_monotonic** loses three model calls to one gap: `mo()`'s
  simplex takes no prior, so the vignette's whole prior section cannot
  be written.
- **brms_families** is a coverage script, not a translation: 28 of 43
  named families are accepted as brms spells them.

## The coexistence story

`dev/brms-vignettes/brms-coexistence.R` is the one script that attaches
brms on purpose.

With **frmtmb first, brms second** (the natural order for a brms user
trying frmtmb), 92 of frmtmb's 134 exports are masked. The consequences
divide in three:

1. **A loud, correct failure.** `frm(bf(y ~ x), ...)` refuses: *"this
   formula was built by brms::bf(): attaching brms after frmtmb masks
   frmtmb's bf(), so a bare bf() call now reaches brms. Call
   frmtmb::bf() explicitly, or attach brms before frmtmb."* This is the
   best message in the audit. Without it the failure would be an
   unreadable structure error.
2. **Silent wrong resolution.** `set_prior()`, `get_prior()`,
   `mixture()`, `custom_family()`, `lf()` and `nlf()` all resolve to
   brms and SUCCEED, returning a brms object. `get_prior()` is the
   worst: a user reads brms's default priors and believes they describe
   the frmtmb fit. Nothing warns. `set_prior()` is caught one step
   later, when `frm(priors =)` refuses the `brmsprior`. (v0.43.0:
   `frm(prior =)` now translates a `brmsprior` row by row, so the
   `set_prior()` collision is absorbed rather than caught late.)
3. **Silent right resolution.** 29 family constructors are masked, so
   `family = student()` hands `frm()` a `brmsfamily`. `frm()` converts
   it. That is the good outcome, but it also means a link or dpar
   difference between the two packages would never be reported.

Every S3 generic dispatches correctly in both attach orders:
`conditional_effects`, `hypothesis`, `pp_check`, `ngrps`, `fixef`,
`ranef` and `VarCorr` all reach the frmtmb method on a `frmtmb_fit` and
a `frmtmb_draws` with brms attached. 29 of 34 draws-surface calls run
clean under a brms attach. The collisions are all on CONSTRUCTORS, not
on methods.

With **brms first, frmtmb second**, the masking reverses and a copied
brms line silently changes meaning in the other direction. At the
measured v0.42.0, `prior()` resolved to brms in both orders, because
frmtmb did not export it; v0.43.0 exports `prior()`, so it masks like
the other constructors, and either package's object now works in
`frm(prior =)` because a `brmsprior` is translated.

## The conditional_effects faceting diagnosis

The mandated deep dive, from the insurance-loss model in
`brms_nonlinear`. The finding is that frmtmb's grid machinery is
CORRECT and two separate things downstream of it are not.

### What brms actually does

The vignette's own lines, verbatim from
`brms/doc/brms_nonlinear.Rmd:215`:

```r
conditions <- data.frame(AY = unique(loss$AY))
rownames(conditions) <- unique(loss$AY)
me_loss <- conditional_effects(
  fit_loss, conditions = conditions,
  re_formula = NULL, method = "predict"
)
plot(me_loss, ncol = 5, points = TRUE)
```

`re_formula` defaults to `NA` in brms too, so the default is not what
produces the per-year curves. The vignette passes `re_formula = NULL`
explicitly. Each row of `conditions` becomes one condition set; with
`re_formula = NULL` the `AY` value in that row selects that year's
group-level effect; the row name becomes the `cond__` label; and
`brms:::plot.brms_conditional_effects` adds
`facet_wrap(facets = "cond__", ncol = ncol)` to a ggplot. Ten rows give
ten facets in a 2 by 5 grid.

### What frmtmb does, measured

The grid layer is right, on both paths.

```
conditions df + re_formula = NULL, fit path, band = "boot"
  1991   1992   1993   1994   1995   1996   1997   1998   1999   2000
3872.3 5205.2 5163.6 5355.7 4684.3 4971.7 5517.9 5994.4 5229.3 5119.1

the same call with re_formula = NA (the default)
  1991   1992   1993   1994   1995   1996   1997   1998   1999   2000
5111.3 5111.3 5111.3 5111.3 5111.3 5111.3 5111.3 5111.3 5111.3 5111.3
```

Ten distinct curves under `NULL`, one repeated curve under `NA`, exactly
as brms documents. The `cond__` column is present, and it carries the
data frame's ROW NAMES, so `rownames(conditions) <- unique(loss$AY)`
does what it does in brms. A grouping factor as the varying column of
`conditions` is accepted. `conditions` as a plain named list works and
correctly produces no `cond__`. A one-sided `re_formula = ~ (1 | AY)`
works. `re.form`, the lme4 spelling, is refused with a message that
names `re_formula`.

So the chain does not break in grid construction, in the group column,
or in the `cond__` labeling. It breaks twice, further down.

### Break 1: the vignette's exact call is refused, and the refusal is wrong

```r
conditional_effects(fit_loss, conditions = conditions,
                    re_formula = NULL, method = "predict")
#> Error: conditional_effects() cannot put a wald band on a nonlinear
#>   predictor: predict() has no standard error for it. Use band =
#>   "boot", which refits instead of differentiating, or display one
#>   nonlinear parameter with dpar = "ult"
```

The guard is `R/conditional-effects.R:994`:

```r
if (!is.null(lp$nl_body) && band != "boot") {
```

`band` defaults to `"wald"`, so the guard fires. But `method` is never
consulted, and `method = "predict"` does not use a Wald band. Follow the
code: the `else` branch at `R/conditional-effects.R:1058` computes
`predict(..., se.fit = TRUE)` and fills `estimate__`, `se__`, `lower__`
and `upper__`; then the `if (method == "predict")` block immediately
below overwrites all four from `ndraws` simulated responses. The Wald
computation is dead work on that path. The guard rejects a call that
would have succeeded.

This is what the maintainer met. The refusal names `band = "boot"`,
which needs `method = "epred"`, so the porter drops `method` too; then
drops `conditions` and `re_formula` to get anything to draw; and lands
on a single population-level curve. The message is accurate about the
Wald band and misleading about this call.

### Break 2: plot() does not facet

`plot.frmtmb_conditional_effects` at `R/conditional-effects.R:1267`:

```r
if (!is.null(df$cond__) && length(unique(df$cond__)) > 1L) {
  for (cv in unique(df$cond__)) {
    ...
    ce_plot_one(sub, cond = cv, points = points)
  }
}
```

It loops and draws one BASE-GRAPHICS panel per condition. Ten conditions
send ten plots to the device in sequence. There is no facet grid, no
`par(mfrow)` is set, and `ncol` is absorbed by `...` and silently
discarded with no warning. On a single-page device the first nine are
overwritten. That is the second half of "a group-level single plot".

The workaround is `par(mfrow = c(2, 5))` before `plot()`.

### The sampling path breaks in a third place, and better

```r
conditional_effects(fl_s, conditions = conditions,
                    re_formula = NULL, method = "predict")
#> Error: conditional_effects() on draws has no method =: the curves ARE
#>   posterior expected-response draws. For predictive bands, quantile
#>   posterior_predict() over your own grid
```

`conditional_effects.frmtmb_draws` has no `method` and no `band`
argument at all. Drop `method = "predict"` and the vignette's display
comes out correctly with no other change, per-year curves and all: the
draws supply the band, so the nonlinear predictor is no obstacle. The
message names the right replacement. This is a good refusal.

Break 2 is unchanged on this path. `ce_grids_build` and `ce_finalize`
are shared by both methods (`R/conditional-effects.R:1005` and `:1135`
for the fit, `:1189` and `:1245` for the draws), and both return the
same class, so both go through the same non-faceting `plot()`.

### The smallest fix, and its location

Two lines, both in `R/conditional-effects.R`, for break 1:

```r
# line 994
if (!is.null(lp$nl_body) && band != "boot" && method != "predict") {

# line 1058
if (band == "boot" || method == "predict") {
```

The second line stops the dead `se.fit` call that the simulation block
overwrites anyway.

Verified by patching the closure in place and comparing against the
unpatched one. On a linear Poisson mixed model the two agree exactly:
identical under `conditions` plus `re_formula = NULL` with
`method = "predict"`, identical with no conditions, and identical on the
default epred path. The patched version is also about 34% faster on
`method = "predict"`, because it no longer runs a standard-error report
whose result is discarded. On the loss model the patched call returns
the ten per-year prediction bands the vignette asks for.

Break 2 is a bigger decision and is NOT a two-line fix, because it is a
choice about the plotting stack. The cheapest honest version is to make
`plot.frmtmb_conditional_effects` accept `ncol` and set
`par(mfrow = c(ceiling(n / ncol), ncol))` around the `cond__` loop,
restoring it on exit. Faceting the way brms does needs ggplot2, which
frmtmb does not depend on.

A third, smaller item found on the way: `re_formula = NULL` with NO
`conditions` conditions on whatever group level the reference rule
lands on. `ce_ref_value` takes the MEAN of a numeric column, so on
`loss` it happens to give exactly 1994, a real level. On any unbalanced
numeric grouping variable whose mean is not a level, the same call ends
in `New levels in grouping factor`. The error is clear, but the vignette
idiom that avoids it is `conditions =`, and nothing says so.

## The ranked edge list

308 of the 567 calls carry an edge. They are not 308 different
problems. Grouped by cause and ranked by frequency times severity,
where severity asks how badly the edge misleads a reader who is
comparing frmtmb's output against the vignette's printed output.

Counts are rows matching the theme, so a row that carries two causes is
counted under both.

| rank | edge | rows | ML / SAMPLE | vignettes | severity |
|---|---|---|---|---|---|
| 1 | prior interface | 49 | 28 / 21 | 8 | high |
| 2 | `loo` family absent on a point fit | 47 | 24 / 23 | 10 | medium |
| 3 | internal or link scale in `summary()` | 33 | 13 / 20 | 7 | **high** |
| 4 | `conditional_effects()`: bands, `method =`, faceting | 30 | 18 / 12 | 7 | high |
| 5 | ordinal thresholds and the `mo()` simplex | 19 | 13 / 6 | 3 | high |
| 6 | `plot()` draws diagnostics, arguments absorbed | 18 | 10 / 8 | 6 | medium |
| 7 | `hypothesis()` has no evidence ratio | 19 | 11 / 8 | 7 | medium |
| 8 | Stan-only surface | 18 | 16 / 2 | 3 | low |
| 9 | `summary()` header and block layout | 14 | 10 / 4 | 7 | medium |
| 10 | multivariate post-processing | 9 | 6 / 3 | 2 | medium |
| 11 | no `plot()` method for `frmtmb_draws` | 2 | 0 / 2 | 2 | **high** |
| - | short chain, not a porting gap | 17 | 2 / 15 | 9 | none |

Ranked by what to fix first, which is not the same order:

### 1. `summary()` prints internal and link-scale parameters (rank 3)

The highest-severity item, and the one a reader meets FIRST. It is not
a missing feature; it is a presentation choice that makes correct
numbers look wrong.

- On an ML fit, random-effect standard deviations print naturally
  (`Std.Dev. 598.51`), but a distributional parameter prints its
  linear-predictor coefficient. `sigma 4.89` is `log(133)` where brms
  reports `139.93`. A constant `zi` prints `-0.37` where brms reports
  `0.41`. Ordinal thresholds print as `tau_raw` in a (first threshold,
  log increments) parameterization. An `mo()` simplex prints as
  unconstrained `zeta` and never as a simplex.
- On a draws object it is worse: EVERYTHING is internal. Variance
  components appear as `theta_1 ... theta_k` on the log or Cholesky
  scale with no names at all. `VarCorr()` is the only route to a named
  natural-scale value, and nothing in `summary()` says so.

Seven of the ten vignettes hit some form of this. The fix is
presentational and does not touch any estimate: give `summary()` a
natural-scale block for distributional parameters and named rows for
variance components on the draws path.

### 2. `mo()`'s simplex has no prior, so it has no posterior (rank 5)

The most consequential single finding after the faceting diagnosis.
brms's implicit `dirichlet(1)` on a monotonic simplex is load-bearing.
frmtmb has no way to supply it: `set_prior()` has no `dirichlet` (the
message lists five densities, none over a simplex) and no `simo` class
(a bare `match.arg` message). Measured on the formula route,
`brms_monotonic`'s `fit1` sampled without that prior gives 14 divergent
transitions, `Rhat` 1.65, and a `zeta1_2` posterior mean near -2e6 with
`n_eff` 3. The fit route survives only because it starts at the ML
mode. `get_prior()` on an `mo()` model does not list a `simo` row, so
the gap is invisible from the discovery call.

A second modeling difference sits next to it: `mo(income) * age` makes
the interaction SHARE the main effect's simplex. brms fits two. The
shape of the monotonic effect therefore cannot vary with age, which is
what that vignette section is about.

### 3. `conditional_effects()` (rank 4)

Four separate things, diagnosed in the section above: the spurious
nonlinear guard under `method = "predict"`, the absent faceting,
`conditional_smooths()` missing, and `surface = TRUE` refused. Only the
first is a two-line fix.

### 4. No `plot()` method for `frmtmb_draws` (rank 11)

Two rows, but the worst message in the audit. `plot(draws)` falls
through to `plot.default` and says *"'x' is a list, but does not have
components 'x' and 'y'"*. It reads as a bug in the user's call.
`mcmc_plot()` is the replacement and nothing names it. A one-line
refusal would fix it. The `frmtmb_fit` methods that refuse
(`plot()` on a `frm_multiple`, `as_draws_array()`, `nchains()`) all do
this well already, so the pattern exists.

### 5. `hypothesis()` has no evidence ratio (rank 7)

The SYNTAX ports completely. brms's one-sided `"x > 0"` is accepted on
both paths, including on a plain `frmtmb_fit`, where it reports
`method = wald` with a one-sided `p` and a half-open interval, and the
print footer says which rows are one-sided. FN-5 is fixed further than
the earlier audit recorded.

What differs is only the REPORT. brms gives `Evid.Ratio`, `Post.Prob`
and `Star`; frmtmb gives a standard error, a `z` and a `p`. On the ML
path a Wald reading is the honest answer, because there is no
posterior. On the SAMPLE path the draws ARE there, so `Evid.Ratio` and
`Post.Prob` could be computed and are not. `plot(hypothesis_object)`
likewise still draws the estimate with an interval where brms shades a
posterior density, even when the draws exist. The one-sided `p` is also
bounded below by `1/ndraws`, which a short chain makes visible.

### 6. The `loo` family (rank 2)

The most FREQUENT edge and among the least severe, because the SAMPLE
path answers it completely. On a `frmtmb_fit`, `loo()`, `LOO()`,
`waic()` and `bayes_R2()` all fail with a bare dispatch message that
names neither `frm_sample()` nor `AIC()`. All ten vignettes hit it. The
cheap fix is not to implement anything: it is to make those dispatch
failures into refusals that name `frm_sample()`.

`add_criterion()` is the one item with no answer on either path. It is
how a brms user attaches ANY criterion to a fit.

### 7. Multivariate post-processing splits (rank 10)

A multivariate fit models well and reports badly. `pp_check()`,
`residuals()` and `fitted()` each refuse with *"is not supported yet
for multivariate fits"*, while `predict(resp = )` works, so the
accessors are split with no rule a reader can guess and `resp =` is
never consulted by the three that refuse. On draws,
`posterior_epred()`, `posterior_predict()` and `bayes_R2()` all return
the FIRST response only where brms returns every response; the second
is simply unreachable through those entry points, and nothing says so.

Worst of the group: **`rescor_matrix()` on a draws object returns
`NULL` in silence.** A silent `NULL` is the only edge in the audit that
neither refuses nor warns nor gives a wrong-looking number. On the ML
path the same fit prints the residual correlation in `summary()` under
a `Residual correlation:` heading, so the quantity exists and the draws
accessor simply drops it.

### 8. Two more silent or unhelpful failures worth a message

- `loo(a, b)` with two draws objects, which is the LITERAL vignette
  line in four vignettes, fails with `'list' object cannot be coerced
  to type 'integer'`. That is an internal message: it neither refuses
  nor names `loo_compare()`, which is the answer.
- `update()` has no `frmtmb_draws` method, so `update(draws, ...)`
  lands in `frm()`'s argument list and reports `unused argument
  (newdata = ...)`. The identical line succeeds on the ML path, which
  makes the failure read as a bug rather than a gap.

### 9. Prior interface (rank 1, owned elsewhere)

The largest single group and the one this audit deliberately does not
judge. See below.

## Follow-up, grouped by the work rather than by the edge

The ranked list above is grouped by CAUSE, which is the right shape for
understanding the surface and the wrong shape for filing tickets.
Several separate-looking edges collapse into one piece of work. Ordered
by vignette lines unblocked per unit of effort.

### 1. One message pass over the draws surface

The cheapest item by a wide margin, because it adds no feature. Five
edges, all on `frmtmb_draws`, all reached by a line the vignettes write
literally, all fixed the same way: a refusal that names the
replacement.

| call | today | should say |
|---|---|---|
| `rescor_matrix(draws)` | silent `NULL` | no draws method; read it from the fit's `summary()` |
| `loo(a, b)` | `'list' object cannot be coerced to type 'integer'` | one `loo()` per fit, then `loo_compare()` |
| `update(draws, ...)` | `unused argument (newdata = ...)` | no draws method; update the fit and re-sample |
| `plot(draws)` | `'x' is a list, but does not have components 'x' and 'y'` | use `mcmc_plot()` |
| `loo()`/`waic()`/`bayes_R2()` on a `frmtmb_fit` | bare dispatch failure | use `frm_sample()`, or `AIC()`/`BIC()` |

The silent `NULL` and the two internal coercion messages are the worst
of the audit: they neither refuse, nor warn, nor produce a
wrong-looking number, so a reader has nothing to act on. The refusal
machinery already exists and is used well elsewhere in the package
(`frm_multiple`'s `plot()`, `as_draws_array()` and `nchains()` all name
their workaround), so this is pattern-matching, not design.

### 2. The `conditional_effects()` nonlinear guard

Two lines, located and verified above. Unblocks the vignette's own
faceting display in two vignettes.

### 3. `summary()` presentation

The highest-severity group and a genuine piece of work, not a message
change. A natural-scale block for distributional parameters on the fit
path, and named natural-scale rows for variance components on the draws
path, so `theta_1` stops standing in for `sd(ult_Intercept)`.

### 4. `plot()` faceting for `cond__`

Honor `ncol` with `par(mfrow = ...)` around the existing loop. Not
brms's ggplot facet grid, but it removes the "one plot" surprise.

### 5. Features, in the order the vignettes ask for them

`add_criterion()` returning AIC/BIC; a `simo` prior class so `mo()`'s
simplex has a posterior; `conditional_smooths()`; the multivariate
accessors (`pp_check`, `residuals`, `fitted`) honoring `resp =`;
per-response `posterior_epred`, `posterior_predict` and `bayes_R2`;
`threshold =` for ordinal parameterizations.

## Which edges the priorcompat sibling already covers

Of the 49 prior-interface rows, these are the distinct items, and all
of them belong to the `wt-priorcompat` lane (which landed with
v0.43.0: `prior()` is exported, `nlpar`/`resp` work, a `brmsprior`
translates; the density/class gaps and `get_prior()` fidelity remain
its recorded follow-ups):

- `prior()`, brms's NSE constructor, is not exported at all.
- `set_prior()` parses five densities. No `dirichlet`, no `horseshoe`.
- `set_prior()` takes five classes. No `simo`, no `sigma` (which is
  `class = "Intercept", dpar = "sigma"`, and the message does not say
  so), and no reachable target for a `cs()` coefficient.
- `get_prior()` returns brms's columns with every row reading `(flat)`,
  plus a `theta` class brms has no name for.
- `prior_summary()` on a plain ML fit returns `NULL`.
- `priors =` penalizes the likelihood rather than defining a posterior,
  and the summary does not say the fit is penalized.
- `frm_sample()`'s default-prior banner omits the monotonic simplex and
  ordinal thresholds, which are exactly the parameters that then fail
  to mix.
- Under a brms attach, `set_prior()` and `get_prior()` silently resolve
  to brms.

**New, and NOT covered by that lane:** everything in ranks 3 through
11 above. In particular the `summary()` scale problem, the absent
`frmtmb_draws` `plot()` method, the `conditional_effects()` faceting,
the shared `mo()` simplex across an interaction, `add_criterion()`, the
multivariate accessor split, the silent `rescor_matrix()` `NULL`, the
two-model `loo()` failure, and the missing draws `update()`.

## What the v0.34 audit recorded that is now fixed

Measured, not assumed. Each was re-run at its original spelling.

| item | v0.34 | now |
|---|---|---|
| FN-1 gaussian default | 13 model calls failed | fixed |
| FN-2 `mm()` | missing function | fixed, with `weights =` and `mmc()` |
| FN-3 `conditional_effects()` on nlpar / `mo()` fits | no plottable predictors | fixed |
| FN-5 one-sided `hypothesis()` | refused | fixed, on `frmtmb_fit` as well as on draws |
| FN-6 bare family constructor | refused | fixed |
| FN-8 multi-parameter nlpar formula | refused | fixed |
| FN-9 `update()` argument names | three failures | fixed |
| FN-10 `lf()` | missing function | fixed, `nlf()` too |
| FN-11 `frm_multiple` post-processing | base-graphics error | now a refusal naming the workaround |
| FN-15 `bayes_R2()` | missing | present on draws, but a multivariate fit returns only the FIRST response where brms returns one row per response |

**Still open:** FN-4 (`surface = TRUE`, now a deliberate refusal), FN-7
(`threshold =`), FN-12 (`pp_check()` on multivariate, and `residuals()`
and `fitted()` too), FN-14 (`conditional_smooths()`), and
`add_criterion()`.

**Not in the earlier audit at all**, and found only because this one
ran both paths and read the output: the multivariate accessor split and
its one-response draws entry points, the silent `rescor_matrix()`
`NULL`, the two-model `loo()` call failing with an internal message, no
`update()` on draws, the `mo()` simplex having no prior and therefore
no posterior, and the whole `summary()` scale family.

## Estimate plausibility

Where a model ports, the ML point estimates sit on the vignettes'
posterior means. Spot checks against the printed vignette output:

| fit | brms | frmtmb ML |
|---|---|---|
| `fit_zinb1` | -1.01, 0.87, -1.36, 0.80 | matches to two decimals |
| `fit_loss` ult / omega / theta | 5273.7, 1.34, 46.07 | 5292.5, 1.337, 45.90 |
| `fit1` (nonlinear) b1 / b2 | 1.90, 0.75 | 2.079, 0.722 |
| `fit_ir3` eta / guess | 1.32, 0.31 | 1.005, 0.382 |
| `model_simple` (phylo) | 38.38, 5.17 | matches |
| `fit_mm` (multi-membership) | 19.0, sd 2.76 | matches |

The systematic divergence is in VARIANCE COMPONENTS at a boundary.
`fit_loss` sd(AY) is 598 against 745.74; `model_fisher` sd(phylo)
collapses to zero against 0.05; `kidney` sd(patient) is 3e-4 against
0.40. In every case brms's half-t prior holds the component away from
zero and maximum likelihood does not. Nothing is wrong, and a reader
comparing the two summaries will notice it first. The SAMPLE path with
the FORMULA route recovers brms's behavior, because that route supplies
the prior; the FIT route does not, and `frm_sample()` says so.



## Files

- `dev/brms-vignette-audit.md` (this file)
- `dev/brms-vignettes/_harness.R`: the `bv()` recorder
- `dev/brms-vignettes/_summarize.R`: the per-edge inventory
- `dev/brms-vignettes/_scoreboard.R`: generates the Headline and Per
  vignette tables in this file, so they cannot go stale
- `dev/brms-vignettes/<name>.R`: one translated vignette each
- `dev/brms-vignettes/brms-coexistence.R`: the one script that attaches
  brms on purpose
- `dev/brms-vignette-port.md`: the earlier mechanical audit, kept for
  the v0.34.0 baseline

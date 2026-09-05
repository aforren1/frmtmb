# Example pages still generating data by hand

The roxygen sweep is deferred. Seven core files were owned by other
lanes during the round that added the agreement tier
(`tests/testthat/test-simulate-density.R`), so no `@examples` block was
touched. This is the list to work from.

## What the sweep is

Replace a hand-written response with `frm_simulate()` on the model the
example is about. The covariates stay hand-written: `frm_simulate()`
draws a RESPONSE from a formula and parameters, so `x <- rnorm(n)` is
not what this is about.

The idiom, from `vignettes/inputs.Rmd`:

```r
set.seed(1)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)), y = 0)
dd$y <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                     newparams = list(Intercept = 0, x = 0.4, sigma = 1,
                                      sd_g__Intercept = 0.5),
                     nsim = 1, seed = 1001)[[1]]
```

Two rules learned the hard way while converting the test fixtures:

- **Give the draw its own seed.** Reusing the seed that made the
  covariates restarts the same random stream, so the residuals come
  back equal to `x` and the example fits a noiseless line. This
  produced convergence warnings in two files before it was caught.
- **Leave a page alone when the model cannot state its own generative
  process.** See "Do not convert" below.

The policy this serves is in `CONTRIBUTING.md`, under "Where test data
comes from".

## Do not convert

These are not oversights. Each one has a reason, and the reason is
worth keeping in the example.

| Page | Why |
| --- | --- |
| `R/families.R:113` (`frmtmb_family`) | The example DEFINES a custom family. A family being written for the first time has no `sim` slot, so `frm_simulate()` refuses it by name. Any example whose point is a user-written density is in this position. |
| `R/interop.R:25` (`check_custom_family`) | Same: the probe response feeds a family that does not exist yet. |
| `R/families.R:4375` (`frmtmb-families`) | Builds three responses on one data frame with `rt()`, `rnbinom()` and an `ifelse(runif(n) < 0.3, ...)` zero-inflation step. Three models, one frame; a conversion would need three `frm_simulate()` calls and would read worse than the thing it replaced. |
| `R/covstruct.R:225` (`frmtmb-student-re`) | Plants a deliberate outlying group (`b[20] <- b[20] + 6`). The outlier is the point, and it is not something the model draws. |
| `R/autocor.R:184` (`frmtmb-autocor`) | `stats::filter(rnorm(5), 0.6, "recursive")`. Convertible in principle: `ar()` needs its `thetaac` entry in the INTERNAL newparams spelling, since a correlation has no natural-scale name. Judge whether the internal spelling reads better than the recursion before converting. |

Two more gaps found while converting the vignettes, which block a
conversion anywhere they appear:

- **`mo()`** has no natural-scale name for its simplex (`zeta1`), so
  `frm_simulate()` refuses the natural spelling and the internal one
  would put raw simplex parameters in an example. This is why
  `vignettes/frmtmb.Rmd`'s `mo()` chunk was left alone.
- **`car()`** and the other structures without a `from_natural` map
  need the internal `theta` spelling. Legal, but opaque in an example.
  This is why the CAR chunk in the same vignette was left alone.

Closing either gap would unblock several pages at once, and is probably
better value than converting around them.

## The pages

The tables below hold 65 `@examples` pages, of which 64 need work.
The 65th, `R/simulate-new.R:447`, is listed for completeness and its
own row says there is nothing to do there. Files marked
**(owned)** were held by another lane and are the reason this list
exists rather than a diff.

### R/confint.R (owned) - 9 pages

| Line | Block |
| --- | --- |
| 340 | `confint.frmtmb_fit` |
| 647 | `confint_varcorr` |
| 914 | `diagnose` |
| 1276 | `anova.frmtmb_fit` |
| 1383 | `drop1.frmtmb_fit` |
| 1532 | `update.frmtmb_fit` |
| 1597 | `profile.frmtmb_fit` |
| 2183 | `hypothesis` |
| 2231 | `variables` |

### R/methods-fit.R - 7 pages

| Line | Block |
| --- | --- |
| 259 | `rescor_matrix` |
| 361 | `vcov.frmtmb_fit` |
| 459 | `coef.frmtmb_fit` |
| 535 | `fixef` |
| 574 | `ranef` |
| 706 | `VarCorr` |
| 804 | `expose_functions` |

### R/predict.R (owned) - 4 pages

| Line | Block |
| --- | --- |
| 991 | `predict.frmtmb_fit` |
| 1677 | `fitted.frmtmb_fit` |
| 2201 | `residuals.frmtmb_fit` |
| 2562 | `simulate.frmtmb_fit` |

### R/sugar.R - 4 pages

| Line | Block |
| --- | --- |
| 35 | `sigma.frmtmb_fit` |
| 116 | `ngrps` |
| 157 | `prior_summary` |
| 203 | `refit` |

### R/families.R (partly owned) - 6 pages

| Line | Block | Note |
| --- | --- | --- |
| 113 | `frmtmb_family` | do not convert (see above) |
| 2724 | `mixture` | convert; the group-level draw is a good demonstration |
| 3083 | `mixture_probs` | convert |
| 3317 | `mixture_mvn` | convert |
| 4057 | `cox_baseline` | `cox()` declares no simulator, so this page stays hand-written |
| 4375 | `frmtmb-families` | do not convert (see above) |

### R/conditional-effects.R (owned) - 3 pages

| Line | Block |
| --- | --- |
| 779 | `conditional_effects` |
| 1570 | `plot.frmtmb_fit` |
| 1628 | `pp_check` |

### R/draws-generics.R - 3 pages

| Line | Block |
| --- | --- |
| 30 | `posterior_summary` |
| 69 | `as_draws` |
| 115 | `draws-dimensions` |

### R/priors.R - 3 pages

| Line | Block |
| --- | --- |
| 202 | `set_prior` |
| 694 | `get_prior` |
| 1366 | `frmtmb-priors` |

### R/influence.R - 2 pages

| Line | Block |
| --- | --- |
| 25 | `influence.frmtmb_fit` |
| 156 | `plot.frmtmb_influence` |

### R/loo.R - 2 pages

| Line | Block |
| --- | --- |
| 47 | `loo` |
| 151 | `bayes_R2` |

### R/multiple.R - 2 pages

| Line | Block |
| --- | --- |
| 48 | `frm_multiple` |
| 226 | `anova.frmtmb_multiple` |

### R/sandwich.R - 2 pages

| Line | Block |
| --- | --- |
| 283 | `cluster_scores` |
| 428 | `vcov_cluster` |

### R/autocor.R - 2 pages

| Line | Block | Note |
| --- | --- | --- |
| 184 | `frmtmb-autocor` | see "Do not convert" |
| 853 | `autocor_matrix` | convert |

### R/structure.R - 2 pages

| Line | Block |
| --- | --- |
| 444 | `latent_probs` |
| 581 | `frmtmb-extension-api` |

### R/fit.R - 2 pages

| Line | Block | Note |
| --- | --- | --- |
| 402 | `frm` | the block also carries a `\dontrun{}` `sleepstudy` example at line 395; only the 402 one generates data |
| 1354 | `frmtmb_control` | |

### One page each

| File | Line | Block |
| --- | --- | --- |
| `R/allfit.R` | 65 | `frm_allfit` |
| `R/bf.R` | 406 | `mvbf` |
| `R/bootstrap.R` | 26 | `frm_bootstrap` |
| `R/covstruct.R` | 225 | `frmtmb-student-re` (do not convert) |
| `R/diagnostics.R` | 39 | `dharma_residuals` |
| `R/interop.R` | 25 | `check_custom_family` (do not convert) |
| `R/interop.R` | 336 | `getME.frmtmb_fit` |
| `R/parse.R` | 491 | `frmtmb-multimembership`; convert. The response draws in the natural spelling, `sd_mmschool1school2__Intercept`, on `(1 \| mm(school1, school2))`, which is the equal-membership structure the page builds by hand. The school factors and `share1` are covariates and stay hand-written. The page draws `y` once and fits it twice, so the weighted term needs no second draw; its own natural name is `sd_mmschool1school2weightscbindshare1share2__Intercept`, which is why the draw should use the unweighted spelling. |
| `R/par-template.R` | 107 | `par_template` |
| `R/sampling-api.R` | 155 | `frmtmb-sampling-api` |
| `R/simulate-new.R` | 447 | `frm_simulate` (its own page; the hand-written frame is only a dummy response column, so nothing to do) |
| `R/utils.R` | 225 | `num_factor` |

## Order to work in

1. The unowned files with many pages and one shared idiom:
   `R/methods-fit.R`, `R/sugar.R`, `R/draws-generics.R`,
   `R/priors.R`. Roughly 17 pages take the same three-line
   replacement.
2. `R/families.R`'s mixture pages, which are the ones where
   `frm_simulate()` shows something the hand-written version cannot:
   the class draw belongs to the group.
3. The owned files, once their lanes land: `R/confint.R`,
   `R/predict.R`, `R/conditional-effects.R`, `R/importance.R`,
   `R/frame.R`, `R/objective.R`, `R/compat.R`.
4. Re-run `R CMD check` after each file. Example code runs during
   check, so a page that stops converging is a check failure, not a
   test failure.

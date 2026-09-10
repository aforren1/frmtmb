# Latent class analysis

`lca(K)` fits the classic latent class measurement model, the one
'poLCA' fits (Linzer and Lewis 2011): `J` polytomous indicator items per
subject, one latent class `c` in `1..K` per subject, and the items
conditionally independent given the class,

## Usage

``` r
lca(K, ncat = NULL, na.rm = TRUE)
```

## Arguments

- K:

  Number of latent classes (at least 2).

- ncat:

  Number of categories per item: a vector of length `J`, a single value
  for equally-sized items, or `NULL` (the default) to infer each item's
  count as its largest observed code, which is what poLCA does.

- na.rm:

  If `TRUE` (default) subjects with any missing item are dropped by
  `na.action`; if `FALSE` they are kept and each missing item's term is
  masked out of that subject's likelihood.

## Value

A `frmtmb_family`.

## Details

\$\$P(y_i) = \sum_c w_c \prod_j \pi\_{j,c,y\_{ij}}.\$\$

The response is a MATRIX, one row per subject and one column per item,
holding integer category codes `1..C_j`. Write it as
`cbind(item1, item2, ...)` on the left of the formula, or attach a
matrix column to the data (`dd$Y <- data.matrix(dd[items])`) and name
it. Items may have different numbers of categories.

The class-membership weights are the `theta1 ... theta{K-1}` dpars,
multinomial logit against class `K`, and the main model formula applies
to every one of them. So `bf(cbind(a, b, c) ~ 1)` is the plain
measurement model and `bf(cbind(a, b, c) ~ age + educ)` is poLCA's
latent class regression: covariates on class membership come from the
ordinary linear-predictor machinery, with the usual `fixef()`,
[`confint()`](https://rdrr.io/r/stats/confint.html) and `hypothesis()`
on top. Individual gating predictors are overridable as
`bf(Y ~ x, theta2 ~ 1)` (all but `theta1`, which the main formula owns).

The item profiles are family extra parameters, held as
reference-category logits, one vector `pi<j>` per item (item `j`'s
`K * (C_j - 1)` free logits, class-major), and reported on the
probability scale by
[`lca_profiles()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_profiles.md).
They take no linear predictor - a covariate acts on class membership,
never on an item's conditional response probability, which is what makes
the classes interpretable as measurement.

## Labeling and starting values

A latent class likelihood is invariant to relabeling the classes and is
genuinely multimodal, so the answer depends on where the optimizer
starts. poLCA uses random restarts (`nrep`). The starting values here
are deterministic instead: subjects are clustered by their RESPONSE
PATTERN, and each cluster's smoothed empirical category proportions,
shrunk halfway toward the pooled proportions, become one class's
starting profile. The clusters are then ordered by the mean of each
subject's item codes, so class 1 is the low-score end where that score
carries information, and re-running the same data gives the same
labeling either way.

**What this replaced, and why.** Until 0.2.2 the rule was to score each
subject by the mean of its item codes and cut the scores into `K`
equal-count slices. A scalar score is blind to any design whose classes
differ in WHICH items they endorse rather than in how many, and such
designs are the ordinary case. On the four-class, ten-item design in the
realistic-scale table of `dev/extension-gaps-plan.md`, where every class
endorses three items strongly and seven weakly, the score explains
**0.08 percent** of its own variance between the true classes and the
four slices are four samples of the same mixture. The fit then reached a
local optimum **243 to 284 log-likelihood units** below
`poLCA(nrep = 10)` on **8 of 200** replicate data sets, and seven of
those eight passed every convergence test this project has: positive
definite Hessian, relative gradient between 1.9e-09 and 2.7e-08. Only
one of the eight would have shown a user anything at all.

**How the rule now shipped was tuned, and what it scored on data it was
not tuned on.** Its one free constant is the shrink weight, and it was
chosen on the SAME 200 replicates the failure above was found on. A
constant chosen that way has to be scored somewhere else before it means
anything, so it was, on a second block of 200 replicates of the same
design drawn from disjoint seeds:


                          tuned on         scored on
                          20260910-        20270401-
                          20261109         20270600
      the 0.2.2 score cut   8 lost           not run
      the rule now shipped  0 lost           0 lost

Both blocks are 200 data sets, scored against the best log-likelihood
anything reached on each one, and on none of the 400 did any start beat
`poLCA(nrep = 10)`. The scripts are `dev/latent-2p4-shrink.R`,
`-shrink2.R` and `-oos.R`, and `dev/latent-2p4-weight-summarize.R` reads
them without pooling the two blocks.

The shrink is not decoration. Clustering with NO shrink fixes all eight
data sets the score cut lost and breaks six it had won, by 79 to 303
units, and 7 more in the second block; a hard partition is a more
confident start than any fitted profile, and that confidence is its own
basin. At the other end `w = 1` is the label-symmetry axis and loses 199
of the first block and 200 of the second, by up to 1558 units. Both
edges of the usable interval are measured in `lca_init_extras()`'s own
notes.

**None of that makes the surface unimodal.** It is multimodal for
everyone: on those same eight data sets poLCA's own single-start EM
landed below its best of ten on 0 to 4 of 10 seeds. But the range
contains its own counterexample and the honest reading is narrower than
"the surface is the cause". On 2 of the 8 (seeds 20261010 and 20261013)
all TEN poLCA single starts found the global optimum while the old
deterministic start landed 274 to 278 units below every one of them, and
no poLCA start ever visited the mode it found. On those two the starting
rule and not the surface is what failed.

So compare starts before reading a solution, exactly as
`poLCA(..., nrep = 10)` does. Perturb the item parameters and refit:

    p0 <- fit$frame$par_template[fit$frame$extra_names]
    refits <- replicate(10, simplify = FALSE,
      frm(bf(Y ~ 1), family = lca(K = 3), data = dd,
          start = lapply(p0, function(v) v + rnorm(length(v)))))
    sapply(refits, logLik)

and keep the best. On the eight data sets above, run against the OLD
start, ten such refits recovered poLCA's optimum on 8 of 8 and took 6 to
17 seconds. Perturbing the OPTIMUM instead does not: two standard errors
around the local optimum, ten refits, still left 72 to 273
log-likelihood units on the table. The thing to scatter is the item
parameters at the start, not the answer the optimizer already found.

**Reproducing a fit made before 0.3.0.** The old rule is gone from the
package, so a published result from 0.2.2 or earlier has to be refit
from that start explicitly. This is it: score each subject by the mean
item code, cut the ranks into `K` equal-count slices, and take each
slice's Laplace-smoothed category proportions with no shrink.

    old_start <- function(Y, K) {
      ncat <- apply(Y, 2, max, na.rm = TRUE)
      sc <- rowMeans(sweep(Y - 1, 2, pmax(ncat - 1, 1), "/"), na.rm = TRUE)
      sl <- pmin(1 + ((rank(sc, ties.method = "first") - 1) * K) %/%
                   nrow(Y), K)
      setNames(lapply(seq_len(ncol(Y)), function(j) {
        unlist(lapply(seq_len(K), function(k) {
          p <- tabulate(Y[sl == k, j], nbins = ncat[j]) + 1
          p <- p / sum(p)
          log(p[-1]) - log(p[1])
        }))
      }), paste0("pi", seq_len(ncol(Y))))
    }
    frm(bf(Y ~ 1), family = lca(K = 3), data = dd,
        start = old_start(dd$Y, 3))

Checked rather than asserted: on three seeds of the design below,
including two the old rule got wrong, that block reproduces the old
starting values to the last bit and the old log-likelihood to a relative
difference of exactly 0 (`dev/latent-2p4-oldstart.R`). It handles binary
and polytomous items alike but not missing responses, which the shipped
rule masks and this block does not.

[`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.html)
is a different check: it re-runs the four optimizers from the SAME
start, so it tests the optimizer, not the surface. Ordering the gating
intercepts through `lower` and `upper` is the way to pin the labeling
itself.

Asking for more classes than the data hold drives item probabilities to
0 and 1. The optimizer then reports singular convergence and the
standard errors come back `NaN`, which is the boundary showing through
rather than a fault; poLCA lands on the same solutions. Compare
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) across `K` and read
[`lca_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_probs.md)'s
entropy before committing to a `K`.

## Recovery at a realistic scale

Measured on TWO blocks of 200 replicate data sets at the design named in
the realistic-scale table of `dev/extension-gaps-plan.md`: `K = 4`,
`n = 2000`, ten binary items whose endorsement probabilities are 0.85
inside a class's own block and 0.2 outside it, and two covariates on
membership. Classes are matched to the truth, and to poLCA, by item
profile before anything is scored. The first block, seeds 20260910 to
20261109, is the one the starting rule's shrink weight was tuned on; the
second, 20270401 to 20270600, is not. Both are reported, and they are
never pooled. `dev/latent-2p4-recheck.R` is the script and
`dev/latent-findings.md` holds the full tables.


                                       tuned on      not tuned on
      reaches poLCA(nrep = 10)'s
        optimum                        200 of 200    200 of 200
      the two log-likelihoods agree
        to, relative                   1.5e-14       4.4e-11
      Wald coverage, 9 gating
        coefficients x 200             94.3 percent  93.4 percent
        Monte Carlo interval           93.3 to 95.4  92.3 to 94.6
      positive definite Hessian        200 of 200    200 of 200
      largest bias on a coefficient    0.012         0.032
      item probabilities, max error    0.098         0.105
      modal-class accuracy, median     89.7 percent  89.6 percent

**Two coefficients under-cover, and there is no remedy for it yet.** At
200 replicates the Monte Carlo half-width on one coefficient is about 3
points, so neither block on its own can separate a 91 percent row from
95. Pooled over all 400 (`dev/latent-2p4-pooled.R`), two of the nine
rows have a Wald interval that excludes 95 and every other row contains
it:


      coefficient      covered      rate    95 percent interval   se/sd
      class3:x2        367 / 400   91.75    89.05 to 94.45        0.940
      class4:x2        363 / 400   90.75    87.91 to 93.59        0.922
      overall         3380 / 3600  93.89    93.11 to 94.67

So it belongs to the DESIGN and not to one draw. It is **two of the
three slopes on the binary covariate, not all three**: `class2:x2` pools
to 94.75 percent and is fine. The shortfall is 3.3 and 4.3 points below
nominal pooled, and 4.5 and 6.0 on the second block alone.

The mechanism is measured, and it rules out the obvious remedy. Writing
`z` for the error over the reported standard error, `sd(z)` is 1.06 and
1.08 on these two, the robust `IQR(z) / 1.349` agrees with it to within
6 percent, and the kurtosis is 2.82 and 3.10. There is no heavy tail and
no curvature: the whole distribution is inflated, which is to say the
standard error the fit reports is uniformly 6 to 8 percent below the
spread the estimator actually has. A profile interval fixes the SHAPE of
a likelihood, so it does not address this, and on the subset where it
can even be computed it does not: 103 of 118 both ways, wider on 118 of
118 by a median ratio of 1.001, and not one interval changed a verdict.
(It reaches only the 59 replicates of 200 where the scored quantity is a
raw parameter rather than a contrast; `hypothesis(method = "profile")`
with a `lincomb` is the only profile route to the rest.)

**So: on a binary gating covariate at this size, treat the Wald interval
as 3 to 4 points narrow, and know that nothing here has been measured to
fix it.** A bootstrap interval estimates the spread instead of reading
it off the Hessian and is the obvious candidate, but it has not been
measured on this design and is not recommended on that basis.

A fit is only as good as the mode it found. Read "Labeling and starting
values" above before quoting any of this.

## Missing item responses

With `na.rm = TRUE` (the default, poLCA's default) a subject with any
missing item is dropped by the usual `na.action`, and the message naming
the dropped rows is the ordinary one. With `na.rm = FALSE` the subject
is kept and only the missing item's factor leaves that subject's
likelihood, which is poLCA's `na.rm = FALSE` behavior; the mask is data,
so the tape never sees a branch. A subject missing EVERY item
contributes a constant and is better dropped.

## What is and is not supported

Post-fit,
[`lca_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_probs.md)
gives posterior class-membership probabilities per subject (with the
relative-entropy classification diagnostic attached),
[`lca_profiles()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_profiles.md)
gives the item profile table, and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws a class per
subject and then its items.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
`predict(type = "response")` and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) are refused: the
response is a matrix of nominal codes, so an "expected item code" would
be an average of arbitrary labels. Read
[`lca_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_probs.md)
and
[`lca_profiles()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_profiles.md)
instead. [`predict()`](https://rdrr.io/r/stats/predict.html) itself
returns the gating linear predictor (`theta1` by default, any `theta`
with `dpar =`), on `newdata` as well.

The gating coefficients are ordinary fixed effects, so `fixef()`,
[`confint()`](https://rdrr.io/r/stats/confint.html) (Wald, profile and
uniroot), `hypothesis()`, `set_prior()`, `lower`/`upper` bounds and
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
all work on them; [`anova()`](https://rdrr.io/r/stats/anova.html)
compares nested gating formulas at one `K`.

Refused in this version: random effects and smooths anywhere in the
model (latent classes plus continuous random effects is the
growth-mixture shape - use `mixture(..., groups = ~g)` for that), `REML`
and `frmtmb_control(profile = TRUE)` (the classes are exchangeable, so
there is no single inner mode), `quadrature`, every addition term
([`weights()`](https://rdrr.io/r/stats/weights.html), `cens()`,
[`trunc()`](https://rdrr.io/r/base/Round.html), `se()`, `mi()`,
`trials()`),
[`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.html), and
`residuals(type = "osa")`.

## References

Linzer, D. A. and Lewis, J. B. (2011). poLCA: An R Package for
Polytomous Variable Latent Class Analysis. *Journal of Statistical
Software*, 42(10), 1-29.

## See also

[`lca_probs()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_probs.md),
[`lca_profiles()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca_profiles.md),
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
for a mixture of continuous responses,
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.html)
for model-based clustering of a numeric matrix.

## Examples

``` r
set.seed(1)
# four binary items measuring two well-separated classes
n <- 300
cl <- rbinom(n, 1, 0.4) + 1
pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
Y <- matrix(0L, n, 4)
for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(x = rnorm(n))
dd$Y <- Y

fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)
lca_profiles(fit)
#> <lca profiles> 2 classes, 4 items
#> 
#> Estimated class sizes (mean prior probability):
#> class1 class2 
#> 0.2863 0.7137 
#> 
#> item1:
#>          cat1   cat2
#> class1 0.9661 0.0339
#> class2 0.1869 0.8131
#> 
#> item2:
#>          cat1   cat2
#> class1 0.8423 0.1577
#> class2 0.2413 0.7587
#> 
#> item3:
#>          cat1   cat2
#> class1 0.7982 0.2018
#> class2 0.2777 0.7223
#> 
#> item4:
#>          cat1   cat2
#> class1 0.9520 0.0480
#> class2 0.1459 0.8541
head(lca_probs(fit))
#>            class1      class2
#> [1,] 0.0093358936 0.990664106
#> [2,] 0.0093358936 0.990664106
#> [3,] 0.0062911864 0.993708814
#> [4,] 0.9926850273 0.007314973
#> [5,] 0.0009150042 0.999084996
#> [6,] 0.0611608386 0.938839161

# latent class regression: covariates gate class membership
frm(bf(Y ~ x), family = lca(K = 2), data = dd)
#> frmtmb fit: Y ~ x 
#> Family: lca(K = 2)   Method: ML 
#>  Links: theta1 = identity
#> 
#> logLik: -697.815  AIC: 1415.63  nobs: 300 
#> 
#> Fixed effects:
#>  theta1:
#> (Intercept)           x 
#>    -0.91700     0.04477 
```

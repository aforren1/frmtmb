# Drift-diffusion models of two-choice response times

## What the model says

A participant sees a letter string and decides whether it is a word. The
drift-diffusion model says that evidence accumulates noisily from a
starting point until it reaches one of two boundaries, and the response
time is when it gets there. Four quantities describe that:

- `mu`, the **drift rate**: how fast evidence arrives, and in which
  direction. This is the one that usually carries the experimental
  effect.
- `bs`, the **boundary separation**: how much evidence is required. Wide
  boundaries mean slow, accurate responses.
- `ndt`, the **non-decision time**: encoding the stimulus and moving the
  finger, which no amount of deciding accounts for.
- `bias`, the **relative start point**: where accumulation begins
  between the boundaries. 0.5 is unbiased.

The names are brms’s, so a model written for
[`brms::wiener()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
reads the same here.

## Simulating an experiment

[`ddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/ddm_simulate.md)
draws response times and boundary choices jointly, the way an experiment
produces them. Here, words and non-words differ in drift rate and in
nothing else, which is the standard claim about a lexical decision task.

``` r

set.seed(20240903)
n <- 1200
word <- rep(c(0, 1), each = n / 2)
dat <- ddm_simulate(n,
                    mu = 0.3 + 1.0 * word,   # words drift faster
                    bs = 1.5,
                    ndt = 0.30,
                    bias = 0.5)
dat$lex <- factor(word, labels = c("nonword", "word"))
head(dat)
#>          rt upper     lex
#> 1 1.3754639     1 nonword
#> 2 1.1195175     1 nonword
#> 3 0.5259899     1 nonword
#> 4 0.8249955     1 nonword
#> 5 0.6424872     0 nonword
#> 6 0.5434753     1 nonword
```

`upper` is the decision: 1 for a response at the upper boundary, 0 for
the lower one.

``` r

with(dat, tapply(rt, list(lex, upper), mean))
#>                 0         1
#> nonword 0.8567651 0.8680013
#> word    0.7121997 0.7176207
with(dat, tapply(rt, lex, length))
#> nonword    word 
#>     600     600
prop.table(table(dat$lex, dat$upper), margin = 1)
#>          
#>                   0         1
#>   nonword 0.4033333 0.5966667
#>   word    0.1366667 0.8633333
```

## The decision indicator

brms spells the boundary a trial ended at as `rt | dec(decision)`.
frmtmb’s addition terms are a closed set and `dec()` is not one of them,
so this family reads the indicator from `vint()` instead:

``` r

# brms
brm(rt | dec(decision) ~ lex, family = wiener(), data = dat)

# here
frm(bf(rt | vint(upper) ~ lex), family = wiener(), data = dat)
```

`upper` must already be coded 0/1. Unlike `dec()`, `vint()` will not
take a factor, and passing anything else is refused rather than silently
coerced.

## Fitting

Drift depends on lexicality; everything else is constant across trials.
`bias = 0.5` fixes the start point on the response scale, which is the
usual choice when nothing in the design would bias it.

``` r

fit <- frm(bf(rt | vint(upper) ~ lex, bias = 0.5),
           family = wiener(), data = dat)
summary(fit)
#> Family: wiener 
#> Formula: rt | vint(upper) ~ lex 
#> Method: ML   nobs: 1200 
#> logLik: -799.08  AIC: 1606.16  BIC: 1626.52 
#> 
#> Coefficients (mu):
#>             Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 0.256576   0.054434  4.7135 2.435e-06
#> lexword     1.046222   0.084275 12.4144 < 2.2e-16
#> 
#> Coefficients (bs):
#>             Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 0.404003   0.014818  27.265 < 2.2e-16
#> 
#> Coefficients (ndt):
#>             Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)  2.04596    0.11509  17.777 < 2.2e-16
#> 
#> Fixed dpar: bias = 0.5
```

Read the coefficients on their own scales. Drift is on the identity
scale, so the `lexword` coefficient is the change in drift rate
directly. Boundary separation is on the log scale. The non-decision time
is on a logit scaled onto `(0, min(rt))`, which is what keeps it inside
the support:

``` r

e <- unlist(fixef(fit))
c(drift_nonword = e[["mu.(Intercept)"]],
  drift_word_gain = e[["mu.lexword"]],
  boundary = exp(e[["bs.(Intercept)"]]),
  ndt = min(dat$rt) / (1 + exp(-e[["ndt.(Intercept)"]])))
#>   drift_nonword drift_word_gain        boundary             ndt 
#>       0.2565761       1.0462219       1.4978084       0.2991607
```

Against the values the data were generated from – drift 0.3 and 1.3,
boundary 1.5, non-decision time 0.30 – that is the recovery this design
supports at 1200 trials.

## Why the non-decision time needs a bounded link

The density is exactly zero for a response time at or below `ndt`. A log
link would let the optimizer propose a non-decision time above the
fastest observed response, where the likelihood is not merely small but
undefined, and it would find out only by stepping there.

So `ndt` gets a logit scaled onto `(0, max_ndt)`, with `max_ndt`
defaulting to the smallest response time in the data. The constraint
becomes structural: no finite linear predictor can violate it.

``` r

lk <- wiener(max_ndt = 0.4)$links[["ndt"]]
lk$linkinv(c(-10, -2, 0, 2, 10))   # always strictly inside (0, 0.4)
#> [1] 1.815915e-05 4.768117e-02 2.000000e-01 3.523188e-01 3.999818e-01
```

Pass `max_ndt` explicitly when you intend to
[`predict()`](https://rdrr.io/r/stats/predict.html) on new data whose
fastest response differs from the training set’s.

## Predictions

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
`predict(type = "response")` give the mean response time, conditional on
the boundary that row ended at, in closed form.

``` r

nd <- data.frame(lex = factor(c("nonword", "word")), upper = c(1, 1))
predict(fit, newdata = nd, type = "response")
#> [1] 0.8532160 0.7309745
```

`vint()` is required on `newdata` as well; there is no default boundary,
because the conditional mean genuinely depends on it.

One caution that looks like a bug and is not. With `bias` fixed at 0.5,
the conditional mean decision time is the **same** at both boundaries,
for any drift rate. That is a theorem about the diffusion, not an
oversight: an unbiased start gives fast and slow responses the same
average time whichever way they go. Free the bias and the two separate:

``` r

set.seed(7)
d2 <- ddm_simulate(900, mu = 0.7, bs = 1.4, ndt = 0.25, bias = 0.35)
f2 <- frm(bf(rt | vint(upper) ~ 1), family = wiener(), data = d2)
c(pred_upper = predict(f2, newdata = data.frame(upper = 1),
                       type = "response"),
  pred_lower = predict(f2, newdata = data.frame(upper = 0),
                       type = "response"),
  obs_upper = mean(d2$rt[d2$upper == 1]),
  obs_lower = mean(d2$rt[d2$upper == 0]))
#> pred_upper pred_lower  obs_upper  obs_lower 
#>  0.7923043  0.6105395  0.8003306  0.6084389
```

## Simulating from a fit

[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws response
times at each row’s own boundary, by inverse transform through the
first-passage quantile function. Comparing those draws with the observed
times is the usual posterior-predictive sanity check.

``` r

set.seed(3)
sims <- as.matrix(simulate(fit, nsim = 100))
c(observed = mean(dat$rt), simulated = mean(sims))
#>  observed simulated 
#> 0.7901746 0.7911800
quantile(dat$rt, c(0.1, 0.5, 0.9))
#>       10%       50%       90% 
#> 0.4413758 0.6824099 1.2794110
quantile(sims,   c(0.1, 0.5, 0.9))
#>       10%       50%       90% 
#> 0.4334025 0.6719358 1.3045922
```

## What this family does not do

- **`cens()` and [`trunc()`](https://rdrr.io/r/base/Round.html) are
  refused.** They need a distribution function, and the family declares
  none. The Wiener first-passage CDF is a third series with its own
  truncation problem, and it is not written here.
- **`residuals(type = "pearson")` and `"deviance"` are refused.** No
  variance function and no unit deviance. `type = "response"` works.
- **The decision is data, not a modelled outcome.** This family models
  the response time given the boundary. It does not jointly model which
  boundary was reached, so it is not a substitute for fitting choice and
  time together.
- **There is no across-trial variability, and no contaminant mixture.**
  Drift, start point and non-decision time are constant within a
  condition, which makes this the pure Wiener model rather than the full
  Ratcliff model. The Ratcliff likelihood integrates over a drift
  distribution at every row, and fast guesses and slow lapses need a
  mixture component that is not here. A random effect over subjects
  models variation BETWEEN participants, which is a different claim from
  variation between trials.

`frm_compat("wiener")` states each of these, along with what was tested
and what merely was not.

``` r

frm_compat("wiener")[, c("feature_b", "status")]
#>             feature_b      status
#> 1                  us       works
#> 2                diag       works
#> 3             homdiag       works
#> 4                  cs       works
#> 5                 ar1 conditional
#> 6              hetar1 conditional
#> 7                  ou conditional
#> 8                toep       works
#> 9             homtoep       works
#> 10              homcs       works
#> 11                exp conditional
#> 12                gau conditional
#> 13                mat conditional
#> 14                 rr conditional
#> 15            equalto conditional
#> 16             gr_cov conditional
#> 17            gr_prec conditional
#> 18             smooth conditional
#> 19                 gp conditional
#> 20               hsgp conditional
#> 21                car conditional
#> 22               spde conditional
#> 23               us_t       works
#> 24             diag_t       works
#> 25          weights()       works
#> 26           trials()    untested
#> 27             cens()     refused
#> 28            trunc()     refused
#> 29               se()     refused
#> 30               mi()     refused
#> 31             vint()       works
#> 32            vreal()    untested
#> 33                s()       works
#> 34               t2()       works
#> 35               mo() conditional
#> 36          mi_pred()       works
#> 37          gp_pred()       works
#> 38          cs_pred()     refused
#> 39               ar()     refused
#> 40               ma()     refused
#> 41             arma()     refused
#> 42             cosy()     refused
#> 43            unstr()     refused
#> 44               REML    untested
#> 45         quadrature       works
#> 46            profile       works
#> 47          autoscale       works
#> 48           sparse_x       works
#> 49              prior       works
#> 50             bounds       works
#> 51            verbose       works
#> 52               mvbf       works
#> 53             rescor     refused
#> 54               |ID|    untested
#> 55                 nl     refused
#> 56            mixture    untested
#> 57        mixture_mvn    untested
#> 58             fitted       works
#> 59            predict       works
#> 60           simulate       works
#> 61          residuals conditional
#> 62      residuals_osa    untested
#> 63            emmeans conditional
#> 64    confint_profile conditional
#> 65 hypothesis_profile conditional
#> 66       bar_crossing     refused
#> 67         call_group       works
#> 68         double_bar       works
#> 69               mm()       works
#> 70              mmc() conditional
```

## The density

The likelihood is the Wiener first-passage time density of Navarro and
Fuss (2009). They give two series for it: one that converges quickly at
small normalized times `u = t / a^2`, one that converges quickly at
large ones, and a criterion for choosing between them per evaluation.

An automatic-differentiation tape cannot make that choice. `u` depends
on the boundary separation, which is a parameter, so a comparison on `u`
is a branch on a parameter – and RTMB advectors carry no comparison
operators, so it does not even tape. Branching on the response alone
would tape and would be wrong, because the same response time falls on
either side of the criterion depending on the boundary separation.

The two BOUNDARIES are the case where that rule works in your favor.
Only the lower-boundary density is written, because reflecting the
process maps an upper crossing onto a lower one with `bias -> 1 - bias`
and `mu -> -mu`, and the reflection is driven by the `vint()` column,
which is data. Arithmetic on a data column tapes.

This package evaluates **both** series at a fixed truncation and blends
their logs with a logistic weight in `log(u)`: Navarro and Fuss’s
criterion with the step smoothed. The weight saturates to exactly 0 and
exactly 1 outside a narrow band, and inside that band both series are
accurate, so any convex combination of them is accurate too.

The cost is that every evaluation pays for both series. What it buys is
a density that holds over the whole range of normalized times rather
than the part a single truncation happens to cover:

``` r

u <- c(0.005, 0.05, 0.4, 2, 8, 30)
ref <- vapply(u, function(uu) {
  log(RWiener::dwiener(uu + 1e-12, 1, 1e-12, 0.45, 1, resp = "lower"))
}, 0)
mine <- vapply(u, function(uu) {
  frmtmb.ddm:::ddm_lpdf_lower(uu, 1, 1, 0.45)
}, 0)
small_only <- vapply(u, function(uu) {
  frmtmb.ddm:::ddm_log_gs(uu, 0.45) - 1 * 1 * 0.45 - 1 * 1 * uu / 2
}, 0)
data.frame(u,
           blended_rel_err = abs(mine - ref) / abs(ref),
           small_time_only_rel_err = abs(small_only - ref) / abs(ref))
#>       u blended_rel_err small_time_only_rel_err
#> 1 5e-03    1.227404e-16            0.000000e+00
#> 2 5e-02    1.005083e-15            1.005083e-15
#> 3 4e-01    0.000000e+00            1.490329e-16
#> 4 2e+00    1.743704e-16            6.800445e-15
#> 5 8e+00    0.000000e+00            5.552216e-03
#> 6 3e+01    0.000000e+00            8.180276e-01
```

The last column is what a fixed truncation of the small-time series
alone would give. It is fine until it is not.

## Reference

Navarro, D. J. and Fuss, I. G. (2009). Fast and accurate calculations
for first-passage times in Wiener diffusion models. *Journal of
Mathematical Psychology*, 53(4), 222-230.

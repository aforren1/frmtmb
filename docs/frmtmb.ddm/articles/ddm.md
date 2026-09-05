# Response time models: two-choice and multi-alternative

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

The boundary a trial ended at is data, and it reaches the density as an
addition term spelled the way brms spells it:

``` r

# brms
brm(rt | dec(decision) ~ lex, family = wiener(), data = dat)

# here
frm(bf(rt | dec(decision) ~ lex), family = wiener(), data = dat)
```

`dec()` takes what brms takes: a factor or a character vector, read on
its levels with the second one as the upper boundary, or a 0/1 column.
The package contributes the term to frmtmb’s addition-term registry when
it loads, so nothing has to be enabled.

`vint()` carries the same thing as a plain 0/1 integer and also works.
It was the only route before frmtmb had a registry to contribute to, and
models written that way keep running:

``` r

frm(bf(rt | vint(upper) ~ lex), family = wiener(), data = dat)
```

Supplying neither is refused. An absent addition term reaches a density
as `NULL`, `NULL` in arithmetic is a zero-length result, and the fit
would otherwise return with a log likelihood summed over no rows.

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
- **The start-point range is not held inside the boundaries by
  construction.** `wiener(variability = "sz")` gives `sz` a logit link,
  so the width is below 1 and the uniform start point stays inside the
  boundaries whenever `bias` is 0.5, which is the usual case. At a
  strongly biased start a wide `sz` can push the range past a boundary.
  Nor is `ndt - st / 2` held above zero. Both are joint constraints on
  two distributional parameters, and a link is a property of one, so
  neither can be made structural from outside frmtmb.

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
#> 45         quadrature     refused
#> 46         importance       works
#> 47            profile       works
#> 48          autoscale       works
#> 49           sparse_x       works
#> 50              prior       works
#> 51             bounds       works
#> 52            verbose       works
#> 53               mvbf       works
#> 54             rescor     refused
#> 55               |ID|    untested
#> 56                 nl     refused
#> 57            mixture       works
#> 58        mixture_mvn    untested
#> 59             fitted       works
#> 60            predict       works
#> 61           simulate       works
#> 62          residuals conditional
#> 63      residuals_osa    untested
#> 64            emmeans conditional
#> 65    confint_profile conditional
#> 66 hypothesis_profile conditional
#> 67       bar_crossing     refused
#> 68         call_group       works
#> 69         double_bar       works
#> 70               mm()       works
#> 71              mmc() conditional
#> 72              dec()       works
```

## The full model: across-trial variability

Everything above holds the four parameters fixed within a condition.
Ratcliff’s full diffusion model does not: it draws the drift rate, the
start point and the non-decision time afresh on every trial, which is
what lets one model account for both the fast errors a start-point range
produces and the slow errors a drift range produces.

`variability` turns those on. Each one it names becomes an ordinary
distributional parameter and takes its own formula like any other:

``` r

set.seed(4)
full <- ddm_simulate(600, mu = 1.2, bs = 1.5, ndt = 0.30,
                     sv = 1.0, st = 0.10)
ffit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
            family = wiener(variability = c("sv", "st")), data = full)
e <- unlist(fixef(ffit))
ub <- min(full$rt)
c(mu  = e[["mu.(Intercept)"]],
  bs  = exp(e[["bs.(Intercept)"]]),
  ndt = ub / (1 + exp(-e[["ndt.(Intercept)"]])),
  sv  = exp(e[["sv.(Intercept)"]]),
  st  = 2 * ub / (1 + exp(-e[["st.(Intercept)"]])))
#>        mu        bs       ndt        sv        st 
#> 1.3188861 1.5183809 0.3159405 1.2936135 0.1380081
```

### Three integrals, three answers

The likelihood is the analytic Wiener density averaged over those
distributions. They are not averaged the same way, because they are not
the same integral.

**The drift rate has a closed form.** The drift enters the log density
only as `-v a w - v^2 t / 2`: the series that carries the shape of the
distribution is the driftless one, and the drift multiplies it by an
exponential-quadratic. Averaging an exponential-quadratic against a
normal completes the square, so

    mean of exp(-v a w - v^2 t / 2)
      = (1 + t sv^2)^(-1/2)
        exp( ((a w)^2 sv^2 - 2 v a w - v^2 t) / (2 (1 + t sv^2)) )

That is exact. There is no node count for `sv` and no accuracy to trade,
and estimating it is free relative to the plain density: over three
seeds a fit with `sv` took 0.97 to 1.07 times a fit without it, and at
the density level both cost 0.0002 s per call. The two quadratures below
are the ones that cost.

**The other two are uniform, and get fixed Gauss-Legendre nodes.** The
node positions and counts are decided when the family object is built.
Nothing may branch on a parameter on an automatic-differentiation tape,
and a node count that moved with `st` would be exactly that; a parameter
only rescales the interval the fixed nodes are mapped onto.

The two defaults differ because the two integrands do. The start-point
integrand is analytic across its whole range. The non-decision-time
integrand is not free to use its whole range: the density is zero for a
non-decision time at or past the response time, so on a fast trial the
range is cut and the integrand turns on sharply at the cut. Measured
against a far finer rule on such a trial:

``` r

lp <- function(nsz, nst) {
  nd <- frmtmb.ddm:::ddm_nodes(c("sz", "st"), c(sz = nsz, st = nst))
  frmtmb.ddm:::ddm_lpdf_var(0.36, 1.0, 1.5, 0.45, 0.30, 1.0, 0.2, 0.2,
                            0, nd, TRUE, 1e-9)
}
best <- lp(31L, 81L)
n <- c(3, 5, 7, 11, 15, 21, 31)
data.frame(nodes = n,
           sz_err = abs(vapply(n, function(k) lp(k, 81L), 0) - best),
           st_err = abs(vapply(n, function(k) lp(31L, k), 0) - best))
#>   nodes       sz_err       st_err
#> 1     3 8.800603e-07 1.190146e-02
#> 2     5 7.849277e-13 1.666151e-03
#> 3     7 1.776357e-15 1.247713e-04
#> 4    11 1.776357e-15 1.746504e-05
#> 5    15 1.998401e-15 3.019358e-07
#> 6    21 1.554312e-15 3.011576e-09
#> 7    31 0.000000e+00 1.085518e-10
```

The start-point column is at machine precision by 7 nodes. The
non-decision-time column takes 21 to reach 1e-9, which is where the
default sits. Lower it with `wiener(nodes = c(st = 11))` if the fit is
slow and the extra digits are not worth the time.

### A contaminant mixture

Fast guesses are not diffusion trials, and a model with no component for
them pays by dragging the non-decision time down. The standard treatment
is a mixture, which needs the Wiener component’s non-decision time to
sit above the fastest response time. That is refused by default, for a
bare model correctly, and `allow_unreachable = TRUE` says the other
component covers those trials:

``` r

frm(bf(rt | dec(decision) ~ 1, bias1 = 0.5),
    family = mixture(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                     lognormal()),
    data = dat)
```

A trial below the non-decision time then gets a log density that
exponentiates to exactly zero, which is the right likelihood for that
component, and differentiates to exactly zero, which a true `-Inf` would
not: `-Inf` gives the correct value and a `NaN` gradient, and one `NaN`
stops the fit.

### Why this is not `quadrature = TRUE`

frmtmb has a `quadrature` argument, and it is not this. It marginalizes
RANDOM EFFECTS by adaptive Gauss-Kronrod, is wired to the random-effect
coefficient vector by name, and refuses a model that has no
random-effect block at all. Across-trial variability shares nothing
between trials, has no level to estimate, and gives every row its own
integral whether or not the model has a grouping factor. It belongs
inside the density, and that is where it is.

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

## More than two alternatives

Everything above is one accumulator running between two absorbing
boundaries. Two is not a simplification there, it is the whole geometry:
the process is one-dimensional and there are two ends, so a
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)
fit can describe a yes/no decision and nothing wider. No choice of
parameters reaches a third response.

The linear ballistic accumulator gets there by changing what races. Each
alternative gets its own accumulator, rising in a straight line from a
start point drawn uniformly on `(0, A)` to a common threshold
`b = A + k`. The rate is drawn once per trial from a normal distribution
rather than varying within the trial, and it is that missing
within-trial noise which buys the closed form: the outcome depends on
one draw per accumulator, so the likelihood is the winner’s density
times the others’ survival probabilities, for any number of
accumulators.

The chunks below run whether or not RWiener is installed, since the
accumulator race needs neither it nor any other reference package.

### An experiment with three responses

A cue makes one of three alternatives more attractive, trial by trial,
and leaves the other two alone. That is a design neither two-choice
family can express, so it is worth building from scratch.
[`lba_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/lba_simulate.md)
takes a matrix of drift means with one row per trial, which is how a
covariate enters the generative model.

``` r

set.seed(20)
N <- 1500
cue <- rnorm(N)
# the cue raises the evidence for alternative 2 and nothing else
V <- cbind(2.4, 1.5 + 1.1 * cue, 0.9)
dat3 <- lba_simulate(N, v = V, A = 0.5, k = 0.4, ndt = 0.2)
dat3$cue <- cue
head(dat3)
#>   choice        rt        cue
#> 1      2 0.4528505  1.1626853
#> 2      1 0.3538470 -0.5859245
#> 3      2 0.3788586  1.7854650
#> 4      1 0.6926590 -1.3325937
#> 5      1 0.4368551 -0.4465668
#> 6      2 0.3145631  0.5696061
table(dat3$choice)
#> 
#>   1   2   3 
#> 904 433 163
```

The choice is data, exactly as the boundary indicator was, and reaches
the family through `vint()`. It names the accumulator that responded,
counting from one.

``` r

fit3 <- frm(bf(rt | vint(choice) ~ cue), family = lba(3), data = dat3)
fixef(fit3)
#> $v1
#> (Intercept)         cue 
#>  2.33038268  0.01772203 
#> 
#> $v2
#> (Intercept)         cue 
#>    1.464786    1.039147 
#> 
#> $v3
#> (Intercept)         cue 
#>   0.7283536  -0.1080327 
#> 
#> $A
#> (Intercept) 
#>    -0.87533 
#> 
#> $k
#> (Intercept) 
#>  -0.8251791 
#> 
#> $ndt
#> (Intercept) 
#>   0.8213659
```

Every drift is a primary parameter, so the main formula went to all
three and each got its own coefficients. That is what makes the answer
readable: the cue slope is near its true 1.1 on `v2` and near zero on
`v1` and `v3`, which is the model recovering that the cue moved one
alternative and not the others.

``` r

vapply(c("v1", "v2", "v3"), function(p) fixef(fit3)[[p]][["cue"]],
       numeric(1))
#>          v1          v2          v3 
#>  0.01772203  1.03914697 -0.10803273
```

Had we known which alternative the cue acts on, we would have said so
and left the other two drifts an intercept each:

``` r

fit3b <- frm(bf(rt | vint(choice) ~ 1, v2 ~ cue), family = lba(3),
             data = dat3)
c(all_drifts = AIC(fit3), v2_only = AIC(fit3b))
#> all_drifts    v2_only 
#>  -742.0252  -743.2697
```

The remaining parameters came back near the values they were drawn from.
`A` and `k` carry log links and `ndt` the bounded link the Wiener family
also uses, so read them through the link:

``` r

e <- fixef(fit3b)
lk <- family(fit3b)$links$ndt
c(A = exp(e$A[[1]]), k = exp(e$k[[1]]), ndt = lk$linkinv(e$ndt[[1]]),
  truth_A = 0.5, truth_k = 0.4, truth_ndt = 0.2)
#>         A         k       ndt   truth_A   truth_k truth_ndt 
#> 0.4095901 0.4447096 0.1906341 0.5000000 0.4000000 0.2000000
```

### Two conventions to know before comparing packages

The model is identified only up to scale: multiply `A`, the threshold,
every drift mean and the drift standard deviation by one constant and
the distribution of `(choice, rt)` does not move. Something has to be
held fixed, and this family fixes the drift standard deviation. It is
the `sd_v` argument rather than a hidden constant, and it is carried on
the family object, so a fit says what it held:

``` r

family(fit3b)[["lba_sd_v"]]
#> [1] 1 1 1
```

The second convention is what happens to an accumulator whose drift
comes out negative, since it would never reach the threshold. Following
`rtdists`, drift rates are truncated at zero by default, which makes
every accumulator arrive eventually and the choice probabilities sum to
one. `lba(3, posdrift = FALSE)` leaves them untruncated, and the
response distribution is then defective by design. The two are different
models rather than a rescaling of each other, so a log-likelihood from
one is not comparable with a log-likelihood from the other, nor with
another package set the other way.

### What it refuses

A choice outside `1..n` is a data error, not a small likelihood, and is
named as such rather than quietly fitted:

``` r

bad <- dat3
bad$choice[7] <- 4L
frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = bad)
#> Error:
#> ! lba(3): the vint() choice indicator names which accumulator responded and must be a whole number from 1 to 3. Saw 4. A factor is not accepted; recode it with as.integer(factor(choice)) and check the level order matches the accumulator numbering.
```

So is omitting the choice altogether, which the family declares through
`required_aterms` so that frmtmb refuses it by name instead of building
a log likelihood over no rows:

``` r

frm(bf(rt ~ 1), family = lba(3), data = dat3)
#> Error:
#> ! lba: the density needs `vint1`, which nothing on this response supplies. Write the addition term: rt | vint(<column>) ~ ...
```

A threshold inside the start-point range needs no refusal at all. The
threshold is `A + k` with a log link on `k`, so it stays above the
start-point range at every value of the linear predictor, and a trial
that begins already finished is not a state the optimizer can reach.

## References

Brown, S. D. and Heathcote, A. (2008). The simplest complete model of
choice response time: Linear ballistic accumulation. *Cognitive
Psychology*, 57(3), 153-178.

Donkin, C., Brown, S. D. and Heathcote, A. (2009). The overconstraint of
response time models: Rethinking the scaling problem. *Psychonomic
Bulletin & Review*, 16(6), 1129-1135.

Navarro, D. J. and Fuss, I. G. (2009). Fast and accurate calculations
for first-passage times in Wiener diffusion models. *Journal of
Mathematical Psychology*, 53(4), 222-230.

Ratcliff, R. and Tuerlinckx, F. (2002). Estimating parameters of the
diffusion model: approaches to dealing with contaminant reaction times
and parameter variability. *Psychonomic Bulletin & Review*, 9(3),
438-481.

# Exact distributional parameters for a custom density

The `dpars` list a
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
log-density is handed carries each distributional parameter on its own
response scale. That is what a density can usually use directly, and it
is not enough in the tail. An inverse link SATURATES:
`stats::plogis(eta)` is exactly one in double precision from
`eta = 36.7368005696771` upward, so a density that forms `1 - mu` by
subtraction gets exactly zero there and returns `NaN` for the value AND
for the gradient, where the truth is an ordinary large negative number.
Below that point the subtraction is not fatal, only wrong: the
complement it forms is off by 1.0e-3 relative at `eta = 30`.

## Usage

``` r
dpar_log(dpars, dpar, link)

dpar_log1m(dpars, dpar, link)

dpar_log_complement(dpars, dpar, link)

dpar_complement(dpars, dpar, link)
```

## Arguments

- dpars:

  The named list of distributional parameters the density was handed.
  Read it by name, never by position.

- dpar:

  The name of the distributional parameter to read, as `dpars` names it:
  `"mu"`, `"zi"`, `"shape"`.

- link:

  That dpar's OWN link, as the link name your family gave in `links = `
  (`"logit"`) or as a link object. It is not optional: the linear
  predictor is on the link scale, so nothing can be recovered from it
  without knowing which link put it there.

## Value

`dpar_log()` and `dpar_log1m()` a numeric or advector vector. The other
two a list of two such vectors: `l` and `l1m` for
`dpar_log_complement()`, which are what the two one-sided functions
return, and `p` and `q` for `dpar_complement()`.

## Details

The linear predictor never saturated.
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) stores it
beside each distributional parameter while the objective is taped, and
these four accessors recover from it the quantity the density needs.
Write a density over them instead of over `log(mu)`, `log(1 - mu)` and
`1 - mu`.

- `dpar_log()`:

  `log(x)`, for a dpar that must stay positive (a `shape`, a `phi`, a
  `sigma`) and for one on the unit interval alike.

- `dpar_log1m()`:

  `log(1 - x)` for a dpar on the unit interval: a `zi` or `hu` gate, a
  probability, a coherence. This is the term that dies first, and it
  stays finite and exactly differentiable at a value the optimizer has
  pushed against 1.

- `dpar_log_complement()`:

  Both of the above at once, as `l` and `l1m`, from ONE log odds, for a
  density that needs the pair. A mixture gate does. Take a one-sided
  function when you need one term. The pair records a second
  [`RTMB::logspace_add()`](https://rdrr.io/pkg/RTMB/man/Distributions.html)
  over the whole response and the tape REPLAYS it on every gradient
  evaluation: at 1000 rows the accessor tapes 14002 AD nodes against
  10002, and one gradient sweep of the `cross_wishart()` density costs
  633 us against 480 us. What you are not paying for is the R call,
  which happens once per fit, when the tape is built.

- `dpar_complement()`:

  `x` as `p` and `1 - x` as `q`, for a density that needs the natural
  scale. You call it for `q`; `p` comes back with it because one
  log-odds gives both, and a pair taken from one source cannot disagree
  with itself.

## How far this reaches

The relative error of the plain arithmetic on a logit gate, against what
these accessors return, at the linear predictors an optimizer can visit.
One process, the log-scale form as the reference:

|                  |                       |                            |
|------------------|-----------------------|----------------------------|
| **eta**          | **`1 - plogis(eta)`** | **`log(1 - plogis(eta))`** |
| 20               | 3.6e-08               | 1.8e-09                    |
| 30               | 1.0e-03               | 3.4e-05                    |
| 36               | 4.3e-02               | 1.2e-03                    |
| 36.7368005696771 | 1 (exactly 0)         | 1.9e-02                    |
| 40               | 1 (exactly 0)         | Inf, the log is -Inf       |
| 700              | 1 (exactly 0)         | Inf, the log is -Inf       |

The plain form does not fail suddenly. It loses digits from `eta = 20`,
is down to three at 30 and to one at 36, and only then becomes `-Inf`
with a `NaN` gradient. A fit that walks out there converges to a floored
likelihood without a word.

## On the tape and off it

The `.eta_` entries exist only while the objective is taped. Your
density reads them exactly once, when the tape is built; what runs on
every gradient sweep after that is the arithmetic that read RECORDED,
not your R code. Anywhere else the dpar values arrive alone, every
accessor falls back to the plain arithmetic, and that is correct because
nothing outside the tape is differentiated.

Anywhere else is a smaller place than it sounds.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) do not evaluate the
density at all: they read `post$mean_fn`, the link, the optimizer's own
value and `sim`. What can reach these accessors off the tape is a
family's OWN numeric helpers, `post$dev_fn` and `post$var_fn` among
them, and
[`check_custom_family()`](https://aforren1.github.io/frmtmb/reference/check_custom_family.md).
Counted on a fitted `cross_wishart()`, whose deviance helper is the one
that calls them, fourteen post-fit entry points reach the arithmetic
exactly once between them, through `residuals(type = "deviance")`.

**Test your density on both paths.** A family author who exercises only
the off-tape path never runs the branch these accessors exist for, and a
density that is wrong on the tape still fits: it converges to the wrong
place, or dies at `NA/NaN gradient evaluation` with nothing naming the
cause. Build the on-tape list by hand to test it, as the example below
does, and put the value under `.eta_<dpar>` yourself for that one
purpose.

## The `.eta_` entries are not the API

`.eta_<dpar>` is a reserved name, and reading it directly is not
supported even though you can see it in `dpars`. It holds the linear
predictor on the LINK scale, so what it means depends on the dpar's
link, and a density that assumes one link reads a different quantity the
moment the family gains a `link_<dpar>` argument. `-log1p(exp(e))` is
`log(1 - x)` on a logit and is nothing at all on an identity link; `e`
is `log(x)` on a log link and is not on a softplus. The accessors take
the link and branch on what it can supply, which is why they take it as
an argument and why it is not optional. Whether an entry is present is
also a property of the phase, not of the model, and folding that branch
in is most of what these do.

## If you need the raw log odds

There is no accessor for it, and there does not need to be:
`dpar_log(d, p, lk) - dpar_log1m(d, p, lk)` reconstructs it. Measured on
a logit at `eta` in `{-700, -40, 0, 40, 700}` the absolute error is 0,
and it peaks at 5.0e-17 near `eta = 0` where the two terms nearly
cancel. A density term linear in the natural parameter cares about that
absolute error, not the relative one.

## See also

[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
for the density these serve,
[frmtmb-links](https://aforren1.github.io/frmtmb/reference/frmtmb-links.md)
for which links carry an exact form and which fall back, and
[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
for the accessors a family uses AFTER a fit, at the estimates

## Examples

``` r
# a gate at plogis(40), which is exactly 1 in double precision
on_tape <- list(zi = stats::plogis(40), .eta_zi = 40)
off_tape <- list(zi = stats::plogis(40))

# the value a density needs, and what the subtraction gives instead
dpar_log1m(on_tape, "zi", "logit")
#> [1] -40
log(1 - off_tape$zi)
#> [1] -Inf

# off the tape the entry is absent and the plain form comes back
dpar_log1m(off_tape, "zi", "logit")
#> [1] -Inf

# both terms from one log odds, for a density that needs the pair
unlist(dpar_log_complement(on_tape, "zi", "logit"))
#>             l           l1m 
#> -4.248354e-18 -4.000000e+01 

# the log of a positive dpar, through that dpar's own link
dpar_log(list(shape = exp(3), .eta_shape = 3), "shape", "log")
#> [1] 3
dpar_log(list(shape = log1p(exp(3)), .eta_shape = 3), "shape",
         "softplus")
#> [1] 1.114678

# a density written over the pair rather than over 1 - mu
d <- list(mu = stats::plogis(40), .eta_mu = 40)
mp <- dpar_complement(d, "mu", "logit")
c(p = mp$p, q = mp$q)
#>            p            q 
#> 1.000000e+00 4.248354e-18 
```

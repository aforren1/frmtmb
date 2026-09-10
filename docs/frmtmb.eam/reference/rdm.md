# The racing diffusion model

A race between `n` independent diffusions, for choices with more than
two alternatives. Each accumulator is a Wiener process with its own
positive drift rate and unit diffusion, starting from a point drawn
uniformly on `(0, A)` and running to a common threshold `b`. The first
accumulator to reach the threshold is the response, and the observed
time is its arrival time plus a non-decision time.

## Usage

``` r
rdm(n, max_ndt = NULL)
```

## Arguments

- n:

  Number of accumulators, so the number of response alternatives. At
  least 2.

- max_ndt:

  Upper bound for the non-decision time, in the units of the response,
  applied to every row. `NULL`, the default, takes the fastest response
  of each row's `ndt_group()`, or of the whole data set when the model
  has no `ndt_group()`. It cannot be combined with `ndt_group()`. Give
  it when a component of a
  [`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
  needs the bound up front.

## Value

A `frmtmb_family`.

## Details

The model is
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)'s
geometry with the ballistic assumption removed: same start-point range,
same threshold, but each accumulator is noisy WITHIN a trial rather than
drawing one rate at the start of it. One accumulator's arrival time is
therefore an inverse Gaussian rather than the reciprocal of a normal,
and the whole race stays closed form.

Use it rather than
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
when a design has more than two response alternatives, and rather than
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)
when the story you want for trial-to-trial spread is within-trial noise
rather than a rate drawn afresh each trial.

## The model

Write the decision time as \\t = rt - ndt\\. One accumulator with drift
\\v\\ travels a distance \\b - z\\ from its start point \\z\\, and \\z\\
is uniform on \\(0, A)\\, so the DISTANCE is uniform on \\(k, k + A)\\
with \\k = b - A\\. For a fixed distance \\c\\ the first-passage time is
inverse Gaussian with density \\c\\\phi((c - vt)/\sqrt{t})\\t^{-3/2}\\,
and averaging that over the distance integrates in closed form. Put

\$\$x_1 = \frac{k - v t}{\sqrt t}, \qquad x_2 = \frac{k + A - v t}{\sqrt
t},\$\$

the standardized positions of the two ends of the distance range at time
\\t\\. Then that accumulator's defective density and survival are

\$\$f(t) = \frac{v\\(\Phi(x_2) - \Phi(x_1)) - (\phi(x_2) -
\phi(x_1))/\sqrt t}{A},\$\$ \$\$S(t) = \frac{\sqrt t\\(x_2\\\Delta\Phi +
\Delta\phi) + A\\\Phi(x_1) - (\Delta E + \Delta\Phi)/(2 v)}{A},\$\$

where \\\Delta\\ is the change of a quantity across the distance range
and \\E(c) = e^{2 v c}\\\Phi(-(v t + c)/\sqrt t)\\ is the reflected term
of the inverse-Gaussian distribution function. A trial on which
accumulator \\j\\ responded at time \\rt\\ contributes

\$\$\log f_j(t) + \sum\_{i \neq j} \log S_i(t),\$\$

the same race
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)
runs, over a different single-accumulator law.

## Parameters

- `v1`, ..., `vn`:

  Drift rates, one per accumulator. LOG link, which is where this family
  parts company with
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md):
  an LBA drift is the MEAN of a normal drawn once per trial and may be
  negative, but a racing-diffusion drift is the rate itself, and an
  accumulator with a rate of zero or less never reaches the threshold at
  all. These are the primary parameters, so the main formula goes to all
  of them and each gets its own coefficients. Give one its own formula
  to move it alone.

- `A`:

  Upper end of the start-point range, so start points are uniform on
  `(0, A)`. Log link.

- `k`:

  Distance from the top of the start-point range to the threshold, so
  that `b = A + k`. Log link.

- `ndt`:

  Non-decision time. Bounded link; see below.

The threshold is `A + k` rather than a free `b`, for the reason
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)
gives: a log link on `k` makes `b > A` structural, so a threshold inside
the start-point range, where a fraction of trials would begin already
finished, is not a state the optimizer can reach.

`A` divides the density, so it cannot be zero. The log link keeps it
positive at every value of its linear predictor, and a race with no
start-point variability is a limit this family approaches rather than a
model it fits.

## Mapping to other software and to the paper

The same model is written three ways. This family's spelling is
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)'s,
because the geometry is the same and moving between the two families
should cost one word.

|             |         |                          |
|-------------|---------|--------------------------|
| this family | EMC2    | Tillman et al. (2020)    |
| `v1`..`vn`  | `v`     | drift rate               |
| `A`         | `A`     | start-point range        |
| `k`         | `B`     | threshold less the range |
| `A + k`     | `B + A` | threshold                |
| `ndt`       | `t0`    | non-decision time        |
| fixed at 1  | `s`     | diffusion coefficient    |

EMC2's `B` is this family's `k` and NOT the threshold: its own sampler
draws the distance to travel as `B + runif(1, 0, A)`, so `B` is the gap
above the start-point range, exactly as `k` is here. Reading it as the
threshold moves every threshold by `A`.

The diffusion coefficient is fixed at one, which is what identifies the
scale: multiplying `A`, `k` and every drift by one constant, and the
diffusion coefficient with them, leaves the distribution of
`(choice, rt)` unchanged. EMC2 divides `A`, `B` and `v` by its `s` and
so fixes the same quantity in the same place.

## Non-decision time, and what its coefficients mean

The density is zero at and below `ndt`, so the likelihood has a hard
edge at the fastest response and a log link would let the optimizer walk
over it.

Without `ndt_group()` the bound is one number, the fastest response in
the data or `max_ndt`, and it stays IN the link: `ndt` is a time, as it
has always been. Write `rt | vint(choice) + ndt_group(subject) ~ ...`
when `ndt` carries a subject deviation, and the bound becomes that
subject's own fastest response. A per-row bound cannot live in a link,
so under a grouping `ndt` is a FRACTION of the row's bound and the
density multiplies: `predict(dpar = "ndt", type = "response")` then
reports the fraction and
[`ndt_time()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_time.md)
reports the time either way. A single bound is the global fastest
response, and a subject whose non-decision time is above it cannot be
represented at any value of the random effect;
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
has the measurement.

## The response, and why not `dec()`

A trial is a `(choice, time)` pair. The time is the response and the
choice is per-row data, which reaches the family through `vint()`:

    frm(bf(rt | vint(choice) ~ cond), family = rdm(3), data = dat)

`vint1` is the ACCUMULATOR INDEX: a whole number in `1..n` naming which
accumulator reached the threshold, counting from one. This is
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)'s
rule, including that a factor is not accepted, so recode it with
`as.integer(factor(choice))` and check that the level order matches the
accumulator numbering.

`dec()` does not apply, for
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)'s
reason: that term carries a 0/1 indicator naming one of two boundaries,
and a race of `n` accumulators needs `1..n`.

## Censoring and truncation

A race has a survivor function in closed form, so `cens()` and
[`trunc()`](https://rdrr.io/r/base/Round.html) both work.

The reason is that the race is over as soon as ANY accumulator finishes,
so the probability that the response time is past `t` is the probability
that none of them has:

\$\$S(t) = \prod_i S_i(t),\$\$

one factor per accumulator, each the start-point integral the likelihood
already forms for the losers of an observed trial. The family declares
that product as its log survivor function and its complement as its
distribution function, so nothing new is derived for censoring: it is
the same `lsurv` the density's loser terms use, multiplied over all `n`
accumulators instead of over `n - 1` of them.

    frm(bf(rt | vint(choice) + cens(censored) ~ cond), family = rdm(3),
        data = dat)

A censored trial has no winner to report, because the race had not
finished when the clock ran out. `vint()` is still required, since a
declaration cannot be conditional on a censoring code, so give such a
row any accumulator index: the likelihood does not read it. The
distribution function is written as `-expm1(log S)`, which keeps its
digits where `1 - S` would lose them, and the log survivor goes to
`frmtmb` on the LOG scale, so a right-censored row stays exact past the
point where `log(1 - F)` is a constant with a zero gradient.

## Accuracy

Every piece is written in the form that keeps its digits rather than the
form the paper prints, and the reasons are
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)'s
reasons.

The normal difference \\\Phi(x_2) - \Phi(x_1)\\ goes through the same
log-space helper the linear ballistic accumulator uses, because the
subtractive spelling returns exactly zero once the lower tail rounds to
one. The density difference \\\phi(x_2) - \phi(x_1)\\ and the
boundary-term difference \\\Delta E\\ both go through `ddm_expdiff()`,
which anchors at the larger of the two logs: written directly, the first
underflows one term to zero while the other is still finite, and the
second overflows.

The survival is written out rather than as `1 - F`. Measured against a
Gauss-Legendre reference that never subtracts, on a grid of 648 points,
the form here holds 1.3e-12 relative down to a survival of 1e-3, 6.4e-11
down to 1e-15 and 4.1e-8 down to 1e-250, and returns no exact zeros.
`1 - EMC2:::pWald()` returns exactly zero on 35 of the 315 points of a
comparable grid, the first at a survival near 1e-13, which would send a
loser's log contribution to `-Inf` on an ordinary row.

The density agrees with `EMC2:::dWald()` to 2.1e-4 relative and with the
same non-subtracting reference to 1.6e-13. The gap is the subtractive
normal difference in EMC2's form, not a disagreement about the model.

Below all of this the density underflows in double precision; the log
density is floored rather than returning `-Inf`, so the optimizer sees a
finite wall instead of a hole. Decision times at or below zero, which
[`predict()`](https://rdrr.io/r/stats/predict.html) on faster new data
can reach, are floored the same way.

## References

Tillman, G., Van Zandt, T. and Logan, G. D. (2020). Sequential sampling
models without random between-trial variability: The racing diffusion
model of speeded decision making. *Psychonomic Bulletin and Review*,
27(5), 911-936.

## See also

[`rdm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm_simulate.md)
to generate from the model,
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md)
for the ballistic race with the same geometry, and
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
for the two-choice diffusion.

## Examples

``` r
set.seed(1)
dat <- rdm_simulate(400, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5,
                    ndt = 0.2)
fit <- frm(bf(rt | vint(choice) ~ 1), family = rdm(3), data = dat)
fixef(fit)
#> $v1
#> (Intercept) 
#>    1.054443 
#> 
#> $v2
#> (Intercept) 
#>   0.5458444 
#> 
#> $v3
#> (Intercept) 
#> -0.09667655 
#> 
#> $A
#> (Intercept) 
#>  -0.8922818 
#> 
#> $k
#> (Intercept) 
#>  -0.7532663 
#> 
#> $ndt
#> (Intercept) 
#>    2.203433 
#> 
```

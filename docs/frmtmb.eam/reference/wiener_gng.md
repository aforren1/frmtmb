# The go/no-go diffusion model

A two-boundary diffusion in which only ONE boundary is observed. The
evidence accumulator runs exactly as in
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
but the experiment records a response only when it reaches the upper
boundary, and only if it does so before a deadline. Reaching the lower
boundary and never reaching either produce the same observation: no
response.

## Usage

``` r
wiener_gng(
  deadline = NULL,
  max_ndt = NULL,
  variability = character(0),
  nodes = c(sz = 7L, st = 21L),
  nogo_nodes = c(sv = 15L, sz = 7L, st = 7L)
)
```

## Arguments

- deadline:

  The response deadline, in the units of the response. One positive
  number when every trial shares it. `NULL`, the default, takes it per
  row from `vreal()`.

- max_ndt:

  Upper bound for the non-decision time. `NULL`, the default, takes it
  from the fastest go response.

- variability:

  Which across-trial variability parameters to estimate: any of `"sv"`
  (drift rate), `"sz"` (start point) and `"st"` (non-decision time),
  exactly as
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  takes them. The default estimates none.

- nodes:

  Gauss-Legendre node counts for the GO branch, which is
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)'s
  density: the same names, the same meaning and the same defaults, so a
  go trial is scored identically by the two families.

- nogo_nodes:

  Quadrature node counts for the NO-GO probability, which is a different
  integral over the same three distributions. `sz` and `st` are
  Gauss-Legendre counts, and `st` defaults lower than the density's
  because the probability's range is not cut by the response time; `sv`
  is a Gauss-Hermite count with no counterpart in
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
  because there the drift integral is closed form. Only the entries for
  the `variability` parameters in use are read. Raise `sv` for a wide
  drift variability; see Across-trial variability for the measured
  table.

## Value

A `frmtmb_family`.

## Details

This is the design Gomez, Ratcliff and Perea (2007) fitted, and the
family EMC2 calls `DDMGNG`. It is the right model whenever one response
alternative has no overt response: go/no-go, lexical decision with a
single key, signal detection with a single report.

Fitting such data with
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
is not an option, because
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
needs to know which boundary each trial ended at and half the trials
here do not say.

## The likelihood

A trial contributes one of two things, and which one is data.

A GO trial, a response at the upper boundary at time `rt`, contributes
the ordinary defective first-passage density of
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
at the upper boundary. It is defective on purpose: it integrates to the
probability of a go response, not to one.

A NO-GO trial contributes the probability that the upper boundary was
not reached before the deadline,

\$\$P(\text{no go}) = 1 - F\_{\text{upper}}(TD),\$\$

which is the lower-boundary mass accumulated by the deadline plus the
mass still diffusing at it. The two add up: integrating the go density
from the non-decision time to the deadline gives
\\F\_{\text{upper}}(TD)\\, and adding the no-go probability gives
exactly one. Nothing has gone missing, and `test-rdm-gng.R` asserts it.

## Parameters

The same four
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
has, with the same names and the same links, so a model moves between
the two families by changing the family and nothing else.

- `mu`:

  Drift rate. Identity link, so it is signed: positive drift favours the
  go boundary.

- `bs`:

  Boundary separation. Log link.

- `ndt`:

  Non-decision time. Bounded link, taken from the GO trials only; see
  below.

- `bias`:

  Relative start point in `(0, 1)`. Logit link. 0.5 is unbiased.

The go boundary is the UPPER one, which is a convention rather than a
restriction: a model in which the observed response is the lower
boundary is this one with the drift negated and the bias reflected.

## The data

which trials, and by when: Two things are per-row data. Which trials
produced a response reaches the family through `dec()`, and the deadline
reaches it either as one number on the family or through `vreal()`:

    frm(bf(rt | dec(responded) ~ cond), family = wiener_gng(deadline = 1.5),
        data = dat)

    frm(bf(rt | dec(responded) + vreal(deadline) ~ cond),
        family = wiener_gng(), data = dat)

`dec()` here means RESPONDED, coded 1, or did not, coded 0. That is the
same 0/1 the term carries for
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
and it means the same thing, because with one observable boundary
"reached the upper boundary" and "responded" are the same event. It
takes a factor, a character vector or a logical the way brms does,
reading the second level as a response, so a plain logical column works
as it reads.

Give the deadline on the family when every trial had the same one, which
is the usual design, and through `vreal()` when it varied. Supplying
neither is refused by name.

A no-go trial has no response time. The response column still needs a
number, because it is the response and frmtmb will not carry an `NA`
through a rowwise family; the likelihood does not read it, and the
deadline is the natural thing to put there. Any positive finite value
gives the same answer.

A go trial whose response time is past its own deadline is refused
rather than fitted, because the model gives that trial no probability at
all: the deadline is what stopped the accumulator.

## Non-decision time

The bound is the fastest GO response, not the fastest row. A no-go row's
response entry is a placeholder, so letting it into the bound would let
a placeholder decide a parameter's range. A model with no go trials at
all is refused: nothing in it identifies the non-decision time.

## Across-trial variability

Ratcliff's `sv`, `sz` and `st` are here, named and linked exactly as
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
names and links them, so a model moves between the two families by
changing the family and nothing else:

    frm(bf(rt | dec(responded) ~ cond, bias = 0.5),
        family = wiener_gng(deadline = 1.5,
                            variability = c("sv", "sz", "st")),
        data = dat)

The GO branch is
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)'s
density and inherits all three unchanged: it calls the same averaged
density at the same nodes, so the two families agree to the last bit on
a go trial. The NO-GO branch is the part that had to be written, and it
is a different integral, because it averages the DISTRIBUTION FUNCTION
rather than the density.

Each of the three enters it its own way.

- `sv`:

  Quadrature, and the only one of the three that costs more here than in
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md).
  The drift enters the density only as `exp(-v a w - v^2 t / 2)`, an
  exponential-quadratic that a normal average integrates in closed form.
  It enters the distribution function through the eigenvalues as well,
  as `1 / (v^2 a^2 + k^2 pi^2)` in every term of the large-time series
  and as the gambler's-ruin probability that series corrects, and
  neither is an exponential-quadratic. So the no-go branch integrates
  the drift by Gauss-Hermite, and `nogo_nodes` carries an `sv` entry
  that
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
  has no use for.

- `sz`:

  The same Gauss-Legendre nodes the density uses, over the same uniform
  start point.

- `st`:

  Shifts the DEADLINE. The density's non-decision-time range is cut at
  the response time, because past that the decision time is negative and
  the density is zero; the no-go probability has no such cut, because a
  non-decision time past the deadline leaves the accumulator no time at
  all and the probability of no response is then exactly one. The whole
  range is averaged.

The identity that ties the two branches together holds WHILE THE
START-POINT RANGE STAYS INSIDE THE BOUNDARIES, which is where the model
is defined: the go density integrated to the deadline plus the no-go
probability is one, to 8.9e-16 in the plain family and to whatever the
quadrature gives once a variability parameter is on. With the default
`nogo_nodes` that is 3.9e-14 at `sv` = 0.6 and 2.0e-10 at `st` = 0.20,
but only 8.1e-05 at `sv` = 2.0 and 9.9e-08 at `st` = 0.45. It is not
"exactly" one, and the size of the gap is the `sv` row of the node table
below.

Outside the boundaries the two branches are no longer the same average
of the same pair, and the mass is not conserved. Only the NO-GO branch
clamps the start point; the go branch is left unclamped deliberately,
because clamping it would break the bit-identity with
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md).
Measured, at a deadline of 1.5:

|                           |          |          |          |
|---------------------------|----------|----------|----------|
|                           | go       | no-go    | total    |
| `sz` = 0.5, `bias` = 0.85 | 0.919017 | 0.051622 | 0.970639 |
| `sz` = 0.9, `bias` = 0.90 | 0.780184 | 0.063186 | 0.843370 |

so up to 16 percent of the pair's mass is missing there. What saves it
is that the region is strictly downhill. Profiling `sz` on 2000 trials
generated at `bias` = 0.85 and a true `sz` of 0.20, where the range
leaves the boundary at `sz` = 0.30, the best point inside is 1867.975 at
`sz` = 0.190 and the best point outside is 1678.640 at `sz` = 0.305: the
interior peak wins by 189 log units and the surface is monotone across
the crossing. The clamp is a barrier an optimizer walks away from, not
an attractor, so a fit does not end up there. Read a fitted `sz` whose
range crosses a boundary as a fit that has run out of model, not as an
estimate.

## How many nodes, and what they cost

Measured against a 200-bit reference, over six parameter settings
spanning deadlines 0.8 to 2.5, drifts 0.5 to 3, boundary separations 0.8
to 2.5 and start points 0.4 to 0.8. Worst relative error, one parameter
live at a time:

|       |            |            |     |             |            |
|-------|------------|------------|-----|-------------|------------|
| nodes | `sz` = 0.1 | `sz` = 0.3 |     | `st` = 0.05 | `st` = 0.3 |
| 3     | 4.5e-06    | 2.1e-03    |     | 3.4e-11     | 1.6e-06    |
| 5     | 1.8e-11    | 6.4e-07    |     | 8.6e-15     | 4.0e-10    |
| 7     | 5.1e-15    | 4.2e-11    |     | 6.3e-15     | 6.1e-14    |
| 9     | 4.0e-15    | 2.1e-15    |     | 4.2e-15     | 2.7e-15    |

Both saturate by seven nodes, and `st` saturates there where the
DENSITY's own `st` integral needs 21. That is not a discrepancy: the
density's range is cut at the response time and its integrand turns on
sharply at the cut, and the probability has no cut at all.

`sv` is the one that does not saturate:

|       |            |            |            |
|-------|------------|------------|------------|
| nodes | `sv` = 0.3 | `sv` = 0.8 | `sv` = 1.5 |
| 7     | 9.0e-07    | 8.9e-04    | 1.8e-02    |
| 11    | 9.5e-12    | 3.8e-06    | 1.0e-02    |
| 15    | 8.4e-15    | 1.3e-06    | 1.4e-03    |
| 21    | 9.5e-15    | 3.3e-09    | 3.5e-05    |
| 31    | 1.2e-14    | 1.2e-12    | 6.6e-07    |
| 41    | 1.4e-14    | 2.8e-14    | 9.1e-08    |

The node count a given accuracy needs rises with the product of the
boundary separation and `sv`, because that product sets how sharp the
transition in the drift is. The default of 15 holds 1e-14 at a narrow
drift variability and 1e-6 at a moderate one; **raise
`nogo_nodes = c(sv = 31)` or higher for a wide one**, and read the table
rather than assuming the default is enough.

The three do not compound. Measured jointly, the error of the full
three-dimensional rule tracks the `sv` error alone to within a factor of
three at every setting, so `sz` and `st` at seven nodes are there to not
be the binding term, and `sv` is the only knob worth turning.

The cost is a product, and it is the reason `nogo_nodes` exists as an
argument separate from `nodes`. On 500 rows with all three live:

|              |      |       |
|--------------|------|-------|
| `nogo_nodes` | grid | wall  |
| 7, 5, 5      | 175  | 124 s |
| 11, 5, 5     | 275  | 185 s |
| 11, 7, 7     | 539  | 610 s |

A model with one variability parameter pays one dimension of that and is
cheap; the three-at-once model is the expensive one, and it is expensive
because a distribution function with no closed form in any of the three
has to be evaluated on a product grid.

EMC2's `DDMGNG` carries the same three parameters, by calling a compiled
distribution function that integrates them numerically. The two agree:
over 108 grid points the no-go probability matches `1 - EMC2:::pDDM()`
to 9.4e-14 with no variability, 2.4e-13 under `sv` and 3.9e-13 under
`sz`, and the go density matches `EMC2:::dDDM()` to 1.2e-15.

`st` needs a shift before the two are comparable, and the difference is
a CONVENTION rather than a defect on either side. EMC2 takes `st0` as a
uniform on `[t0, t0 + st0]`; this family takes `st` centred on `ndt`,
following brms. Compare EMC2 at `t0` with this family at
`ndt = t0 + st0 / 2` and the agreement is 1.5e-13 for the probability
and 1.6e-15 for the density. Compare them without the shift and they
differ by 14 to 45 percent, in the density as much as in the
probability.

## What a go/no-go design can and cannot identify

`sv` is weakly identified here, and it is the design rather than the
likelihood. With only one boundary observed, drift variability trades
against the drift and the boundary separation along a ridge.

Measured on 3000 simulated trials with a true `sv` of 0.6, this family
returns 0.001 while
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
on the same generative parameters, seeing BOTH boundaries, returns
0.612. Profiled at the true values of everything else the likelihood
does peak at the truth, but the fit reaches a HIGHER log likelihood at
`sv` near zero by moving the drift from 1.20 to 1.09 and the separation
from 1.40 to 1.33. Nothing is wrong with the surface; the sample simply
prefers that corner of the ridge.

How much of that is the ONE observable boundary is less settled than the
contrast suggests, and the honest statement is narrower. A grid search
on two-boundary data from the same generative parameters, with `ndt` and
`bias` held at the truth, finds the same corner there too: `sv` near
zero beats the truth by 2.4 log units, against 1.8 on the go/no-go data.
So the ridge is a property of the drift-diffusion likelihood and not
only of this design; what the go/no-go design does is remove the
boundary-proportion information that would otherwise help pin `sv`, and
the evidence for that is the pair of fits above rather than a study.
Either way the practical advice is the same.

So read a small fitted `sv` here as "the data did not pin it" rather
than as "there is no drift variability", and prefer to fix it, or to
estimate it from a two-choice condition of the same experiment, over
reading it off a go/no-go block. `sz` and `st` are better behaved,
because both change the SHAPE of the go response time distribution
rather than trading against its location.

## Censoring

A trial whose clock was stopped before it responded is a RIGHT-CENSORED
observation, and the probability of it is the no-go probability at the
time the clock stopped. This family declares that as its log survivor
function, so `cens()` works:

    frm(bf(rt | dec(responded) + cens(stopped) ~ cond),
        family = wiener_gng(deadline = 1.5), data = dat)

Left censoring, interval censoring and
[`trunc()`](https://rdrr.io/r/base/Round.html) are refused, and not for
want of an integral. The likelihood here is a defective density plus a
point mass at "no response", and a truncation window on the response
scale renormalizes the density while saying nothing about the mass, so
the two halves of every row would be divided by different things. Right
censoring has no such problem, because it replaces a whole row rather
than reweighting it. So this family declares an `lccdf` and,
deliberately, no `lcdf`.

The refusal you will see is frmtmb's own and does not name this family:

    cens()/trunc() need a family with a CDF (currently: gaussian,
    lognormal, poisson, exponential, weibull, inverse.gaussian, cox) ...

Read it as a decision rather than as an omission. It arrives from frame
assembly, which runs before any family-supplied check, so this family
has no seam that does not require claiming a CDF it does not have.

## Accuracy

The no-go branch needs the Wiener defective distribution function, which
this package did not have:
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
declares no `lcdf`, and says so in its own compatibility table. It is
written for this family in `R/wiener-cdf.R`, as two series blended in
`log(u)` the way the density's two are.

The blend is centred much lower than the density's, at `u = 0.02` rather
than 0.35, and that is the substantive choice here. The small-time route
reaches the no-go probability as `1 - F_upper`, which cancels precisely
where a no-go trial is surprising; the large-time route computes it
directly, as a gambler's-ruin probability plus a correction, and never
subtracts. Handing over as early as the large-time route is accurate is
what keeps the likelihood honest on the rows that pull hardest on it.

Measured against a 260-bit reference over 1200 points spanning `t` in
0.05 to 15, drift in -2 to 5, boundary separation 0.8 to 4 and relative
start point 0.25 to 0.9, the blended probability holds **2.2e-12**
relative, with one row worse than 1e-12 and none worse than 1e-9.

The centre was 0.06 at 0.3.0 and is 0.02 because of what a wider grid
showed. At `t = 2.5`, drift 5, separation 4 and start point 0.90, where
the true no-go probability is 8.2e-16, the 0.06 blend was 7.75e-05
relative. Neither series was at fault: the large-time route alone was
6.7e-15 there, and the small-time route was 12.1 RELATIVE, because
`1 - F_upper` had nothing left to subtract from. A weight of about 3e-05
on a hopeless number is what the error was. A smooth blend has no safe
side unless both branches degrade gently, and this one has a branch that
does not, so the weight has to saturate before that branch collapses.

The two series are independent derivations and agree with each other to
4.8e-78 at 260 bits, which is the check that the derivation is right
rather than merely stable. Against
[`WienR::pWDM()`](https://rdrr.io/pkg/WienR/man/WienerCDF.html), which
is what EMC2 calls, agreement is 3.3e-13 wherever the no-go probability
is above 0.01. Below that the reference is the weaker of the two, though
by how much depends on what it is asked for. At one grid point with a
no-go probability of 1.19e-13, `1 - WienR::pWDM()` is 4.4 percent wrong
at WienR's DEFAULT precision and at every setting down to
`precision = 1e-12`, 0.13 percent at 1e-14 and 0.033 percent at 1e-16.
What matters for the comparison with EMC2 is the first row of that list:
`EMC2:::pDDM` calls `pWDM` with `precision = 0.005`, looser than any of
them, so EMC2 in practice sits at the 4.4 percent end. This family's
large-time route is 3.0e-15 there.

## References

Gomez, P., Ratcliff, R. and Perea, M. (2007). A model of the go/no-go
task. *Journal of Experimental Psychology: General*, 136(3), 389-413.

## See also

[`wiener_gng_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng_simulate.md)
to generate from the model, and
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
for the two-choice diffusion this is the one-boundary version of.

## Examples

``` r
set.seed(1)
dat <- wiener_gng_simulate(600, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
mean(dat$responded)
#> [1] 0.7516667
fit <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
           family = wiener_gng(deadline = 1.5), data = dat)
fixef(fit)
#> $mu
#> (Intercept) 
#>   0.8922841 
#> 
#> $bs
#> (Intercept) 
#>    0.341704 
#> 
#> $ndt
#> (Intercept) 
#>    1.706139 
#> 
#> $bias
#> (Intercept) 
#>           0 
#> 
```

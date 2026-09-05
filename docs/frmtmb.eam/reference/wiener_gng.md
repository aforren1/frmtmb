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
wiener_gng(deadline = NULL, max_ndt = NULL)
```

## Arguments

- deadline:

  The response deadline, in the units of the response. One positive
  number when every trial shares it. `NULL`, the default, takes it per
  row from `vreal()`.

- max_ndt:

  Upper bound for the non-decision time. `NULL`, the default, takes it
  from the fastest go response.

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

## Across-trial variability, and why there is none here

[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
offers Ratcliff's `sv`, `sz` and `st0`. This family offers none of the
three, and the reason is the no-go branch rather than an oversight.

The go branch would inherit all three for free, because it is
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)'s
density. The no-go branch would not. The drift enters the DENSITY as an
exponential-quadratic, which is what makes averaging over a normal drift
exact and free; it enters the DISTRIBUTION FUNCTION through the
eigenvalues of both series, where it is not, so `sv` would need a
quadrature of its own rather than a completed square. `sz` and `st0`
could be reached with the existing Gauss-Legendre nodes, but shipping
two of three would make `variability =` mean something different here
than it does on
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
which is worse than not having it.

EMC2's `DDMGNG` does carry all three, by calling a compiled distribution
function that integrates them numerically. If you need them, that is
where they are.

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
rather than merely stable. Against `WienR::pWDM()`, which is what EMC2
calls, agreement is 3.3e-13 wherever the no-go probability is above
0.01. Below that the reference is the weaker of the two, though by how
much depends on what it is asked for. At one grid point with a no-go
probability of 1.19e-13, `1 - WienR::pWDM()` is 4.4 percent wrong at
WienR's DEFAULT precision and at every setting down to
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
#>     1.70614 
#> 
#> $bias
#> (Intercept) 
#>           0 
#> 
```

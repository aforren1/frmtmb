# The Wiener first-passage log density, for another package's family

The log density of the two-boundary Wiener diffusion's first passage
time, as
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
evaluates it: Navarro and Fuss's (2009) two series combined with a
smooth weight, so that the result differentiates exactly on an 'RTMB'
tape rather than branching on a parameter. This is the exported seam
another extension package builds a joint choice-and-response-time family
on, and it is the whole of what this package promises outside itself.

## Usage

``` r
wiener_lpdf(dt, drift, bs, bias, upper)
```

## Arguments

- dt:

  Decision time: the response time less the non-decision time. Not the
  response time.

- drift:

  Drift rate, on the real line. Positive drift moves the accumulator
  toward the upper boundary.

- bs:

  Boundary separation, positive.

- bias:

  Relative start point, in `(0, 1)`. `0.5` starts halfway.

- upper:

  Which boundary the response landed on, `1` for the upper and `0` for
  the lower. Data, not a parameter: the reflection that turns the
  lower-boundary series into the upper-boundary one is arithmetic on
  this rather than a branch.

## Value

A numeric or advector vector, the recycled length of the arguments.

## Tape safety

Every argument may be an 'RTMB' advector except `upper`, which selects a
boundary and must be data. The function contains no comparison and no
branch on a parameter, so it tapes; the price is that both series are
evaluated on every call, which is what makes the density correct over
the whole range of normalized times rather than over the part one
truncation happens to cover.

## What a decision time of zero or less gives

`-Inf`, not `NaN`. The density is zero there, so the log density is
`-Inf`, and that is a value a mixture's log-sum-exp can use. A family
calling this is responsible for keeping a fit away from that region;
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
does it with a bounded link on the non-decision time.

## References

Navarro, D. J. and Fuss, I. G. (2009). Fast and accurate calculations
for first-passage times in Wiener diffusion models. *Journal of
Mathematical Psychology* 53, 222-230.

## See also

[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
the family this density serves inside this package.

## Examples

``` r
# the density over the upper boundary integrates to the probability
# of ever reaching it, which the diffusion has in closed form
v <- 1.2; a <- 1.5; w <- 0.3
t <- seq(1e-5, 12, length.out = 200000)
p_up <- sum(exp(wiener_lpdf(t, v, a, w, 1))) * diff(t)[1]
p_lo <- (exp(-2 * v * a) - exp(-2 * v * a * w)) / (exp(-2 * v * a) - 1)
c(quadrature = p_up, analytic = 1 - p_lo)
#> quadrature   analytic 
#>  0.6789561  0.6789561 

# and it differentiates on an RTMB tape, which is the point
tp <- RTMB::MakeTape(function(v) sum(wiener_lpdf(c(0.4, 0.9), v,
                                                 1.5, 0.3, 1)), 1)
tp$jacobian(1.2)
#>      [,1]
#> [1,] 0.54
```

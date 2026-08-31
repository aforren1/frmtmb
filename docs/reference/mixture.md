# Finite mixture families

`mixture(fam1, fam2, ...)` builds a K-component mixture: each component
keeps its own distributional parameters, suffixed by the component index
(`mu1`, `sigma1`, `mu2`, ...), and the mixing proportions come from
`theta1 ... theta{K-1}` (multinomial-logit against the last component,
each with its own linear predictor - so mixing weights may depend on
covariates). The main model formula applies to every component mean;
override per component with `bf(y ~ x, mu2 ~ 1)`.

## Usage

``` r
mixture(...)
```

## Arguments

- ...:

  Two or more component families.

## Value

A `frmtmb_family`.

## Details

The likelihood is a parameter-branch-free logsumexp, so Laplace
machinery is untouched; the usual finite-mixture ML caveats apply
instead: the likelihood is invariant to component relabeling, so the
component means are initialized on spread-out response quantiles, and
multimodality is real (compare starts, or order the intercepts through
`lower`/`upper`). Component families with extra parameters (ordinal) are
not supported.

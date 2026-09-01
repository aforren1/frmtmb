# Likelihood-ratio tests between nested frmtmb fits

ML fits compare freely. REML fits compare only with each other, and only
when their fixed-effect designs span the same column space: a REML
likelihood is a likelihood for the error contrasts of that design, so
two of them are on a common scale exactly when the design is the same.
That covers the usual REML use - testing variance-component structures
with the fixed effects held fixed - and refuses the rest with the reason
(glmmTMB#776).

## Usage

``` r
# S3 method for class 'frmtmb_fit'
anova(object, ..., refit = FALSE)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Further `frmtmb_fit` objects, nested with `object`.

- refit:

  If `TRUE`, refit every REML fit in the comparison with ML and compare
  those, with a message naming what was refit. The refits reuse the
  assembled design and warm-start at the REML estimates. `FALSE` (the
  default) keeps the REML fits and refuses the comparisons a restricted
  likelihood cannot support.

## Value

An `anova` table.

## Details

`refit = TRUE` is the lme4 convenience for the refused case: every REML
fit in the comparison is refit with `REML = FALSE` and the ML fits are
compared instead. lme4 does this silently by default; here it is opt-in
and the message names the models that were refit.

When the smaller model removes a variance component, the null value sits
on the boundary of the parameter space and the usual chi-square
reference is wrong: the asymptotic null is a mixture (for one component,
half a point mass at zero and half a chi-square with one df), so the
reported p-value is conservative - up to a factor of two for a single
component. lme4 and glmmTMB report the same naive p-value; halve it for
the one-component case, or use
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
for a simulation-based reference.

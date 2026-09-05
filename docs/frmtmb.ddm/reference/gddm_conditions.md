# Build a condition index for [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)

One solve of the generalized drift-diffusion likelihood serves every
trial that shares a parameter vector, and the family finds those trials
through an index it is given rather than by comparing parameter values,
which it cannot do on a tape. This builds the index: the distinct
combinations of the variables you name, numbered.

## Usage

``` r
gddm_conditions(data, ...)
```

## Arguments

- data:

  A data frame.

- ...:

  Bare variable names, or a one-sided formula naming them.

## Value

An integer vector, one entry per row of `data`.

## Details

Name every variable that appears on the right-hand side of any formula
in the model, and every covariate a drift term reads. Naming more than
that is safe and only costs solves; naming fewer is wrong, and wrong in
a way nothing downstream can detect.

## See also

[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)

## Examples

``` r
d <- data.frame(coh = c(0, 0, 0.5, 0.5), block = c(1, 2, 1, 2))
gddm_conditions(d, coh, block)
#> [1] 1 2 3 4
gddm_conditions(d, ~ coh)
#> [1] 1 1 2 2
```

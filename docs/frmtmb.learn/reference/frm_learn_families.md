# The families this package supplies, and what each one is

A reference table rather than a computation. It carries the name each
family has here, the name hBayesDM gives the same model where it has
one, the task, the options per trial, the addition terms and the
parameters with their links.

## Usage

``` r
frm_learn_families()
```

## Value

A data frame, one row per family.

## Details

The parameter names are this package's, not hBayesDM's, and the
`hbayesdm_pars` column is the map, in the same order as `pars`. A cell
that is a name is a RENAME; a cell that is an expression is a TRANSFORM,
and three of the six rows carry one. Two of those three involve an
estimated parameter (`-tau * alpha` for
[`prl_fictitious()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/prl_fictitious.md)'s
`bias`, `pi / tau1` for
[`ts_par7()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/ts_par7.md)'s
`pers`), so no constant relates them and a value carried across
unchanged fits a different model without complaint.
[`?prl_fictitious`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/prl_fictitious.md)
and
[`?ts_par7`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/ts_par7.md)
derive both. One convention across six families is worth more than six
inherited ones: `alpha` is always a learning rate on `(0, 1)` and `tau`
is always the softmax sensitivity on `(0, Inf)`, so what a formula does
to `alpha` means the same thing whichever family it is written against.
hBayesDM spells the two-armed bandit's learning rate `A`, the reversal
models' sensitivity `beta`, and the Iowa gambling task's UTILITY
EXPONENT `alpha`.

ONE ROW HAS NO COUNTERPART.
[`rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.md)
is a delta rule feeding a drift-diffusion choice rule, and hBayesDM
ships the two halves separately (`bandit2arm_delta` and `choiceRT_ddm`)
rather than the join. Its `hbayesdm` cell says `none` rather than naming
a model it is not, and
[`?rlddm`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.md)
says which of that package's parameters map onto which of its own.
[`igt_orl()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_orl.md)'s
cell is a rename with one transform: hBayesDM bounds `K` to `(0, 5)`
through a probit where this family puts it on `(0, Inf)` with a log
link, and the two are the same number on the same scale wherever
hBayesDM's bound is not binding.

Two cautions the map does not remove. Core refuses a distributional
parameter name containing a dot or an underscore, which is why the
dual-rate family spells its rates `Arew` and `Apun` rather than anything
more readable. And `lambda` is used by the literature for three
unrelated quantities: the arm decay in the Kalman filter, the
eligibility trace in the two-step model, and loss aversion in the Iowa
gambling task. Each family's help says which.

## Examples

``` r
frm_learn_families()[, c("family", "hbayesdm", "pars")]
#>                      family                   hbayesdm
#> 1          bandit2arm_delta           bandit2arm_delta
#> 2           bandit2arm_dual prl_rp (split = 'outcome')
#> 3            prl_fictitious             prl_fictitious
#> 4 bandit4arm2_kalman_filter  bandit4arm2_kalman_filter
#> 5                   ts_par7                    ts_par7
#> 6             igt_pvl_delta              igt_pvl_delta
#> 7                   igt_orl                    igt_orl
#> 8                     rlddm                       none
#>                                          pars
#> 1                                  alpha, tau
#> 2                             Arew, Apun, tau
#> 3                            alpha, bias, tau
#> 4    tau, lambda, center, mu0, sigma0, sigmaD
#> 5 w, alpha1, tau1, alpha2, tau2, lambda, pers
#> 6                   alpha, shape, lambda, tau
#> 7                 Arew, Apun, k, betaF, betaP
#> 8                 alpha, drift, bs, ndt, bias
```

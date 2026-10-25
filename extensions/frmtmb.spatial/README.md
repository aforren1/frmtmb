# frmtmb.spatial

Spatial random fields for [frmtmb](https://github.com/aforren1/frmtmb)
models. The package registers its fields with frmtmb through the public
covariance-structure seam, `frmtmb::frmtmb_register_covstruct()`.

* `car(M, gr = g, type = )`: brms's conditional autoregressive field
  over an adjacency matrix, in its four types.
* `spde(fem, gr = node)`: a Matern field over a finite-element mesh,
  from the mesh's three matrices (`fmesher::fm_fem()`).
* `sar(M, type = )`: brms's simultaneous autoregressive term, in its lag
  and error forms, registered through the residual correlation seam,
  `frmtmb::frmtmb_register_autocor()`.

## Installation

frmtmb.spatial lives in the `extensions/` directory of the frmtmb
repository. Install frmtmb first, then this package:

```r
remotes::install_github("aforren1/frmtmb")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.spatial")
```

## Use

```r
library(frmtmb.spatial)
fit <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")),
           family = gaussian(), data = d, data2 = list(W = W))
confint_varcorr(fit)
```

Without the package loaded, `frm()` refuses `car()`, `spde()` and
`sar()` and names this package.
`vignette("spatial", package = "frmtmb.spatial")` covers all three.

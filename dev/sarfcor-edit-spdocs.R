setwd("C:/Users/adf44/source/r/frmtmb-wt-sarfcor/extensions/frmtmb.spatial")
rep1 <- function(p, old, new) {
  txt <- paste(readLines(p), collapse = "\n")
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop(p, " matches: ", n, " for ", substr(old, 1, 60))
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], p)
}
rep1("vignettes/spatial.Rmd", "frmtmb.spatial adds two spatial random fields to `frm()`: `car()`,
brms's conditional autoregressive field over an adjacency matrix, and
`spde()`, a Matern field over a finite-element mesh. Load it with
`library(frmtmb.spatial)`; without it, `frm()` refuses the terms below
and names this package. The fields are registered with core through
`frmtmb::frmtmb_register_covstruct()`, the seam any package can use.",
"frmtmb.spatial adds two spatial random fields to `frm()`: `car()`,
brms's conditional autoregressive field over an adjacency matrix, and
`spde()`, a Matern field over a finite-element mesh. It also adds
`sar()`, brms's simultaneous autoregressive term, which correlates the
residuals rather than adding a field. Load it with
`library(frmtmb.spatial)`; without it, `frm()` refuses the terms below
and names this package. The fields are registered with core through
`frmtmb::frmtmb_register_covstruct()`, and `sar()` through
`frmtmb::frmtmb_register_autocor()`, the seams any package can use.")
txt <- paste(readLines("vignettes/spatial.Rmd"), collapse = "\n")
txt <- paste0(txt, "

## Simultaneous autoregressive terms

`sar(M, type = )` makes the residuals of a `gaussian()` or `student()`
response one spatially autoregressive draw over a weights matrix `M`,
as in brms 2.23.0. With `A = I - rho W`, `type = \"lag\"` (the default)
is SAR of the response, `A y - mu ~ N(0, sigma^2 I)`, so the mean of
the response is `A^-1 mu`, which `fitted()` reports; `type = \"error\"`
is SAR of the residuals, `A (y - mu) ~ N(0, sigma^2 I)`, with mean
`mu`. `summary()` reports `rho` as brms names it, `lagsar` or
`errorsar`.

`M` covers every row of the data, as a matrix, a sparse `Matrix`, or
spdep's `nb` or `listw` object, which is read as brms reads it, without
spdep. Unlike `car()`, `sar()` takes no grouping variable: each row is
a location. `rho` is estimated on `(1 / min eigenvalue, 1 / max
eigenvalue)` of `M`, where `I - rho M` has a positive determinant;
brms takes the eigenvalue range itself, which for a row-standardized
`M` cuts off the negative values below the smallest eigenvalue.

```{r sar}
set.seed(2)
g <- expand.grid(r = 1:6, c = 1:6)
n <- nrow(g)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- A / rowSums(A)
ds <- data.frame(x = rnorm(n))
ds$y <- as.vector(solve(diag(n) - 0.5 * W, 1 + 0.5 * ds$x + rnorm(n)))
fit_sar <- frm(bf(y ~ x + sar(W)), family = gaussian(), data = ds,
               data2 = list(W = W))
confint_varcorr(fit_sar)
```

`sigma` and `nu` must be constant, as in brms. The likelihood is one
joint density over all the rows, so `weights()`, `cens()`, `trunc()`,
`se()` and `mi()` are refused, and so is `newdata` in `fitted()`,
`predict()` and `simulate()`: `M` describes the fitted rows only.
brms refuses the same calls unless its `newdata2` holds a matrix over
the new rows.
")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], "vignettes/spatial.Rmd")
rep1("README.md", "* `spde(fem, gr = node)`: a Matern field over a finite-element mesh,
  from the mesh's three matrices (`fmesher::fm_fem()`).",
"* `spde(fem, gr = node)`: a Matern field over a finite-element mesh,
  from the mesh's three matrices (`fmesher::fm_fem()`).
* `sar(M, type = )`: brms's simultaneous autoregressive term, in its lag
  and error forms, registered through the residual correlation seam,
  `frmtmb::frmtmb_register_autocor()`.")
rep1("README.md", "Without the package loaded, `frm()` refuses `car()` and `spde()` and
names this package. `vignette(\"spatial\", package = \"frmtmb.spatial\")`
covers both.",
"Without the package loaded, `frm()` refuses `car()`, `spde()` and
`sar()` and names this package.
`vignette(\"spatial\", package = \"frmtmb.spatial\")` covers all three.")

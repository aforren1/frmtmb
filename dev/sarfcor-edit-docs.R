# One-off doc edits for lane sarfcor (README, vignettes).
setwd("C:/Users/adf44/source/r/frmtmb-wt-sarfcor")
rep1 <- function(p, old, new) {
  txt <- paste(readLines(p), collapse = "\n")
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop(p, " matches: ", n, " for ", substr(old, 1, 60))
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], p)
}
rep1("README.md", "  REML. brms's default `cov = FALSE` form of `ar()`, `ma()` and
  `arma()` is brms's own likelihood, a regression on earlier
  residuals.",
"  REML. brms's default `cov = FALSE` form of `ar()`, `ma()` and
  `arma()` is brms's own likelihood, a regression on earlier
  residuals. brms's `fcor(M)`, a fixed residual covariance over the
  observations, is in frmtmb too, and its `sar()` in frmtmb.spatial.")
rep1("vignettes/brms-migration.Rmd",
"| `ar(t, g)`, `ma()`, `arma()`, `cosy()`, `unstr()` | same | gaussian/student only; both `cov = FALSE` (brms's default) and `cov = TRUE` (below) |",
"| `ar(t, g)`, `ma()`, `arma()`, `cosy()`, `unstr()` | same | gaussian/student only; both `cov = FALSE` (brms's default) and `cov = TRUE` (below) |
| `fcor(M)` | same | gaussian/student only; `M` in `data2`, `sigma` scales it, and a predicted `sigma` gives `D M D`, as in brms; no `newdata`, as brms without `newdata2` |
| `sar(M, type = \"lag\")`, `sar(M, type = \"error\")` | same, with `library(frmtmb.spatial)` | gaussian/student only, constant `sigma`, as in brms; `rho` on `(1 / min eigenvalue, 1 / max eigenvalue)`, where brms takes the eigenvalue range itself |")
rep1("vignettes/compatibility.Rmd",
"`ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` replace the
response's per-row density with one joint density per group.",
"`ar()`, `ma()`, `arma()`, `cosy()`, `unstr()` and `fcor()` replace the
response's per-row density with one joint density per group (for
`fcor()`, one group holding every row).")

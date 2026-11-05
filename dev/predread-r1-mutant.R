# The test of fix round 2's R1 against the two earlier readings of me_band_end():
# fix round 1's body, and the one-line change with the direction still off the
# first predictor. Each must fail it.
#   Rscript dev/predread-r1-mutant.R <library>
args <- commandArgs(TRUE)
.libPaths(c(args[1L], "C:/Users/adf44/source/r/rellib-r16", .libPaths()))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")
muts <- list(
  round1 = function(mm, se_coef, p) {
    gr <- abs(mm$grad)
    extra <- ifelse(gr > 0 & mm$s2 > 0, (as.numeric(se_coef) / gr)^2 / mm$s2, 0)
    z <- stats::qnorm(p) * sqrt(1 + extra) * ifelse(mm$grad < 0, -1, 1)
    v <- drop(mm$val(z))
    v[!is.finite(v)] <- NA_real_
    v
  },
  sign_of_lead = function(mm, se_coef, p) {
    extra <- ifelse(mm$var > 0, as.numeric(se_coef)^2 / mm$var, 0)
    z <- stats::qnorm(p) * sqrt(1 + extra) * ifelse(mm$grad < 0, -1, 1)
    v <- drop(mm$val(z))
    v[!is.finite(v)] <- NA_real_
    v
  })
muts$current <- get("me_band_end", ns)
for (m in names(muts)) {
  f <- muts[[m]]
  environment(f) <- ns
  assignInNamespace("me_band_end", f, "frmtmb")
  env <- new.env(parent = ns)
  sys.source("tests/testthat/helper-warnings.R", envir = env)
  r <- as.data.frame(test_file(
    "tests/testthat/test-me-newdata.R", env = env, reporter = "silent",
    desc = "fitted()'s interval holds when two predictors read the term"))
  cat(sprintf("%-14s tests %d failed %d\n", m, sum(r$nb), sum(r$failed)))
}

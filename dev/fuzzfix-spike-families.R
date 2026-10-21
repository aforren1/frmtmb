# For each family the unbounded-weights refusal names, fit an
# observation-level intercept with constant weights 1.2 (sum(w) > n)
# and 0.9 (sum(w) < n). The first must run away (no maximum), the
# second must fit. On the lane build the first must be refused.
#   Rscript dev/fuzzfix-spike-families.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261020)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
eta <- 0.3 + 0.4 * d$x + rnorm(n, 0, 0.3)
pos <- exp(eta) * rgamma(n, 4, 4)
fams <- list(
  gaussian = list(gaussian(), eta + rnorm(n, 0, 0.5)),
  student = list(student(), eta + 0.5 * rt(n, 5)),
  lognormal = list(lognormal(), pos),
  Gamma = list(Gamma(link = "log"), pos),
  weibull = list(weibull(), pos),
  Beta = list(Beta(), plogis(eta + rnorm(n, 0, 0.5))),
  inverse.gaussian = list(inverse.gaussian(link = "log"), pos),
  skew_normal = list(skew_normal(), eta + rnorm(n, 0, 0.5)),
  shifted_lognormal = list(shifted_lognormal(), pos + 0.2),
  asym_laplace = list(asym_laplace(), eta + rnorm(n, 0, 0.5)),
  huber = list(huber(), eta + rnorm(n, 0, 0.5)),
  von_mises = list(von_mises(), 2 * atan(0.5 * (eta - 0.3) +
                                           rnorm(n, 0, 0.3))),
  frechet = list(frechet(), pos),
  gen_extreme_value = list(gen_extreme_value(), eta + rnorm(n, 0, 0.5)),
  # bounded: a probability mass never exceeds one
  poisson = list(poisson(), rpois(n, exp(eta))),
  exponential = list(exponential(), pos)
)
for (nm in names(fams)) {
  for (wv in c(1.2, 0.9)) {
    dd <- d; dd$y <- fams[[nm]][[2]]; dd$w <- wv
    warns <- character()
    f <- withCallingHandlers(
      tryCatch(frm(bf(y | weights(w) ~ 1 + x + (1 | obs)),
                   family = fams[[nm]][[1]], data = dd),
               error = function(e) e),
      warning = function(x) {
        warns <<- c(warns, conditionMessage(x))
        invokeRestart("muffleWarning")
      }, message = function(m) invokeRestart("muffleMessage"))
    res <- if (inherits(f, "error")) {
      paste("ERROR", substr(gsub("\n", " ", conditionMessage(f)), 1, 70))
    } else {
      sprintf("logLik %.6g betad %s", as.numeric(logLik(f)),
              paste(format(f$estimates$betad, digits = 3),
                    collapse = ","))
    }
    cat(sprintf("FAM %-18s w=%.1f %s | warn %d\n", nm, wv, res,
                length(warns)))
  }
}
cat("END\n")

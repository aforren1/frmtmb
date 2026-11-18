# Review of lane nbrobust, check 6: the families' pointwise log density
# off the tape (the route log_lik() and the bin-2 adapters take) against
# brms 2.23.0's own log_lik_<family>() on a brms-style prep, as
# tests.log_lik.R (rows log_lik:285, 288, 295, 414-436) builds it, at
# brms's shapes (rgamma(4)) and at poisson-limit shapes (exp(15..35)).
# Also one fit's rows at its MLE against brms's formulas.
#   Rscript dev/reviews/nbrobust-review-parity.R <base|lib>
args <- commandArgs(TRUE)
arm <- args[1]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages({
  library(frmtmb)
})
cat("lib:", dirname(system.file(package = "frmtmb")), " brms",
    as.character(utils::packageVersion("brms")), "\n")
ll_frm <- function(fam, y, dp) fam$lpdf(rep(y, length(dp$mu)), dp, list())
cmp <- function(label, a, b) {
  ok <- is.finite(a) & is.finite(b)
  rel <- max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300))
  cat(sprintf(paste("%-34s n %3d nonfinite frm %d brms %d max rel %9.3g",
                    "expect_equal %s\n"),
              label, length(a), sum(!is.finite(a)), sum(!is.finite(b)), rel,
              isTRUE(all.equal(a, b))))
}
for (regime in c("brms", "limit")) {
  set.seed(if (regime == "brms") 31 else 32)
  ns <- 50
  nobs <- 8
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  shape <- if (regime == "brms") stats::rgamma(ns, 4) else
    exp(stats::runif(ns, 15, 35))
  prep$dpars <- list(eta = matrix(stats::rnorm(ns * nobs), ncol = nobs),
                     shape = shape, zi = stats::rbeta(ns, 1, 1))
  prep$dpars$hu <- prep$dpars$zi
  prep$dpars$sigma <- 1 / prep$dpars$shape
  prep$dpars$mu <- exp(prep$dpars$eta)
  prep$data <- list(Y = c(0, 1, 3, 7, 12, 0, 2, 40))
  for (i in seq_len(nobs)) {
    y <- prep$data$Y[i]
    mu <- prep$dpars$mu[, i]
    tag <- sprintf("%-5s i %d y %2d", regime, i, y)
    cmp(paste(tag, "negbinomial"),
        ll_frm(negbinomial(), y, list(mu = mu, shape = shape)),
        brms:::log_lik_negbinomial(i, prep = prep))
    cmp(paste(tag, "negbinomial2"),
        ll_frm(brmsfamily("negbinomial2"), y,
               list(mu = mu, sigma = prep$dpars$sigma)),
        brms:::log_lik_negbinomial2(i, prep = prep))
    cmp(paste(tag, "geometric"),
        ll_frm(geometric(), y, list(mu = mu)),
        brms:::log_lik_geometric(i, prep = prep))
    cmp(paste(tag, "zero_inflated_negbinomial"),
        ll_frm(zero_inflated_negbinomial(), y,
               list(mu = mu, shape = shape, zi = prep$dpars$zi)),
        brms:::log_lik_zero_inflated_negbinomial(i, prep = prep))
    cmp(paste(tag, "hurdle_negbinomial"),
        ll_frm(hurdle_negbinomial(), y,
               list(mu = mu, shape = shape, hu = prep$dpars$hu)),
        brms:::log_lik_hurdle_negbinomial(i, prep = prep))
  }
}

# Lane nbrobust: the margin of test-numerical-robustness.R's "reads the
# poisson limit" assertion, over 21 ulp perturbations of the linear
# predictor (eta * (1 + 2e-15 * k)), as error / (eps * sum(y + mu + 1));
# the test fails at 64.
#   Rscript dev/nbrobust-gradmargin.R
.libPaths(c("C:/Users/adf44/source/r/wt-nbrobust-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), " R:", R.home(),
    " CORETYPE:", Sys.getenv("OPENBLAS_CORETYPE"), "\n")
dpars_at <- function(fam, etas) {
  out <- list()
  for (nm in names(etas)) {
    out[[nm]] <- fam$links[[nm]]$linkinv(etas[[nm]])
    out[[paste0(".eta_", nm)]] <- etas[[nm]]
  }
  out
}
set.seed(2)
x <- stats::rnorm(200)
y <- stats::rpois(200, exp(1 + 0.4 * x))
worst <- 0
for (k in 0:20) {
  eta <- (1.0544184774 + 0.4005673777 * x) * (1 + 2e-15 * k)
  mu <- exp(eta)
  for (fam in list(negbinomial(), zero_inflated_negbinomial())) {
    base <- list(mu = eta, shape = 0)
    if (fam$family != "negbinomial") base$zi <- -40
    f <- function(p) {
      e <- base
      e$shape <- p[1] + 0 * base$shape
      -sum(fam$lpdf(y, dpars_at(fam, e), list()))
    }
    tp <- RTMB::MakeTape(f, 20)
    for (ls in c(20, 24.4791351853, 30, 35, 40)) {
      n <- exp(ls)
      r <- mu / n
      rise <- vapply(y, function(kk) {
        j <- seq_len(kk) - 1
        sum(j / (n + j))
      }, 0)
      tail <- r^2 / 2 - 2 * r^3 / 3 + 3 * r^4 / 4 - 4 * r^5 / 5
      exact <- -sum(-rise + y * r / (1 + r) - n * tail)
      g <- as.numeric(tp$jacobian(ls))
      u <- abs(g - exact) / (.Machine$double.eps * sum(y + mu + 1))
      worst <- max(worst, u)
      cat(sprintf("k %2d %-25s ls %7.4f grad % .6e exact % .6e units %.3f\n",
                  k, fam$family, ls, g, exact, u))
    }
  }
}
cat(sprintf("RESULT worst units %.3f of 64\n", worst))

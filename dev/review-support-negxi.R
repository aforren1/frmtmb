# Reviewer probe (lane support): GEV fits with xi near -1, where the
# upper end of the support sits on the largest row. Does the doubles
# reader (frmtmb.sample draws_row_loglik()) still sum to logLik()?
# Also RTMB::logspace_add on doubles at -Inf.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
print(RTMB::logspace_add(c(-Inf, -Inf, 0), c(-Inf, 1, -Inf)))
rgev <- function(n, mu, sigma, xi) mu + sigma * ((-log(runif(n)))^(-xi) - 1) / xi
for (xi in c(-0.6, -0.8, -0.95)) for (seed in 1:5) {
  set.seed(seed)
  d <- data.frame(y = rgev(60, 0, 1, xi))
  fit <- tryCatch(suppressWarnings(frm(y ~ 1, family = gen_extreme_value(),
                                       data = d)),
                  error = function(e) conditionMessage(e))
  if (is.character(fit)) {
    cat(sprintf("xi %5.2f seed %d ERROR %s\n", xi, seed, substr(fit, 1, 100)))
    next
  }
  ll <- as.numeric(logLik(fit))
  rows <- frmtmb.sample:::draws_row_loglik(fit, NULL)
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  ub <- dp$mu - dp$sigma / dp$xi
  tmax <- 1 + dp$xi[1] * (max(d$y) - dp$mu[1]) / dp$sigma[1]
  cat(sprintf(paste("xi %5.2f seed %d  xi_hat %.4f  t at max(y) %.3g  logLik",
                    "%.10g  rows not finite %d  sum - logLik %.3g\n"),
              xi, seed, dp$xi[1], tmax, ll, sum(!is.finite(rows)),
              sum(rows) - ll))
}

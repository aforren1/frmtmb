## Why fit$opt$nonfinite_trials is 0 on two fits whose base run warned
## NA/NaN 7 and 4 times: list every run_optimizer() call and its count.
ARM <- "lane"
source("dev/famlink-rev-common.R")
calls <- list()
trace("run_optimizer", where = asNamespace("frmtmb"), print = FALSE,
      exit = quote(calls[[length(calls) + 1L]] <<- list(n = returnValue()$nonfinite_trials,
                                                        obj = returnValue()$objective,
                                                        start = par)))
mk <- function(n = 200) { set.seed(20260916); data.frame(x = runif(n, -2, 2), g = factor(rep(1:20, length.out = n))) }
for (nm in c("bern_identity_wall", "gamma_inverse_cross")) {
  calls <- list()
  d <- mk()
  if (nm == "bern_identity_wall") { d$y <- rbinom(200, 1, pmin(0.98, pmax(0.02, 0.5 + 0.45 * d$x))); fam <- bernoulli("identity") }
  else { mu <- 1 / pmax(0.05, 0.6 + 0.35 * d$x); d$y <- rgamma(200, 3, 3 / mu); fam <- Gamma("inverse") }
  nanw <- 0L
  fit <- withCallingHandlers(frm(y ~ x, d, family = fam, verbose = TRUE),
    warning = function(w) { if (grepl("NA/NaN", conditionMessage(w))) nanw <<- nanw + 1L; invokeRestart("muffleWarning") },
    message = function(m) { cat("  |", conditionMessage(m)); invokeRestart("muffleMessage") })
  cat(nm, ": run_optimizer calls", length(calls), "; counts", vapply(calls, function(z) z$n %||% NA_integer_, 1L),
      "; fit$opt$nonfinite_trials", fit$opt$nonfinite_trials, "; NaN warnings", nanw, "\n")
}

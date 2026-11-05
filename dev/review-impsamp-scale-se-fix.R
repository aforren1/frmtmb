# SEs and consumers on scaled or shifted columns: does the tape behind the
# reported fit give the Laplace SEs back? Prints each case with the message of
# any refusal in full. Usage: REV_ARM=after Rscript dev/review-impsamp-scale-se.R
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("mutant")
suppressPackageStartupMessages(library(frmtmb))
options(width = 180)
mk <- function(fam, seed = 1L) {
  set.seed(seed)
  G <- 15L
  g <- factor(rep(seq_len(G), each = 10))
  x <- rnorm(150)
  u <- rnorm(G, 0, 0.8)
  y <- if (fam == "poisson") rpois(150, exp(0.5 + 0.5 * x + u[g])) else
    rbinom(150, 1, plogis(0.2 + 0.7 * x + u[g]))
  data.frame(y, x, g)
}
famf <- list(poisson = poisson(), bernoulli = bernoulli())
for (fam in names(famf)) {
  d0 <- mk(fam)
  base <- suppressWarnings(frm(bf(y ~ xx + (1 | g)), family = famf[[fam]],
                               data = transform(d0, xx = x), importance = 300L))
  sb <- sqrt(diag(frmtmb:::sdr_of(base)$cov.fixed))
  for (cs in list(c(1e3, 0), c(1e4, 0), c(1e5, 0), c(1, 2000), c(100, 1e4),
                  c(10, 1000))) {
    d <- transform(d0, xx = cs[1] * x + cs[2])
    w <- character(0)
    t0 <- Sys.time()
    fit <- withCallingHandlers(
      tryCatch(frm(bf(y ~ xx + (1 | g)), family = famf[[fam]], data = d,
                   importance = 300L), error = function(e) e),
      warning = function(x2) {
        w <<- c(w, conditionMessage(x2))
        invokeRestart("muffleWarning")
      })
    tt <- as.numeric(Sys.time() - t0, units = "secs")
    cat(sprintf("\n== %s scale %g shift %g [%.1fs]\n", fam, cs[1], cs[2], tt))
    if (inherits(fit, "error")) {
      cat("REFUSAL:", conditionMessage(fit), "\n")
      next
    }
    s <- tryCatch(sqrt(diag(frmtmb:::sdr_of(fit)$cov.fixed)),
                  error = function(e) rep(NA, 3))
    # the slope's SE scales by 1/cs[1]; the intercept mixes under a shift, so
    # compare the slope and theta only
    cat("rounds", fit$importance$rounds, "grads", format(fit$importance$grads, digits = 3),
        "moves", format(fit$importance$moves, digits = 3),
        "conv", fit$opt$convergence, "warns", length(w), "\n")
    if (length(w)) cat("  warn:", substr(w[1], 1, 200), "\n")
    cat("SE ratio to standard column: slope",
        format(s[2] * cs[1] / sb[2], digits = 5), " theta",
        format(s[3] / sb[3], digits = 5), "\n")
    print(tryCatch(summary(fit)$fixed[, 1:2], error = function(e) conditionMessage(e)))
    prof <- suppressWarnings(tryCatch(confint(fit, parm = "theta_1",
                                              method = "profile"),
                                      error = function(e) conditionMessage(e)))
    print(prof)
  }
  print(suppressWarnings(confint(base, parm = "theta_1", method = "profile")))
}

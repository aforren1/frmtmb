# Fix round, MAJOR 2: a one-block quadrature runaway that lane wcheck's
# weights check leaves alone (a family outside its table). hurdle_lognormal
# with every y > 0, weights 1.5 and `(1 | obs)`: the positive part is a
# lognormal, whose weighted marginal likelihood rises without bound as
# sigma falls. Exact marginal log-likelihood of the positive rows (the
# hurdle part adds n log(1 - hu), which does not depend on sigma):
#   w [log N(log y; m, sig) - log y] integrated over u ~ N(0, s^2)
#   = sum( -w/2 log(2 pi sig^2) + 1/2 log(2 pi sig^2 / w)
#          + log N(log y - m; 0, sqrt(s^2 + sig^2 / w)) - w log y )
# Printed at the estimate and with log sigma moved down by 1 and 3.
#   Rscript dev/quadseam-runaway2.R <lib|base|wcheck>
# log dev/quadseam-log/runaway2-<tag>.txt
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
wc <- "C:/Users/adf44/source/r/wt-wcheck-lib"
.libPaths(switch(a[1], base = c(base, ul), wcheck = c(wc, base, ul),
                 c(a[1], base, ul)))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
exact_ll <- function(ly, m, w, s, sig) {
  v <- s^2 + sig^2 / w
  sum(-w / 2 * log(2 * pi * sig^2) + 0.5 * log(2 * pi * sig^2 / w) +
        stats::dnorm(ly - m, 0, sqrt(v), log = TRUE) - w * ly)
}
set.seed(20261112)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
d$y <- exp(1 + 0.5 * d$x + rnorm(n, 0, 0.7))
d$w <- 1.5
for (quad in c(TRUE, FALSE)) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(
    frm(bf(y | weights(w) ~ x + (1 | obs)), family = hurdle_lognormal(),
        data = d, quadrature = quad),
    warning = function(x) {
      w <<- c(w, paste(class(x)[1], substr(gsub("[[:space:]]+", " ",
                                                conditionMessage(x)), 1, 150)))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t0
  lab <- if (quad) "quadrature" else "Laplace"
  if (inherits(f, "error")) {
    cat(sprintf("%-10s ERROR: %s\n", lab,
                substr(gsub("[[:space:]]+", " ", conditionMessage(f)), 1, 200)))
    next
  }
  b <- f$estimates$beta
  ls <- f$estimates$betad[grep("sigma", names(f$estimates$betad))][1]
  s <- exp(f$estimates$theta[1])
  qs <- f$opt[["quad_settle"]]
  cat(sprintf("%-10s logLik %.4f log sigma %.3f conv %d secs %.1f settle %s\n",
              lab, as.numeric(logLik(f)), ls, f$opt$convergence, secs,
              if (is.null(qs)) "-" else paste(format(qs, digits = 3),
                                              collapse = "/")))
  for (x in w) cat("    ", x, "\n")
  m <- b[1] + b[2] * d$x
  ex <- vapply(c(0, -1, -3), function(dl) {
    exact_ll(log(d$y), m, d$w, s, exp(ls + dl))
  }, 0)
  cat(sprintf("    exact positive part at est %.4f  log sigma -1: %.4f  -3: %.4f\n",
              ex[1], ex[2], ex[3]))
}
cat("END\n")

# Fix round, MAJOR 2 with the coordinator's addendum: the one-block
# quadrature runaway that lane wcheck's weights check skips by design, a
# gaussian response with mi() and 3 rows missing (wcheck's review,
# frmtmb-wt-wcheck/dev/reviews/wcheck-review-misc.R part "quad"; same
# data). Weights 1.5 and `(1 | obs)`: the weighted marginal likelihood of
# the 77 observed rows rises without bound as sigma falls. Its exact
# value, from the model alone:
#   int N(u; 0, s^2) N(y; m + u, sig)^w du
#   = -w/2 log(2 pi sig^2) + 1/2 log(2 pi sig^2 / w)
#     + log N(y - m; 0, sqrt(s^2 + sig^2 / w))
# printed at the estimate and with log sigma moved down by 1 and 3.
#   Rscript dev/quadseam-runaway3.R <lib|base|wcheck>
# log dev/quadseam-log/runaway3-<tag>.txt
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
wc <- "C:/Users/adf44/source/r/wt-wcheck-lib"
.libPaths(switch(a[1], base = c(base, ul), wcheck = c(wc, base, ul),
                 c(a[1], base, ul)))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
exact_ll <- function(y, m, w, s, sig) {
  v <- s^2 + sig^2 / w
  sum(-w / 2 * log(2 * pi * sig^2) + 0.5 * log(2 * pi * sig^2 / w) +
        stats::dnorm(y - m, 0, sqrt(v), log = TRUE))
}
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                g = factor(rep(1:16, each = 5)), t = rep(1:5, 16))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$w <- 1.5
d$ym <- d$y
d$ym[c(5, 17, 33)] <- NA
ok <- !is.na(d$ym)
for (quad in c(TRUE, FALSE)) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(
    frm(bf(ym | weights(w) + mi() ~ x + (1 | obs)), family = gaussian(),
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
  ls <- f$estimates$betad[1]
  s <- exp(f$estimates$theta[1])
  qs <- f$opt[["quad_settle"]]
  cat(sprintf("%-10s logLik %.4f log sigma %.3f conv %d secs %.1f settle %s\n",
              lab, as.numeric(logLik(f)), ls, f$opt$convergence, secs,
              if (is.null(qs)) "-" else paste(format(qs, digits = 3),
                                              collapse = "/")))
  for (x in w) cat("    ", x, "\n")
  m <- (b[1] + b[2] * d$x)[ok]
  ex <- vapply(c(0, -1, -3), function(dl) {
    exact_ll(d$ym[ok], m, d$w[ok], s, exp(ls + dl))
  }, 0)
  cat(sprintf("    exact (observed rows) at est %.4f  log sigma -1: %.4f  -3: %.4f\n",
              ex[1], ex[2], ex[3]))
}
cat("END\n")

# Items 3466 / 3614: quadrature = TRUE stops a runaway at a spurious
# maximum. Gaussian rows with an observation-level effect, where the
# exact marginal likelihood has a closed form:
#   uncensored, weight w: int N(u; 0, s^2) N(y; m + u, sig)^w du
#     = w log N(y; m, sig) ... written out in exact_ll() below;
#   right censored at c: 1 - pnorm((c - m) / sqrt(s^2 + sig^2)).
# Arms: the rel073 wcens design (80 rows, 16 censored, uncensored excess
# E), and a two-response mvbf() the weights check does not reach. For
# each quadrature fit: the reported logLik, the exact one at the
# estimate, and the exact one with log sigma moved down by 1 and 3 (the
# rest held), which keeps rising on a runaway.
#   Rscript dev/quadseam-runaway.R <lib|base>
# log dev/quadseam-log/runaway-<tag>.txt
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
exact_ll <- function(y, m, w, cens, s, sig) {
  unc <- cens == "none"
  v <- s^2 + sig^2 / w[unc]
  # int N(u;0,s^2) prod N(y;m+u,sig)^w du
  ll_u <- -w[unc] / 2 * log(2 * pi * sig^2) + 0.5 * log(2 * pi * sig^2 /
            w[unc]) + stats::dnorm(y[unc] - m[unc], 0, sqrt(v), log = TRUE)
  # right censored: P(Y > c) = pnorm((m - c) / sd)
  ll_c <- stats::pnorm((m[!unc] - y[!unc]) / sqrt(s^2 + sig^2),
                       log.p = TRUE)
  sum(ll_u) + sum(w[!unc] * ll_c)
}
run <- function(label, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, paste(class(x)[1], conditionMessage(x)))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) {
    cat(sprintf("%-34s ERROR: %s\n", label,
                substr(gsub("[[:space:]]+", " ", r), 1, 150)))
    return(invisible(NULL))
  }
  cat(sprintf("%-34s logLik %12.4f  log sigma %8.3f  conv %d  warn %d\n",
              label, as.numeric(logLik(r)), r$estimates$betad[1],
              r$opt$convergence, length(w)))
  for (x in w) cat("    ", substr(gsub("[[:space:]]+", " ", x), 1, 160), "\n")
  invisible(r)
}
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
q <- quantile(d$y, 0.8)
d$cc <- ifelse(d$y > q, "right", "none")
d$yc <- pmin(d$y, q)
nu <- sum(d$cc == "none")
for (E in c(0, 2, 8, 32)) {
  d$w <- ifelse(d$cc == "none", 1 + E / nu, 1)
  f <- run(sprintf("cens gaussian E %2g quadrature", E),
           frm(bf(yc | weights(w) + cens(cc) ~ x + (1 | obs)),
               family = gaussian(), data = d, quadrature = TRUE))
  if (is.null(f)) next
  b <- f$estimates$beta
  m <- b[1] + b[2] * d$x
  s <- exp(f$estimates$theta[1])
  ls <- f$estimates$betad[1]
  ex <- vapply(c(0, -1, -3), function(dl) {
    exact_ll(d$yc, m, d$w, d$cc, s, exp(ls + dl))
  }, 0)
  cat(sprintf("    exact at est %.4f  log sigma -1: %.4f  -3: %.4f\n",
              ex[1], ex[2], ex[3]))
}
# two gaussian responses, each with its own (1 | obs), weights summing
# above n: the weights check returns for more than one response
set.seed(20261111)
d2 <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
d2$y1 <- 1 + 0.5 * d2$x + rnorm(n)
d2$y2 <- -1 + 0.3 * d2$x + rnorm(n)
d2$w <- 1.5
f <- run("mvbf two gaussian w 1.5 quadrature",
         frm(mvbf(bf(y1 | weights(w) ~ x + (1 | obs)),
                  bf(y2 | weights(w) ~ x + (1 | obs))),
             family = gaussian(), data = d2, quadrature = TRUE))
if (!is.null(f)) {
  # the exact marginal log-likelihood at the estimate, and with both log
  # sigmas moved down (rest held); positions by the outer vector's order
  b <- f$estimates$beta
  ls <- f$estimates$betad
  s <- exp(f$estimates$theta)
  cat("    names beta:", names(b), " betad:", names(ls), "\n")
  none <- rep("none", n)
  ww <- rep(1.5, n)
  ex <- vapply(c(0, -1, -3), function(dl) {
    exact_ll(d2$y1, b[1] + b[2] * d2$x, ww, none, s[1], exp(ls[1] + dl)) +
      exact_ll(d2$y2, b[3] + b[4] * d2$x, ww, none, s[2], exp(ls[2] + dl))
  }, 0)
  cat(sprintf("    exact at est %.4f  log sigmas -1: %.4f  -3: %.4f\n",
              ex[1], ex[2], ex[3]))
}
f <- run("mvbf two gaussian w 1.5 Laplace",
         frm(mvbf(bf(y1 | weights(w) ~ x + (1 | obs)),
                  bf(y2 | weights(w) ~ x + (1 | obs))),
             family = gaussian(), data = d2, quadrature = FALSE))
cat("END\n")

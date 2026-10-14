# Lane postfix, fix round 2: a copy of the review's
# dev/postfix-rev2-nlnew.R whose candidate step skips the syntax example
# of the new message ("<name> = <value>"). The review's header:
# Reviewer of lane postfix, round 2: nl_scan_accept() on new nonlinear
# models that fail at the default start. For each model: the default
# fit, classified as rescued (and how far from the best of a hand start
# and a multistart) or refused; for a refusal, the candidate `start` the
# message names, passed back.
# Usage: Rscript dev/postfix-rev2-nlnew.R <lib or "base">
a <- commandArgs(TRUE)
lib <- if (length(a) && a[1] != "base") a[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
run <- function(expr) {
  msgs <- character(); warns <- character()
  v <- tryCatch(withCallingHandlers(expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
  list(v = v, msgs = msgs, warns = warns)
}
ll <- function(r) if (inherits(r$v, "error")) NA_real_ else
  as.numeric(logLik(r$v))
sefin <- function(r) {
  if (inherits(r$v, "error")) return(NA)
  s <- tryCatch(fixef(r$v)[, "Est.Error"], error = function(e) NA)
  all(is.finite(s))
}
desc <- function(r) {
  if (inherits(r$v, "error")) return("ERROR")
  sprintf("conv %d, SEs finite %s, %d warnings", r$v$opt$convergence,
          sefin(r), length(unique(r$warns)))
}
report <- function(lbl, f, data, family, hand, multi, prior = NULL) {
  fit1 <- function(st = NULL) {
    run(frm(f, data = data, family = family, prior = prior,
            start = if (!is.null(st)) list(beta = st)))
  }
  d0 <- fit1()
  lls <- vapply(c(list(hand), multi), function(s) ll(fit1(s)), 0)
  best <- max(lls, na.rm = TRUE)
  cat(sprintf("\n== %s\n", lbl))
  cat("default:", desc(d0), sprintf(" logLik %.6f  gap %.4e\n", ll(d0),
                                    ll(d0) - best))
  if (!inherits(d0$v, "error")) {
    cat("  est:", signif(fixef(d0$v)[, "Estimate"], 6), "\n")
    for (w in unique(d0$warns)) cat("  [warn]", substr(w, 1, 200), "\n")
    for (m in grep("scan", d0$msgs, value = TRUE))
      cat("  [msg]", substr(m, 1, 120), "\n")
  } else {
    em <- conditionMessage(d0$v)
    tail <- regmatches(em, regexpr("A scan of the nonlinear.*$", em))
    cat("  [error tail]", if (length(tail)) tail else "(no scan note)", "\n")
    m <- regmatches(em, regexpr("start = list[(]beta = c[(].*[)][)]", em))
    if (length(m) && !grepl("<", m, fixed = TRUE)) {
      st <- eval(parse(text = sub("^start = list[(]beta = ", "",
                                  sub("[)]$", "", m))))
      c1 <- fit1(st)
      cat("  candidate start passed back:", desc(c1),
          sprintf(" logLik %.6f  gap %.4e\n", ll(c1), ll(c1) - best))
      if (inherits(c1$v, "error"))
        cat("    ", substr(conditionMessage(c1$v), 1, 200), "\n")
      for (w in unique(c1$warns)) cat("    [warn]", substr(w, 1, 160), "\n")
    }
  }
  cat(sprintf("hand %.6f; best of %d starts %.6f (%d finite)\n", lls[1],
              length(lls), best, sum(is.finite(lls))))
  invisible(list(ll = ll(d0), best = best, err = inherits(d0$v, "error"),
                 conv = if (!inherits(d0$v, "error")) d0$v$opt$convergence,
                 se = sefin(d0)))
}

## N1 sinusoid, true period 3 (between grid values 2 and 5)
set.seed(201)
x <- seq(0, 12, length.out = 150)
d1 <- data.frame(x = x, y = 2 * sin(6.283185307179586 * x / 3) +
                   rnorm(150, 0, 0.5))
f1 <- bf(y ~ A * sin(6.283185307179586 * x / P), A + P ~ 1, nl = TRUE)
ms1 <- lapply(1:20, function(i) c(runif(1, 0.5, 4), runif(1, 0.5, 10)))
report("N1 gaussian sinusoid, A 2, P 3", f1, d1, gaussian(), c(2, 3), ms1)

## N2 biexponential decay, Gamma(identity), four intercepts
set.seed(202)
x2 <- runif(200, 0, 15)
mu2 <- 50 * exp(-x2 / 0.5) + 10 * exp(-x2 / 5)
d2 <- data.frame(x = x2, y = rgamma(200, shape = 20, rate = 20 / mu2))
f2 <- bf(y ~ a1 * exp(-x / t1) + a2 * exp(-x / t2), a1 + t1 + a2 + t2 ~ 1,
         nl = TRUE)
ms2 <- lapply(1:15, function(i) c(runif(1, 10, 100), runif(1, 0.1, 2),
                                  runif(1, 1, 30), runif(1, 2, 20)))
report("N2 Gamma(identity) biexponential", f2, d2, Gamma("identity"),
       c(50, 0.5, 10, 5), ms2)

## N3 Gompertz, Gamma(identity), located prior with a bound on k
set.seed(203)
x3 <- runif(150, 0, 30)
mu3 <- 100 * exp(-5 * exp(-0.3 * x3))
d3 <- data.frame(x = x3, y = rgamma(150, shape = 30, rate = 30 / mu3))
f3 <- bf(y ~ A * exp(-B * exp(-k * x)), A + B + k ~ 1, nl = TRUE)
ms3 <- lapply(1:15, function(i) c(runif(1, 20, 300), runif(1, 0.5, 20),
                                  runif(1, 0.05, 1)))
report("N3 Gamma(identity) Gompertz, prior normal(0.5, 1) on k, lb 0",
       f3, d3, Gamma("identity"), c(100, 5, 0.3), ms3,
       prior = prior(normal(0.5, 1), nlpar = "k", lb = 0))

## N4 Michaelis-Menten, poisson(identity), Vm ~ 1 + (1 | g)
set.seed(204)
G <- 25
g <- factor(rep(seq_len(G), each = 16))
u <- rnorm(G, 0, 3)
x4 <- rep(c(0.5, 1, 2, 4, 8, 12, 16, 20), length.out = G * 16)
mu4 <- (30 + u[g]) * x4 / (2 + x4)
d4 <- data.frame(x = x4, g = g, y = rpois(G * 16, mu4))
f4 <- bf(y ~ Vm * x / (K + x), Vm ~ 1 + (1 | g), K ~ 1, nl = TRUE)
report("N4 poisson(identity) Michaelis-Menten, Vm ~ 1 + (1 | g)", f4, d4,
       poisson("identity"), c(30, 2), list(c(10, 1), c(50, 5), c(20, 0.5)))

## N5 power law with a negative exponent, Gamma(identity)
set.seed(205)
x5 <- runif(150, 0.1, 10)
d5 <- data.frame(x = x5, y = rgamma(150, shape = 10,
                                    rate = 10 / (3 * x5^-0.7)))
f5 <- bf(y ~ a * x^b, a + b ~ 1, nl = TRUE)
report("N5 Gamma(identity) a x^b, b = -0.7", f5, d5, Gamma("identity"),
       c(3, -0.7), lapply(1:10, function(i) c(runif(1, 0.5, 10),
                                              runif(1, -2, 1))))

## N6 tiny scale, poisson(identity) a exp(b x), a = 1e-3, b = 2
set.seed(206)
x6 <- runif(200, 0, 5)
d6 <- data.frame(x = x6, y = rpois(200, 1e-3 * exp(2 * x6)))
f6 <- bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE)
report("N6 poisson(identity) a exp(b x), a = 1e-3, b = 2", f6, d6,
       poisson("identity"), c(1e-3, 2),
       lapply(1:10, function(i) c(10^runif(1, -4, -1), runif(1, 0.5, 3))))

## Stress: the sinusoid over periods and seeds. Count the wrong fits that
## pass the rule: no error, conv 0, every SE finite, gap < -0.01.
cat("\n== stress: sinusoid, P in 1.3, 3, 4, 7, 9 x seeds 1:6\n")
res <- NULL
for (P0 in c(1.3, 3, 4, 7, 9)) for (s in 1:6) {
  set.seed(1000 + s)
  xs <- seq(0, 15, length.out = 120)
  ds <- data.frame(x = xs, y = 2 * sin(6.283185307179586 * xs / P0) +
                     rnorm(120, 0, 0.6))
  r <- run(frm(f1, data = ds))
  starts <- c(list(c(2, P0)), lapply(1:12, function(i) {
    c(runif(1, 0.5, 4), runif(1, 0.5, 12))
  }))
  best <- max(vapply(starts, function(st) {
    ll(run(frm(f1, data = ds, start = list(beta = st))))
  }, 0), na.rm = TRUE)
  out <- if (inherits(r$v, "error")) "error" else {
    gap <- ll(r) - best
    if (r$v$opt$convergence == 0 && isTRUE(sefin(r)) && gap < -0.01) {
      "WRONG-PASSED"
    } else if (gap >= -0.01) "rescued" else "wrong-flagged"
  }
  gapv <- if (inherits(r$v, "error")) NA else ll(r) - best
  res <- rbind(res, data.frame(P = P0, seed = 1000 + s, out = out,
                               gap = gapv, nwarn = length(unique(r$warns))))
}
print(res)
print(table(res$out))

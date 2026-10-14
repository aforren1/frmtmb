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
    m <- if (grepl("A candidate to try: start = ", em, fixed = TRUE)) trimws(sub("^.*A candidate to try: ", "", em)) else character(0)
    if (length(m)) {
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

## The candidate start of a refused scan, passed back (M2, M4, M7, M9 of
## dev/postfix-rev-nlscan.R and dev/postfix-rev-nlscan7.R, same seeds)
set.seed(102)
x2 <- rep(c(0.02, 0.06, 0.11, 0.22, 0.56, 1.1), each = 8)
mu2 <- 200 * x2 / (0.1 + x2)
d2 <- data.frame(x = x2, y = rgamma(length(x2), shape = 50, rate = 50 / mu2))
f2 <- bf(y ~ Vm * x / (K + x), Vm + K ~ 1, nl = TRUE)
report("M2 Michaelis-Menten", f2, d2, Gamma("identity"), c(200, 0.1),
       lapply(1:10, function(i) c(runif(1, 20, 500), runif(1, 0.01, 1))))
set.seed(104)
x4 <- runif(120, 0, 1000)
mu4 <- 2e5 * exp(-x4 / 300)
d4 <- data.frame(x = x4, y = rgamma(120, shape = 20, rate = 20 / mu4))
f4 <- bf(y ~ a * exp(-x / tau), a + tau ~ 1, nl = TRUE)
report("M4 a exp(-x / tau)", f4, d4, Gamma("identity"), c(2e5, 300),
       lapply(1:10, function(i) c(10^runif(1, 4, 6), 10^runif(1, 1.5, 3.5))))
set.seed(107)
x7 <- seq(0, 10, length.out = 120)
d7 <- data.frame(x = x7, y = 3 * sin(2 * pi * x7 / 2.5) + rnorm(120, 0, 0.5))
f7 <- bf(y ~ A * sin(6.283185307179586 * x / P), A + P ~ 1, nl = TRUE)
report("M7 sinusoid", f7, d7, gaussian(), c(3, 2.5),
       lapply(1:20, function(i) c(runif(1, 0.5, 5), runif(1, 0.5, 10))))
set.seed(109)
x9 <- runif(300, 0, 3)
d9 <- data.frame(x = x9)
d9$y <- rgamma(300, shape = 10, rate = 10 / (1 + 2 + 3 * x9 + 0.5 + 1 + 0.2))
f9 <- bf(y ~ p1 + p2 + p3 * x + p4 + p5 + p6, p1 + p2 + p3 + p4 + p5 + p6 ~ 1,
         nl = TRUE)
report("M9 ridge", f9, d9, Gamma("identity"), c(1, 2, 3, 0.5, 1, 0.2),
       list(c(1, 1, 1, 1, 1, 1)))

# Reviewer of lane postfix: nl_intercept_scan() outside the lane's six
# models. For each model: the default fit (with the scan's objective
# evaluations counted and timed), the fit from a start placed by hand
# near the truth, and the best of a multistart, to see whether the scan
# lands on a worse local optimum.
# Usage: Rscript dev/postfix-rev-nlscan.R <lib or "base">
a <- commandArgs(TRUE)
lib <- if (length(a) && a[1] != "base") a[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
scan_n <- 0L; scan_t <- 0
if (exists("nl_intercept_scan", ns)) {
  trace("nl_intercept_scan", where = ns, print = FALSE,
        tracer = quote({
          .t0 <- Sys.time()
          obj$fn <- local({
            g <- obj$fn
            function(...) {
              assign("scan_n", get("scan_n", globalenv()) + 1L,
                     envir = globalenv())
              g(...)
            }
          })
        }),
        exit = quote(assign("scan_t", get("scan_t", globalenv()) +
                              as.numeric(Sys.time() - .t0, units = "secs"),
                            envir = globalenv())))
}
run <- function(expr) {
  msgs <- character(); warns <- character()
  t0 <- Sys.time()
  v <- tryCatch(withCallingHandlers(expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
  list(v = v, msgs = msgs, warns = warns,
       t = as.numeric(Sys.time() - t0, units = "secs"))
}
ll <- function(r) {
  if (inherits(r$v, "error")) return(NA_real_)
  as.numeric(logLik(r$v))
}
report <- function(lbl, f, data, family, hand, multi, prior = NULL) {
  scan_n <<- 0L; scan_t <<- 0
  d0 <- run(frm(f, data = data, family = family, prior = prior))
  sn <- scan_n; st <- scan_t
  h <- run(frm(f, data = data, family = family, prior = prior,
               start = list(beta = hand)))
  ms <- vapply(multi, function(s) {
    ll(run(frm(f, data = data, family = family, prior = prior,
               start = list(beta = s))))
  }, 0)
  best <- max(c(ms, ll(h)), na.rm = TRUE)
  cat(sprintf("\n== %s\n", lbl))
  if (inherits(d0$v, "error")) {
    cat("default: ERROR", substr(conditionMessage(d0$v), 1, 160), "\n")
  } else {
    cat(sprintf("default: conv %d  logLik %.6f  est %s\n",
                d0$v$opt$convergence, ll(d0),
                paste(signif(fixef(d0$v)[, "Estimate"], 6), collapse = " ")))
  }
  cat(sprintf("  time %.2fs; scan fn evals %d, scan time %.2fs\n",
              d0$t, sn, st))
  for (m in grep("scan", d0$msgs, value = TRUE)) cat("  [msg]", substr(m, 1, 300), "\n")
  for (w in unique(d0$warns)) cat("  [warn]", substr(w, 1, 160), "\n")
  cat(sprintf("hand: logLik %.6f (time %.2fs); multistart best %.6f of %d finite\n",
              ll(h), h$t, best, sum(is.finite(ms))))
  cat(sprintf("default - best = %.4e\n", ll(d0) - best))
}

set.seed(101)
## M1 poisson identity, a * exp(b x) with b negative
x <- runif(150, 0, 5)
d1 <- data.frame(x = x, y = rpois(150, 20 * exp(-0.5 * x)))
f1 <- bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE)
ms1 <- lapply(1:10, function(i) c(runif(1, 1, 50), runif(1, -2, 0.5)))
report("M1 poisson(identity) a exp(b x), b = -0.5", f1, d1,
       poisson("identity"), c(20, -0.5), ms1)

## M1b same, b bounded above by a located prior (b placed at 0)
report("M1b with prior normal(0, 5) on b, ub = 0", f1, d1,
       poisson("identity"), c(20, -0.5), ms1[1:4],
       prior = prior(normal(0, 5), nlpar = "b", ub = 0))
## M1c bound on a, location at the bound
report("M1c with prior normal(0, 100) on a, lb = 0", f1, d1,
       poisson("identity"), c(20, -0.5), ms1[1:4],
       prior = prior(normal(0, 100), nlpar = "a", lb = 0))

## M2 Michaelis-Menten, Gamma identity, Vm = 200, K = 0.1
set.seed(102)
x2 <- rep(c(0.02, 0.06, 0.11, 0.22, 0.56, 1.1), each = 8)
mu2 <- 200 * x2 / (0.1 + x2)
d2 <- data.frame(x = x2, y = rgamma(length(x2), shape = 50, rate = 50 / mu2))
f2 <- bf(y ~ Vm * x / (K + x), Vm + K ~ 1, nl = TRUE)
ms2 <- lapply(1:10, function(i) c(runif(1, 20, 500), runif(1, 0.01, 1)))
report("M2 Gamma(identity) Michaelis-Menten", f2, d2, Gamma("identity"),
       c(200, 0.1), ms2)

## M4 large scale: a = 2e5, tau = 300
set.seed(104)
x4 <- runif(120, 0, 1000)
mu4 <- 2e5 * exp(-x4 / 300)
d4 <- data.frame(x = x4, y = rgamma(120, shape = 20, rate = 20 / mu4))
f4 <- bf(y ~ a * exp(-x / tau), a + tau ~ 1, nl = TRUE)
ms4 <- lapply(1:10, function(i) c(10^runif(1, 4, 6), 10^runif(1, 1.5, 3.5)))
report("M4 Gamma(identity) a exp(-x / tau), a = 2e5, tau = 300", f4, d4,
       Gamma("identity"), c(2e5, 300), ms4)

## M6 Emax with negative emax and x = 0 in the data
set.seed(106)
x6 <- rep(c(0, 0.5, 1, 2, 4, 8, 16), each = 10)
d6 <- data.frame(x = x6, y = 10 - 8 * x6 / (2 + x6) + rnorm(70, 0, 0.5))
f6 <- bf(y ~ e0 + emax * x / (ed50 + x), e0 + emax + ed50 ~ 1, nl = TRUE)
ms6 <- lapply(1:10, function(i) c(runif(1, 0, 20), runif(1, -20, 20),
                                  runif(1, 0.1, 10)))
report("M6 gaussian Emax, emax = -8, ed50 = 2, x = 0 present", f6, d6,
       gaussian(), c(10, -8, 2), ms6)

## M7 multimodal: sinusoid with period in the denominator
set.seed(107)
x7 <- seq(0, 10, length.out = 120)
d7 <- data.frame(x = x7, y = 3 * sin(2 * pi * x7 / 2.5) + rnorm(120, 0, 0.5))
f7 <- bf(y ~ A * sin(2 * pi * x / P), A + P ~ 1, nl = TRUE)
ms7 <- lapply(1:20, function(i) c(runif(1, 0.5, 5), runif(1, 0.5, 10)))
report("M7 gaussian sinusoid, P = 2.5 (multimodal in P)", f7, d7,
       gaussian(), c(3, 2.5), ms7)

## M8 nonlinear with a random intercept, Gamma identity: scan cost with
## a Laplace inner problem
set.seed(108)
G <- 30
g <- factor(rep(seq_len(G), each = 20))
u <- rnorm(G, 0, 0.5)
x8 <- runif(G * 20, 0, 2)
mu8 <- 5 + u[g] + 2 * x8
d8 <- data.frame(x = x8, g = g, y = rgamma(G * 20, shape = 10, rate = 10 / mu8))
f8 <- bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE)
report("M8 Gamma(identity) a + b x, a ~ 1 + (1 | g)", f8, d8,
       Gamma("identity"), c(5, 2), list(c(3, 1), c(8, 3)))
## M9 many nonlinear intercepts: six parameters
set.seed(109)
x9 <- runif(300, 0, 3)
d9 <- data.frame(x = x9)
d9$y <- rgamma(300, shape = 10, rate = 10 / (1 + 2 + 3 * x9 + 0.5 + 1 + 0.2))
f9 <- bf(y ~ p1 + p2 + p3 * x + p4 + p5 + p6, p1 + p2 + p3 + p4 + p5 + p6 ~ 1,
         nl = TRUE)
report("M9 Gamma(identity), six nonlinear intercepts (ridge)", f9, d9,
       Gamma("identity"), c(1, 2, 3, 0.5, 1, 0.2), list(c(1, 1, 1, 1, 1, 1)))

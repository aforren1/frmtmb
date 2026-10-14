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
## M7 multimodal: sinusoid with period in the denominator
set.seed(107)
x7 <- seq(0, 10, length.out = 120)
d7 <- data.frame(x = x7, y = 3 * sin(2 * pi * x7 / 2.5) + rnorm(120, 0, 0.5))
f7 <- bf(y ~ A * sin(6.283185307179586 * x / P), A + P ~ 1, nl = TRUE)
ms7 <- lapply(1:20, function(i) c(runif(1, 0.5, 5), runif(1, 0.5, 10)))
report("M7 gaussian sinusoid, P = 2.5 (multimodal in P)", f7, d7,
       gaussian(), c(3, 2.5), ms7)


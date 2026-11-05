# Profile of the latent-residual inner problem, planned before the build was
# trusted: what grows with n, and whether the chain coordinates pay.
#
# Counts first, which do not depend on load: the inner dimension, the stored
# entries of the inner Hessian, and the tape's size (values and operators).
# Then clocks, by the lane rules: arms interleaved in one process, blocks grown
# past 1.2 s, the minimum of several rounds, and a CONTROL built from the same
# code that must read 1.0.
#
# Arms. The shipped chains, then two references: `dense_toep` on toep(ft + 0 | g),
# a dense k x k density per group with k - 1 free lags, what brms's coordinates
# cost; and `ar1_re` on ar1(ft + 0 | g), the same AR(1) latent model through
# frmtmb's own sparse random-effect structure, so that its fitted logLik against
# ar(t, g, cov = TRUE) is an identity worth printing. Seed 167.
args <- commandArgs(TRUE)
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

mk <- function(G, k) {
  set.seed(167)
  d <- expand.grid(t = seq_len(k), g = factor(seq_len(G)))
  n <- nrow(d)
  d$x <- rnorm(n)
  e <- as.vector(replicate(G, arima.sim(list(ar = 0.6), k, sd = 0.5)))
  d$y <- rpois(n, exp(0.5 + 0.3 * d$x + e))
  d$ft <- factor(d$t)
  d
}
# nodes on the tape, the load-independent measure of what an evaluation costs
tape_size <- function(obj) {
  f <- obj$env$ADFun
  v <- tryCatch(RTMB:::getTapeInfo(f), error = function(e) NULL)
  if (is.null(v)) NA_character_ else paste(names(v), v, collapse = " ")
}
block_time <- function(f, target = 1.2) {
  reps <- 1L
  repeat {
    t0 <- proc.time()[[3]]
    for (i in seq_len(reps)) f()
    el <- proc.time()[[3]] - t0
    if (el >= target) return(el / reps)
    reps <- reps * 2L
  }
}
report <- function(lab, fit, t_fit) {
  force(t_fit)
  obj <- fit$obj
  par <- fit$opt$par
  H <- obj$env$spHess(obj$env$last.par.best, random = TRUE)
  fn <- function() obj$fn(par)
  gr <- function() obj$gr(par)
  tf <- min(replicate(3L, block_time(fn)))
  tg <- min(replicate(3L, block_time(gr)))
  ctl <- min(replicate(3L, block_time(fn)))
  cat(sprintf(paste0("%-10s inner %6d nnz %9d fit %7.2fs fn %8.2fms gr %8.2fms ",
                     "ctrl/fn %.2f logLik %.6f\n"),
              lab, nrow(H), length(H@x), t_fit, tf * 1e3, tg * 1e3, ctl / tf,
              as.numeric(logLik(fit))))
}
shapes <- list(c(200L, 5L), c(50L, 20L), c(1600L, 5L), c(400L, 20L),
               c(20L, 100L), c(4L, 500L), c(1L, 2000L))
for (s in shapes) {
  G <- s[[1]]
  k <- s[[2]]
  d <- mk(G, k)
  cat(sprintf("\n=== G = %d, k = %d, n = %d\n", G, k, nrow(d)))
  for (tm in c("ar(t, g, cov = TRUE)", "ma(t, g, cov = TRUE)",
               "arma(t, g, cov = TRUE)", "ar(t, g, p = 2)", "cosy(t, g)")) {
    t0 <- proc.time()[[3]]
    fit <- tryCatch(suppressWarnings(frm(bf(as.formula(paste("y ~ x +", tm))),
                                         family = poisson(), data = d)),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      cat(tm, "ERROR", conditionMessage(fit), "\n")
      next
    }
    report(substr(gsub("[(, =)]", "", tm), 1L, 10L), fit,
           proc.time()[[3]] - t0)
  }
  if (k <= 100L) {
    t0 <- proc.time()[[3]]
    fit <- tryCatch(suppressWarnings(frm(bf(y ~ x + toep(ft + 0 | g)),
                                         family = poisson(), data = d)),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      cat("toep ERROR", conditionMessage(fit), "
")
    } else {
      report("dense_toep", fit, proc.time()[[3]] - t0)
    }
    t0 <- proc.time()[[3]]
    fit <- tryCatch(suppressWarnings(frm(bf(y ~ x + ar1(ft + 0 | g)),
                                         family = poisson(), data = d)),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      cat("ar1_re ERROR", conditionMessage(fit), "\n")
    } else {
      report("ar1_re", fit, proc.time()[[3]] - t0)
    }
  } else {
    cat("references skipped above k = 100: a dense stored Hessian is",
        format(G * k * k), "entries\n")
  }
}

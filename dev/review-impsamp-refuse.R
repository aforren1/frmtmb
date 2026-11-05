# Item 4: what the refusal costs a user at its worst, and what its message says.
# (a) the wrapper from the lane's own test, displacing every round, so no count
# settles: rows n just under imp_refuse_cells / 8192, so the ladder reaches 8192.
# (b) natural firing: a covariate on a large scale and shift.
# (c) message truth: top < 2 * draws, one doubling tried.
# Usage: Rscript dev/review-impsamp-refuse.R a|b|c
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
peak_mb <- function() {
  cmd <- sprintf("(Get-Process -Id %d).PeakWorkingSet64 / 1MB", Sys.getpid())
  as.numeric(system2("powershell", c("-NoProfile", "-Command", shQuote(cmd)),
                     stdout = TRUE))
}
displaced <- function(below, code) {
  real_plan <- frmtmb:::imp_plan
  real_opt <- frmtmb:::optimize_obj
  last_n <- 0L
  testthat::with_mocked_bindings(
    code,
    imp_plan = function(lap_obj, frame, lay, par, n_draw, seed) {
      last_n <<- n_draw
      real_plan(lap_obj, frame, lay, par, n_draw, seed)
    },
    optimize_obj = function(obj, control, bounds = list(lower = -Inf, upper = Inf),
                            par_units = NULL, verbose = 0L, start_par = obj$par,
                            tally = NULL) {
      o <- real_opt(obj, control, bounds, par_units, verbose, start_par, tally)
      if (last_n > 0L && last_n < below) o$par[[1L]] <- o$par[[1L]] + 1
      o
    },
    .package = "frmtmb")
}
mode <- commandArgs(TRUE)[[1]]
mk <- function(G, n, seed = 1L) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), each = n))
  x <- rnorm(G * n)
  d <- data.frame(y = rpois(G * n, exp(0.3 + 0.4 * x + rnorm(G, 0, 0.5)[g])), x, g)
  d
}
if (mode == "a") {
  d <- mk(100L, 24L)
  cat("rows", nrow(d), "top", min(8192, floor(2e7 / nrow(d))), "\n")
  f <- bf(y ~ x + (1 | g))
  t0 <- system.time(ok <- frm(f, family = poisson(), data = d, importance = 512L))
  cat("the ordinary fit at 512:", t0[["elapsed"]], "s, rounds", ok$importance$rounds,
      " peak Mb", peak_mb(), "\n")
  t1 <- system.time(msg <- tryCatch(displaced(1e9, frm(f, family = poisson(), data = d,
                                                       importance = 512L)),
                                    error = conditionMessage))[["elapsed"]]
  cat("refusal after", t1, "s, ", t1 / t0[["elapsed"]], "x the fit; peak Mb", peak_mb(), "\n")
  cat(msg, "\n")
}
if (mode == "b") {
  set.seed(1)
  g <- factor(rep(seq_len(15L), each = 10))
  x <- rnorm(150)
  y <- rbinom(150, 1, plogis(0.2 + 0.7 * x + rnorm(15, 0, 0.8)[g]))
  d <- data.frame(y, xx = 100 * x + 1e4, g, xc = x)
  for (N in c(64L, 300L, 1200L)) {
    t <- system.time(m <- tryCatch(suppressWarnings(frm(bf(y ~ xx + (1 | g)),
      family = bernoulli(), data = d, importance = N)), error = conditionMessage))
    cat("N", N, t[["elapsed"]], "s:", if (is.character(m)) m else
      paste("fit, rounds", m$importance$rounds, "grad", m$importance$grad,
            "moves", paste(format(m$importance$moves, digits = 3), collapse = " ")), "\n")
  }
  t <- system.time(cc <- frm(bf(y ~ xc + (1 | g)), family = bernoulli(), data = d,
                             importance = 300L))
  cat("the same data, column centered and scaled, N 300:", t[["elapsed"]], "s rounds",
      cc$importance$rounds, "\n")
  cat("Laplace's own fit at xx:\n")
  print(suppressWarnings(frm(bf(y ~ xx + (1 | g)), family = bernoulli(), data = d))$opt$par)
}
if (mode == "c") {
  # n large enough that one doubling is all the memory bound allows
  d <- mk(200L, 100L)
  f <- bf(y ~ x + (1 | g))
  cat("rows", nrow(d), "top", min(8192, floor(2e7 / nrow(d))), "\n")
  t <- system.time(msg <- tryCatch(displaced(1e9, frm(f, family = poisson(), data = d,
                                                      importance = 500L)),
                                   error = conditionMessage))[["elapsed"]]
  cat("draws 500:", t, "s peak", peak_mb(), "\n", msg, "\n\n")
  t <- system.time(msg <- tryCatch(displaced(1e9, frm(f, family = poisson(), data = d,
                                                      importance = 1500L)),
                                   error = conditionMessage))[["elapsed"]]
  cat("draws 1500:", t, "s peak", peak_mb(), "\n", msg, "\n")
}

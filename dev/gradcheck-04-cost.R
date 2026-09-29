# What a curvature-aware criterion COSTS, measured against the fit that
# would pay for it, and against one restart, which is what the current
# criterion buys with the same gradient.
#
#   Rscript dev/gradcheck-04-cost.R base

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "base"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

# proc.time() ticks at 10 ms here, so every arm is a block grown past
# 1.5 s and reported per call, with a minimum over 3 rounds.
block <- function(fun, secs = 1.5, rounds = 3L) {
  best <- Inf
  for (k in seq_len(rounds)) {
    t1 <- proc.time()[["elapsed"]]
    reps <- 0L
    repeat {
      fun()
      reps <- reps + 1L
      if (proc.time()[["elapsed"]] - t1 > secs) break
    }
    best <- min(best, (proc.time()[["elapsed"]] - t1) / reps)
  }
  best
}

probe <- function(lbl, fit, refit) {
  p <- fit$opt$par
  np <- length(p)
  t_gr <- block(function() fit$obj$gr(p))
  he_ok <- !inherits(tryCatch(fit$obj$he(p), error = function(e) e),
                     "error")
  t_he <- if (he_ok) block(function() fit$obj$he(p)) else NA_real_
  t_oh <- block(function() {
    stats::optimHess(p, function(q) fit$obj$fn(q),
                     function(q) drop(fit$obj$gr(q)))
  })
  # one restart is what the CURRENT criterion spends when it trips
  t_fit <- block(function() suppressWarnings(refit()), secs = 3)
  cat(sprintf("%-30s np %3d  gr %9.5f  he %9.5f  optimHess %9.5f  fit %9.4f\n",
              lbl, np, t_gr, t_he, t_oh, t_fit))
  cat(sprintf("%-30s     optimHess / gr %7.1f   optimHess / whole fit %7.4f\n",
              "", t_oh / t_gr, t_oh / t_fit))
  invisible(NULL)
}

set.seed(301)
n <- 500
d1 <- data.frame(x = rnorm(n), z = rnorm(n))
d1$y <- rnorm(n, 3 + 2 * d1$x - 1.5 * d1$z, 1)
probe("gaussian GLM n=500 np=4",
      frm(bf(y ~ x + z), family = gaussian(), data = d1),
      function() frm(bf(y ~ x + z), family = gaussian(), data = d1))

set.seed(302)
n <- 20000
d2 <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
e2 <- 0.8 * d2$x1 - 0.5 * d2$x2
d2$yo <- cut(e2 + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
             labels = FALSE)
probe("cumulative n=20000 np=5",
      frm(bf(yo ~ x1 + x2), family = cumulative(), data = d2),
      function() frm(bf(yo ~ x1 + x2), family = cumulative(), data = d2))

set.seed(303)
n <- 2000
d3 <- data.frame(x = rnorm(n), g = factor(rep(1:50, each = 40)))
d3$y <- rpois(n, exp(0.4 + 0.5 * d3$x + rnorm(50, 0, 0.6)[d3$g]))
probe("poisson GLMM n=2000 q=50",
      frm(bf(y ~ x + (1 | g)), family = poisson(), data = d3),
      function() frm(bf(y ~ x + (1 | g)), family = poisson(), data = d3))

# the case a Hessian on the warning path could actually hurt: many
# fixed effects AND a Laplace inner problem
set.seed(304)
n <- 8000
d4 <- data.frame(g = factor(rep(1:100, each = 80)))
d4$fx <- factor(rep_len(seq_len(60), n))
d4$x <- rnorm(n)
d4$y <- rnorm(n, 1 + 0.4 * d4$x + rnorm(100, 0, 0.6)[d4$g] +
               rnorm(60, 0, 0.5)[d4$fx], 1)
probe("gaussian LMM np=63",
      frm(bf(y ~ x + fx + (1 | g)), family = gaussian(), data = d4),
      function() frm(bf(y ~ x + fx + (1 | g)), family = gaussian(),
                     data = d4))

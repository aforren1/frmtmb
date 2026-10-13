# Reviewer of lane seflags, claim 3's cost on large ordinary identified
# models. Per design: the fit with check_se = "ignore" and "warning",
# interleaved, minimum of `reps`; the check's own objective and gradient
# evaluations on a fit made with the check off (load-independent); the
# fit's own evaluation count for scale; the number of
# se_ridge_curvature() calls.
#   Rscript dev/seflags-rev-cost.R <lane|base> <design> [reps]
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-common.R")
a <- commandArgs(trailingOnly = TRUE)
dz <- a[2]
reps <- if (length(a) > 2) as.integer(a[3]) else 2L
n_ridge <- 0L
if (exists("se_ridge_curvature", ns)) {
  trace("se_ridge_curvature", quote(n_ridge <<- n_ridge + 1L),
        print = FALSE, where = ns)
}
set.seed(1)
s <- switch(dz,
  glmm_1e5 = {
    n <- 1e5
    d <- data.frame(x = rnorm(n), z = rnorm(n),
                    g = factor(sample(500, n, TRUE)))
    d$y <- rbinom(n, 1, plogis(-0.5 + 0.5 * d$x + 0.3 * d$z +
                                 rnorm(500, 0, 0.7)[d$g]))
    list(f = y ~ x + z + (1 | g), fam = bernoulli(), d = d)
  },
  glmm_1e5_year = {
    n <- 1e5
    d <- data.frame(year = 2000 + sample(0:20, n, TRUE), z = rnorm(n),
                    g = factor(sample(500, n, TRUE)))
    d$y <- rbinom(n, 1, plogis(-0.5 + 0.05 * (d$year - 2010) + 0.3 * d$z +
                                 rnorm(500, 0, 0.7)[d$g]))
    list(f = y ~ year + z + (1 | g), fam = bernoulli(), d = d)
  },
  smooth_mid = {
    n <- 5000
    d <- data.frame(x = runif(n), g = factor(sample(50, n, TRUE)))
    d$y <- sin(8 * d$x) + rnorm(50, 0, 0.3)[d$g] + rnorm(n, 0, 0.5)
    list(f = y ~ s(x, k = 40) + (1 | g), fam = gaussian(), d = d)
  },
  smooth_big = {
    n <- 2e4
    d <- data.frame(x = runif(n), g = factor(sample(100, n, TRUE)))
    d$y <- sin(8 * d$x) + rnorm(100, 0, 0.3)[d$g] + rnorm(n, 0, 0.5)
    list(f = y ~ s(x, k = 100) + (1 | g), fam = gaussian(), d = d)
  },
  lmm_400 = {
    n <- 1e4
    d <- data.frame(f = factor(rep(1:400, each = 25)),
                    g = factor(sample(50, n, TRUE)), x = rnorm(n))
    d$y <- rnorm(400)[d$f] + 0.3 * d$x + rnorm(50, 0, 0.5)[d$g] + rnorm(n)
    list(f = y ~ f + x + (1 | g), fam = gaussian(), d = d)
  },
  glm_400 = {
    n <- 2e4
    d <- data.frame(f = factor(rep(1:400, each = 50)), x = rnorm(n))
    d$y <- rpois(n, exp(0.5 + rnorm(400, 0, 0.3)[d$f] + 0.3 * d$x))
    list(f = y ~ f + x, fam = poisson(), d = d)
  })
fitit <- function(ch) suppressMessages(suppressWarnings(
  frm(s$f, family = s$fam, data = s$d,
      control = frmtmb_control(check_se = ch))))
t <- matrix(NA_real_, reps, 2)
for (r in seq_len(reps)) {
  t[r, 1] <- system.time(f1 <- fitit("ignore"))[["elapsed"]]
  n_ridge <- 0L
  t[r, 2] <- system.time(f2 <- fitit("warning"))[["elapsed"]]
}
nr_fit <- n_ridge
f <- fitit("ignore")
n_fn <- 0L; n_gr <- 0L
fn0 <- f$obj$fn; gr0 <- f$obj$gr
f$obj$fn <- function(...) { n_fn <<- n_fn + 1L; fn0(...) }
f$obj$gr <- function(...) { n_gr <<- n_gr + 1L; gr0(...) }
n_ridge <- 0L
tc <- system.time(invisible(suppressWarnings(suppressMessages(
  ns$se_check(f, frmtmb_control(check_se = "warning"))))))[["elapsed"]]
cat(sprintf(paste0("%-14s %-4s ignore %.2f s | warning %.2f s | ratio %.3f |",
                   " se_check alone %.2f s | check evals fn %d gr %d |",
                   " fit evals %s | ridge calls %d | lost %s\n"),
            dz, rev_arm, min(t[, 1]), min(t[, 2]),
            min(t[, 2]) / min(t[, 1]), tc, n_fn, n_gr,
            paste(f$opt$evals, collapse = "/"), n_ridge, lost_of(f2)))
cat("times ignore:", t[, 1], "| warning:", t[, 2], "\n")

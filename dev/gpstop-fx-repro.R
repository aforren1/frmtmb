# Lane gpstop, fix round: the review's rows (dev/reviews/2026-11-11-gpstop.md)
# on the installed build, arms lane and norelin in one process, through the
# reviewer's helpers (dev/reviews/gpstop-review-lib.R). Each row also prints
# the inner Hessian's smallest eigenvalue at the fit's optimum. Data:
# gpstop-review-sweep2.R mk(), set.seed(81100 + 1000 fam + 100 seed + n +
# 10 body). Args: rows as body:fam:n:seed (comma).
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
rv_instrument()
bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
fams <- c("gaussian", "poisson", "bernoulli", "student")
spell <- function(body) {
  switch(body,
    exp_a = bf(y ~ a * exp(b * z), a ~ gp(x), b ~ 1, nl = TRUE),
    exp_b = bf(y ~ a * exp(b * z), a ~ 1, b ~ gp(x), nl = TRUE),
    logis_a = bf(y ~ a / (1 + exp(-b * z)), a ~ gp(x), b ~ 1, nl = TRUE),
    mm_V = bf(y ~ V * w / (K + w), V ~ gp(x), K ~ 1, nl = TRUE),
    mm_K = bf(y ~ V * w / (K + w), V ~ 1, K ~ gp(x), nl = TRUE),
    pow_a = bf(y ~ a * w^b, a ~ gp(x), b ~ 1, nl = TRUE),
    pow_b = bf(y ~ a * w^b, a ~ 1, b ~ gp(x), nl = TRUE))
}
mk <- function(body, fam, n, seed) {
  set.seed(81100 + 1000 * match(fam, fams) + 100 * seed + n +
             10 * match(body, bodies))
  x <- sort(round(runif(n, 0, 10), 3))
  z <- round(rnorm(n), 3)
  w <- round(runif(n, 0.2, 5), 3)
  f <- if (seed %% 2) sin(x) else 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  eta <- switch(body,
    exp_a = (2 + fs) * exp(0.4 * z),
    exp_b = 2 * exp((0.3 + 0.2 * fs) * z),
    logis_a = (3 + fs) / (1 + exp(-1.5 * z)),
    mm_V = (3 + fs) * w / (1 + w),
    mm_K = 3 * w / (1 + 0.6 * fs + w),
    pow_a = (2 + fs) * w^0.5,
    pow_b = w^(0.5 + 0.3 * fs))
  y <- switch(fam,
    gaussian = eta + 0.1 * rnorm(n),
    student = eta + 0.1 * rt(n, 3),
    poisson = rpois(n, exp(0.6 * eta)),
    bernoulli = rbinom(n, 1, plogis(1.5 * (eta - median(eta)) / sd(eta))))
  data.frame(x, z, w, y)
}
famf <- function(fam) switch(fam, gaussian = gaussian(), poisson = poisson(),
                             student = student(), bernoulli = bernoulli())
inner_ev <- function(fit) {
  env <- fit$obj$env
  fit$obj$fn(fit$opt$par)
  H <- as.matrix(env$spHess(env$last.par, random = TRUE))
  min(eigen(0.5 * (H + t(H)), symmetric = TRUE, only.values = TRUE)$values)
}
rows <- strsplit(a[1], ",")[[1]]
for (rw in rows) {
  p <- strsplit(rw, ":")[[1]]
  d <- mk(p[1], p[2], as.integer(p[3]), as.integer(p[4]))
  for (arm in c("lane", "norelin")) {
    q <- rv_run(rw, arm, spell(p[1]), data = d, family = famf(p[2]))
    if (!inherits(q$fit, "error")) {
      cat(sprintf("INNER %s %s ev_min %.3e | warnings: %s\n", rw, arm,
                  inner_ev(q$fit), paste(q$warn, collapse = " || ")))
    }
  }
}
rv_arm("lane")
cat("DONE\n")

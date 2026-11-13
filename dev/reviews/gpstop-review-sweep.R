# Review of lane gpstop: fresh nonlinear bodies with an exact (and some
# approximate) gp() in a body parameter, arms lane and norelin interleaved
# in one process, order alternating by row. Every row prints ll, code,
# warnings, max |grad|, the outer Hessian's smallest eigenvalue, theta.
# Args: families (comma list of gaussian, poisson, bernoulli), ns (comma),
#   seeds (comma), extra ("approx", "big" or "none").
# Data: set.seed(81100 + 1000 * fam index + 100 * seed + n + 10 * body index)
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
rv_instrument()
fams <- strsplit(a[1], ",")[[1]]
ns <- as.integer(strsplit(a[2], ",")[[1]])
seeds <- as.integer(strsplit(a[3], ",")[[1]])
extra <- a[4]
bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
spell <- function(body, k = NULL) {
  g <- if (is.null(k)) quote(gp(x)) else bquote(gp(x, k = .(k)))
  switch(body,
    exp_a = bf(y ~ a * exp(b * z), as.formula(bquote(a ~ .(g))), b ~ 1,
               nl = TRUE),
    exp_b = bf(y ~ a * exp(b * z), a ~ 1, as.formula(bquote(b ~ .(g))),
               nl = TRUE),
    logis_a = bf(y ~ a / (1 + exp(-b * z)), as.formula(bquote(a ~ .(g))),
                 b ~ 1, nl = TRUE),
    mm_V = bf(y ~ V * w / (K + w), as.formula(bquote(V ~ .(g))), K ~ 1,
              nl = TRUE),
    mm_K = bf(y ~ V * w / (K + w), V ~ 1, as.formula(bquote(K ~ .(g))),
              nl = TRUE),
    pow_a = bf(y ~ a * w^b, as.formula(bquote(a ~ .(g))), b ~ 1, nl = TRUE),
    pow_b = bf(y ~ a * w^b, a ~ 1, as.formula(bquote(b ~ .(g))), nl = TRUE))
}
mk <- function(body, fam, n, seed) {
  set.seed(81100 + 1000 * match(fam, c("gaussian", "poisson", "bernoulli")) +
             100 * seed + n + 10 * match(body, bodies))
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
    poisson = rpois(n, exp(0.6 * eta)),
    bernoulli = rbinom(n, 1, plogis(1.5 * (eta - median(eta)) / sd(eta))))
  data.frame(x, z, w, y)
}
famf <- function(fam) switch(fam, gaussian = gaussian(), poisson = poisson(),
                             bernoulli = bernoulli())
rows <- list()
for (fam in fams) for (n in ns) for (seed in seeds) for (body in bodies) {
  rows[[length(rows) + 1L]] <- list(body = body, fam = fam, n = n,
                                    seed = seed, k = NULL)
}
if (extra == "approx") {
  for (body in bodies) rows[[length(rows) + 1L]] <-
    list(body = body, fam = "gaussian", n = 120L, seed = 2L, k = 12L)
}
i <- 0L
for (r in rows) {
  i <- i + 1L
  d <- mk(r$body, r$fam, r$n, r$seed)
  tag <- sprintf("%s_%s_n%d_s%d%s", r$body, r$fam, r$n, r$seed,
                 if (is.null(r$k)) "" else paste0("_k", r$k))
  arms <- if (i %% 2) c("lane", "norelin") else c("norelin", "lane")
  for (arm in arms) {
    rv_run(tag, arm, spell(r$body, r$k), data = d, family = famf(r$fam))
  }
}
cat("DONE\n")

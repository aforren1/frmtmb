# Lane illse2, item 3: the relative asymmetry of the Hessian in the
# design's coordinates, max |H - H'| / 2 over max |H|, at the fit's
# point, on the designs the route takes (seeds 1 to 3): gaussian,
# bernoulli, poisson and cumulative raw polynomials of degree 8 (7 for
# cumulative), with and without (1 | g), and the mixed-scale binomial of
# the illse review. The guard that refuses the route above a bound must
# sit far above all of these.
#   ILLSE2_ARM=lane Rscript dev/illse2-asym2.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
asym_of <- function(fit) {
  p <- fit$opt$par
  M <- ns$reparam_basis(fit$frame, names(p), est = fit$estimates)
  if (is.null(M)) return(NA_real_)
  h <- ns$reparam_hessian(fit$obj, p, M)
  if (is.null(h)) return(NA_real_)
  max(h$E) / max(abs(h$H))
}
poly_d <- function(seed, n, re, deg) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(1:20, length.out = n)))
  for (k in seq_len(deg)) d[[paste0("x", k)]] <- d$x^k
  d$eta <- sin(3 * d$x) + if (re) stats::rnorm(20, 0, 0.5)[d$g] else 0
  d
}
fo_of <- function(deg, re) {
  stats::as.formula(paste("y ~", paste0("x", seq_len(deg), collapse = " + "),
                          if (re) "+ (1 | g)"))
}
out <- list()
for (seed in 1:3) for (re in c(FALSE, TRUE)) {
  d <- poly_d(seed, 300, re, 8)
  d$y <- d$eta + stats::rnorm(300, 0, 0.3)
  out[[length(out) + 1]] <- c("gaussian", seed, re,
    asym_of(suppressWarnings(frm(fo_of(8, re), data = d))))
  d$y <- stats::rbinom(300, 1, stats::plogis(d$eta))
  out[[length(out) + 1]] <- c("bernoulli", seed, re,
    asym_of(suppressWarnings(frm(fo_of(8, re), data = d,
                                 family = bernoulli()))))
  d$y <- stats::rpois(300, exp(d$eta))
  out[[length(out) + 1]] <- c("poisson", seed, re,
    asym_of(suppressWarnings(frm(fo_of(8, re), data = d,
                                 family = poisson()))))
  d$y <- cut(d$eta + stats::rlogis(300), c(-Inf, -0.5, 0.5, 1.5, Inf),
             labels = FALSE)
  out[[length(out) + 1]] <- c("cumulative", seed, re,
    asym_of(suppressWarnings(frm(fo_of(7, re), data = d,
                                 family = cumulative()))))
}
for (seed in 201:203) {
  set.seed(seed)
  n <- 300
  d <- data.frame(u = 100 + stats::runif(n), z = 1e-3 * stats::rnorm(n),
                  w = 1e4 * stats::rnorm(n))
  d$u2 <- d$u^2
  d$N <- 20L
  eta <- -0.3 + 3 * (d$u - 100.5) - 2 * (d$u - 100.5)^2 + 300 * d$z +
    5e-5 * d$w
  d$y <- stats::rbinom(n, d$N, stats::plogis(eta))
  out[[length(out) + 1]] <- c("binomial mixed", seed, FALSE,
    asym_of(suppressWarnings(frm(y | trials(N) ~ u + u2 + z + w, data = d,
                                 family = binomial()))))
}
for (o in out) {
  cat(sprintf("%-15s seed %s re %-5s | asym %s\n", o[1], o[2], o[3],
              format(as.numeric(o[4]), digits = 3)))
}
v <- as.numeric(vapply(out, `[`, "", 4))
cat(sprintf("largest %.2e over %d fits\n", max(v, na.rm = TRUE), sum(!is.na(v))))

# Review of lane illse2, item 2: the difference scheme of the Hessian in
# the design's coordinates on non-gaussian fits, where a central
# difference along a design column is not exact. Raw polynomial on
# [1, 2] beside two badly scaled columns (sd 1e4 and 1e-4), no random
# effects, so the exact covariance at the fit's point is
# (X' W X)^-1 formed and inverted in 256 bits.
#   ILLSE2_ARM=lane|base I2_SCHEME=E|A Rscript dev/illse2-rev-hess.R \
#     <poisson|bernoulli> <deg> <n> <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
fam <- a[1]
deg <- as.integer(a[2])
n <- as.numeric(a[3])
seeds <- eval(parse(text = a[4]))
scheme <- Sys.getenv("I2_SCHEME", "E")
if (arm == "lane" && scheme == "A") {
  # 0.71.0's scheme on the lane's code: Richardson on every column
  rh <- ns$reparam_hessian
  environment(rh) <- ns
  wrap <- function(obj, p, M, step = 1e-3, rich_all = FALSE) {
    rh(obj, p, M, step = step, rich_all = TRUE)
  }
  environment(wrap) <- environment()
  utils::assignInNamespace("reparam_hessian", wrap, "frmtmb")
}
cat(sprintf("scheme %s | %s deg %d n %g\n", if (arm == "base") "A (0.71.0)"
            else scheme, fam, deg, n))
for (seed in seeds) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, 1, 2), w = stats::rnorm(n) * 1e4,
                  v = stats::rnorm(n) * 1e-4)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  eta <- sin(3 * d$x) + 0.2 * d$w / 1e4 - 0.2 * d$v / 1e-4
  d$y <- switch(fam,
                bernoulli = stats::rbinom(n, 1, stats::plogis(eta)),
                poisson = stats::rpois(n, exp(eta)))
  rhs <- paste(c(xn, "w", "v"), collapse = " + ")
  fo <- stats::as.formula(paste("y ~", rhs))
  famo <- switch(fam, bernoulli = bernoulli(), poisson = poisson())
  t0 <- proc.time()[["elapsed"]]
  r <- cap(frm(fo, data = d, family = famo))
  tt <- proc.time()[["elapsed"]] - t0
  if (is_err(r$v)) { cat(seed, "ERROR", r$v, "\n"); next }
  f <- r$v
  b <- f$opt$par
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
  e <- drop(X %*% b)
  w <- switch(fam, poisson = exp(e), bernoulli = stats::plogis(e) *
                (1 - stats::plogis(e)))
  Xm <- mM(X)
  H <- crossprod(Xm, Xm * m(w))
  Ve <- mp_inv(H)
  V <- suppressWarnings(vcov(f))
  dg <- sqrt(num(diag(Ve)))
  E <- if (all(is.finite(V))) max(abs(num(mM(V) - Ve)) / outer(dg, dg)) else
    NaN
  sr <- if (all(is.finite(V))) sqrt(diag(V)) / dg else NaN
  # the fit's shortfall: the exact Newton gain at its point
  g <- drop(crossprod(X, switch(fam, poisson = exp(e) - d$y,
                                bernoulli = stats::plogis(e) - d$y)))
  gm <- col1(m(g))
  gain <- num(sum(gm * (Ve %*% gm))) / 2
  ev <- eigen(cov2cor(num(Ve)), only.values = TRUE)$values
  cat(sprintf(paste0("seed %d | code %d | kappa %.1e | lost %s | root %s | ",
                     "cor %.2e | SE dev %.2e | gain %.1e | %.1f s | warn %s\n"),
              seed, f$opt$convergence, max(ev) / min(ev), lost_of(f),
              !is.null(ns$sdr_of(f)$cov_fixed_root), E,
              max(abs(sr - 1)), gain, tt, wtag(r$w)))
}

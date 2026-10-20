# Review of lane adhess: copy of dev/illse2-acc.R (lane illse2), REV_ARM arms, plus the SE error.
# coordinates under a given difference scheme. Every entry of
# vcov(full = TRUE) against the exact covariance at the fit's own point
# (256-bit Hessian, inverted in 256 bits), largest error over
# sqrt(Ve_ii Ve_jj), and the fit's shortfall (Newton gain at its point).
#   ILLSE2_ARM=lane I2_STEP=1e-3 I2_RICH=1 Rscript dev/illse2-acc.R \
#     <family gaussian|bernoulli|poisson> <re 0|1> <degs> <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-i2common.R")
options(frmtmb.i2_step = as.numeric(Sys.getenv("I2_STEP", "1e-3")))
ri <- Sys.getenv("I2_RICH", "1")
options(frmtmb.i2_rich = if (ri == "single") "single" else ri == "1")
a <- commandArgs(TRUE)
fam <- a[1]
re <- as.integer(a[2]) == 1L
degs <- eval(parse(text = a[3]))
seeds <- eval(parse(text = a[4]))
cat(sprintf("scheme: step %g, richardson %s\n", getOption("frmtmb.i2_step"),
            getOption("frmtmb.i2_rich")))
for (deg in degs) for (seed in seeds) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(1:20, length.out = n)))
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  u <- if (re) stats::rnorm(20, 0, 0.5)[d$g] else 0
  eta <- sin(3 * d$x) + u
  d$y <- switch(fam,
                gaussian = eta + stats::rnorm(n, 0, 0.3),
                bernoulli = stats::rbinom(n, 1, stats::plogis(eta)),
                poisson = stats::rpois(n, exp(eta)))
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                                if (re) "+ (1 | g)"))
  famo <- switch(fam, gaussian = gaussian(), bernoulli = bernoulli(),
                 poisson = poisson())
  t0 <- proc.time()[["elapsed"]]
  r <- cap(frm(fo, data = d, family = famo))
  tt <- proc.time()[["elapsed"]] - t0
  if (is_err(r$v)) { cat(deg, seed, "ERROR", r$v, "\n"); next }
  f <- r$v
  p <- f$opt$par
  X <- stats::model.matrix(stats::as.formula(paste("~",
                                    paste(xn, collapse = " + "))), d)
  Xm <- mM(X); ym <- m(d$y)
  nll <- if (re) {
    if (fam != "gaussian") stop("exact reference with (1 | g): gaussian")
    lmm_nll(X, d$y, d$g)
  } else switch(fam,
    gaussian = gauss_nll(X, d$y),
    bernoulli = function(th) {
      e <- Xm %*% col1(th)
      sum(log(1 + exp(e))) - sum(ym * e)
    },
    poisson = function(th) {
      e <- Xm %*% col1(th)
      sum(exp(e)) - sum(ym * e)
    })
  He <- mp_hessian(nll, p)
  Ve <- mp_inv(He)
  ge <- mp_grad(nll, p)
  gain <- num(sum(ge * (Ve %*% col1(ge)))) / 2
  V <- suppressWarnings(vcov(f, full = TRUE))
  dg <- sqrt(num(diag(Ve)))
  E <- if (all(is.finite(V))) {
    max(abs(num(mM(V) - Ve)) / outer(dg, dg))
  } else NaN
  Ese <- max(abs(sqrt(pmax(diag(V), 0)) / dg - 1))
  cat(sprintf(paste0("%s re %d deg %d seed %d | code %d | lost %s | ",
                     "root %s | se %.2e | cor %.2e | gain %.1e | %.1f s | warn %s\n"),
              fam, re, deg, seed, f$opt$convergence, lost_of(f),
              !is.null(ns$sdr_of(f)$cov_fixed_root), Ese, E, gain, tt,
              wtag(r$w)))
}

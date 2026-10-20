# Review of lane adhess, item 2: on illse2's raw polynomials without
# random effects, how accurate the plain inverse of the AD Hessian is
# (what adhess's route would give if the QR route did not take the fit),
# against the 256-bit covariance at the fit's own point.
#   REV_ARM=lane Rscript dev/adhess-rev-adinv.R <family> <degs> <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-i2common.R")
a <- commandArgs(TRUE)
fam <- a[1]
degs <- eval(parse(text = a[2]))
seeds <- eval(parse(text = a[3]))
for (deg in degs) for (seed in seeds) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2))
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  eta <- sin(3 * d$x)
  d$y <- switch(fam, gaussian = eta + stats::rnorm(n, 0, 0.3),
                poisson = stats::rpois(n, exp(eta)))
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + ")))
  f <- suppressWarnings(frm(fo, data = d, family = switch(fam,
    gaussian = gaussian(), poisson = poisson())))
  p <- f$opt$par
  X <- stats::model.matrix(stats::as.formula(paste("~",
                                    paste(xn, collapse = " + "))), d)
  Xm <- mM(X); ym <- m(d$y)
  nll <- switch(fam, gaussian = gauss_nll(X, d$y),
    poisson = function(th) {
      e <- Xm %*% col1(th)
      sum(exp(e)) - sum(ym * e)
    })
  Ve <- mp_inv(mp_hessian(nll, p))
  dg <- sqrt(num(diag(Ve)))
  st <- ns$obj_state_save(f$obj)
  Hx <- f$obj$he(p)
  ns$obj_state_restore(f$obj, st)
  Hs <- (Hx + t(Hx)) / 2
  D <- sqrt(diag(Hs))
  ev <- eigen(Hs / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  Va <- tryCatch(solve(Hs), error = function(e) NULL)
  ea <- if (is.null(Va)) c(NA, NA) else
    c(max(abs(sqrt(pmax(diag(Va), 0)) / dg - 1)),
      max(abs(num(mM(Va) - Ve)) / outer(dg, dg)))
  Vf <- suppressWarnings(vcov(f, full = TRUE))
  ef <- max(abs(num(mM(Vf) - Ve)) / outer(dg, dg))
  cat(sprintf(paste0("%s deg %d seed %d | kappa(unit diag) %.1e | ",
                     "solve(AD H): se %.2e cor %.2e | vcov (route) cor %.2e\n"),
              fam, deg, seed, max(ev) / min(ev), ea[1], ea[2], ef))
}

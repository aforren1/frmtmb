# Lane illse2, item 3: nonlinear formulas on ill-conditioned designs.
# Case "poly": the degree-8 raw polynomial written as an nl body, each
# coefficient its own parameter (test-se-reparam.R); the reference is
# lm()'s fit, the same likelihood. Case "rat": a degree-6 raw polynomial
# over (1 + dd x), nonlinear in dd. Case "exp": y = a exp(b x) with x on
# c0 + [0, 1]. The reference is the exact Hessian at the fit's own
# point in 256 bits (and nls() for the log-likelihood).
#   ILLSE2_ARM=lane|base Rscript dev/illse2-p3nl.R <case> <seeds> [c0]
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
case <- a[1]
seeds <- eval(parse(text = a[2]))
c0 <- if (length(a) > 2) as.numeric(a[3]) else 20
for (seed in seeds) {
  set.seed(seed)
  start <- NULL
  if (case == "poly") {
    d <- data.frame(x = stats::runif(200, 1, 2))
    d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
    f <- bf(y ~ b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 + b5 * x^5 +
              b6 * x^6 + b7 * x^7 + b8 * x^8,
            b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 ~ 1, nl = TRUE)
    X <- outer(d$x, 0:8, `^`)
    nll <- gauss_nll(X, d$y)
    lmf <- stats::lm(d$y ~ X - 1)
    ref_ll <- as.numeric(stats::logLik(lmf))
  } else if (case == "rat") {
    d <- data.frame(x = stats::runif(200, 1, 2))
    d$y <- sin(3 * d$x) / (1 + 0.3 * d$x) + stats::rnorm(200, 0, 0.3)
    f <- bf(y ~ (b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 +
                   b5 * x^5 + b6 * x^6) / (1 + dd * x),
            b0 + b1 + b2 + b3 + b4 + b5 + b6 + dd ~ 1, nl = TRUE)
    X <- outer(d$x, 0:6, `^`)
    Xm <- mM(X); xm <- m(d$x); ym <- m(d$y)
    nll <- function(th) {
      mu <- (Xm %*% col1(th[1:7])) / (1 + th[8] * xm)
      r <- ym - mu
      s2 <- exp(2 * th[9])
      (length(d$y) * log(s2) + sum(r^2) / s2) / 2
    }
    b_lm <- unname(stats::coef(stats::lm(d$y ~ X - 1)))
    nf <- tryCatch(stats::nls(y ~ (X %*% b) / (1 + dd * x),
                              data = c(as.list(d), list(X = X)),
                              start = list(b = b_lm, dd = 0)),
                   error = function(e) e)
    ref_ll <- if (inherits(nf, "error")) NA else
      as.numeric(stats::logLik(nf))
  } else if (case == "sin") {
    # a degree-6 raw polynomial in x beside A sin(w z): the nonlinear
    # part is identified, the polynomial is ill-conditioned
    d <- data.frame(x = stats::runif(300, 1, 2), z = stats::runif(300, 0, 10))
    d$y <- sin(3 * d$x) + 1.5 * sin(1.3 * d$z) + stats::rnorm(300, 0, 0.3)
    f <- bf(y ~ b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 + b5 * x^5 +
              b6 * x^6 + A * sin(w * z),
            b0 + b1 + b2 + b3 + b4 + b5 + b6 + A + w ~ 1, nl = TRUE)
    X <- outer(d$x, 0:6, `^`)
    Xm <- mM(X); zm <- m(d$z); ym <- m(d$y)
    nll <- function(th) {
      mu <- (Xm %*% col1(th[1:7])) + th[8] * sin(th[9] * zm)
      r <- ym - mu
      s2 <- exp(2 * th[10])
      (length(d$y) * log(s2) + sum(r^2) / s2) / 2
    }
    b_lm <- unname(stats::coef(stats::lm(d$y ~ X - 1)))
    nf <- tryCatch(stats::nls(y ~ X %*% b + A * sin(w * z),
                              data = c(as.list(d), list(X = X)),
                              start = list(b = b_lm, A = 1, w = 1.3)),
                   error = function(e) e)
    ref_ll <- if (inherits(nf, "error")) NA else
      as.numeric(stats::logLik(nf))
    start <- list(beta = c(b_lm, 1, 1.3))
  } else {
    d <- data.frame(x = c0 + stats::runif(200))
    d$y <- 2 * exp(0.8 * (d$x - c0)) + stats::rnorm(200, 0, 0.3)
    f <- bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE)
    xm <- m(d$x); ym <- m(d$y)
    nll <- function(th) {
      r <- ym - th[1] * exp(th[2] * xm)
      s2 <- exp(2 * th[3])
      (length(d$y) * log(s2) + sum(r^2) / s2) / 2
    }
    st <- list(a = 2 * exp(-0.8 * c0), b = 0.8)
    nf <- tryCatch(stats::nls(y ~ a * exp(b * x), data = d, start = st),
                   error = function(e) e)
    ref_ll <- if (inherits(nf, "error")) NA else
      as.numeric(stats::logLik(nf))
    start <- list(beta = c(2 * exp(-0.8 * c0), 0.8))
  }
  r <- cap(if (is.null(start)) frm(f, data = d) else
    frm(f, data = d, start = start))
  if (is_err(r$v)) { cat("seed", seed, "ERROR", r$v, "\n"); next }
  fit <- r$v
  p <- fit$opt$par
  # the objective and the written-out nll agree up to a constant
  pp <- p
  pp[length(p)] <- pp[length(p)] + 0.01
  pp[1] <- pp[1] + 1e-4
  cdiff <- (as.numeric(fit$obj$fn(p)) - num(nll(m(p)))) -
    (as.numeric(fit$obj$fn(pp)) - num(nll(m(pp))))
  fit$obj$fn(p)
  He <- mp_hessian(nll, p)
  Ve <- mp_inv(He)
  V <- suppressWarnings(vcov(fit, full = TRUE))
  dg <- sqrt(num(diag(Ve)))
  E <- if (all(is.finite(V))) {
    max(abs(num(mM(V) - Ve)) / outer(dg, dg))
  } else NaN
  # kappa at unit diagonal from the exact Hessian: the largest eigenvalue
  # of the scaled matrix and of its exact inverse, both well conditioned
  D <- sqrt(num(diag(He)))
  Sd <- num(He) / outer(D, D)
  Si <- num(Ve) * outer(D, D)
  kap <- max(eigen(Sd, symmetric = TRUE, only.values = TRUE)$values) *
    max(eigen(Si, symmetric = TRUE, only.values = TRUE)$values)
  cat(sprintf(paste0("%s seed %d | code %d | logLik - ref %+.2e | ",
                     "kappa %.2e | lost %s | root %s | cor err %.2e | ",
                     "SE/exact %s | warn %s | const %.1e\n"),
              case, seed, fit$opt$convergence,
              as.numeric(logLik(fit)) - ref_ll, kap, lost_of(fit),
              !is.null(ns$sdr_of(fit)$cov_fixed_root), E,
              ratio_str(sqrt(diag(V)) / dg), wtag(r$w), cdiff))
}

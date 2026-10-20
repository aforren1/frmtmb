# Lane illse2, item 3: REML on a raw polynomial with (1 | g), gaussian.
#
# Under REML beta is integrated with the group effects, so frmtmb's
# coefficient covariance is the beta block of the joint covariance its
# joint precision carries: the conditional covariance at the variance
# parameters, sigma^2 (X' V^-1 X)^-1 (lme4's vcov()), plus what the
# variance parameters' uncertainty adds, J V_theta J', with J the
# derivative of the GLS estimate in (log sigma, log sd) and V_theta the
# inverse Hessian of the restricted log-likelihood there. Both parts are
# computed exactly in 256 bits at frmtmb's own REML estimates. lme4's
# REML vcov() is compared with the conditional part at lme4's estimates.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-p3reml.R <deg> <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
deg <- as.integer(a[1])
seeds <- eval(parse(text = a[2]))
xn <- paste0("x", seq_len(deg))
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                              "+ (1 | g)"))
fo0 <- stats::as.formula(paste("~", paste(xn, collapse = " + ")))
# log determinant of an mpfr matrix by elimination with partial pivoting
mp_logdet <- function(A) {
  p <- nrow(A)
  out <- m(0)
  for (k in seq_len(p)) {
    piv <- k - 1L + which.max(num(abs(A[k:p, k])))
    if (piv != k) {
      t1 <- A[k, ]; A[k, ] <- A[piv, ]; A[piv, ] <- t1
    }
    out <- out + log(abs(A[k, k]))
    if (k < p) for (i in (k + 1L):p) {
      A[i, ] <- A[i, ] - (A[i, k] / A[k, k]) * A[k, ]
    }
  }
  out
}
reml_parts <- function(X, y, g) {
  Xm <- mM(X); ym <- m(y)
  XX <- crossprod(Xm); Xy <- crossprod(Xm, ym); yy <- sum(ym^2)
  S <- lapply(split(seq_along(y), g), function(ix) {
    list(s = col1(colSums(Xm[ix, , drop = FALSE])), sy = sum(ym[ix]),
         n = length(ix))
  })
  n <- length(y)
  # th = (log sigma, log sd): A = X' V^-1 X sigma^2, bvec, c0 = y' V^-1 y
  # sigma^2, and log|V|
  pieces <- function(th) {
    s2 <- exp(2 * th[1]); t2 <- exp(2 * th[2])
    A <- XX; bv <- Xy; c0 <- yy; ld <- m(0)
    for (st in S) {
      cg <- t2 / (s2 + st$n * t2)
      A <- A - cg * (st$s %*% t(st$s))
      bv <- bv - cg * st$sy * st$s
      c0 <- c0 - cg * st$sy^2
      ld <- ld + (st$n - 1) * log(s2) + log(s2 + st$n * t2)
    }
    list(A = A, bv = bv, c0 = c0, ld = ld, s2 = s2)
  }
  beta_hat <- function(th) {
    pc <- pieces(th)
    mp_inv(pc$A) %*% pc$bv
  }
  # negative restricted log-likelihood, up to a constant
  nreml <- function(th) {
    pc <- pieces(th)
    Ai <- mp_inv(pc$A)
    q <- pc$c0 - sum(pc$bv * (Ai %*% pc$bv))
    p <- ncol(X)
    (pc$ld + mp_logdet(pc$A) - p * log(pc$s2) + q / pc$s2) / 2
  }
  cond_cov <- function(th) {
    pc <- pieces(th)
    pc$s2 * mp_inv(pc$A)
  }
  list(beta_hat = beta_hat, nreml = nreml, cond_cov = cond_cov)
}
for (seed in seeds) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.5)[d$g] +
    stats::rnorm(200, 0, 0.3)
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  X <- stats::model.matrix(fo0, d)
  rp <- reml_parts(X, d$y, d$g)
  r <- cap(frm(fo, data = d, REML = TRUE))
  if (is_err(r$v)) { cat("seed", seed, "ERROR", r$v, "\n"); next }
  f <- r$v
  th <- unname(f$opt$par)
  stopifnot(identical(names(f$opt$par), c("betad", "theta")))
  # the REML maximum: frmtmb's stop against the exact optimum in th
  opt <- stats::optim(th, function(t) num(rp$nreml(m(t))), method = "BFGS",
                      control = list(reltol = 1e-14))
  short <- num(rp$nreml(m(th))) - opt$value
  Vc <- rp$cond_cov(m(th))
  Ht <- mp_hessian(rp$nreml, th)
  Vt <- mp_inv(Ht)
  # J: d beta_hat / d th by central differences in 256 bits
  h <- 1e-20
  J <- mM(matrix(0, ncol(X), 2))
  for (j in 1:2) {
    e <- m(c(0, 0)); e[j] <- m(h)
    J[, j] <- (rp$beta_hat(m(th) + e) - rp$beta_hat(m(th) - e)) / (2 * h)
  }
  Ve <- Vc + J %*% Vt %*% t(J)
  vc <- cap(vcov(f))
  V <- vc$v
  dg <- sqrt(num(diag(Ve)))
  cat(sprintf(paste0("seed %d deg %d | code %d | REML short %.1e | lost %s | ",
                     "whitened %s | warn %s | vcov warn %s\n"),
              seed, deg, f$opt$convergence, short, lost_of(f),
              !is.null(f$beta_whiten), wtag(r$w), wtag(vc$w)))
  if (all(is.finite(V))) {
    E <- abs(num(mM(V) - Ve)) / outer(dg, dg)
    cat(sprintf("  frmtmb vs exact (its definition): SE %s | cor %.2e\n",
                ratio_str(sqrt(diag(V)) / dg), max(E)))
    cat(sprintf("  share of the variance theta adds: %.1e to %.1e\n",
                min(1 - num(diag(Vc)) / num(diag(Ve))),
                max(1 - num(diag(Vc)) / num(diag(Ve)))))
  } else cat("  frmtmb SEs:", paste(signif(sqrt(diag(V)), 3), collapse = " "),
             "\n")
  # fitted values at x = 1, 1.5, 2 through the factor
  nd <- data.frame(x = c(1, 1.5, 2), g = d$g[1])
  for (k in seq_len(deg)) nd[[xn[k]]] <- nd$x^k
  A <- stats::model.matrix(fo0, nd)
  fp <- cap(fitted(f, newdata = nd, re_formula = NA))
  vpe <- num(rowSums((mM(A) %*% Ve) * mM(A)))
  if (!is_err(fp$v)) {
    cat(sprintf("  fitted(re_formula = NA) SE / exact %s | warn %s\n",
                ratio_str(fp$v[, "Est.Error"] / sqrt(vpe)), wtag(fp$w)))
  }
  l <- cap(lme4::lmer(fo, data = d, REML = TRUE))
  if (!is_err(l$v)) {
    Vl <- as.matrix(stats::vcov(l$v))
    thl <- log(c(sigma(l$v), attr(lme4::VarCorr(l$v)$g, "stddev")))
    Vcl <- rp$cond_cov(m(thl))
    cat(sprintf(paste0("  lme4: REML short %.1e | SE / exact conditional ",
                       "at its point %s | warn %s\n"),
                num(rp$nreml(m(thl))) - opt$value,
                ratio_str(sqrt(diag(Vl)) / sqrt(num(diag(Vcl)))),
                wtag(l$w)))
  } else cat("  lme4 ERROR", substr(l$v, 1, 80), "\n")
}

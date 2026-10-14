# Reviewer of lane illse, attack 1: a binomial with mixed-scale
# covariates (a covariate near 1000 and its square, one on 1e-3, one on
# 1e4) and a cumulative logit with a degree-7 raw polynomial. Exact
# covariance at the fit's own point (256-bit Hessian of the
# log-likelihood written out here), every vcov() entry, and the
# consumers against the exact delta-method answers.
#   REV_ARM=lane|base Rscript dev/illse-rev-glm.R <binom|cum> <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
a <- commandArgs(TRUE)
fam <- a[1]
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 201:203
plogis_m <- function(z) 1 / (1 + exp(-z))
for (seed in seeds) {
  set.seed(seed)
  n <- 300
  if (fam == "binom") {
    d <- data.frame(u = 100 + stats::runif(n), z = 1e-3 * stats::rnorm(n),
                    w = 1e4 * stats::rnorm(n))
    d$u2 <- d$u^2
    d$N <- 20L
    eta <- -0.3 + 3 * (d$u - 100.5) - 2 * (d$u - 100.5)^2 +
      300 * d$z + 5e-5 * d$w
    d$y <- stats::rbinom(n, d$N, stats::plogis(eta))
    fo <- y | trials(N) ~ u + u2 + z + w
    fo0 <- ~ u + u2 + z + w
    r <- cap(frm(fo, data = d, family = binomial()))
    nd <- data.frame(u = 100 + c(0.1, 0.5, 0.9), z = c(-1e-3, 0, 1e-3),
                     w = c(-1e4, 0, 1e4), N = 20L)
    nd$u2 <- nd$u^2
  } else {
    d <- data.frame(x = stats::runif(n, 1, 2))
    d$y <- factor(cut(2 * sin(3 * d$x) + stats::rlogis(n),
                      c(-Inf, -1, 0, 1, Inf), labels = FALSE),
                  ordered = TRUE)
    for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
    fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7
    fo0 <- ~ 0 + x1 + x2 + x3 + x4 + x5 + x6 + x7
    r <- cap(frm(fo, data = d, family = cumulative()))
    nd <- data.frame(x = c(1.2, 1.5, 1.8))
    for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
  }
  if (inherits(r$v, "rev_err")) { cat("seed", seed, "ERROR", r$v, "\n"); next }
  fit <- r$v
  p <- fit$opt$par
  pn <- names(p)
  pb <- which(pn == "beta")
  X <- stats::model.matrix(fo0, d)
  Xf <- as.matrix(fit$frame$linpreds[[1]]$X)
  if (!identical(dim(X), dim(Xf))) cat("  frame X columns:", colnames(Xf), "
")
  stopifnot(identical(dim(X), dim(Xf)), max(abs(X - Xf)) == 0)
  Xm <- mM(X)
  if (fam == "binom") {
    yy <- m(d$y); NN <- m(d$N)
    nll <- function(th) {
      e <- Xm %*% mpfr2array(th[pb], dim = c(length(pb), 1L))
      -sum(yy * e - NN * log(1 + exp(e)))
    }
  } else {
    pt <- which(pn == "tau_raw")
    K <- nlevels(d$y)
    yi <- as.integer(d$y)
    nll <- function(th) {
      e <- as.vector(Xm %*% mpfr2array(th[pb], dim = c(length(pb), 1L)))
      raw <- th[pt]
      tau <- c(raw[1], raw[1] + cumsum(exp(raw[-1])))
      out <- m(0)
      for (k in seq_len(K)) {
        ix <- which(yi == k)
        if (!length(ix)) next
        hi <- if (k < K) plogis_m(tau[k] - e[ix]) else m(rep(1, length(ix)))
        lo <- if (k > 1) plogis_m(tau[k - 1] - e[ix]) else m(rep(0, length(ix)))
        out <- out - sum(log(hi - lo))
      }
      out
    }
  }
  c0 <- as.numeric(fit$obj$fn(p)) - num(nll(m(p)))
  pp <- p + 1e-6 * abs(p)
  c1 <- as.numeric(fit$obj$fn(pp)) - num(nll(m(pp)))
  fit$obj$fn(p)
  H <- mp_hessian(nll, p)
  Ve <- mp_inv(H)
  gx <- mp_grad(nll, p)
  ng <- num(sum(gx * (Ve %*% mpfr2array(gx, dim = c(length(p), 1L))))) / 2
  V <- cap(vcov(fit, full = TRUE))
  e <- cov_err(V$v, Ve)
  D <- sqrt(abs(num(diag(H))))
  kap <- {
    Hs <- num(H) / outer(D, D)
    ev <- eigen(Hs, symmetric = TRUE, only.values = TRUE)$values
    max(ev) / min(ev)
  }
  cat(sprintf(paste0("%s seed %d | code %d | fit warn %s | lost %d | ",
                     "obj const diff %.1e | Newton gain %.1e | kappa %.1e |",
                     " root %s lossy %s | vcov warn %s\n"),
              fam, seed, fit$opt$convergence, wtag(r$w),
              length(ns$sdr_of(fit)$se_lost), abs(c0 - c1), ng, kap,
              !is.null(ns$sdr_of(fit)$cov_fixed_root),
              isTRUE(attr(ns$sdr_of(fit)$cov_fixed_root, "lossy")),
              wtag(V$w)))
  cat(sprintf("  vcov(full) every entry: cor %.2e | rel %.2e | diag %.2e\n",
              e[["cor"]], e[["rel"]], e[["diag"]]))
  se_ex <- sqrt(num(diag(Ve)))
  out <- list()
  s <- cap(summary(fit))
  if (!inherits(s$v, "rev_err")) {
    out$summary_fixed <- s$v$fixed[, "Est.Error"] /
      se_ex[c(if (fam == "cum") pt, pb)][seq_len(nrow(s$v$fixed))]
    if (fam == "cum") {
      # thresholds are tau, not tau_raw: delta through the map
      Jt <- function(th) {
        raw <- th[pt]; c(raw[1], raw[1] + cumsum(exp(raw[-1])))
      }
      vt <- vapply(seq_along(pt), function(k) {
        gk <- mp_grad(function(th) Jt(th)[k], p)
        num(sum(gk * (Ve %*% mpfr2array(gk, dim = c(length(p), 1L)))))
      }, 0)
      out$summary_fixed <- s$v$fixed[, "Est.Error"] /
        c(sqrt(vt), se_ex[pb])
    }
  }
  out$summary_w <- wtag(s$w)
  ci <- cap(confint(fit))
  if (!inherits(ci$v, "rev_err")) {
    hw <- (ci$v[, 2] - ci$v[, 1]) / 2
    ref <- if (fam == "cum") c(sqrt(vt), se_ex[pb]) else se_ex[pb]
    out$confint <- hw[seq_along(ref)] / (stats::qnorm(0.975) * ref)
  }
  out$confint_w <- wtag(ci$w)
  A <- stats::model.matrix(fo0, nd)
  v_eta <- mp_rowvar(A, Ve[pb, pb])
  lp <- cap(frm_linpred(fit, newdata = nd, se.fit = TRUE))
  if (!inherits(lp$v, "rev_err")) out$linpred <- lp$v$se.fit / sqrt(v_eta)
  out$linpred_w <- wtag(lp$w)
  fl <- cap(fitted(fit, newdata = nd, scale = "linear"))
  if (!inherits(fl$v, "rev_err")) out$fitted_linear <- fl$v[, "Est.Error"] / sqrt(v_eta)
  out$fitted_linear_w <- wtag(fl$w)
  # the response scale
  fr <- cap(fitted(fit, newdata = nd))
  if (fam == "binom") {
    pr <- function(th, i) {
      ei <- sum(m(A[i, ]) * th[pb]); 20 * plogis_m(ei)
    }
    rows <- seq_len(nrow(nd))
    vr <- vapply(rows, function(i) {
      gi <- mp_grad(function(th) pr(th, i), p)
      num(sum(gi * (Ve %*% mpfr2array(gi, dim = c(length(p), 1L)))))
    }, 0)
    if (!inherits(fr$v, "rev_err")) out$fitted_resp <- fr$v[, "Est.Error"] / sqrt(vr)
  } else {
    Am <- mM(A)
    pk <- function(th, i, k) {
      ei <- sum(Am[i, ] * th[pb])
      raw <- th[pt]; tau <- c(raw[1], raw[1] + cumsum(exp(raw[-1])))
      hi <- if (k < 4) plogis_m(tau[k] - ei) else m(1)
      lo <- if (k > 1) plogis_m(tau[k - 1] - ei) else m(0)
      hi - lo
    }
    vr <- matrix(NA_real_, nrow(nd), 4)
    for (i in seq_len(nrow(nd))) for (k in 1:4) {
      gi <- mp_grad(function(th) pk(th, i, k), p)
      vr[i, k] <- num(sum(gi * (Ve %*% mpfr2array(gi, dim = c(length(p), 1L)))))
    }
    if (!inherits(fr$v, "rev_err")) out$fitted_cat <- fr$v[, "Est.Error", ] / sqrt(vr)
  }
  out$fitted_resp_w <- wtag(fr$w)
  hs <- if (fam == "binom") {
    paste0("Intercept + ", 100.5, " * u + ", 100.5^2, " * u2 = 0")
  } else {
    paste0(paste0(1.5^(1:7), " * x", 1:7, collapse = " + "), " = 0")
  }
  ha <- if (fam == "binom") c(1, 100.5, 100.5^2, 0, 0) else 1.5^(1:7)
  vh <- mp_rowvar(matrix(ha, 1L), Ve[pb, pb])
  h <- cap(hypothesis(fit, hs, class = "b"))
  if (!inherits(h$v, "rev_err")) out$hypothesis <- h$v$hypothesis$Est.Error / sqrt(vh)
  out$hypothesis_w <- wtag(h$w)
  if (fam == "binom") {
    em <- cap(summary(emmeans::emmeans(fit, ~ u + u2 + z + w,
                                       at = as.list(nd[2, c("u", "u2", "z", "w")]))))
    if (!inherits(em$v, "rev_err")) out$emmeans_link <- em$v$SE / sqrt(v_eta[2])
    out$emmeans_w <- wtag(em$w)
  }
  mf <- cap(marginaleffects::predictions(fit, newdata = nd))
  if (!inherits(mf$v, "rev_err")) {
    est_f <- if (fam == "binom") fr$v[, "Estimate"] else as.vector(fr$v[, "Estimate", ])
    se_f <- if (fam == "binom") fr$v[, "Estimate"] * 0 + sqrt(vr) else as.vector(sqrt(vr))
    est_m <- mf$v$estimate
    j <- vapply(est_m, function(v) which.min(abs(est_f - v)), 1L)
    cat("  mfx est match rel", sprintf("%.1e", max(abs(est_m / est_f[j] - 1))), "
")
    out$marginaleffects <- mf$v$std.error / se_f[j]
  }
  out$marginaleffects_w <- wtag(mf$w)
  mm <- cap(frm_multiple(fo, data = list(d, d),
                         family = if (fam == "binom") binomial() else cumulative()))
  if (!inherits(mm$v, "rev_err")) {
    hm <- cap(hypothesis(mm$v, hs, class = "b"))
    if (!inherits(hm$v, "rev_err")) out$multiple_hyp <- hm$v$hypothesis$Est.Error / sqrt(vh)
    out$multiple_hyp_w <- wtag(hm$w)
  } else out$multiple <- mm$v
  for (nm in names(out)) {
    v <- out[[nm]]
    cat(sprintf("  %-18s %s\n", nm,
                if (is.character(v)) substr(v[1], 1, 100) else
                  if (length(v) == 1L) sprintf("%.9g", v) else ratio_str(as.vector(v))))
  }
}

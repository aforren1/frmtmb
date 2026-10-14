# Reviewer of lane illse, attack 2: what a caller's own a' vcov() a and
# marginaleffects give on fits whose covariance is held as a factor but
# not marked lossy (kappa below se_product_kappa = 1e12), and where the
# marginaleffects error comes from. Gaussian raw polynomials of degree 5
# and 6 on [1, 2] (exact: gaussian Hessian in 256 bits) and the binomial
# of dev/illse-rev-glm.R.
#   REV_ARM=lane|base Rscript dev/illse-rev-mfx.R <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
a <- commandArgs(TRUE)
seeds <- if (length(a)) eval(parse(text = a[1])) else 301:303
gauss_case <- function(seed, deg) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + ")))
  nd <- data.frame(x = c(1.1, 1.5, 1.9))
  for (k in seq_len(deg)) nd[[xn[k]]] <- nd$x^k
  list(d = d, fo = fo, fo0 = stats::as.formula(paste("~", paste(xn, collapse = " + "))),
       nd = nd, fam = gaussian(), lab = paste0("gauss deg ", deg))
}
binom_case <- function(seed) {
  set.seed(seed)
  n <- 300
  d <- data.frame(u = 100 + stats::runif(n), z = 1e-3 * stats::rnorm(n),
                  w = 1e4 * stats::rnorm(n))
  d$u2 <- d$u^2
  d$N <- 20L
  eta <- -0.3 + 3 * (d$u - 100.5) - 2 * (d$u - 100.5)^2 + 300 * d$z +
    5e-5 * d$w
  d$y <- stats::rbinom(n, d$N, stats::plogis(eta))
  nd <- data.frame(u = 100 + c(0.1, 0.5, 0.9), z = c(-1e-3, 0, 1e-3),
                   w = c(-1e4, 0, 1e4), N = 20L)
  nd$u2 <- nd$u^2
  list(d = d, fo = y | trials(N) ~ u + u2 + z + w, fo0 = ~ u + u2 + z + w,
       nd = nd, fam = binomial(), lab = "binom mixed")
}
for (seed in seeds) {
  cases <- list(gauss_case(seed, 5), gauss_case(seed, 6), binom_case(seed))
  for (cs in cases) {
    r <- cap(frm(cs$fo, data = cs$d, family = cs$fam))
    fit <- r$v
    p <- fit$opt$par; pn <- names(p); pb <- which(pn == "beta")
    X <- stats::model.matrix(cs$fo0, cs$d)
    Xf <- as.matrix(fit$frame$linpreds[[1]]$X)
    stopifnot(identical(dim(X), dim(Xf)))
    Xm <- mM(X)
    if (cs$lab == "binom mixed") {
      yy <- m(cs$d$y); NN <- m(cs$d$N)
      nll <- function(th) {
        e <- Xm %*% mpfr2array(th[pb], dim = c(length(pb), 1L))
        -sum(yy * e - NN * log(1 + exp(e)))
      }
      H <- mp_hessian(nll, p)
    } else {
      # exact gaussian Hessian in (beta, log sigma)
      ps <- which(pn == "betad")
      rr <- m(cs$d$y) - Xm %*% mpfr2array(m(p[pb]), dim = c(length(pb), 1L))
      s2 <- exp(2 * m(p[ps]))
      H <- mM(matrix(0, length(p), length(p)))
      H[pb, pb] <- crossprod(Xm) / s2
      hb <- 2 * crossprod(Xm, rr) / s2
      H[pb, ps] <- hb; H[ps, pb] <- hb
      H[ps, ps] <- 2 * sum(rr^2) / s2
    }
    Ve <- mp_inv(H)
    A <- stats::model.matrix(cs$fo0, cs$nd)
    ve <- mp_rowvar(A, Ve[pb, pb])
    root <- ns$sdr_of(fit)$cov_fixed_root
    V <- suppressWarnings(vcov(fit, full = TRUE))
    Vb <- V[pb, pb]
    own <- rowSums((A %*% Vb) * A)
    canc <- num(rowSums((mM(abs(A)) %*% abs(Ve[pb, pb])) * mM(abs(A)))) / ve
    fl <- cap(fitted(fit, newdata = cs$nd, scale = "linear"))
    mf <- cap(marginaleffects::predictions(fit, newdata = cs$nd, type = "link"))
    mfr <- if (!inherits(mf$v, "rev_err")) mf$v$std.error / sqrt(ve) else NA
    # marginaleffects' own Jacobian with frmtmb's V, against the exact J
    # (A on the link scale) with the same V
    gv <- cap(ns$get_vcov.frmtmb_fit(fit))
    cat(sprintf(paste0("%-13s s%d | code %d | root %s lossy %s | ",
                       "cancellation |a||V||a|/a'Va %.1e-%.1e | ",
                       "fitted %s | own a' vcov() a %s | ",
                       "marginaleffects link SE %s | warn mfx %s, ",
                       "get_vcov finite %s\n"),
                cs$lab, seed, fit$opt$convergence, !is.null(root),
                isTRUE(attr(root, "lossy")), min(canc), max(canc),
                if (inherits(fl$v, "rev_err")) "ERR" else
                  ratio_str(fl$v[, "Est.Error"] / sqrt(ve)),
                ratio_str(sqrt(pmax(own, 0)) / sqrt(ve)),
                ratio_str(mfr), wtag(mf$w),
                all(is.finite(gv$v))))
  }
}

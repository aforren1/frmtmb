# Reviewer of lane illse, attack 2: marginaleffects on the response scale
# of a GLM whose covariance is held as a factor but not marked lossy.
# Bernoulli and poisson, raw polynomial through I(x^k) on [1, 2], degree
# 4 and 5. Exact: 256-bit Hessian of the log-likelihood at the fit's
# point; predictions() and slopes(x) against the delta method in 256
# bits. Also emmeans(type = "response") and fitted().
#   REV_ARM=lane|base Rscript dev/illse-rev-mfx2.R <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
a <- commandArgs(TRUE)
seeds <- if (length(a)) eval(parse(text = a[1])) else 401:403
xs <- c(1.2, 1.5, 1.8)
for (seed in seeds) for (fam in c("bernoulli", "poisson")) for (deg in 4:5) {
  set.seed(seed)
  n <- 400
  d <- data.frame(x = stats::runif(n, 1, 2))
  eta <- if (fam == "bernoulli") 2 * sin(3 * d$x) else 1 + sin(3 * d$x)
  d$y <- if (fam == "bernoulli") stats::rbinom(n, 1, stats::plogis(eta)) else
    stats::rpois(n, exp(eta))
  rhs <- paste(c("x", if (deg > 1) paste0("I(x^", 2:deg, ")")), collapse = " + ")
  fo <- stats::as.formula(paste("y ~", rhs))
  r <- cap(frm(fo, data = d, family = if (fam == "bernoulli") bernoulli() else poisson()))
  fit <- r$v
  if (inherits(fit, "rev_err")) { cat(fam, deg, seed, "ERROR", fit, "\n"); next }
  p <- fit$opt$par; pb <- which(names(p) == "beta")
  X <- outer(d$x, 0:deg, `^`)
  Xf <- as.matrix(fit$frame$linpreds[[1]]$X)
  stopifnot(identical(dim(X), dim(Xf)), max(abs(X - Xf)) == 0)
  Xm <- mM(X); yy <- m(d$y)
  nll <- if (fam == "bernoulli") function(th) {
    e <- Xm %*% mpfr2array(th, dim = c(length(th), 1L))
    -sum(yy * e - log(1 + exp(e)))
  } else function(th) {
    e <- Xm %*% mpfr2array(th, dim = c(length(th), 1L))
    -sum(yy * e - exp(e))
  }
  H <- mp_hessian(nll, p)
  Ve <- mp_inv(H)
  mu_m <- function(e) if (fam == "bernoulli") 1 / (1 + exp(-e)) else exp(e)
  dv <- function(fun) vapply(xs, function(x0) {
    gi <- mp_grad(function(th) fun(th, x0), p)
    num(sum(gi * (Ve %*% mpfr2array(gi, dim = c(length(p), 1L)))))
  }, 0)
  pred_f <- function(th, x0) mu_m(sum(m(x0^(0:deg)) * th))
  slope_f <- function(th, x0) {
    e <- sum(m(x0^(0:deg)) * th)
    de <- sum(m(c(0, (1:deg) * x0^(0:(deg - 1)))) * th)
    mu <- mu_m(e)
    (if (fam == "bernoulli") mu * (1 - mu) else mu) * de
  }
  v_pred <- dv(pred_f)
  v_slope <- dv(slope_f)
  nd <- data.frame(x = xs)
  root <- ns$sdr_of(fit)$cov_fixed_root
  ft <- cap(fitted(fit, newdata = nd))
  mp1 <- cap(marginaleffects::predictions(fit, newdata = nd))
  ms <- cap(marginaleffects::slopes(fit, newdata = nd, variables = "x"))
  em <- cap(summary(emmeans::emmeans(fit, ~ x, at = list(x = 1.5),
                                     type = "response")))
  D <- sqrt(num(diag(H)))
  ev <- eigen(num(H) / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  rs <- function(v, ref) if (inherits(v, "rev_err")) "ERR" else ratio_str(v / sqrt(ref))
  cat(sprintf(paste0("%-9s deg %d s%d | code %d | kappa %.1e | root %s lossy %s | ",
                     "fitted %s | mfx pred %s (warn %s) | mfx slope %s ",
                     "(warn %s) | emmeans resp %s\n"),
              fam, deg, seed, fit$opt$convergence, max(ev) / min(ev),
              !is.null(root), isTRUE(attr(root, "lossy")),
              rs(ft$v[, "Est.Error"], v_pred),
              if (inherits(mp1$v, "rev_err")) "ERR" else rs(mp1$v$std.error, v_pred),
              wtag(mp1$w),
              if (inherits(ms$v, "rev_err")) "ERR" else rs(ms$v$std.error, v_slope),
              wtag(ms$w),
              if (inherits(em$v, "rev_err")) "ERR" else rs(em$v$SE, v_pred[2])))
}

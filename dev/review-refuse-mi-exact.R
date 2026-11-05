# Item 2/3: mi() against exact references, all 13 families. For each family:
# (a) an unread response: complete-case fit, joint Laplace fit, gaps in logLik,
#     and shift in SE of EVERY parameter, for several n and missing shares;
# (b) a response read by an mi() predictor: frmtmb's Laplace marginal at the
#     estimate against integrate() over each missing value of the exact joint.
# Usage: Rscript dev/review-refuse-mi-exact.R <family>. Seed 9200 + family.
fam <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/review-refuse-lib.R")
fams <- c("gaussian", "student", "skew_normal", "beta", "Gamma", "weibull",
          "exponential", "frechet", "inverse.gaussian", "lognormal",
          "exgaussian", "gen_extreme_value", "von_mises")
fi <- match(fam, fams)
ff <- switch(fam,
  beta = frmtmb::Beta(),
  Gamma = stats::Gamma(link = "log"),
  inverse.gaussian = stats::inverse.gaussian(link = "log"),
  gaussian = stats::gaussian(),
  get(fam, asNamespace("frmtmb"))())
gen <- function(n) {
  x <- rnorm(n)
  # a design of mine: stronger slope, a shape of its own, not the lane's
  eta <- 0.3 + 0.7 * x
  y <- switch(fam,
    gaussian = rnorm(n, eta, 0.8),
    student = eta + 0.8 * rt(n, 3),
    skew_normal = eta + 0.5 * rnorm(n) + 1.2 * abs(rnorm(n)),
    beta = rbeta(n, plogis(eta) * 3, (1 - plogis(eta)) * 3),
    Gamma = rgamma(n, 1.3, 1.3 / exp(eta)),
    weibull = rweibull(n, 1.4, exp(eta)),
    exponential = rexp(n, 1 / exp(eta)),
    frechet = 1 / rweibull(n, 2.2, exp(-eta)),
    inverse.gaussian = {
      m <- exp(eta); v <- rnorm(n)^2; lam <- 2
      yy <- m + m^2 * v / (2 * lam) -
        m / (2 * lam) * sqrt(4 * m * lam * v + m^2 * v^2)
      ifelse(runif(n) < m / (m + yy), yy, m^2 / yy)
    },
    lognormal = rlnorm(n, eta, 0.9),
    exgaussian = eta + rnorm(n, 0, 0.3) + rexp(n, 1 / 1.2),
    gen_extreme_value = eta + 0.8 * (rexp(n)^(-0.2) - 1) / 0.2,
    von_mises = ((2 * atan(0.4 * eta) + rnorm(n, 0, 0.7) + pi) %% (2 * pi)) - pi)
  data.frame(x = x, y = y)
}
set.seed(9200 + fi)
for (cfg in list(c(200, 15), c(200, 60), c(1000, 75))) {
  n <- cfg[1]; m <- cfg[2]
  d <- gen(n)
  d$y[sample(n, m)] <- NA
  a <- tr(frm(bf(y | mi() ~ x), family = ff, data = d))
  b <- tr(frm(bf(y ~ x), family = ff, data = d[!is.na(d$y), ]))
  if (is_err(a) || is_err(b)) {
    pr("A %-18s n=%d m=%d ERROR joint: %s; cc: %s", fam, n, m, msg(a), msg(b))
    next
  }
  gap <- as.numeric(logLik(a)) - as.numeric(logLik(b))
  ea <- c(a$estimates$beta, a$estimates$betad)
  eb <- c(b$estimates$beta, b$estimates$betad)
  se <- tr(sqrt(diag(vcov(b, include_dpars = TRUE))))
  H <- stats::optimHess(b$opt$par, b$obj$fn, b$obj$gr)
  sb <- sqrt(diag(solve(H)))
  shifts <- (unlist(ea) - unlist(eb)) / sb[seq_along(eb)]
  pr("V %-18s n=%4d m=%3d joint=[%s] cc=[%s]", fam, n, m,
     paste(sprintf("%.4f", unlist(ea)), collapse = ","),
     paste(sprintf("%.4f", unlist(eb)), collapse = ","))
  pr("A %-18s n=%4d m=%3d conv=%d/%d gap/miss=%+.5f gap=%+.3f shifts(SE): %s", fam, n, m,
     a$opt$convergence, b$opt$convergence, gap / m, gap,
     paste(sprintf("%s=%+.4f", names(eb), shifts), collapse = " "))
}

# (b) read by a predictor: frmtmb's Laplace against integrate(). z ~ mi(y) + x
set.seed(9300 + fi)
n <- 80
d <- gen(n)
d$z <- 0.5 + 0.6 * d$y + 0.3 * d$x + rnorm(n, 0, 0.7)
mrows <- sample(n, 8)
d$y[mrows] <- NA
fit <- tr(frm(bf(y | mi() ~ x) + ff + bf(z ~ mi(y) + x) + gaussian(), data = d))
if (is_err(fit)) {
  pr("B %-18s ERROR %s", fam, msg(fit))
} else {
  lp <- fit$frame$linpreds
  env <- fit$obj$env
  # the exact marginal at the estimate: the observed rows' densities plus
  # one integrate() per missing row, over the family's support on the value
  # scale, of family density x z density. Evaluated by tape: per-row deltas
  # of the joint in u, holding every other value at its mode, are exact,
  # because the joint is additive in the missing values.
  pv <- env$last.par.best
  mi_i <- which(names(pv) == "miss")
  base <- env$f(pv)
  chart <- fit$frame$mi_map$y$chart %||% "identity"
  err <- numeric(length(mi_i))
  for (k in seq_along(mi_i)) {
    u0 <- pv[mi_i[k]]
    # the integrand in u, relative to the mode, on the tape itself: the
    # chart's Jacobian is in it, so the integral over u is chart free
    h1 <- function(u) {
      vapply(u, function(ui) {
        p <- pv; p[mi_i[k]] <- ui
        v <- env$f(p)
        exp(base - v)
      }, 0)
    }
    # second derivative at the mode by finite differences for the step
    eps <- 1e-3
    d2 <- (env$f(replace(pv, mi_i[k], u0 + eps)) +
             env$f(replace(pv, mi_i[k], u0 - eps)) - 2 * base) / eps^2
    w <- 1 / sqrt(max(d2, 1e-6))
    lo <- u0 - 40 * w; hi <- u0 + 40 * w
    I <- tryCatch(integrate(h1, lo, hi, subdivisions = 2000L,
                            rel.tol = 1e-9, stop.on.error = FALSE)$value,
                  error = function(e) NA_real_)
    lap <- sqrt(2 * pi) * w
    err[k] <- log(lap) - log(I)
  }
  pr("B %-18s chart=%-8s per missing, log Laplace - log exact: mean %+.5f range [%+.5f, %+.5f]; total %+.4f; conv %d",
     fam, chart, mean(err), min(err), max(err), sum(err), fit$opt$convergence)
  # coverage of the readers
  ft <- tr(fitted(fit, resp = "z"))
  pr("B %-18s fitted(resp=z) finite: %s; predictor reads value in support: %s",
     fam, if (is_err(ft)) msg(ft) else all(is.finite(ft[, 1])),
     if (fam %in% c("Gamma", "weibull", "exponential", "frechet",
                    "inverse.gaussian", "lognormal")) {
       all(frmtmb:::mi_values(fit, "y") > 0)
     } else if (fam == "beta") {
       all(frmtmb:::mi_values(fit, "y") > 0 & frmtmb:::mi_values(fit, "y") < 1)
     } else if (fam == "von_mises") {
       all(abs(frmtmb:::mi_values(fit, "y")) < pi)
     } else NA)
}

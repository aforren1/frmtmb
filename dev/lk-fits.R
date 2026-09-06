suppressMessages(library(frmtmb))
options(digits = 10)
set.seed(404)
n <- 400
d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
d$nt <- sample(4:10, n, TRUE)
tight <- frmtmb_control(optCtrl = list(rel.tol = 1e-14, x.tol = 1e-12,
                                        eval.max = 5000, iter.max = 5000))

cat("== binomial GLM agreement, frmtmb vs stats::glm ==\n")
for (lk in c("probit", "cauchit", "cloglog", "logit")) {
  p <- switch(lk,
    probit = stats::pnorm(0.3 + 0.8 * d$x),
    cauchit = stats::pcauchy(0.3 + 0.8 * d$x),
    cloglog = 1 - exp(-exp(0.3 + 0.8 * d$x)),
    stats::plogis(0.3 + 0.8 * d$x))
  d$y <- stats::rbinom(n, d$nt, p)
  g <- stats::glm(cbind(y, nt - y) ~ x + z, data = d,
                  family = stats::binomial(link = lk),
                  control = stats::glm.control(epsilon = 1e-14, maxit = 200))
  for (tg in c(FALSE, TRUE)) {
    fit <- suppressWarnings(frm(bf(y | trials(nt) ~ x + z),
                                family = binomial(link = lk), data = d,
                                control = if (tg) tight else frmtmb_control()))
    cat(sprintf("%-9s %-8s max|coef diff| = %.3g   logLik diff = %.3g\n",
                lk, if (tg) "tight" else "default",
                max(abs(fixef(fit)$mu - stats::coef(g))),
                abs(as.numeric(stats::logLik(fit)) -
                      as.numeric(stats::logLik(g)))))
  }
}

cat("\n== poisson sqrt against stats::glm ==\n")
d$yp <- stats::rpois(n, (1.2 + 0.4 * d$x)^2)
gp <- stats::glm(yp ~ x, data = d, family = stats::poisson(link = "sqrt"),
                 control = stats::glm.control(epsilon = 1e-14, maxit = 200))
for (tg in c(FALSE, TRUE)) {
  fp <- suppressWarnings(frm(bf(yp ~ x), family = poisson(link = "sqrt"),
                             data = d,
                             control = if (tg) tight else frmtmb_control()))
  cat(sprintf("sqrt      %-8s max|coef diff| = %.3g   logLik diff = %.3g\n",
              if (tg) "tight" else "default",
              max(abs(fixef(fp)$mu - stats::coef(gp))),
              abs(as.numeric(stats::logLik(fp)) -
                    as.numeric(stats::logLik(gp)))))
}

cat("\n== inverse.gaussian 1/mu^2 against stats::glm ==\n")
mu <- exp(0.4 + 0.2 * d$x)
d$yi <- stats::rgamma(n, 8, 8 / mu)
gi <- stats::glm(yi ~ x, data = d,
                 family = stats::inverse.gaussian(link = "1/mu^2"),
                 start = c(1 / mean(d$yi)^2, 0),
                 control = stats::glm.control(epsilon = 1e-14, maxit = 200))
fi <- suppressWarnings(frm(bf(yi ~ x),
                           family = inverse.gaussian(link = "1/mu^2"),
                           data = d, control = tight))
cat(sprintf("1/mu^2    max|coef diff| = %.3g\n",
            max(abs(fixef(fi)$mu - stats::coef(gi)))))

cat("\n== links stats::make.link rejects: reachable via frmtmb families ==\n")
d$yb <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.8 * d$x))
for (lk in c("probit", "probit_approx", "cauchit", "softit")) {
  f <- suppressWarnings(frm(bf(yb ~ x), family = bernoulli(link = lk),
                            data = d))
  cat(sprintf("bernoulli(%-14s) logLik %.5f coefs %s\n", lk,
              as.numeric(stats::logLik(f)),
              paste(round(fixef(f)$mu, 4), collapse = " ")))
}
cat("\nstats::binomial(link = 'softit'): ",
    tryCatch({ binomial(link = "softit"); "accepted" },
             error = function(e) paste("rejected -", conditionMessage(e))),
    "\n", sep = "")

d$yw <- stats::rweibull(n, 2, exp(0.5 + 0.3 * d$x))
for (lk in c("softplus", "squareplus")) {
  f <- suppressWarnings(frm(bf(yw ~ x), family = weibull(link = lk),
                            data = d))
  cat(sprintf("weibull(%-10s) logLik %.5f coefs %s\n", lk,
              as.numeric(stats::logLik(f)),
              paste(round(fixef(f)$mu, 4), collapse = " ")))
}
d$yn <- stats::rnbinom(n, size = 3, mu = exp(0.6 + 0.4 * d$x))
for (lk in c("softplus", "squareplus", "sqrt")) {
  f <- suppressWarnings(frm(bf(yn ~ x), family = negbinomial(link = lk),
                            data = d))
  cat(sprintf("negbinomial(%-10s) logLik %.5f coefs %s\n", lk,
              as.numeric(stats::logLik(f)),
              paste(round(fixef(f)$mu, 4), collapse = " ")))
}
d$ybe <- stats::rbeta(n, stats::pnorm(0.2 + 0.5 * d$x) * 8,
                      (1 - stats::pnorm(0.2 + 0.5 * d$x)) * 8)
for (lk in c("probit", "cauchit", "softit")) {
  f <- suppressWarnings(frm(bf(ybe ~ x), family = Beta(link = lk), data = d))
  cat(sprintf("Beta(%-8s) logLik %.5f coefs %s\n", lk,
              as.numeric(stats::logLik(f)),
              paste(round(fixef(f)$mu, 4), collapse = " ")))
}

cat("\n== predict(se.fit = TRUE) against the delta method ==\n")
for (lk in c("probit", "cauchit", "probit_approx", "softit")) {
  f <- suppressWarnings(frm(bf(yb ~ x), family = bernoulli(link = lk),
                            data = d))
  nd <- data.frame(x = c(-1.5, -0.4, 0.6, 2))
  pl <- stats::predict(f, newdata = nd, type = "link", se.fit = TRUE)
  pr <- stats::predict(f, newdata = nd, type = "response", se.fit = TRUE)
  lo <- frmtmb:::frmtmb_links[[lk]]
  dn <- vapply(pl$fit, function(z) numDeriv::grad(lo$linkinv, z), 0)
  cat(sprintf("%-14s max rel diff %.3g\n", lk,
              max(abs(pr$se.fit - abs(dn) * pl$se.fit) /
                    pmax(1e-12, abs(pr$se.fit)))))
}

cat("\n== conditional_effects bands ==\n")
for (lk in c("probit", "cauchit", "softit", "probit_approx")) {
  f <- suppressWarnings(frm(bf(yb ~ x), family = bernoulli(link = lk),
                            data = d))
  r <- conditional_effects(f)[[1]]
  cat(sprintf("%-14s rows %d finite %s inside (0,1) %s\n", lk, nrow(r),
              all(is.finite(r$estimate__)),
              all(r$lower__ > 0 & r$upper__ < 1)))
}

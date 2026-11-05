# Short check, items 6 and 7. von_mises on the line against brms's window;
# frm_sample(formula, laplace = TRUE) on an unread row. Seed 11400.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
quiet <- function(e) suppressMessages(suppressWarnings(e))
set.seed(11400)
n <- 150; x <- rnorm(n)
rvm <- function(mu, k) {
  vapply(mu, function(m) { repeat { w <- runif(1, -pi, pi); if (runif(1) < exp(k * (cos(w) - 1))) break }; ((w + m + pi) %% (2 * pi)) - pi }, 0)
}
for (cfg in list(c(0, 3), c(3.05, 3), c(3.14, 20))) {
  th <- cfg[1] + 0.05 * x
  y0 <- rvm(th, cfg[2])
  d <- data.frame(x = x, y = y0)
  # z reads the value linearly, so the window and the line are different models
  d$z <- 0.5 + 0.5 * y0 + rnorm(n, 0, 0.5)
  rows <- seq(5, n, by = 10); d$y[rows] <- NA
  fit <- tr(frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y)) + gaussian(), data = d))
  if (is_err(fit)) { pr("mean %.2f kappa %g: %s", cfg[1], cfg[2], msg(fit, 150)); next }
  ch <- fit$frame$mi_map$y$chart
  v <- frmtmb:::mi_values(fit, "y")[rows]
  pr("mean %.2f kappa %g: chart %s; latent range %.3f %.3f; outside (-pi, pi): %d of %d; true y at those rows range %.3f %.3f",
     cfg[1], cfg[2], ch, min(v), max(v), sum(abs(v) >= pi), length(v), min(y0[rows]), max(y0[rows]))
  mu <- fitted(fit, resp = "y", dpar = "mu")
  pr("   reported mean direction range %.3f %.3f (brms's interval is (-pi, pi))", min(mu[, 1]), max(mu[, 1]))
  ey <- fit$estimates
  ey_mu <- 2 * atan(fixef(fit)["y_Intercept", 1])
  pr("   fixef mean direction at x = 0: %.4f", ey_mu)
  # the same point on the two routes: the ML objective and the sampling one
  res <- tr(mi_sampling(fit))
  pr("   mi_sampling: %s; chart %s; start u range %s", msg(res, 100), if (is_err(res)) "-" else res$frame$mi_map$y$chart,
     if (is_err(res)) "-" else paste(round(range(res$estimates$miss), 2), collapse = " "))
  if (!is_err(res)) {
    pv1 <- fit$obj$env$last.par.best; pv2 <- res$obj$env$last.par
    i1 <- which(names(pv1) == "miss"); i2 <- which(names(pv2) == "miss")
    # one value inside the window on both: the density equal but the Jacobian
    w <- ((v + pi) %% (2 * pi)) - pi
    pv1[i1] <- w
    p <- (w + pi) / (2 * pi); pv2[i2] <- qlogis(p); pv2[-i2] <- pv1[-i1]
    jac <- sum(log(2 * pi) + log(p) + log1p(-p))
    pr("   one point in the window: -f on the line %.8f; -f on the window %.8f; difference minus the Jacobian %.2e", -fit$obj$env$f(pv1), -res$obj$env$f(pv2), -res$obj$env$f(pv2) + fit$obj$env$f(pv1) - jac)
    # the point the line's mode reports, when it leaves the window: z reads
    # v, brms reads a value in the window
    if (any(abs(v) >= pi)) {
      pr("   line's mode leaves the window on %d rows: density the two read differ", sum(abs(v) >= pi))
    }
  }
  # brms at frm's point, through the instrument
  spec <- list(bf = "mvbf(bf(y | mi() ~ x, family = von_mises()), bf(z ~ mi(y)), rescor = FALSE)", family = NA, data = d)
  r <- tr(meas(spec))
  if (is_err(r)) pr("   L1: %s", msg(r, 160)) else pr("   L1: verdict %s const %s pts [%s] grad %s err %s", r$verdict, format(r$l1_const %||% NA, digits = 3),
    paste(format(unlist(r$l1_extra_const) %||% NA, digits = 3), collapse = ","), format(r$l1_grad %||% NA, digits = 3), substr(r$l1_error %||% "", 1, 110))
}
# an unread response on the line: zero Laplace
d <- data.frame(x = x, y = rvm(3.05 + 0.05 * x, 3)); d$y[seq(5, n, by = 10)] <- NA
a <- frm(bf(y | mi() ~ x), family = von_mises(), data = d); b <- frm(bf(y ~ x), family = von_mises(), data = d[-seq(5, n, by = 10), ])
pr("von_mises unread: logLik diff %.2e", as.numeric(logLik(a)) - as.numeric(logLik(b)))
sa <- tr(quiet(frm_sample(a, chains = 2, iter = 600, seed = 1, refresh = 0)))
pr("frm_sample(fit), von_mises near the end: %s", msg(sa, 150))
if (!is_err(sa)) {
  dm <- as.data.frame(posterior::as_draws_df(sa)); nn <- grep("^Ymi", names(dm), value = TRUE)
  pr("  Ymi range %.3f %.3f", min(dm[nn]), max(dm[nn]))
}
# exact against the line, a read von_mises, per value by integrate over the
# window: brms's model, a value in (-pi, pi), read linearly by z
set.seed(11401)
th <- 3.05 + 0.05 * x; y0 <- rvm(th, 3)
d <- data.frame(x = x, y = y0); d$z <- 0.5 + 0.5 * y0 + rnorm(n, 0, 1.5); rows <- seq(5, n, by = 10); d$y[rows] <- NA
fit <- frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y)) + gaussian(), data = d)
env <- fit$obj$env; pv <- env$last.par.best; ii <- which(names(pv) == "miss"); b0 <- env$f(pv)
err <- win <- numeric(length(ii))
for (k in seq_along(ii)) {
  h <- function(v) vapply(v, function(s) exp(b0 - env$f(replace(pv, ii[k], s))), 0)
  I <- integrate(h, -pi, pi, subdivisions = 2000L, stop.on.error = FALSE)$value
  H <- env$spHess(pv, random = TRUE)
  err[k] <- 0.5 * log(2 * pi) - 0.5 * log(H[k, k]) - log(I)
  win[k] <- pv[ii[k]]
}
pr("read von_mises near pi, weak reader, Laplace on the line minus exact on brms's window: mean %+.4f range [%+.4f, %+.4f]; latent range %.3f %.3f", mean(err), min(err), max(err), min(win), max(win))

# item 7: laplace = TRUE from a formula on an unread row
set.seed(10600)
n <- 200; x <- rnorm(n)
d <- data.frame(x = x, y = 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2))
rows <- sample(n, 60); d$y[rows] <- NA; dc <- d[-rows, ]
pf <- set_prior("normal(0, 5)", class = "b") + set_prior("normal(0, 5)", class = "Intercept") + set_prior("exponential(1)", class = "sigma") + set_prior("exponential(0.5)", class = "beta")
s1 <- tr(quiet(frm_sample(bf(y | mi() ~ x), family = exgaussian(), data = d, prior = pf, laplace = TRUE, chains = 4, iter = 1500, seed = 1, refresh = 0, cores = 1)))
s2 <- tr(quiet(frm_sample(bf(y ~ x), family = exgaussian(), data = dc, prior = pf, laplace = FALSE, chains = 4, iter = 1500, seed = 1, refresh = 0, cores = 1)))
s3 <- tr(quiet(frm_sample(frm(bf(y | mi() ~ x), family = exgaussian(), data = d), prior = pf, laplace = TRUE, chains = 4, iter = 1500, seed = 1, refresh = 0, cores = 1)))
sm <- function(lab, s) {
  if (is_err(s)) { pr("%-34s %s", lab, msg(s, 160)); return(invisible()) }
  dd <- as.data.frame(posterior::as_draws_df(s))
  pr("%-34s sigma %.3f (sd %.3f); beta %.3f (sd %.3f); b_x %.3f", lab, mean(dd$sigma), sd(dd$sigma), mean(dd$beta), sd(dd$beta), mean(dd$b_x))
}
sm("frm_sample(formula, laplace = TRUE)", s1)
sm("observed rows, full sampling", s2)
sm("frm_sample(fit, laplace = TRUE)", s3)
s4 <- tr(quiet(frm_sample(bf(y | mi() ~ x), family = student(), data = transform(d, y = 0.3 + 0.5 * x + rt(n, 6) + ifelse(is.na(y), NA, 0)), laplace = TRUE, chains = 2, iter = 600, seed = 1, refresh = 0, cores = 1)))
pr("student formula laplace: %s", msg(s4, 100))
g <- tr(grep("laplace", readLines("extensions/frmtmb.sample/man/frm_sample.Rd"), value = TRUE))
pr("doc lines on laplace = %d", length(g))

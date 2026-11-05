# Items 2, 6: the cens x mixture gradient against finite differences; a Newton
# polish; the ordered default's 2.9e-3; the no_mode chart on a positive mixture
# under frm_sample, where the placeholder is 0. Seed 10300.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(10200); n <- 200; x <- rnorm(n)
set.seed(10202)
k <- rbinom(n, 1, 0.35)
y <- ifelse(k == 1, 2.5 + 0.5 * x, -2 - 0.3 * x) + rnorm(n)
ce <- ifelse(y > 3, 1L, ifelse(y < -3, -1L, 0L))
d <- data.frame(x = x, y = pmax(pmin(y, 3), -3), ce = ce)
pr("censored rows: %d", sum(ce != 0))
chk <- function(lab, form, fam) {
  a <- frm(form, family = fam, data = d)
  p <- a$opt$par
  g <- a$obj$gr(p)
  fd <- vapply(seq_along(p), function(i) {
    h <- 1e-5 * max(1, abs(p[i]))
    (a$obj$fn(replace(p, i, p[i] + h)) - a$obj$fn(replace(p, i, p[i] - h))) / (2 * h)
  }, 0)
  H <- stats::optimHess(p, a$obj$fn, a$obj$gr)
  p2 <- p
  for (it in 1:6) p2 <- p2 - as.numeric(solve(H, as.numeric(a$obj$gr(p2))))
  pr("%-26s max|gr| %.2e; AD vs FD max diff %.2e; after 6 Newton steps max|gr| %.2e, objective gain %.2e, largest step %.2e",
     lab, max(abs(g)), max(abs(g - fd)), max(abs(a$obj$gr(p2))), a$obj$fn(p) - a$obj$fn(p2), max(abs(p2 - p)))
}
mx <- mixture(gaussian(), gaussian(), order = "none")
chk("none, no cens", bf(y ~ 1, mu1 ~ x, mu2 ~ x), mx)
chk("none, cens", bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x), mx)
chk("ordered, no cens", bf(y ~ x), mixture(gaussian(), gaussian()))
chk("ordered, cens", bf(y | cens(ce) ~ x), mixture(gaussian(), gaussian()))

fe <- frm(bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x), family = mx, data = d)
cf <- fixef(fe); est <- fe$estimates
mu1 <- cf['mu1_Intercept', 1] + cf['mu1_x', 1] * x; mu2 <- cf['mu2_Intercept', 1] + cf['mu2_x', 1] * x
s1 <- exp(est$betad[['sigma1_(Intercept)']]); s2 <- exp(est$betad[['sigma2_(Intercept)']])
w1 <- plogis(est$betad[['theta1_(Intercept)']]); w2 <- 1 - w1
ld <- ifelse(ce == 0, log(w1 * dnorm(d$y, mu1, s1) + w2 * dnorm(d$y, mu2, s2)), ifelse(ce == 1, log(w1 * pnorm(d$y, mu1, s1, FALSE) + w2 * pnorm(d$y, mu2, s2, FALSE)), log(w1 * pnorm(d$y, mu1, s1) + w2 * pnorm(d$y, mu2, s2))))
pr('cens mix: logLik - hand, right and left rows by R pnorm %.3e', as.numeric(logLik(fe)) - sum(ld))
show <- function(lab, r) { if (is_err(r)) { pr('%-34s ERROR %s', lab, msg(r, 170)); return(invisible()) }; pr('%-34s verdict=%s const=%s pts=[%s] grad=%s l1[%s]', lab, r$verdict, format(r$l1_const %||% NA, digits = 3), paste(format(unlist(r$l1_extra_const) %||% NA, digits = 3), collapse = ','), format(r$l1_grad %||% NA, digits = 3), substr(r$l1_error %||% r$frm_error %||% '', 1, 110)) }
set.seed(10301); n0 <- 200; x0 <- rnorm(n0); dm <- data.frame(x = x0, ce = 0L, cc = sample(c(0L, 1L, 1L, -1L, 0L, 0L), n0, TRUE))
dm$yc <- rpois(n0, exp(ifelse(rbinom(n0, 1, 0.4) == 1, 0.2, 1.8) + 0.3 * x0)); dm$cp <- ifelse(dm$yc >= 8, 1L, ifelse(dm$yc <= 0, -1L, 0L)); dm$cp[dm$yc <= 0] <- 0L
show('cens mix poisson poisson', tr(meas(list(bf = 'bf(yc | cens(cp) ~ 1, mu1 ~ x, mu2 ~ x)', family = "mixture(poisson(), poisson(), order = 'none')", data = dm))))
dm$yw <- rweibull(n0, 1.5, exp(ifelse(rbinom(n0, 1, 0.4) == 1, -0.5, 1) + 0.3 * x0)); dm$cw <- ifelse(dm$yw > 3, 1L, 0L)
show('cens mix weibull exponential', tr(meas(list(bf = 'bf(yw | cens(cw) ~ 1, mu1 ~ x, mu2 ~ x)', family = "mixture(weibull(), exponential(), order = 'none')", data = dm))))
dm$yv <- rbeta(n0, 2, 5); dm$cv <- ifelse(dm$yv > 0.6, 1L, 0L)
show('cens mix beta beta', tr(meas(list(bf = 'bf(yv | cens(cv) ~ 1, mu1 ~ x, mu2 ~ x)', family = "mixture(Beta(), Beta(), order = 'none')", data = dm))))

# positive mixtures on the sampling route
set.seed(10300)
n <- 120; x <- rnorm(n)
kk <- rbinom(n, 1, 0.4)
d2 <- data.frame(x = x, y = ifelse(kk == 1, rlnorm(n, 1, 0.4), rlnorm(n, -1, 0.4)))
d2$y[seq(5, n, 15)] <- NA
fam <- mixture(lognormal(), lognormal())
o <- tr(frm(bf(y | mi() ~ x), family = fam, data = d2, dry_run = "objective"))
pr("mixture(lognormal, lognormal) mi, objective: %s", msg(o, 200))
b <- tr(brms::make_stancode(brms::bf(y | mi() ~ x), family = brms::mixture(lognormal(), lognormal()), data = d2))
pr("  brms stancode: %s; Ymi bound: %s", msg(b, 100), if (is_err(b)) "" else paste(grep("Ymi;", strsplit(b, "\n")[[1]], value = TRUE), collapse = " "))
s <- tr(suppressWarnings(frm_sample(bf(y | mi() ~ x), family = fam, data = d2, chains = 1, iter = 200, seed = 3, refresh = 0)))
pr("frm_sample mixture(lognormal, lognormal) mi: %s", msg(s, 240))
fam2 <- mixture(stats::Gamma(link = "log"), stats::Gamma(link = "log"))
s <- tr(suppressWarnings(frm_sample(bf(y | mi() ~ x), family = fam2, data = d2, chains = 1, iter = 200, seed = 3, refresh = 0)))
pr("frm_sample mixture(Gamma, Gamma) mi: %s", msg(s, 240))
d3 <- data.frame(x = x, y = plogis(rnorm(n))); d3$y[seq(5, n, 15)] <- NA
s <- tr(suppressWarnings(frm_sample(bf(y | mi() ~ x), family = mixture(Beta(), Beta()), data = d3, chains = 1, iter = 200, seed = 3, refresh = 0)))
pr("frm_sample mixture(Beta, Beta) mi: %s", msg(s, 240))
if (!is_err(s)) {
  mm <- as.data.frame(posterior::as_draws_df(s))
  ms <- grep("^miss", names(mm), value = TRUE)
  pr("  miss draws outside (0, 1): %.3f of %d", mean(unlist(mm[ms]) <= 0 | unlist(mm[ms]) >= 1), length(ms) * nrow(mm))
}
d4 <- data.frame(x = x, y = x + rexp(n) - rexp(n)); d4$y[seq(5, n, 15)] <- NA
s <- tr(suppressWarnings(frm_sample(bf(y | mi() ~ x), family = asym_laplace(), data = d4, chains = 2, iter = 600, seed = 3, refresh = 0)))
pr("frm_sample asym_laplace mi: %s", msg(s, 240))
d5 <- data.frame(x = x, y = c(-1, 1)[kk + 1L] * 2 + rnorm(n, 0, 0.5)); d5$y[seq(5, n, 15)] <- NA
s <- tr(suppressWarnings(frm_sample(bf(y | mi() ~ x), family = mixture(gaussian(), gaussian()), data = d5, chains = 2, iter = 600, seed = 3, refresh = 0)))
pr("frm_sample mixture(gaussian, gaussian) mi: %s", msg(s, 240))
if (!is_err(s)) {
  mm <- as.data.frame(posterior::as_draws_df(s))
  ms <- grep("^miss", names(mm), value = TRUE)
  pr("  a missing value's draws: two modes, share below 0: %.3f (want about share of component 1)", mean(mm[[ms[1]]] < 0))
}
# the frm() refusal's advice: "drop them to fit by maximum likelihood" works
a <- tr(frm(bf(y ~ x), family = asym_laplace(), data = d4[!is.na(d4$y), ]))
pr("asym_laplace, rows dropped: %s %s", msg(a, 120), paste(substr(attr(a, "w"), 1, 80), collapse = "; "))
# mi x asym_laplace as a PREDICTOR-read case with a kink: brief's "mode at the kink"
d6 <- data.frame(x = x, w = rnorm(n)); d6$y <- 1 + 0.8 * d6$w + rexp(n) - rexp(n); d6$w[seq(5, n, 15)] <- NA
a <- tr(frm(bf(y ~ mi(w)) + asym_laplace() + bf(w | mi() ~ x) + gaussian(), data = d6))
pr("asym_laplace response reads a gaussian mi() predictor: %s %s", msg(a, 180), paste(substr(attr(a, "w"), 1, 100), collapse = "; "))
if (!is_err(a)) {
  env <- a$obj$env; pv <- env$last.par.best; ii <- which(names(pv) == "miss")
  # exact per-row marginal of a kinked inner problem by integrate
  b0 <- env$f(pv); err <- numeric(length(ii))
  for (k in seq_along(ii)) {
    u0 <- pv[ii[k]]
    h1 <- function(u) vapply(u, function(ui) { p <- pv; p[ii[k]] <- ui; exp(b0 - env$f(p)) }, 0)
    I <- integrate(h1, u0 - 15, u0 + 15, subdivisions = 2000L, stop.on.error = FALSE, rel.tol = 1e-8)$value
    # the Laplace width the tape's Hessian reads at the mode, per row
    he <- env$spHess(pv, random = TRUE)
    w <- 1 / sqrt(he[k, k])
    err[k] <- log(sqrt(2 * pi) * w) - log(I)
  }
  pr("  log Laplace - log exact per missing value: mean %+.4f range [%+.4f, %+.4f]; max |inner grad| %.3f", mean(err), min(err), max(err), max(abs(env$f(pv, order = 1)[ii])))
}

# Short check, item 5, second part. The limit's fixed tolerances: a strong AR
# that is well identified on a long series, and a small sderr that is small only
# in the data's units. M5: nlf(zi ~ ...) and nlf(shape ~ ...) by the Laplace value
# on brms's own program at chosen points, no translator.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
capt <- function(expr) {
  msgs <- character(0)
  wr <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e), message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }, warning = function(w) {
    wr <<- c(wr, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(v = v, m = msgs, w = wr)
}
cat("== a strong AR, long series, 3 seeds\n")
for (sd in 1:3) for (n in c(2000, 8000)) {
  set.seed(sd)
  ar <- .9995; s <- .05
  e <- numeric(n)
  e[1] <- rnorm(1, 0, s / sqrt(1 - ar^2))
  for (i in 2:n) e[i] <- ar * e[i - 1] + rnorm(1, 0, s)
  d <- data.frame(t = seq_len(n), g = factor(1), x = rnorm(n))
  d$y <- rpois(n, exp(1 + .3 * d$x + e))
  r <- capt(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
  f <- r$v
  tt <- f$estimates$theta
  cy <- capt(frm(y ~ x + cosy(t, g), family = poisson(), data = d))$v
  sm <- summary(f)$cor_pars
  vv <- sqrt(diag(as.matrix(vcov(f, full = TRUE))))
  cat(sprintf("seed %d n %d: ar %.5f sderr %.4f (truth .9995, .05) SE log sderr %.3f chart %.1f code %d  ll %.1f  plain ll %.1f  Est.Error shown %s  says cosy %s\n",
              sd, n, tt[2] / sqrt(1 + tt[2]^2), exp(tt[1]), vv["theta_1"], vv["theta_2"], f$opt$convergence, as.numeric(logLik(f)),
              as.numeric(logLik(capt(frm(y ~ x, family = poisson(), data = d))$v)), format(sm[1, 2]), any(grepl("read as cosy", r$m))))
}
cat("\n== a small sderr in small units: Gamma(log) has none; a skew_normal in 1e-4 units\n")
for (u in c(1, 1e-2, 1e-4)) {
  set.seed(7)
  dd <- expand.grid(t = 1:8, g = factor(1:25))
  dd <- dd[order(dd$g, dd$t), ]
  dd$x <- rnorm(nrow(dd))
  e <- as.vector(replicate(25, arima.sim(list(ar = .6), 8, sd = .5)))
  dd$y <- u * (.3 * dd$x + e + rnorm(nrow(dd), 0, .3))
  r <- capt(frm(y ~ x + ar(t, g, cov = TRUE), family = skew_normal(), data = dd))
  f <- r$v
  if (inherits(f, "error")) { cat("u", u, "ERR", substr(conditionMessage(f), 1, 120), "\n"); next }
  tt <- f$estimates$theta
  sm <- summary(f)$cor_pars
  cat(sprintf("units %g: sderr/units %.3f ar %.3f  Est.Error sderr %s  msgs: %s  warns: %s\n", u, exp(tt[1]) / u, tt[2] / sqrt(1 + tt[2]^2), format(sm["sderr", 2]),
              substr(paste(r$m, collapse = " | "), 1, 160), substr(paste(r$w, collapse = " | "), 1, 120)))
}
cat("\n== M5, by the Laplace value on brms's program\n")
suppressPackageStartupMessages({library(rstan)})
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
set.seed(167)
d <- expand.grid(t = 1:4, g = factor(1:8))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$y <- rnbinom(n, size = 3, mu = exp(1 + .3 * d$x))
d$z <- d$y * rbinom(n, 1, .75)
lap <- function(fb, fam_b, fr, fam_f, a0) {
  pri <- c(brms::set_prior("", class = "sderr"), brms::set_prior("", class = "Intercept"))
  code <- brms::stancode(fb, d, family = fam_b, prior = pri)
  sdat <- brms::standata(fb, d, family = fam_b, prior = pri)
  m <- hel$brms_stan_model(code)
  sf <- rstan::sampling(m, data = sdat, chains = 0, refresh = 0)
  f <- capt(frm(fr, family = fam_f, data = d))$v
  env <- f$obj$env
  print(names(f$opt$par))
  vals <- list()
  for (p in list(list(b0 = .5, b = .2, s = .6, ar = .4, a = a0), list(b0 = 1, b = -.1, s = .3, ar = -.5, a = a0 + .4))) {
    pl <- list(b = array(p$b, 1), Intercept = p$b0 + p$b * mean(d$x), ar = array(p$ar, 1), zerr = rep(0, n), sderr = p$s, b_a = array(p$a, 1))
    u <- rstan::unconstrain_pars(sf, pl)
    pl1 <- pl; pl1$zerr <- rep(1, n)
    iz <- which(rstan::unconstrain_pars(sf, pl1) != u)
    gr <- function(z) { u[iz] <- z; as.numeric(rstan::grad_log_prob(sf, u, adjust_transform = FALSE)[iz]) }
    lp <- function(z) { u[iz] <- z; as.numeric(rstan::log_prob(sf, u, adjust_transform = FALSE)) }
    z <- rep(0, n)
    hs <- function(z, h = 1e-5) { H <- matrix(0, n, n); for (i in 1:n) { e <- rep(0, n); e[i] <- h; H[i, ] <- (gr(z + e) - gr(z - e)) / (2 * h) }; (H + t(H)) / 2 }
    for (it in 1:100) {
      g <- gr(z)
      if (max(abs(g)) < 1e-9) break
      st <- solve(-hs(z), g)
      a <- 1
      while (a > 1e-6 && lp(z + a * st) < lp(z)) a <- a / 2
      z <- z + a * st
    }
    lb <- lp(z) + n / 2 * log(2 * pi) - .5 * as.numeric(determinant(-hs(z))$modulus)
    par <- f$opt$par
    nm <- names(par)
    ib <- which(nm == "beta")
    # beta holds mu's (Intercept, x), then a, in the frame's order
    lab <- unlist(lapply(f$frame$linpreds, function(l) colnames(l$X)))
    par[nm == "beta"][1:2] <- c(p$b0, p$b)
    par[setdiff(seq_along(par), which(nm == "theta"))][3] <- p$a
    par[nm == "theta"] <- c(log(p$s), p$ar / sqrt(1 - p$ar^2))
    env$last.par.best[env$random] <- 0
    lf <- -as.numeric(f$obj$fn(par))
    cat(sprintf("brms %.8f frm %.8f diff %.2e\n", lb, lf, lf - lb))
  }
}
lap(brms::bf(y ~ x + ar(t, g, cov = TRUE), brms::nlf(shape ~ exp(a)), a ~ 1), brms::brmsfamily("negbinomial"),
    frmtmb::bf(y ~ x + ar(t, g, cov = TRUE)) + frmtmb::nlf(shape ~ exp(a)) + frmtmb::lf(a ~ 1), negbinomial(), 1)
lap(brms::bf(z ~ x + ar(t, g, cov = TRUE), brms::nlf(zi ~ exp(a) / (1 + exp(a))), a ~ 1), brms::brmsfamily("zero_inflated_poisson"),
    frmtmb::bf(z ~ x + ar(t, g, cov = TRUE)) + frmtmb::nlf(zi ~ exp(a) / (1 + exp(a))) + frmtmb::lf(a ~ 1), zero_inflated_poisson(), -1)
r <- capt(frm(frmtmb::bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(t, g, cov = TRUE)), family = poisson(), data = d))
cat("nl mu with autocor: ", if (inherits(r$v, "error")) substr(conditionMessage(r$v), 1, 200) else "fits", "\n")
r <- capt(frm(frmtmb::bf(y ~ x + ar(t, g, cov = TRUE)) + frmtmb::nlf(mu ~ a) + frmtmb::lf(a ~ 1), family = poisson(), data = d))
cat("nlf(mu): ", if (inherits(r$v, "error")) substr(conditionMessage(r$v), 1, 200) else "fits", "\n")

# Item 1: the chart identities by hand, reading neither the translator nor
# frmtmb's chart functions. For one family per chart: brms's log_prob at
# hand-mapped parameters against R's own density, and frmtmb's taped joint
# against R's density plus a hand-written Jacobian. Seed 9101.
source("dev/review-refuse-lib.R")
set.seed(9101)
n <- 120
x <- rnorm(n)
mk <- function(y) { d <- data.frame(x = x, y = y); d$y[c(5, 17, 33, 60, 88, 101)] <- NA; d }
cases <- list(
  Gamma = list(
    y = rgamma(n, 2.5, 2.5 / exp(0.3 + 0.5 * x)),
    ffrm = stats::Gamma(link = "log"), fbr = brms::brmsfamily("gamma", "log"),
    dp = "shape",
    ll = function(y, eta, p) dgamma(y, p, p / exp(eta), log = TRUE),
    lj = function(y) log(y), u = function(y) log(y)),
  beta = list(
    y = rbeta(n, plogis(0.2 + 0.6 * x) * 6, (1 - plogis(0.2 + 0.6 * x)) * 6),
    ffrm = frmtmb::Beta(), fbr = brms::Beta(), dp = "phi",
    ll = function(y, eta, p) {
      m <- plogis(eta); dbeta(y, m * p, (1 - m) * p, log = TRUE)
    },
    lj = function(y) log(y) + log1p(-y), u = function(y) qlogis(y)),
  von_mises = list(
    y = ((2 * atan(0.3 + 0.5 * x) + rnorm(n, 0, 0.5) + pi) %% (2 * pi)) - pi,
    ffrm = frmtmb::von_mises(), fbr = brms::von_mises(), dp = "kappa",
    ll = function(y, eta, p) {
      m <- 2 * atan(eta)
      p * cos(y - m) - log(2 * pi) - log(besselI(p, 0))
    },
    # y = 2 pi plogis(u) - pi: dy/du = 2 pi p (1 - p)
    lj = function(y) {
      p <- (y + pi) / (2 * pi); log(2 * pi) + log(p) + log1p(-p)
    },
    u = function(y) qlogis((y + pi) / (2 * pi))))

for (nm in names(cases)) {
  cs <- cases[[nm]]
  d <- mk(cs$y)
  fit <- frm(bf(y | mi() ~ x), family = cs$ffrm, data = d)
  env <- fit$obj$env
  bb <- brms::bf(y | mi() ~ x)
  prior <- E$brms_flat_prior(bb, data = d, family = cs$fbr)
  code <- brms::make_stancode(bb, data = d, family = cs$fbr, prior = prior)
  sdat <- E$brms_standata(bb, data = d, family = cs$fbr, prior = prior)
  sf <- suppressMessages(rstan::sampling(E$brms_stan_model(code), data = sdat,
                                         chains = 0))
  rows <- which(is.na(d$y))
  pts <- list(c(0.3, 0.5, 1.5), c(-0.2, 0.8, 3), c(0.6, 0.1, 0.8))
  for (i in seq_along(pts)) {
    b0 <- pts[[i]][1]; b1 <- pts[[i]][2]; dpv <- pts[[i]][3]
    set.seed(9101 + i)
    ymis <- switch(nm,
      Gamma = rgamma(length(rows), 2, 1),
      beta = rbeta(length(rows), 2, 2),
      von_mises = runif(length(rows), -3, 3))
    yall <- d$y; yall[rows] <- ymis
    eta <- b0 + b1 * d$x
    hand <- sum(cs$ll(yall, eta, dpv))
    pars <- list(b = array(b1, 1L), Intercept = b0 + b1 * mean(d$x),
                 Ymi = array(ymis, length(ymis)))
    pars[[cs$dp]] <- dpv
    up <- rstan::unconstrain_pars(sf, pars)
    lp <- rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
    # the taped joint, from the brms-free side: beta, then log dpar, then u
    p2 <- env$last.par.best
    p2[names(p2) == "beta"] <- c(b0, b1)
    p2[names(p2) == "betad"] <- log(dpv)
    p2[names(p2) == "miss"] <- cs$u(ymis)
    ours <- -env$f(p2)
    jac <- sum(cs$lj(ymis))
    pr("%-10s pt %d: brms lp - R density %10.3e; taped joint - (R density + hand Jacobian) %10.3e [Jacobian %.4f]",
       nm, i, lp - hand, ours - hand - jac, jac)
  }
  # the mode: frmtmb's inner solution must solve Stan's own transformed
  # gradient, by Stan's adjust_transform = TRUE on Ymi, at frmtmb's estimate
  est <- fit$estimates
  ymode <- switch(nm, Gamma = exp(est$miss), beta = plogis(est$miss),
                  von_mises = 2 * pi * plogis(est$miss) - pi)
  pars <- list(b = array(est$beta[[2]], 1L),
               Intercept = est$beta[[1]] + est$beta[[2]] * mean(d$x),
               Ymi = array(ymode, length(ymode)))
  pars[[cs$dp]] <- exp(est$betad[[1]])
  up <- rstan::unconstrain_pars(sf, pars)
  g <- rstan::grad_log_prob(sf, up, adjust_transform = TRUE)
  pr("%-10s Stan-transformed inner gradient at frmtmb's mode, max |g[Ymi]| %.3e",
     nm, {p3 <- pars; p3$Ymi <- p3$Ymi * 0.9
       max(abs(g[rstan::unconstrain_pars(sf, p3) != up]))})
}

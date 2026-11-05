# Short check 2, item 3: writes review-latres-fresh.R from review-latres-m1b.R,
# with a set of designs neither side used to calibrate the warning.
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-m1b.R")
txt <- paste(src, collapse = "\n")
new <- '
D0 <- list(
  f_pois05 = list(fam = quote(poisson()), dist = "pois", b0 = -.7, b1 = .4, s = .7, ar = .3),
  f_pois_ar8 = list(fam = quote(poisson()), dist = "pois", b0 = .7, b1 = .4, s = .3, ar = .8),
  f_pois_short = list(fam = quote(poisson()), dist = "pois", G = 60, K = 4, b0 = .5, b1 = .4, s = .6, ar = .5),
  f_pois_long = list(fam = quote(poisson()), dist = "pois", G = 6, K = 40, b0 = 0, b1 = .4, s = .5, ar = .7),
  f_bin2 = list(fam = quote(binomial()), dist = "binom", tr = 2, b0 = 0, b1 = .5, s = .7, ar = .5),
  f_bin6 = list(fam = quote(binomial()), dist = "binom", tr = 6, b0 = -1, b1 = .5, s = .6, ar = .5),
  f_bin30 = list(fam = quote(binomial()), dist = "binom", tr = 30, b0 = -1.5, b1 = .4, s = .5, ar = .6),
  f_probit = list(fam = quote(bernoulli("probit")), dist = "bern", b0 = 0, b1 = .4, s = .6, ar = .6),
  f_cloglog = list(fam = quote(bernoulli("cloglog")), dist = "bern", b0 = -.5, b1 = .4, s = .8, ar = .5),
  f_geom = list(fam = quote(geometric()), dist = "geom", b0 = 1, b1 = .4, s = .5, ar = .6),
  f_nb1 = list(fam = quote(negbinomial()), dist = "nb", shape = 1, b0 = 1, b1 = .4, s = .5, ar = .6),
  f_nb8 = list(fam = quote(negbinomial()), dist = "nb", shape = 8, b0 = 2, b1 = .4, s = .8, ar = .4),
  f_gam1 = list(fam = quote(Gamma(link = "log")), dist = "gamma", shape = 1, b0 = 0, b1 = .4, s = .5, ar = .6),
  f_gam2s1 = list(fam = quote(Gamma(link = "log")), dist = "gamma", shape = 2, b0 = 0, b1 = .4, s = 1, ar = .5),
  f_gam20 = list(fam = quote(Gamma(link = "log")), dist = "gamma", shape = 20, b0 = 0, b1 = .4, s = .3, ar = .7),
  f_wei2 = list(fam = quote(weibull()), dist = "wei", shape = 2, b0 = 0, b1 = .4, s = .5, ar = .5),
  f_exp = list(fam = quote(exponential()), dist = "expo", b0 = 0, b1 = .4, s = .7, ar = .5),
  f_beta2 = list(fam = quote(Beta()), dist = "beta", phi = 2, b0 = 0, b1 = .4, s = .5, ar = .5),
  f_beta20 = list(fam = quote(Beta()), dist = "beta", phi = 20, b0 = 0, b1 = .4, s = .4, ar = .7),
  f_zinb = list(fam = quote(zero_inflated_negbinomial()), dist = "zinb", zi = .2, shape = 3, b0 = 1.2, b1 = .4, s = .5, ar = .5),
  f_hnb = list(fam = quote(hurdle_negbinomial()), dist = "hnb", hu = .3, shape = 3, b0 = 1, b1 = .4, s = .5, ar = .5),
  f_cumpro = list(fam = quote(cumulative("probit")), dist = "cum", pfun = pnorm, rfun = rnorm, thr = c(-1, 0, 1), b0 = 0, b1 = .4, s = .6, ar = .5),
  f_ln1 = list(fam = quote(lognormal()), dist = "ln", sigma = 1, b0 = 0, b1 = .4, s = .5, ar = .5))
'
a <- regexpr("D0 <- list\\(", txt)
b <- regexpr("D <- D0\\[\\[design\\]\\]", txt)
txt <- paste0(substr(txt, 1, a - 1), new, substr(txt, b, nchar(txt)))
txt <- sub("G <- 25\nK <- 8", "G <- D$G %||% 25\nK <- D$K %||% 8", txt, fixed = TRUE)
# densities
txt <- sub("    zip = if (y == 0)", paste0(
  "    pois = dpois(y, mu), bern = dbinom(y, 1, mu), geom = dnbinom(y, size = 1, mu = mu),\n",
  "    gamma = dgamma(y, p$shape, p$shape / mu), wei = dweibull(y, p$shape, mu / gamma(1 + 1 / p$shape)), expo = dexp(y, 1 / mu),\n",
  "    zinb = if (y == 0) p$zi + (1 - p$zi) * dnbinom(0, size = p$shape, mu = mu) else (1 - p$zi) * dnbinom(y, size = p$shape, mu = mu),\n",
  "    hnb = if (y == 0) p$hu + 0 * mu else (1 - p$hu) * dnbinom(y, size = p$shape, mu = mu) / (1 - dnbinom(0, size = p$shape, mu = mu)),\n",
  "    zip = if (y == 0)"), txt, fixed = TRUE)
txt <- sub("plogis(th[y + 1] - eta) - plogis(th[y] - eta)", "(D$pfun %||% plogis)(th[y + 1] - eta) - (D$pfun %||% plogis)(th[y] - eta)", txt, fixed = TRUE)
# simulation
txt <- sub("    zip = rpois(n, mu) * rbinom(n, 1, 1 - D$zi),", paste0(
  "    pois = rpois(n, mu), bern = rbinom(n, 1, mu), geom = rnbinom(n, size = 1, mu = mu),\n",
  "    gamma = rgamma(n, D$shape, D$shape / mu), wei = rweibull(n, D$shape, mu / gamma(1 + 1 / D$shape)), expo = rexp(n, 1 / mu),\n",
  "    zinb = rnbinom(n, size = D$shape, mu = mu) * rbinom(n, 1, 1 - D$zi),\n",
  "    hnb = { pos <- rbinom(n, 1, 1 - D$hu); pmax(qnbinom(runif(n, dnbinom(0, size = D$shape, mu = mu), 1), size = D$shape, mu = mu), 1) * pos },\n",
  "    zip = rpois(n, mu) * rbinom(n, 1, 1 - D$zi),"), txt, fixed = TRUE)
txt <- sub("cum = findInterval(eta + rlogis(n), D$thr) + 1L,", "cum = findInterval(eta + (D$rfun %||% rlogis)(n), D$thr) + 1L,", txt, fixed = TRUE)
# parameters
txt <- sub("    zip = , zib = list(zi = plogis(p[1])), hp = list(hu = plogis(p[1])),", paste0(
  "    pois = , bern = , geom = , expo = list(), gamma = , wei = list(shape = exp(p[1])),\n",
  "    zinb = list(zi = plogis(p[1]), shape = exp(p[2])), hnb = list(hu = plogis(p[1]), shape = exp(p[2])),\n",
  "    zip = , zib = list(zi = plogis(p[1])), hp = list(hu = plogis(p[1])),"), txt, fixed = TRUE)
txt <- sub("    zip = , zib = qlogis(D$zi), hp = qlogis(D$hu),", paste0(
  "    pois = , bern = , geom = , expo = numeric(0), gamma = , wei = log(D$shape),\n",
  "    zinb = c(qlogis(D$zi), log(D$shape)), hnb = c(qlogis(D$hu), log(D$shape)),\n",
  "    zip = , zib = qlogis(D$zi), hp = qlogis(D$hu),"), txt, fixed = TRUE)
writeLines(txt, "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-fresh.R")

# The reviewer's own L1 specs, instrument format. Seeds in the 9xxx range,
# distinct from the lane's.
base_dat <- function(seed, n = 200) {
  d <- eq_dat(seed = seed, n = n, ngrp = 10)
  set.seed(seed + 1L)
  d$yn <- 1 + 0.5 * d$x + rnorm(n, 0, 0.5)
  d$pos <- exp(0.2 + 0.3 * d$x + rnorm(n, 0, 0.4))
  d$cnt0 <- rpois(n, exp(0.5 + 0.4 * d$x))
  d$zi0 <- ifelse(runif(n) < 0.3, 0L, d$cnt0)
  d$hl <- ifelse(runif(n) < 0.25, 0, d$pos)
  d$bt <- plogis(0.3 * d$x + rnorm(n, 0, 0.7))
  d$zo <- ifelse(runif(n) < 0.1, 0, ifelse(runif(n) < 0.1, 1, d$bt))
  d$zb <- ifelse(runif(n) < 0.2, 0, d$bt)
  d$tr <- 12L
  d$bb <- rbinom(n, 12, plogis(0.2 + 0.3 * d$x))
  d$ang <- atan2(sin(0.5 * d$x + rnorm(n, 0, 0.3)), cos(0.5 * d$x + rnorm(n, 0, 0.3)))
  d$ce <- ifelse(d$pos > quantile(d$pos, 0.85), "right", "none")
  d$posc <- pmin(d$pos, quantile(d$pos, 0.85))
  d$ng <- rnbinom(n, mu = exp(0.5 + 0.4 * d$x), size = 2)
  d
}
ord_bd <- function(seed) {
  d <- eq_dat("ordinal", seed = seed)
  d
}
l1_specs <- list()
add <- function(nm, f) l1_specs[[nm]] <<- f
body <- "a + b * x"
ln <- function(y, fam, extra = "", aterm = "", d, frm_args = NULL, seed) {
  list(bf = sprintf("bf(%s%s ~ b0 + b1 * x, b0 + b1 ~ 1%s, nl = TRUE)", y, aterm, extra), family = fam,
       data = d, frm_args = frm_args)
}
# item 7: compat rows with dpars other than mu, zi, truncation, censoring
add("r_zip_zi", function() ln("zi0", "zero_inflated_poisson()", ", zi ~ z", d = base_dat(9001)))
add("r_zinb_zi", function() ln("zi0", "zero_inflated_negbinomial()", ", zi ~ z, shape ~ z", d = base_dat(9002)))
add("r_hurdle_ln", function() ln("hl", "hurdle_lognormal()", ", hu ~ z, sigma ~ z", d = base_dat(9003)))
add("r_zoib", function() ln("zo", "zero_one_inflated_beta()", ", phi ~ z, zoi ~ z, coi ~ 1", d = base_dat(9004)))
add("r_zib", function() ln("zb", "zero_inflated_beta()", ", phi ~ z, zi ~ z", d = base_dat(9005)))
add("r_student_sig", function() ln("yn", "student()", ", sigma ~ z", d = base_dat(9006)))
add("r_trunc_gauss", function() ln("yn", "gaussian()", aterm = " | trunc(lb = -5, ub = 6)", d = base_dat(9007)))
add("r_trunc_pois", function() ln("cnt0", "poisson()", aterm = " | trunc(ub = 25)", d = base_dat(9008)))
add("r_cens_weibull", function() ln("posc", "weibull()", aterm = " | cens(ce)", d = base_dat(9009)))
add("r_gev", function() ln("yn", "gen_extreme_value()", ", sigma ~ z", d = base_dat(9010)))
add("r_betabin", function() ln("bb", "beta_binomial()", aterm = " | trials(tr)", extra = ", phi ~ z", d = base_dat(9011)))
add("r_skew", function() ln("yn", "skew_normal()", ", alpha ~ z", d = base_dat(9012)))
add("r_vonmises", function() ln("ang", "von_mises()", ", kappa ~ z", d = base_dat(9013)))
add("r_exgauss", function() ln("pos", "exgaussian()", ", beta ~ z", d = base_dat(9014)))
add("r_shifted", function() ln("pos", "shifted_lognormal()", ", sigma ~ z", d = base_dat(9015)))
add("r_compois", function() ln("cnt0", "brmsfamily('com_poisson')", d = base_dat(9016)))
add("r_dw", function() ln("cnt0", "discrete_weibull()", ", shape ~ z", d = base_dat(9017)))
add("r_hgamma", function() ln("hl", "hurdle_gamma()", ", hu ~ z, shape ~ z", d = base_dat(9018)))
add("r_xbeta", function() ln("bt", "xbeta()", ", phi ~ z", d = base_dat(9019)))
add("r_hnb", function() ln("zi0", "hurdle_negbinomial()", ", hu ~ z", d = base_dat(9020)))
add("r_frechet", function() ln("pos", "frechet()", ", nu ~ z", d = base_dat(9021)))
add("r_weights", function() { d <- base_dat(9022); d$w <- runif(nrow(d), 0.5, 2); ln("yn", "gaussian()", aterm = " | weights(w)", d = d) })
add("r_nonlin_gamma_log", function() list(bf = "bf(pos ~ a * exp(b * x), a + b ~ 1, shape ~ z, nl = TRUE)", family = "Gamma('log')", data = base_dat(9023),
    frm_args = list(start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.2)))))
# item 3: multi-mu on the reviewer's specs
mm <- function(seed) {
  set.seed(seed); n <- 200
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3), ntr = 10L)
  d$cat <- factor(sapply(seq_len(n), function(i) sample(c("a", "b", "c"), 1, prob = c(1, exp(0.8 * d$x[i]), exp(-0.5 * d$x[i])))))
  m <- t(vapply(seq_len(n), function(i) stats::rmultinom(1, 10, c(1, exp(0.8 * d$x[i]), exp(-0.5 * d$x[i])))[, 1], numeric(3)))
  colnames(m) <- c("a", "b", "c")
  d$cnt <- m
  s <- (m + 0.5) / 11.5
  colnames(s) <- c("a", "b", "c")
  d$sim <- s
  d$yg <- ifelse(runif(n) < 0.4, -1, 1) + 0.5 * d$x + rnorm(n, 0, 0.5)
  d$yc <- rpois(n, exp(0.5 + 0.2 * d$x))
  d
}
add("r_cat_perloc", function() list(bf = "bf(cat ~ 1) + nlf(mub ~ ab + bb * x) + nlf(muc ~ ac + bc * x) + lf(ab + bb + ac + bc ~ 1)", family = "categorical()", data = mm(9101)))
add("r_mn_perloc", function() list(bf = "bf(cnt | trials(ntr) ~ 1) + nlf(mub ~ ab + bb * x) + nlf(muc ~ ac + bc * x) + lf(ab + bb + ac + bc ~ 1)", family = "multinomial()", data = mm(9102)))
add("r_dir_perloc", function() list(bf = "bf(sim ~ 1) + nlf(mub ~ ab + bb * x) + nlf(muc ~ ac + bc * x) + lf(ab + bb + ac + bc ~ 1)", family = "dirichlet()", data = mm(9103)))
add("r_ln_perloc", function() list(bf = "bf(sim ~ 1) + nlf(mub ~ ab + bb * x) + nlf(muc ~ ac + bc * x) + lf(ab + bb + ac + bc ~ 1)", family = "logistic_normal()", data = mm(9104)))
add("r_mix_perloc", function() list(bf = "bf(yg ~ 1) + nlf(mu1 ~ a1 + b1 * x) + nlf(mu2 ~ a2 + b2 * x) + lf(a1 + b1 + a2 + b2 ~ 1)", family = "mixture(gaussian(), gaussian())", data = mm(9105),
    frm_args = list(start = list(beta = c("a1_(Intercept)" = -1, "b1_(Intercept)" = 0.3, "a2_(Intercept)" = 1, "b2_(Intercept)" = 0.3)))))
add("r_mix_theta", function() list(bf = "bf(yg ~ b0 + b1 * x, theta1 ~ z, b0 + b1 ~ 1, nl = TRUE)", family = "mixture(gaussian(), gaussian())", data = mm(9106)))
add("r_mix_sigma", function() list(bf = "bf(yg ~ b0 + b1 * x, sigma2 ~ z, b0 + b1 ~ 1, nl = TRUE)", family = "mixture(gaussian(), gaussian())", data = mm(9107)))
add("r_mix_pois", function() list(bf = "bf(yc ~ b0 + b1 * x, theta1 ~ 1, b0 + b1 ~ 1, nl = TRUE)", family = "mixture(poisson(), poisson())", data = mm(9108)))
add("r_cat_refcat", function() list(bf = "bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", family = "categorical(refcat = 'c')", data = mm(9109)))
add("r_cat_prior", function() list(bf = "bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", family = "categorical()", data = mm(9110),
    prior = "set_prior('normal(0, 1)', nlpar = 'b1')"))
# item 1: a body with its own intercept, a body from default starts
add("r_ord_ridge", function() list(bf = "bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", family = "cumulative()", data = ord_bd(9201)))
add("r_ord_default", function() list(bf = "bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)", family = "cumulative()", data = ord_bd(9202)))
# item 5
add("r_powzero", function() { d <- base_dat(9301); d$xz <- abs(d$x); d$xz[1:5] <- 0; d$yp <- 2 * d$xz^0.5 + rnorm(nrow(d), 0, 0.3)
  list(bf = "bf(yp ~ a * xz^b, a + b ~ 1, nl = TRUE)", family = "gaussian()", data = d,
       frm_args = list(start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 1)))) })
# item 6
add("r_stan_funs", function() { d <- base_dat(9401); d$p <- plogis(d$x)
  list(bf = "bf(yn ~ a + log1p_exp(b * x) + c * inv_cloglog(z) + e * Phi(x) + f * log_sum_exp(x, z), a + b + c + e + f ~ 1, nl = TRUE)", family = "gaussian()", data = d,
       frm_args = list(start = list(beta = c("b_(Intercept)" = 0.2)))) })
add("r_student_ctl", function() list(bf = "bf(yn ~ x, sigma ~ z)", family = "student()", data = base_dat(9006)))
add("r_student_real", function() { d <- base_dat(9006); set.seed(90061); d$yn <- 1 + 0.5 * d$x + rt(nrow(d), 4) * 0.5; ln("yn", "student()", ", sigma ~ z", d = d) })
add("r_stan_funs2", function() { d <- base_dat(9401)
  list(bf = "bf(yn ~ a0 + log1p_exp(b0 * x) + c0 * inv_cloglog(z) + e0 * Phi(x) + g0 * log_sum_exp(x, z), a0 + b0 + c0 + e0 + g0 ~ 1, nl = TRUE)", family = "gaussian()", data = d,
       frm_args = list(start = list(beta = c("b0_(Intercept)" = 0.2)))) })
add("r_hurdle_disc", function() { d <- eq_dat("ordinal", seed = 4952); d$y <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
  list(bf = "bf(y ~ x, disc ~ 0 + z, hu ~ z)", family = "hurdle_cumulative()", data = d) })
add("r_ext_wiener", function() { set.seed(9501); n <- 200; d <- data.frame(x = rnorm(n)); d$rt <- 0.3 + rlnorm(n, -0.5, 0.4); d$dc <- rbinom(n, 1, plogis(0.8 * d$x))
  list(bf = "bf(rt | dec(dc) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", family = "brmsfamily('wiener')", data = d) })
add("r_ext_cox", function() { set.seed(9502); n <- 200; d <- data.frame(x = rnorm(n)); d$time <- rexp(n, exp(0.5 * d$x)); d$cc <- rbinom(n, 1, 0.2)
  list(bf = "bf(time | cens(cc) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", family = "brmsfamily('cox')", data = d) })
add("r_ext_cox_lin", function() { set.seed(9502); n <- 200; d <- data.frame(x = rnorm(n)); d$time <- rexp(n, exp(0.5 * d$x)); d$cc <- rbinom(n, 1, 0.2)
  list(bf = "bf(time | cens(cc) ~ x)", family = "brmsfamily('cox')", data = d) })
bb_od <- function(seed, n = 200, size = 10) { set.seed(seed); d <- data.frame(x = round(rnorm(n), 3), ntr = as.integer(size)); p <- plogis(0.3 + 0.5 * d$x); d$y <- rbinom(n, d$ntr, rbeta(n, p * 8, (1 - p) * 8)); d }
add("r_bb_cens", function() { d <- bb_od(9601); d$ce <- ifelse(d$y >= 8, "right", ifelse(d$y <= 1, "left", "none")); d$yc <- pmin(pmax(d$y, 1L), 8L)
  list(bf = "bf(yc | trials(ntr) + cens(ce) ~ x)", family = "beta_binomial()", data = d) })
add("r_bb_cens_body", function() { d <- bb_od(9602); d$ce <- ifelse(d$y >= 8, "right", ifelse(d$y <= 1, "left", "none")); d$yc <- pmin(pmax(d$y, 1L), 8L)
  list(bf = "bf(yc | trials(ntr) + cens(ce) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)", family = "beta_binomial()", data = d) })
add("r_bb_trunc", function() { d <- bb_od(9603); d <- d[d$y >= 1 & d$y <= 8, ]; list(bf = "bf(y | trials(ntr) + trunc(lb = 1, ub = 8) ~ x)", family = "beta_binomial()", data = d) })
add("r_bb_nood", function() { set.seed(9); n <- 300; d <- data.frame(x = round(rnorm(n), 3), ntr = 10L); d$y <- rbinom(n, 10, plogis(0.3 + 0.5 * d$x)); list(bf = "bf(y | trials(ntr) ~ x)", family = "beta_binomial()", data = d) })
add("r_bb_join", function() { d <- bb_od(9604); d$ntr <- rep(c(24L, 25L), length.out = nrow(d)); set.seed(96041); p <- plogis(0.3 + 0.5 * d$x); d$y <- rbinom(nrow(d), d$ntr, rbeta(nrow(d), p * 8, (1 - p) * 8)); list(bf = "bf(y | trials(ntr) ~ x)", family = "beta_binomial()", data = d) })
add("r_zibb_cens", function() { d <- bb_od(9605); d$y[runif(nrow(d)) < 0.2] <- 0L; d$ce <- ifelse(d$y >= 8, "right", "none"); d$yc <- pmin(d$y, 8L); list(bf = "bf(yc | trials(ntr) + cens(ce) ~ x)", family = "zero_inflated_beta_binomial()", data = d) })

# Specs for dev/nlfix-l1.R, in the instrument's format (dev/equiv-specs.R).
# Evaluated in a child of the instrument's environment, so eq_dat() is in
# scope. Each seed is the number after the name's family prefix below.
nl_dat <- function(seed, n = 200, ngrp = 10) {
  d <- eq_dat(seed = seed, n = n, ngrp = ngrp)
  d$y <- 2 * exp(0.3 * d$x) + rnorm(nrow(d), 0, 0.5)
  d
}
ord_dat <- function(seed, hurdle = FALSE) {
  d <- eq_dat("ordinal", seed = seed)
  if (hurdle) {
    d$y <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)),
                   levels = 0:4)
  }
  d
}
l1_specs <- list()
add <- function(nm, f) l1_specs[[nm]] <<- f

# D1: every ordinal family, plain, with thres(), disc, a nonlinear nlpar,
# and a random effect in an nlpar
ord_fams <- c("cumulative", "sratio", "cratio", "acat")
for (i in seq_along(ord_fams)) local({
  fam <- ord_fams[i]
  k <- i
  add(paste0("ord_", fam, "_plain"), function() list(
    bf = "bf(y ~ b1 * x, b1 ~ 1, nl = TRUE)",
    family = paste0(fam, "()"), data = ord_dat(4900 + k)))
  add(paste0("ord_", fam, "_thres"), function() list(
    bf = "bf(y | thres(3) ~ b1 * x, b1 ~ 1, nl = TRUE)",
    family = paste0(fam, "()"), data = ord_dat(4910 + k)))
  add(paste0("ord_", fam, "_disc"), function() list(
    bf = "bf(y ~ b1 * x, b1 ~ 1, disc ~ 0 + z, nl = TRUE)",
    family = paste0(fam, "()"), data = ord_dat(4920 + k)))
  add(paste0("ord_", fam, "_nonlin"), function() list(
    bf = "bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE)",
    family = paste0(fam, "()"), data = ord_dat(4930 + k),
    frm_args = list(start = list(beta = c("b1_(Intercept)" = 0.5,
                                          "b2_(Intercept)" = 0.5)))))
  add(paste0("ord_", fam, "_re"), function() list(
    bf = "bf(y ~ b1 * x, b1 ~ 1 + (1 | g), nl = TRUE)",
    family = paste0(fam, "()"), data = ord_dat(4940 + k)))
})
add("ord_hurdle_plain", function() list(
  bf = "bf(y ~ b1 * x, b1 ~ 1, nl = TRUE)",
  family = "hurdle_cumulative()", data = ord_dat(4951, hurdle = TRUE)))
add("ord_hurdle_disc", function() list(
  bf = "bf(y ~ b1 * x, b1 ~ 1, disc ~ 0 + z, hu ~ z, nl = TRUE)",
  family = "hurdle_cumulative()", data = ord_dat(4952, hurdle = TRUE)))
add("ord_cumulative_probit", function() list(
  bf = "bf(y ~ b1 * x, b1 ~ 1, nl = TRUE)",
  family = "cumulative(link = 'probit')", data = ord_dat(4953)))
add("ord_equidistant", function() list(
  bf = "bf(y ~ b1 * x, b1 ~ 1, nl = TRUE)",
  family = "cumulative(threshold = 'equidistant')", data = ord_dat(4954)))

# D7: every core autocor term under an nl body
ac_terms <- list(
  ar_cov = "ar(time = t, gr = g, cov = TRUE)",
  ar = "ar(time = t, gr = g, cov = FALSE)",
  ar2_cov = "ar(time = t, gr = g, p = 2, cov = TRUE)",
  ma_cov = "ma(time = t, gr = g, cov = TRUE)",
  ma = "ma(time = t, gr = g, cov = FALSE)",
  arma_cov = "arma(time = t, gr = g, cov = TRUE)",
  arma = "arma(time = t, gr = g, cov = FALSE)",
  cosy = "cosy(time = t, gr = g)")
for (i in seq_along(ac_terms)) local({
  nm <- names(ac_terms)[i]
  tm <- ac_terms[[i]]
  k <- i
  add(paste0("ac_", nm), function() list(
    bf = paste0("bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, ",
                "autocor = ~ ", tm, ")"),
    data = nl_dat(1580 + k),
    frm_args = list(start = list(beta = c("b1_(Intercept)" = 1,
                                          "b2_(Intercept)" = 0.1)))))
})
add("ac_unstr", function() list(
  bf = paste0("bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, ",
              "autocor = ~ unstr(time = t, gr = g))"),
  data = nl_dat(1590, n = 40, ngrp = 10),
  frm_args = list(start = list(beta = c("b1_(Intercept)" = 1,
                                        "b2_(Intercept)" = 0.1)))))
add("ac_fcor", function() {
  d <- nl_dat(158)
  list(bf = paste0("bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, ",
                   "autocor = ~ fcor(M))"),
       data = d, data2 = list(M = diag(nrow(d))))
})
add("ac_fcor_dense", function() {
  d <- nl_dat(159)
  n <- nrow(d)
  M <- 0.3^abs(outer(seq_len(n), seq_len(n), "-"))
  list(bf = paste0("bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, ",
                   "autocor = ~ fcor(M))"),
       data = d, data2 = list(M = M))
})
add("ac_student", function() {
  list(bf = paste0("bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, ",
                   "autocor = ~ ar(time = t, gr = g, cov = TRUE))"),
       family = "student()", data = nl_dat(1591),
       frm_args = list(start = list(beta = c("b1_(Intercept)" = 1,
                                             "b2_(Intercept)" = 0.1))))
})

# item 3: the default links, from the default starts
pos_dat <- function(seed) {
  d <- eq_dat(seed = seed)
  d$y <- exp(0.4 + 0.3 * d$x + rnorm(nrow(d), 0, 0.4))
  d
}
add("start_gamma", function() list(
  bf = "bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "Gamma()", data = pos_dat(3001)))
add("start_ig", function() list(
  bf = "bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "inverse.gaussian()", data = pos_dat(3002)))

# item 4: the vignette's bodies
add("fun_guess33", function() {
  set.seed(4001)
  d <- data.frame(ability = rnorm(200))
  d$y <- rbinom(200, 1, 0.33 + 0.67 * plogis(1.2 * d$ability))
  list(bf = "bf(y ~ 0.33 + 0.67 * inv_logit(eta), eta ~ ability, nl = TRUE)",
       family = "bernoulli('identity')", data = d)
})
add("fun_guess", function() {
  set.seed(4002)
  d <- data.frame(ability = rnorm(300))
  d$y <- rbinom(300, 1, 0.25 + 0.75 * plogis(1.2 * d$ability))
  list(bf = paste0("bf(y ~ guess + (1 - guess) * inv_logit(eta), ",
                   "eta ~ 0 + ability, guess ~ 1, nl = TRUE)"),
       family = "bernoulli('identity')", data = d,
       frm_args = list(start = list(beta = c("guess_(Intercept)" = 0.3))))
})

# multi-mu families, nl = TRUE: brms copies the body to every location
# predictor and the nl parameters are shared
mm_dat <- function(seed) {
  set.seed(seed)
  n <- 200
  d <- data.frame(x = round(rnorm(n), 3), ntr = 10L)
  d$cat <- factor(sapply(seq_len(n), function(i) {
    sample(c("a", "b", "c"), 1, prob = c(1, exp(0.4 * d$x[i]), 1))
  }))
  m <- t(vapply(seq_len(n), function(i) {
    stats::rmultinom(1, 10, c(1, exp(0.4 * d$x[i]), 1))[, 1]
  }, numeric(3)))
  colnames(m) <- c("a", "b", "c")
  d$cnt <- m
  s <- (m + 0.5) / 11.5
  colnames(s) <- c("a", "b", "c")
  d$sim <- s
  d$yg <- ifelse(runif(n) < 0.4, -1, 1) + 0.5 * d$x + rnorm(n, 0, 0.5)
  d
}
add("mm_categorical", function() list(
  bf = "bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "categorical()", data = mm_dat(5001)))
add("mm_multinomial", function() list(
  bf = "bf(cnt | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "multinomial()", data = mm_dat(5002)))
add("mm_dirichlet", function() list(
  bf = "bf(sim ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "dirichlet()", data = mm_dat(5003)))
add("mm_dirmult", function() list(
  bf = "bf(cnt | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "dirichlet_multinomial()", data = mm_dat(5004)))
add("mm_logistic_normal", function() list(
  bf = "bf(sim ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "logistic_normal()", data = mm_dat(5005)))
add("mm_mixture", function() list(
  bf = "bf(yg ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "mixture(gaussian(), gaussian())", data = mm_dat(5006)))
add("mv_nl", function() {
  d <- nl_dat(5007)
  d$z2 <- 0.3 + 0.5 * d$x + rnorm(nrow(d))
  list(bf = paste0("mvbf(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), ",
                   "bf(z2 ~ x), rescor = FALSE)"), data = d,
       frm_args = list(start = list(beta = c("y_b1_(Intercept)" = 1,
                                             "y_b2_(Intercept)" = 0.1))))
})

# second batch: the rest of the D1 variants and the audit
add("ord_cosy_data", function() NULL)
l1_specs[["ord_cosy_data"]] <- NULL
add("ac_cosy2", function() {
  d <- nl_dat(1592)
  u <- rnorm(nlevels(d$g), 0, 0.5)
  d$y <- d$y + u[as.integer(d$g)]
  list(bf = paste0("bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE, ",
                   "autocor = ~ cosy(time = t, gr = g))"),
       data = d,
       frm_args = list(start = list(beta = c("b1_(Intercept)" = 1,
                                             "b2_(Intercept)" = 0.1))))
})
add("ord_grouped", function() {
  d <- ord_dat(4955)
  d$gg <- factor(ifelse(as.integer(d$g) <= 5, "p", "q"))
  list(bf = "bf(y | thres(gr = gg) ~ b1 * x, b1 ~ 1, nl = TRUE)",
       family = "cumulative()", data = d)
})
add("ord_sum_to_zero", function() list(
  bf = "bf(y ~ b1 * x, b1 ~ 1, nl = TRUE)",
  family = "cumulative(threshold = 'sum_to_zero')", data = ord_dat(4956)))
add("ord_mixture", function() list(
  bf = "bf(y ~ b1 * x, b1 ~ 1, nl = TRUE)",
  family = "mixture(cumulative(), cumulative())", data = ord_dat(4957)))
add("ord_mv", function() {
  d <- ord_dat(4958)
  d$z2 <- 0.2 + 0.5 * d$x + rnorm(nrow(d))
  list(bf = paste0("mvbf(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE, ",
                   "family = cumulative()), bf(z2 ~ x, family = gaussian()), ",
                   "rescor = FALSE)"),
       family = NULL, data = d)
})
add("fam_asym_laplace", function() {
  d <- nl_dat(6001)
  list(bf = "bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
       family = "asym_laplace()", data = d)
})

# third batch: the methods at the same parameter values (L2)
for (nm_i in setdiff(names(l1_specs), c("ac_ar2_cov", "ac_cosy"))) local({
  nm <- nm_i
  base <- l1_specs[[nm]]
  l1_specs[[paste0(nm, "_m")]] <<- function() {
    sp <- base()
    sp$method <- c("fitted", "simulate", "predict")
    sp
  }
})

# fourth batch
ord_exp <- function(seed) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n),
                     c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  d
}
for (nmf in ord_fams) local({
  fam <- nmf
  add(paste0("ord_", fam, "_exp"), function() list(
    bf = "bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE)",
    family = paste0(fam, "()"), data = ord_exp(4960),
    frm_args = list(start = list(beta = c("b1_(Intercept)" = 1,
                                          "b2_(Intercept)" = 0.3)))))
})
add("mm_dirichlet2", function() list(
  bf = "bf(sim ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "dirichlet2()", data = mm_dat(5008)))

# fifth batch: controls for the simulate shapes, on linear formulas
add("ctl_hurdle_m", function() list(
  bf = "bf(y ~ x)", family = "hurdle_cumulative()",
  data = ord_dat(4951, hurdle = TRUE), method = c("fitted", "simulate")))
add("ctl_multinomial_m", function() list(
  bf = "bf(cnt | trials(ntr) ~ x)", family = "multinomial()",
  data = mm_dat(5002), method = c("fitted", "simulate")))
add("ctl_dirichlet_m", function() list(
  bf = "bf(sim ~ x)", family = "dirichlet()",
  data = mm_dat(5003), method = c("fitted", "simulate")))
add("ctl_mv_m", function() {
  d <- nl_dat(5007)
  d$z2 <- 0.3 + 0.5 * d$x + rnorm(nrow(d))
  list(bf = "mvbf(bf(y ~ x), bf(z2 ~ x), rescor = FALSE)", data = d,
       method = c("fitted", "simulate"))
})
add("mm_logistic_normal_fit", function() {
  sp <- l1_specs[["mm_logistic_normal"]]()
  sp$method <- "fitted"
  sp
})
add("ord_hurdle_disc_again", function() l1_specs[["ord_hurdle_disc"]]())

# sixth batch
add("ctl_hurdle_disc", function() list(
  bf = "bf(y ~ x, disc ~ 0 + z, hu ~ z)",
  family = "hurdle_cumulative()", data = ord_dat(4952, hurdle = TRUE)))
add("ord_hurdle_hu", function() list(
  bf = "bf(y ~ b1 * x, b1 ~ 1, hu ~ z, nl = TRUE)",
  family = "hurdle_cumulative()", data = ord_dat(4952, hurdle = TRUE)))
add("mm_dirichlet2b", function() list(
  bf = "bf(sim ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  family = "brmsfamily('dirichlet2')", data = mm_dat(5008)))

# fix round 1: from the DEFAULT start, the bodies the review found on a saddle
sad_dat <- function(hurdle = FALSE) {
  set.seed(9202)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n),
                     c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  if (hurdle) {
    set.seed(1)
    d$y <- ordered(ifelse(runif(n) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
  }
  d
}
for (nmf in c(ord_fams, 'hurdle_cumulative')) local({
  fam <- nmf
  add(paste0('sad_', fam), function() list(
    bf = 'bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)',
    family = paste0(fam, '()'), data = sad_dat(fam == 'hurdle_cumulative')))
})
add('sad_probit', function() list(
  bf = 'bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)',
  family = "cumulative('probit')", data = sad_dat()))
add('sad_sat', function() list(
  bf = 'bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)',
  family = 'cumulative()', data = sad_dat()))
add('sad_gaussian', function() {
  set.seed(9203)
  d <- data.frame(x = round(rnorm(250), 3))
  d$y <- 1 + 0.8 * d$x + rnorm(250, 0, 0.5)
  list(bf = 'bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE)', data = d)
})
add('sad_poisson', function() {
  set.seed(9203)
  d <- data.frame(x = round(rnorm(250), 3))
  d$y <- rpois(250, exp(0.5 + 0.5 * d$x))
  list(bf = 'bf(y ~ a + b1 * b2 * x, a + b1 + b2 ~ 1, nl = TRUE)',
       family = 'poisson()', data = d)
})
add('powzero', function() {
  set.seed(9302)
  n <- 150
  d <- data.frame(x = c(0, 0, runif(n - 2, 0, 3)))
  d$y <- 2 * d$x^0.5 + rnorm(n, 0, 0.3)
  list(bf = 'bf(y ~ a * x^b, a + b ~ 1, nl = TRUE)', data = d,
       frm_args = list(start = list(beta = c(1, 0.5))))
})
add('powzero_fun', function() {
  set.seed(9302)
  n <- 150
  d <- data.frame(x = c(0, 0, runif(n - 2, 0, 3)))
  d$y <- 2 * d$x^0.5 + rnorm(n, 0, 0.3)
  list(bf = 'bf(y ~ a * pow(x, b), a + b ~ 1, nl = TRUE)', data = d,
       frm_args = list(start = list(beta = c(1, 0.5))))
})

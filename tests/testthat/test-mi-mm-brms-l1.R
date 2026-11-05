# Gated: log_prob of brms 2.23.0's program against frmtmb's density, for mi()
# with a value no density reads (the complete-case fit) and with one a density
# reads (the joint at its conditional modes), on each chart; for mm() with a
# key and with a known covariance; for se() on skew_normal; for cens() on a
# mixture. Compiles through FRMTMB_STAN_CACHE.

l1_mi_data <- function(gen, seed, n = 150, nmiss = 12) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  d$y <- gen(n, 0.5 + 0.4 * d$x, d$x)
  d$full <- d$y
  d$y[sample(n, nmiss)] <- NA
  d
}

l1_read_check <- function(fam_frm, fam_brm, gen, seed) {
  d <- l1_mi_data(gen, seed)
  d$z <- rnorm(nrow(d))
  d$resp <- d$z + 0.5 * d$full + rnorm(nrow(d), 0, 0.5)
  fit <- allow_warnings(
    frm(mvbf(bf(resp ~ mi(y) + z), bf(y | mi() ~ x, family = fam_frm),
             rescor = FALSE), data = d),
    "Laplace approximation of the")
  fb <- brms::bf(resp ~ mi(y) + z) + brms::bf(y | mi() ~ x, family = fam_brm) +
    brms::set_rescor(FALSE)
  brms_lp_check(fb, NULL, d, fit, joint = TRUE)
}

gen_gamma <- function(n, eta, x) rgamma(n, 3, 3 / exp(eta))
gen_lognormal <- function(n, eta, x) rlnorm(n, eta, 0.6)
gen_beta <- function(n, eta, x) rbeta(n, plogis(eta) * 8, (1 - plogis(eta)) * 8)
gen_vm <- function(n, eta, x) {
  ((2 * atan(0.5 * x) + rnorm(n, 0, 0.5) + pi) %% (2 * pi)) - pi
}
gen_exg <- function(n, eta, x) eta + rnorm(n, 0, 0.5) + rexp(n, 1)
gen_vm_end <- function(n, eta, x) {
  ((3.05 + 0.05 * x + rnorm(n, 0, 0.6) + pi) %% (2 * pi)) - pi
}

test_that("a value no density reads: frm() is brms on the observed rows", {
  skip_unless_brms_fit()
  for (cs in list(list(stats::Gamma(link = "log"),
                       brms::brmsfamily("gamma", "log"), gen_gamma, 171),
                  list(Beta(), brms::brmsfamily("beta"), gen_beta, 173),
                  list(von_mises(), brms::brmsfamily("von_mises"), gen_vm,
                       174),
                  list(exgaussian(), brms::brmsfamily("exgaussian"), gen_exg,
                       175))) {
    d <- l1_mi_data(cs[[3L]], cs[[4L]])
    fit <- frm(bf(y | mi() ~ x), family = cs[[1L]], data = d)
    expect_null(fit$estimates$miss)
    brms_lp_check(brms::bf(y | mi() ~ x), cs[[2L]], d, fit)
  }
})

test_that("brms's log_prob minus the density at Ymi is the logLik, at any Ymi", {
  skip_unless_brms_fit()
  d <- l1_mi_data(gen_gamma, 171)
  fit <- frm(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d)
  fe <- frmtmb:::fixef_by_dpar(fit)
  dens <- function(v, f, rows) {
    mu <- exp(fe$mu[[1L]] + fe$mu[[2L]] * d$x[rows])
    sh <- exp(fe$shape[[1L]])
    stats::dgamma(v, sh, sh / mu, log = TRUE)
  }
  ymi <- cbind(rep(0.5, 12), rep(1.7, 12), seq(0.2, 5, length.out = 12))
  gaps <- brms_mi_unread_gap(brms::bf(y | mi() ~ x),
                             brms::brmsfamily("gamma", "log"), d, fit, ymi,
                             dens)
  expect_lt(max(abs(gaps)) / abs(as.numeric(logLik(fit))), 1e-10)
})

test_that("a read value under a log chart is brms's model", {
  skip_unless_brms_fit()
  l1_read_check(stats::Gamma(link = "log"), brms::brmsfamily("gamma", "log"),
                gen_gamma, 171)
  l1_read_check(lognormal(), brms::brmsfamily("lognormal"), gen_lognormal, 172)
})

test_that("a read value under the logit chart is brms's model", {
  skip_unless_brms_fit()
  l1_read_check(Beta(), brms::brmsfamily("beta"), gen_beta, 173)
})

test_that("a read value on the window of von_mises is brms's model", {
  skip_unless_brms_fit()
  l1_read_check(von_mises(), brms::brmsfamily("von_mises"), gen_vm, 174)
  # near the end of the window, where a value on the line left brms's support
  # and Stan could not read it
  l1_read_check(von_mises(), brms::brmsfamily("von_mises"), gen_vm_end, 11400)
})

test_that("a read value on the line is brms's model", {
  skip_unless_brms_fit()
  l1_read_check(exgaussian(), brms::brmsfamily("exgaussian"), gen_exg, 175)
})

mm_l1_data <- function(seed = 70, n = 160) {
  set.seed(seed)
  lv <- letters[1:8]
  d <- data.frame(x = rnorm(n), x1 = rnorm(n), x2 = rnorm(n),
                  g1 = factor(sample(lv, n, TRUE), levels = lv),
                  g2 = factor(sample(lv, n, TRUE), levels = lv))
  u <- stats::setNames(rnorm(8, 0, 0.8), lv)
  d$y <- 1 + 0.5 * d$x + 0.5 * (u[as.character(d$g1)] + u[as.character(d$g2)]) +
    rnorm(n, 0, 0.5)
  d
}

test_that("mm() with an |ID| key is brms's model", {
  skip_unless_brms_fit()
  d <- mm_l1_data()
  fit <- frm(bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2))),
             data = d)
  brms_lp_check(brms::bf(y ~ x + (1 | p | mm(g1, g2)),
                         sigma ~ 1 + (1 | p | mm(g1, g2))),
                stats::gaussian(), d, fit, joint = TRUE)
})

test_that("mm(cov = A) is brms's model", {
  skip_unless_brms_fit()
  d <- mm_l1_data(seed = 87)
  lv <- letters[1:8]
  A <- 0.5^abs(outer(1:8, 1:8, "-"))
  dimnames(A) <- list(lv, lv)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2, cov = A))), data = d,
             data2 = list(A = A))
  brms_lp_check(brms::bf(y ~ x + (1 | mm(g1, g2, cov = A))),
                stats::gaussian(), d, fit, joint = TRUE,
                data2 = list(A = A))
})

test_that("se() on skew_normal is brms's model", {
  skip_unless_brms_fit()
  set.seed(4040)
  n <- 160
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n) + abs(rnorm(n))
  d$s <- 0.2 + 0.3 * abs(d$z)
  fit <- frm(bf(y | se(s) ~ x), family = skew_normal(), data = d)
  brms_lp_check(brms::bf(y | se(s) ~ x), brms::brmsfamily("skew_normal"), d,
                fit)
})

test_that("cens() on a mixture is brms's model", {
  skip_unless_brms_fit()
  set.seed(7)
  n <- 200
  d <- data.frame(x = rnorm(n))
  k <- rbinom(n, 1, 0.4)
  d$y <- ifelse(k == 1, 3 + 0.5 * d$x, -3 - 0.5 * d$x) + rnorm(n)
  d$ce <- 0L
  d$ce[d$y > 3.5] <- 1L
  d$ce[d$y < -3.5] <- -1L
  fit <- frm(bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x),
             family = mixture(gaussian(), gaussian()), data = d)
  # brms orders same-family components unless told not to. The bound on the
  # gradient is the optimizer's stop here, 1e-3 of an unpolished fit
  brms_lp_check(brms::bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x),
                brms::mixture(stats::gaussian(), stats::gaussian(),
                              order = "none"), d, fit, tol_grad = 5e-3)
})

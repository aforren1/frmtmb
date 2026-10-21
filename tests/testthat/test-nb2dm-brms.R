# negbinomial2 and dirichlet_multinomial against the Stan program brms
# 2.23.0 generates for them, at frmtmb's estimates: brms_lp_check() in
# helper-brms.R (checks A, B and, with a group term, C; see the head of
# test-brms-likelihood.R). With flat priors brms's log density is the
# log-likelihood, so the admitted constant is 0.
#
# brms writes negbinomial2 as neg_binomial_2_log_glm_lpmf(Y | Xc,
# Intercept, b, inv(sigma)) when sigma has no formula, and as
# neg_binomial_2_log_lpmf(Y | mu, inv(sigma)) when it has one, so both
# spellings are checked. dirichlet_multinomial is brms's own
# dirichlet_multinomial_logit2_lpmf().
#
# Stan compiles here, so the file is opt-in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")

nbdm_tight <- frmtmb_control(grad_tol = 1e-6, restarts = 3)

nbdm_nb_data <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                   g = factor(rep(1:15, length.out = n)))
  dd$e <- exp(stats::runif(n))
  u <- stats::rnorm(15, 0, 0.4)[dd$g]
  dd$y <- stats::rnbinom(n, mu = exp(0.8 + 0.5 * dd$x + u),
                         size = 1 / exp(-0.8 + 0.3 * dd$z))
  dd
}

nbdm_dm_data <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                   g = factor(rep(1:15, length.out = n)))
  E <- cbind(0, 0.3 + 0.6 * dd$x + stats::rnorm(15, 0, 0.4)[dd$g],
             -0.4 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  phi <- exp(1.6 + 0.3 * dd$z)
  dd$N <- sample(4:25, n, TRUE)
  Yc <- t(vapply(seq_len(n), function(i) {
    p <- stats::rgamma(3, P[i, ] * phi[i])
    as.numeric(stats::rmultinom(1, dd$N[i], p / sum(p)))
  }, numeric(3)))
  colnames(Yc) <- c("y1", "y2", "y3")
  dd$Yc <- Yc
  dd
}

test_that("negbinomial2 is brms's likelihood, in both of brms's spellings", {
  skip_unless_brms_fit()
  dd <- nbdm_nb_data(61)
  fit <- frm(y ~ x, family = brmsfamily("negbinomial2"), data = dd,
             control = nbdm_tight)
  brms_lp_check(brms::bf(y ~ x), brms::brmsfamily("negbinomial2"), dd, fit)
  f2 <- frm(bf(y ~ x, sigma ~ z), family = brmsfamily("negbinomial2"),
            data = dd, control = nbdm_tight)
  brms_lp_check(brms::bf(y ~ x, sigma ~ z),
                brms::brmsfamily("negbinomial2"), dd, f2)
})

test_that("negbinomial2 under rate() and other links is brms's likelihood", {
  skip_unless_brms_fit()
  dd <- nbdm_nb_data(62)
  fit <- frm(y | rate(e) ~ x, family = brmsfamily("negbinomial2"), data = dd,
             control = nbdm_tight)
  brms_lp_check(brms::bf(y | rate(e) ~ x), brms::brmsfamily("negbinomial2"),
                dd, fit)
  f2 <- frm(bf(y ~ x, sigma ~ z),
            family = brmsfamily("negbinomial2", "softplus",
                                link_sigma = "softplus"),
            data = dd, control = nbdm_tight)
  brms_lp_check(brms::bf(y ~ x, sigma ~ z),
                brms::brmsfamily("negbinomial2", "softplus",
                                 link_sigma = "softplus"), dd, f2)
})

test_that("negbinomial2 with (1 | g), on the joint density", {
  skip_unless_brms_fit()
  dd <- nbdm_nb_data(63)
  fit <- frm(y ~ x + (1 | g), family = brmsfamily("negbinomial2"), data = dd,
             control = nbdm_tight)
  brms_lp_check(brms::bf(y ~ x + (1 | g)), brms::brmsfamily("negbinomial2"),
                dd, fit, joint = TRUE)
})

test_that("dirichlet_multinomial is brms's likelihood, refcat and phi ~ z", {
  skip_unless_brms_fit()
  dd <- nbdm_dm_data(64)
  fit <- frm(Yc | trials(N) ~ x, family = dirichlet_multinomial(), data = dd,
             control = nbdm_tight)
  brms_lp_check(brms::bf(Yc | trials(N) ~ x), brms::dirichlet_multinomial(),
                dd, fit)
  f2 <- frm(bf(Yc | trials(N) ~ x, phi ~ z),
            family = dirichlet_multinomial(refcat = "y2"), data = dd,
            control = nbdm_tight)
  brms_lp_check(brms::bf(Yc | trials(N) ~ x, phi ~ z),
                brms::dirichlet_multinomial(refcat = "y2"), dd, f2)
})

test_that("dirichlet_multinomial with (1 | ID | g), on the joint density", {
  # (1 | ID | g) and not (1 | g): the joint-density translator misreads
  # uncorrelated per-category blocks, multinomial()'s too, so check C
  # cannot speak for that spelling (dev/nb2dm-findings.md); frm_sample()
  # of (1 | g) agrees with brms's posterior there instead
  skip_unless_brms_fit()
  dd <- nbdm_dm_data(65)
  fit <- frm(Yc | trials(N) ~ x + (1 | ID | g),
             family = dirichlet_multinomial(), data = dd,
             control = nbdm_tight)
  brms_lp_check(brms::bf(Yc | trials(N) ~ x + (1 | ID | g)),
                brms::dirichlet_multinomial(), dd, fit, joint = TRUE)
})

test_that("dirichlet_multinomial inside mvbf() is brms's likelihood", {
  skip_unless_brms_fit()
  dd <- nbdm_dm_data(66)
  dd$w <- 0.5 * dd$x + stats::rnorm(nrow(dd))
  fit <- frm(mvbf(bf(Yc | trials(N) ~ x, family = dirichlet_multinomial()),
                  bf(w ~ x, family = gaussian())), data = dd,
             control = nbdm_tight)
  brms_lp_check(brms::bf(Yc | trials(N) ~ x,
                         family = brms::dirichlet_multinomial()) +
                  brms::bf(w ~ x, family = gaussian()) +
                  brms::set_rescor(FALSE),
                list(brms::dirichlet_multinomial(), gaussian()), dd, fit)
})

# The simplex families' likelihood against the Stan program brms 2.23.0
# generates for them, at frmtmb's estimates: brms_lp_check() in
# helper-brms.R (checks A, B and, with group effects, C; see the head of
# test-brms-likelihood.R). With flat priors brms's log density is the
# log-likelihood itself, so the admitted constant is 0. brms declares
# the logistic-normal correlations as their Cholesky factor Llncor,
# which stan_pars_from_fit() reads off frmtmb's lncor chart.
#
# Stan compiles here, so the file is opt-in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")

lp_tight_sx <- frmtmb_control(grad_tol = 1e-6, restarts = 3)

sx_brms_data <- function(seed, n = 300, K = 3, ng = 0) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  E <- cbind(0, sapply(seq_len(K - 1L), function(k) {
    (-1)^k * (0.2 + 0.3 * k) + (0.5 - 0.3 * k) * dd$x
  }))
  if (ng > 0) {
    dd$g <- factor(rep(seq_len(ng), length.out = n))
    U <- matrix(stats::rnorm(ng * (K - 1L), 0, 0.4), ng)
    E[, -1L] <- E[, -1L] + U[as.integer(dd$g), ]
  }
  P <- exp(E - apply(E, 1L, max))
  P <- P / rowSums(P)
  A <- P * exp(2 + 0.4 * dd$z)
  G <- matrix(stats::rgamma(length(A), shape = A), n)
  Y <- G / rowSums(G)
  colnames(Y) <- letters[seq_len(K)]
  dd$Y <- Y
  dd
}

test_that("row 16h: dirichlet with phi ~ z is brms's likelihood", {
  skip_unless_brms_fit()
  dd <- sx_brms_data(31)
  fit <- frm(bf(Y ~ x, phi ~ z), family = dirichlet(), data = dd,
             control = lp_tight_sx)
  brms_lp_check(brms::bf(Y ~ x, phi ~ z), brms::dirichlet(), dd, fit)
})

test_that("row 16i: dirichlet with refcat and a category's own formula", {
  skip_unless_brms_fit()
  dd <- sx_brms_data(32, K = 4)
  fit <- frm(bf(Y ~ x, mua ~ 1), family = dirichlet(refcat = "c"),
             data = dd, control = lp_tight_sx)
  brms_lp_check(brms::bf(Y ~ x, mua ~ 1), brms::dirichlet(refcat = "c"),
                dd, fit)
})

test_that("row 16j: dirichlet with (1 | ID | g), on the joint density", {
  skip_unless_brms_fit()
  dd <- sx_brms_data(33, n = 400, ng = 15)
  fit <- frm(bf(Y ~ x + (1 | ID | g)), family = dirichlet(), data = dd,
             control = lp_tight_sx)
  brms_lp_check(brms::bf(Y ~ x + (1 | ID | g)), brms::dirichlet(), dd, fit,
                joint = TRUE)
})

test_that("row 16k: logistic_normal with sigma<k> ~ z and three lncor", {
  skip_unless_brms_fit()
  dd <- sx_brms_data(34, K = 4)
  fit <- frm(bf(Y ~ x, sigmab ~ z), family = logistic_normal(), data = dd,
             control = lp_tight_sx)
  expect_identical(rownames(summary(fit)$spec_pars),
                   c("sigmac", "sigmad", "lncor(b,c)", "lncor(b,d)",
                     "lncor(c,d)"))
  brms_lp_check(brms::bf(Y ~ x, sigmab ~ z), brms::logistic_normal(), dd,
                fit)
})

test_that("row 16l: logistic_normal with refcat and (1 | ID | g)", {
  skip_unless_brms_fit()
  dd <- sx_brms_data(35, n = 400, ng = 15)
  fit <- frm(bf(Y ~ x + (1 | ID | g)),
             family = logistic_normal(refcat = "b"), data = dd,
             control = lp_tight_sx)
  brms_lp_check(brms::bf(Y ~ x + (1 | ID | g)),
                brms::logistic_normal(refcat = "b"), dd, fit, joint = TRUE)
})

test_that("row 16m: brmsfamily(\"dirichlet2\") is brms's likelihood", {
  skip_unless_brms_fit()
  dd <- sx_brms_data(36)
  fit <- frm(bf(Y ~ x), family = brmsfamily("dirichlet2"), data = dd,
             control = lp_tight_sx)
  brms_lp_check(brms::bf(Y ~ x), brms::brmsfamily("dirichlet2"), dd, fit)
})

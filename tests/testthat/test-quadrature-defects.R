# The quadrature defect cluster the grammar fuzzer and the
# compatibility-registry probes turned up at v0.23.0
# (dev/fuzz-findings.md N1-N3, dev/test-backlog.md).
#
# Every one of them traces back to the same property of TMBad's
# marginal_gk transform: it rescales each integrand ONCE, at whatever
# parameter values the template holds when MakeADFun tapes it, and
# bakes that (mu, sigma) pair in as a constant. The fix is to fit the
# plain Laplace objective first and tape the marginalized one at its
# optimum, which also gives back the conditional modes that the
# marginalized objective no longer carries.

#' @srrstats {G5.3} Return objects that should contain no missing or
#'   undefined values are explicitly tested for their absence. These
#'   tests sweep the whole post-fit surface after a quadrature fit,
#'   asserting `anyNA()` is false for the conditional modes, `fitted()`,
#'   `ranef()`, `predict()` on new data, and `residuals()`, and that the
#'   log-likelihood is finite. The same assertion appears for
#'   `predict(se.fit = TRUE)` in `test-autoscale.R`,
#'   `test-gp-multidim.R`, `test-review-v25.R`, and
#'   `test-method-residue.R`, for one-step-ahead residuals across
#'   families in `test-osa-inference.R`, and for Cook's distance in
#'   `test-review-fixes.R`. Where `NA` is the correct answer, the tests
#'   are paired: the value must be `NA` in exactly the expected positions
#'   and non-`NA` everywhere else. The fuzz tier applies the same rule to
#'   every generated model, treating a non-finite estimate or
#'   log-likelihood as a defect unless the fit also warned.
#' @noRd
NULL

quad_nested_data <- function(seed = 4, ng = 20, nt = 5) {
  set.seed(seed)
  n <- ng * nt
  d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                  gb = factor(rep(c("i", "ii"), length.out = n)),
                  x = rnorm(n))
  d$ga <- d$g
  d$eta <- 0.5 + 0.4 * d$x + rnorm(ng, 0, 0.4)[d$g] +
    rnorm(ng * 2, 0, 0.3)[as.integer(d$ga) * 2 + as.integer(d$gb) - 2]
  d
}

test_that("quadrature reports every conditional mode, not just the first", {
  set.seed(4)
  ng <- 20; nt <- 5; n <- ng * nt
  d <- data.frame(g = factor(rep(seq_len(ng), each = nt)), x = rnorm(n))
  d$y <- 1 + 0.4 * d$x + rnorm(ng, 0, 0.5)[d$g] + rnorm(n)

  fq <- frm(bf(y ~ 1 + x + (1 | g)) + gaussian(), data = d,
            quadrature = TRUE)
  fl <- frm(bf(y ~ 1 + x + (1 | g)) + gaussian(), data = d)

  # parList() under integrate= used to leave 19 of 20 modes NA and slide
  # an outer value into the twentieth
  expect_false(anyNA(fq$estimates$b))
  expect_equal(length(fq$estimates$b), ng)
  expect_false(anyNA(fitted(fq)))
  expect_false(anyNA(unlist(ranef(fq))))
  expect_false(anyNA(predict(fq, newdata = d)))
  expect_false(anyNA(residuals(fq)))

  # Gauss-Kronrod and Laplace agree exactly for a gaussian response, so
  # the modes must agree too
  expect_loglik_equal(fq, fl, tol = 1e-6)
  expect_vector_equal(fq$estimates$b, fl$estimates$b, tol = 1e-5)
})

test_that("quadrature modes match glmer(nAGQ = 25) where Laplace does not", {
  skip_if_not_installed("lme4")
  set.seed(601)
  ng <- 100; per <- 4
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 1.2)
  dd <- data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.7 * x + u[g])),
                   x = x, g = g)

  fq <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
            quadrature = TRUE)
  fl <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd)
  ref <- lme4::glmer(y ~ x + (1 | g), dd, family = binomial, nAGQ = 25)

  expect_false(anyNA(fq$estimates$b))
  expect_vector_equal(fq$estimates$b, lme4::ranef(ref)$g[[1]], tol = 1e-3)
  # the modes are recomputed at the quadrature optimum, not copied from
  # the Laplace pre-fit: those two differ by more than 0.1 here
  expect_gt(max(abs(fq$estimates$b - fl$estimates$b)), 0.05)
})

test_that("quadrature survives non-gaussian families and nested blocks", {
  # every one of these satisfied the documented scalar-intercept guard
  # and still died at a bare "NA/NaN gradient evaluation": the frozen
  # rescaling taped at the cold start overflowed
  d <- quad_nested_data()
  n <- nrow(d)
  cases <- list(
    list(fam = poisson(), y = function() rpois(n, exp(d$eta))),
    list(fam = Gamma(link = "log"),
         y = function() rgamma(n, shape = 2, scale = exp(d$eta) / 2)),
    list(fam = Beta(), y = function() {
      p <- plogis(d$eta); rbeta(n, p * 5, (1 - p) * 5)
    })
  )
  for (cs in cases) {
    for (re in c("(1 | g)", "(1 | ga/gb)")) {
      set.seed(11)
      dd <- d
      dd$y <- cs$y()
      fo <- stats::as.formula(paste("y ~ 1 + x +", re))
      fit <- frm(bf(fo) + cs$fam, data = dd, quadrature = TRUE)
      lab <- paste(cs$fam$family, re)
      expect_s3_class(fit, "frmtmb_fit")
      expect_true(is.finite(as.numeric(logLik(fit))), label = lab)
      expect_false(anyNA(fit$estimates$b))
      expect_false(anyNA(fitted(fit)))
      expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-2)
    }
  }
})

test_that("quadrature over a nested block agrees with the Laplace fit", {
  # gaussian: Laplace is exact, so the two marginalizations must land on
  # the same optimum even with two scalar blocks in the integral
  d <- quad_nested_data()
  set.seed(11)
  d$y <- d$eta + rnorm(nrow(d))
  fq <- frm(bf(y ~ 1 + x + (1 | ga/gb)) + gaussian(), data = d,
            quadrature = TRUE)
  fl <- frm(bf(y ~ 1 + x + (1 | ga/gb)) + gaussian(), data = d)
  expect_loglik_equal(fq, fl, tol = 1e-4)
  expect_vector_equal(fq$estimates$b, fl$estimates$b, tol = 1e-4)
})

test_that("quadrature refuses trunc() instead of returning logLik = +Inf", {
  set.seed(4)
  ng <- 20; nt <- 5; n <- ng * nt
  d <- data.frame(g = factor(rep(seq_len(ng), each = nt)), x = rnorm(n))
  d$y <- 1 + 0.4 * d$x + rnorm(ng, 0, 0.5)[d$g] + rnorm(n)
  d <- d[d$y > -0.6, ]

  expect_error(
    frm(bf(y | trunc(lb = -0.6) ~ 1 + x + (1 | g)) + gaussian(),
        data = d, quadrature = TRUE),
    "trunc\\(\\)")
  expect_equal(frm_compat("quadrature", "trunc()")$status, "refused")

  # the same model under Laplace keeps a finite objective, which is why
  # the refusal names it as the remedy
  fl <- frm(bf(y | trunc(lb = -0.6) ~ 1 + x + (1 | g)) + gaussian(),
            data = d)
  expect_true(is.finite(as.numeric(logLik(fl))))
})

test_that("mixture() refuses REML and profile", {
  set.seed(3)
  ng <- 25; nt <- 6; n <- ng * nt
  dm <- data.frame(g = factor(rep(seq_len(ng), each = nt)), x = rnorm(n))
  u <- rnorm(ng, 0, 0.5)
  cls <- rbinom(n, 1, 0.4)
  dm$y <- ifelse(cls == 1, rnorm(n, 3 + u[dm$g], 0.7),
                 rnorm(n, u[dm$g], 0.7))
  fo <- bf(y ~ 1 + x + (1 | g), mu2 ~ 1)
  mixfam <- mixture(gaussian(), gaussian())

  # the components are exchangeable, so the objective REML and profile
  # would expand about a single inner mode is multimodal in the fixed
  # effects. Both used to die at "NA/NaN gradient evaluation" or report
  # a gradient near 1e9 with no guard at all.
  expect_error(frm(fo + mixfam, data = dm, REML = TRUE), "mixture\\(\\)")
  expect_error(frm(fo + mixfam, data = dm,
                   control = frmtmb_control(profile = TRUE)),
               "mixture\\(\\)")
  expect_equal(frm_compat("REML", "mixture")$status, "refused")
  expect_equal(frm_compat("profile", "mixture")$status, "refused")

  # quadrature marginalizes the random effects, not the coefficients,
  # so it stays allowed (test-v19.R pins down its exactness)
  fq <- frm(fo + mixfam, data = dm, quadrature = TRUE)
  expect_s3_class(fq, "frmtmb_fit")
  expect_lt(max(abs(fq$obj$gr(fq$opt$par))), 1e-2)
  expect_false(anyNA(fq$estimates$b))
})

test_that("mixture_mvn() refuses REML and profile too", {
  set.seed(5)
  n <- 200
  z <- rbinom(n, 1, 0.4)
  dz <- data.frame(x = rnorm(n))
  dz$Y <- cbind(rnorm(n, ifelse(z == 1, 3, 0)),
                rnorm(n, ifelse(z == 1, -2, 0)))
  expect_s3_class(frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dz),
                  "frmtmb_fit")
  # each family states the refusal in its own name now: this used to be
  # one message in fit.R listing mixture(), mixture_mvn() and lca()
  expect_error(frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dz,
                   REML = TRUE), "mixture_mvn\\(\\)")
  expect_error(frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dz,
                   control = frmtmb_control(profile = TRUE)),
               "mixture_mvn\\(\\)")
})

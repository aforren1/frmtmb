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
#'   `ranef()`, `frm_linpred()` on new data, and `residuals()`, and that the
#'   log-likelihood is finite. The same assertion appears for
#'   `frm_linpred(se.fit = TRUE)` in `test-autoscale.R`,
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

# Muffle exactly the warning that a quadrature fit's Est.Error leaves out
# the group-effect term; any other warning still reaches the test.
quad_quiet <- function(expr) {
  withCallingHandlers(expr, frmtmb_modes_conditional_se = function(w) {
    invokeRestart("muffleWarning")
  })
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
  # fitted() and residuals() say their Est.Error omits the group-effect
  # term (lane wt-predfix); only the estimates are asserted here
  expect_false(anyNA(quad_quiet(fitted(fq))[, "Estimate"]))
  expect_false(anyNA(unlist(ranef(fq))))
  expect_false(anyNA(frm_linpred(fq, newdata = d)))
  expect_false(anyNA(quad_quiet(residuals(fq))[, "Estimate"]))

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

  fq <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd,
            quadrature = TRUE)
  fl <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd)
  ref <- lme4::glmer(y ~ x + (1 | g), dd, family = binomial, nAGQ = 25)

  expect_false(anyNA(fq$estimates$b))
  expect_vector_equal(fq$estimates$b, lme4::ranef(ref)$g[[1]], tol = 1e-3)
  # the modes are recomputed at the quadrature optimum, not copied from
  # the Laplace pre-fit: those two differ by more than 0.1 here
  expect_gt(max(abs(fq$estimates$b - fl$estimates$b)), 0.05)
})

test_that("quadrature survives non-gaussian families", {
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
    # nested blocks take the iterated rule (test-quadrature-nested.R)
    for (re in c("(1 | g)", "(1 | ga/gb)")) {
      set.seed(11)
      dd <- d
      dd$y <- cs$y()
      fo <- stats::as.formula(paste("y ~ 1 + x +", re))
      # a beta fit can put the sd at zero: on the one-block rule that is
      # a boundary fit, as under Laplace; the nested rule keeps the
      # "flat" verdict, since the boundary check cannot run on that
      # objective (lane setier)
      fit <- allow_boundary(allow_warnings(
        frm(bf(fo) + cs$fam, data = dd, quadrature = TRUE),
        "Standard errors are not available"))
      lab <- paste(cs$fam$family, re)
      expect_s3_class(fit, "frmtmb_fit")
      expect_true(is.finite(as.numeric(logLik(fit))), label = lab)
      expect_false(anyNA(fit$estimates$b))
      expect_false(anyNA(quad_quiet(fitted(fit))[, "Estimate"]))
      expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-2)
    }
  }
})

test_that("quadrature over a nested gaussian block is the exact Laplace fit", {
  # the Gauss-Kronrod rule integrated only the ga effects of the pair and
  # left the ga:gb ones to Laplace; 0.73.0 refused the model. The nested
  # rule integrates both, and for a gaussian response the Laplace
  # integral is exact, so the two fits must agree
  d <- quad_nested_data()
  set.seed(11)
  d$y <- d$eta + rnorm(nrow(d))
  fq <- frm(bf(y ~ 1 + x + (1 | ga/gb)) + gaussian(), data = d,
            quadrature = TRUE)
  fl <- frm(bf(y ~ 1 + x + (1 | ga/gb)) + gaussian(), data = d)
  expect_length(fq$obj$env$random, 0L)
  ll_q <- as.numeric(logLik(fq))
  ll_l <- as.numeric(logLik(fl))
  expect_lt(abs(ll_q - ll_l) / abs(ll_l), 1e-7)
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

# The rule is anchored where the fit was when it was taped, so every
# path that moves the estimate afterwards needs a rule taped again there
# (quad_settle()).

test_that("the stationary-point escape runs before the rule is anchored", {
  # a skew normal started at alpha's declared stationary point: the
  # Laplace warm start stalls there, the rule was anchored at the stall,
  # and the escape on that frozen rule landed 1.0e-2 short of the fit
  # from the default start, against 1e-8 under Laplace
  # (dev/quadseam-punch1c.R, wt-skewinit's acceptance seed 1)
  skip_on_cran()
  set.seed(1)
  n <- 200
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  dd <- data.frame(y = y, xs = xs, g = factor(rep(1:20, length.out = n)))
  m3 <- mean((y - mean(y))^3) / stats::sd(y)^3
  bad <- list(betad = c(log(stats::sd(y)), 2 * sign(m3) + 0.5 * m3))
  fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
  f0 <- allow_boundary(frm(fo, family = skew_normal(), data = dd,
                           quadrature = TRUE))
  f1 <- allow_boundary(frm(fo, family = skew_normal(), data = dd,
                           quadrature = TRUE, start = bad))
  gain <- f1$opt[["stationary_escape"]][["gain"]]
  expect_true(isTRUE(gain > 1))
  gap <- abs(as.numeric(logLik(f1)) - as.numeric(logLik(f0)))
  expect_lt(gap / gain, 1e-6)
})

test_that("quadrature does not stop a runaway at a false maximum", {
  # weights of 1.5 with an observation-level effect: the integral rises
  # without bound as sigma falls. A response with mi() is outside the
  # weights check by design, so the fit reaches the rule. The rule frozen
  # at the Laplace start stopped at log sigma -4.06 with code 0 and no
  # warning, where the exact integral rises by 38.5 per unit of log sigma
  # (dev/quadseam-runaway3.R, the data of lane wcheck's review)
  skip_on_cran()
  set.seed(20261111)
  n <- 80
  d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  d$y[c(5, 17, 33)] <- NA
  d$w <- 1.5
  unsettled <- 0L
  f <- withCallingHandlers(
    allow_warnings(
      frm(bf(y | weights(w) + mi() ~ x + (1 | obs)), family = gaussian(),
          data = d, quadrature = TRUE),
      c("gives every observation its own random effect",
        "Optimizer did not report convergence",
        "Large maximum absolute gradient")),
    frmtmb_quad_unsettled = function(w) {
      unsettled <<- unsettled + 1L
      invokeRestart("muffleWarning")
    })
  expect_identical(unsettled, 1L)
})

test_that("a quadrature fit at its maximum settles in silence", {
  # Three fits that are at their maximum must give no warning of any
  # class. The rule taped at the estimate reads a gradient of calibration
  # noise: at grad_tol = 1e-5 the two binomial fits spent 10 rounds that
  # did not move and then warned that the estimate was not the maximum,
  # where it is the exact maximum to 1e-9 (review of lane quadseam,
  # MAJOR 1). The third is a Bernoulli fit at the default tolerance.
  # Every class is counted, not only frmtmb_quad_unsettled, because the
  # convergence check reads the tape the loop returns, and a gradient of
  # noise there warned too
  count_quad <- function(expr) {
    n <- character()
    v <- withCallingHandlers(expr, warning = function(w) {
      n <<- c(n, class(w)[1L])
      invokeRestart("muffleWarning")
    })
    list(fit = v, n = n)
  }
  set.seed(1)
  d1 <- data.frame(g = factor(rep(1:10, each = 5)), x = rnorm(50))
  u <- rnorm(10, 0, 0.4)
  d1$k <- rbinom(50, 500, plogis(0.3 + 0.5 * d1$x + u[d1$g]))
  d1$nk <- 500 - d1$k
  set.seed(4)
  d2 <- data.frame(g = factor(rep(1:6, each = 4)), x = rnorm(24))
  u <- rnorm(6, 0, 0.5)
  d2$z <- c(-3, -2, -1, 1, 2, 3)[d2$g]
  d2$k <- rbinom(24, 100, plogis(1.5 * d2$z + 0.2 * d2$x + u[d2$g]))
  d2$nk <- 100 - d2$k
  ctl <- frmtmb_control(grad_tol = 1e-5)
  r1 <- count_quad(frm(cbind(k, nk) ~ x + (1 | g), family = binomial(),
                       data = d1, quadrature = TRUE, control = ctl))
  r2 <- count_quad(frm(cbind(k, nk) ~ z + x + (1 | g), family = binomial(),
                       data = d2, quadrature = TRUE, control = ctl))
  set.seed(11)
  d3 <- data.frame(g = factor(rep(1:40, each = 5)), x = rnorm(200))
  d3$y <- rbinom(200, 1, plogis(-0.3 + 0.7 * d3$x + rnorm(40)[d3$g]))
  r3 <- count_quad(frm(y ~ x + (1 | g), family = bernoulli(), data = d3,
                       quadrature = TRUE))
  expect_identical(c(r1$n, r2$n, r3$n), character())
})

test_that("a quadrature sd at zero is a boundary fit, as under Laplace", {
  # it had a finite Wald standard error, 0.1355 on an sd of 4.8e-5, and
  # an interval reaching -0.27 (dev/quadseam-boundary.R, seed 3 of
  # dev/fuzzfix-quad-boundary.R)
  set.seed(3)
  d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
  d$y <- rbinom(120, 1, plogis(0.2 + 0.5 * d$x))
  for (quad in c(FALSE, TRUE)) {
    f <- allow_boundary(frm(y ~ x + (1 | g), family = bernoulli(),
                            data = d, quadrature = quad), require = TRUE)
    sdr <- sdr_of(f)
    expect_identical(unname(sdr$se_lost), "boundary")
    v <- diag(sdr$cov.fixed)
    expect_true(is.nan(v[[3L]]))
    expect_true(all(is.finite(v[1:2]) & v[1:2] > 0))
    expect_true(is.nan(VarCorr(f)[["g"]]$sd[1L, "Q2.5"]))
  }
  # an sd inside its range keeps its standard error and says nothing
  set.seed(11)
  d2 <- data.frame(g = factor(rep(1:40, each = 5)), x = rnorm(200))
  d2$y <- rbinom(200, 1, plogis(-0.3 + 0.7 * d2$x + rnorm(40)[d2$g]))
  bm <- 0L
  fq <- withCallingHandlers(
    frm(y ~ x + (1 | g), family = bernoulli(), data = d2,
        quadrature = TRUE),
    frmtmb_boundary_fit = function(m) {
      bm <<- bm + 1L
      invokeRestart("muffleMessage")
    })
  expect_identical(bm, 0L)
  expect_true(all(is.finite(diag(sdr_of(fq)$cov.fixed))))
})

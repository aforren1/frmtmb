# predict() draws the parameters from their asymptotic normal law. Where
# an estimate sits at or near the end of its range, that law has a huge
# variance and puts mass on values the data rule out. The draws now read
# each side of such a parameter off the likelihood: a flat side is not
# drawn into, a steep one is narrowed, the rest are drawn given it, and a
# side the normal law describes is left alone. Each case was seen to
# fail on 0.81.0 (dev/ppdraw-findings.md has the output). The assertions
# are about the size of the defect, not exact values.

pde_zib_data <- function(zi, seed = 7) {
  set.seed(seed)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
  d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
  if (zi > 0) d$y[stats::runif(n) < zi] <- 0L
  d
}

# predict()'s draws and the frmtmb_draw_held messages they raised
pde_draw <- function(expr) {
  said <- character(0)
  v <- withCallingHandlers(expr, frmtmb_draw_held = function(m) {
    said <<- c(said, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(p = v, said = said)
}

# The mean of predict()'s draws, in units of the Monte Carlo standard
# error of the draws' mean (each replicate is one parameter draw, so the
# replicates' row means are independent), against two references:
# fitted(), and the expected response at each replicate's own parameters.
# The second takes the same parameter draws (the same seed, and the
# drawer is the first thing predict() takes from the stream), so it
# differs from predict() by the observation noise alone, whatever the
# spread of the draw law; fitted() is the plug-in, which a wide law moves
# away from the draws' mean.
pde_mean_z <- function(fit, ndraws, seed = 1) {
  set.seed(seed)
  r <- pde_draw(predict(fit, ndraws = ndraws, summary = FALSE))
  set.seed(seed)
  dr <- pde_draw(predict_par_drawer(fit, TRUE, ndraws))$p
  em <- vapply(seq_len(ndraws), function(s) {
    mean(fitted_point(dr(s), NULL, NULL, "response", NULL, NULL, FALSE,
                      NULL))
  }, 0)
  m <- rowMeans(r$p)
  fe <- fitted(fit)[, "Estimate"]
  list(fitted = (mean(m) - mean(fe)) /
         (stats::sd(m) / sqrt(ndraws)),
       draws = (mean(m) - mean(em)) / (stats::sd(m - em) / sqrt(ndraws)),
       zero = mean(r$p == 0), said = r$said)
}

test_that("a zero inflation at its boundary is not drawn toward 1", {
  # 0.81.0 drew the zi intercept (-22.2, standard error 4927) from its
  # normal law, so half the replicates had zi near 1: predict() averaged
  # 3.04 where fitted() gave 5.74, 19 Monte Carlo standard errors apart,
  # and half the draws were 0 where brms has 0.0097
  fit <- frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
             data = pde_zib_data(0))
  z <- pde_mean_z(fit, 400L)
  expect_length(z$said, 1L)
  expect_match(z$said, "flat toward one end of the range of `zi_(Intercept)`",
               fixed = TRUE)
  expect_lt(abs(z$fitted), 6)
  expect_lt(abs(z$draws), 4)
  expect_lt(z$zero, 0.05)
  # said once per fit, not once per call
  expect_no_message(predict(fit, ndraws = 10), class = "frmtmb_draw_held")
  # the absent case: a zero inflation inside its range is drawn from its
  # normal law, and nothing is said, at 0.3 and at 0.02, where its
  # likelihood is skewed
  for (zi in c(0.3, 0.02)) {
    f2 <- frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
              data = pde_zib_data(zi))
    z <- pde_mean_z(f2, 400L)
    expect_length(z$said, 0L)
    expect_lt(abs(z$draws), 4)
  }
  # propagate_error = FALSE draws nothing and says nothing
  f3 <- frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
            data = pde_zib_data(0))
  expect_no_message(predict(f3, ndraws = 10, propagate_error = FALSE),
                    class = "frmtmb_draw_held")
})

test_that("a predictor on a zero inflation or hurdle at its end is narrowed", {
  # with zi ~ x or hu ~ x, the slope's standard error is as huge as the
  # intercept's (4710), and both its sides are steep: 0.81.0's draws were
  # 25.9 and 26.2 Monte Carlo standard errors from fitted(), and holding
  # the intercept alone made them 57.7 and 118.1
  d <- pde_zib_data(0)
  f1 <- frm(bf(y | trials(N) ~ x, zi ~ x), family = zero_inflated_binomial(),
            data = d)
  set.seed(12)
  dh <- data.frame(x = stats::rnorm(300))
  dh$y <- stats::rpois(300, exp(2.5 + 0.2 * dh$x))
  f4 <- frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(), data = dh)
  for (f in list(f1, f4)) {
    z <- pde_mean_z(f, 400L)
    expect_lt(abs(z$fitted), 4.5)
    expect_lt(abs(z$draws), 4)
    # the slope is narrowed here, or held where a standard-error check
    # already calls it undefined at that limit (reason "limit"); either
    # way the draws stop short of what the data rule out
    lost <- sdr_of(f)$se_lost
    slope <- grep("_x$", outer_par_names(f), value = TRUE)
    if (!isTRUE(lost[slope] == "limit")) {
      expect_match(z$said, "_x` are narrowed", fixed = TRUE)
    }
  }
})

test_that("a zero inflation just inside its range is not drawn toward 1", {
  # zero_inflated_poisson on poisson data, seed 3: zi at -7.0 with
  # standard error 18.5, its flat side changing by 0.0017, so no side is
  # flat. 0.81.0 averaged 1.98 against fitted() 3.08 and drew 40% zeros;
  # brms 2.23.0 averages 3.063 with 8.0% zeros
  set.seed(3)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.3 * d$x))
  fit <- frm(bf(y ~ x), family = zero_inflated_poisson(), data = d)
  z <- pde_mean_z(fit, 400L)
  expect_match(z$said, "`zi_(Intercept)` are narrowed", fixed = TRUE)
  expect_false(grepl("flat toward", z$said, fixed = TRUE))
  expect_lt(abs(z$fitted), 4.5)
  expect_lt(z$zero, 0.15)
})

test_that("a group sd near zero does not blow up a new level's draw", {
  # gaussian (1 | g), 6 groups of 20, seed 5: the group log sd is -3.04
  # with standard error 6.84, and lme4 does not call the fit singular.
  # 0.81.0 drew a new level's prediction with a standard deviation of
  # 2.2e7, against 0.995 at the estimates
  set.seed(5)
  d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, 0.1)[d$g] + stats::rnorm(120)
  fit <- frm(bf(y ~ x + (1 | g)), data = d)
  nd <- data.frame(x = 0, g = factor("new"))
  set.seed(1)
  r <- pde_draw(predict(fit, newdata = nd, ndraws = 4000, summary = FALSE,
                        allow_new_levels = TRUE))
  set.seed(1)
  p0 <- predict(fit, newdata = nd, ndraws = 4000, summary = FALSE,
                allow_new_levels = TRUE, propagate_error = FALSE)
  expect_lt(stats::sd(r$p) / stats::sd(p0), 2)
  expect_false(grepl("flat toward", paste(r$said, collapse = ""),
                     fixed = TRUE))
})

test_that("a residual sd that se(sigma = TRUE) leaves at 0 is not drawn up", {
  # 0.81.0 drew log sigma (-10.7, standard error 3598) from its normal
  # law: the draws' standard deviation was 1.6e141 on a row whose known
  # standard error is 0.05
  set.seed(7)
  d <- data.frame(x = stats::rnorm(200), s = stats::runif(200, 0.2, 1))
  d$y <- 1 + d$x + stats::rnorm(200, 0, d$s)
  fit <- frm(bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(), data = d)
  nd <- data.frame(x = c(0, 0), s = c(0.05, 5))
  set.seed(7)
  r <- pde_draw(allow_warnings(
    predict(fit, newdata = nd, summary = FALSE, ndraws = 400),
    c("NAs produced", "simulated parameter draw was not finite")))
  expect_match(r$said, "`sigma_(Intercept)`", fixed = TRUE, all = FALSE)
  sim <- as.matrix(simulate(fit, nsim = 400, newdata = nd, seed = 7))
  # the draws add the uncertainty in the coefficients and in sigma's
  # identified side to the noise that simulate() draws at the estimates
  # (brms 2.23.0: 0.108 on the first row against simulate()'s 0.047),
  # and nothing near 1e141
  ratio <- apply(r$p, 2L, stats::sd, na.rm = TRUE) /
    apply(sim, 1L, stats::sd)
  expect_true(all(ratio > 2 / 3 & ratio < 3),
              info = paste(ratio, collapse = " "))
})

test_that("an estimate near the end of its range keeps its flat side", {
  # the absent case of the flat rule: a random slope's log sd at -1.96,
  # which lme4 does not call singular, changes the log-likelihood by
  # 0.024 two standard errors toward zero. That side keeps the normal
  # law; only a steep side may be narrowed
  set.seed(3)
  G <- 15
  m <- 8
  d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                  h = factor(rep(1:5, times = G * m / 5)),
                  x = stats::rnorm(G * m))
  ug <- stats::rnorm(G, 0, 0.8)
  sg <- stats::rnorm(G, 0, 0.4)
  uh <- stats::rnorm(5, 0, 0.6)
  d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + uh[d$h] +
    stats::rnorm(G * m)
  fit <- frm(bf(y ~ x + (1 + x | g) + (1 | h)), data = d)
  set.seed(4)
  r <- pde_draw(predict(fit, newdata = d[1:3, ], ndraws = 50))
  expect_false(grepl("flat toward", paste(r$said, collapse = ""),
                     fixed = TRUE))
  law <- fit$cache[["draw_edge_law"]]
  if (!isFALSE(law)) {
    j <- match("theta_2", outer_par_names(fit))
    expect_identical(law$lo[[j]], 1)
    expect_length(law$edge, 0L)
  }
})

test_that("a curvature that rounds below zero raises no warning", {
  # (1 + x | g) with the correlation held at its limit: its row of the
  # outer Hessian reads +-1e-9, and sqrt() of the negative ones leaked an
  # unclassed "NaNs produced" from predict() on 3 of the first 3 slope
  # scales here (6 of 16 over 0.5 to 2). testthat counts any escaped
  # warning, so the loop asserts none escapes, whatever the sign lands on
  set.seed(5)
  dd <- data.frame(g = factor(rep(1:20, each = 10)), x = stats::rnorm(200))
  a <- stats::rnorm(20, 0, 0.7)[dd$g]
  b <- stats::rnorm(20, 0, 0.05)[dd$g] * dd$x
  e <- stats::rnorm(200)
  nd <- data.frame(x = 2, g = factor("new"))
  for (lam in c(0.5, 0.6, 0.7, 1)) {
    dd$y <- 1 + 0.5 * dd$x + a + lam * b + e
    fit <- frm(bf(y ~ x + (1 + x | g)), data = dd)
    set.seed(1)
    r <- pde_draw(predict(fit, newdata = nd, ndraws = 100,
                          allow_new_levels = TRUE))
    expect_true(all(is.finite(r$p)))
  }
})

test_that("a fit the normal law describes keeps its draws bit for bit", {
  # sleepstudy-like: every side falls as the normal law says, so the
  # side law is not built and the draws are the plain normal ones
  set.seed(11)
  d <- data.frame(g = factor(rep(1:18, each = 10)), x = rep(0:9, 18))
  d$y <- 250 + 10 * d$x + stats::rnorm(18, 0, 25)[d$g] +
    stats::rnorm(18, 0, 6)[d$g] * d$x + stats::rnorm(180, 0, 25)
  fit <- frm(bf(y ~ x + (x | g)), data = d)
  set.seed(2)
  r <- pde_draw(predict_par_drawer(fit, TRUE, 50L))
  expect_length(r$said, 0L)
  expect_false(fit$cache[["draw_edge_law"]])
  # the same normal draws by hand
  ds <- fit_draw_space(fit)
  set.seed(2)
  Z <- matrix(stats::rnorm(nrow(ds$V) * 50L), nrow(ds$V), 50L)
  D <- fit_outer_vector(fit, ds$map) + crossprod(chol(ds$V), Z)
  expect_identical(attr(r$p, "delta")(7), D[, 7] -
                     fit_outer_vector(fit, ds$map))
})

test_that("a 'limit' parameter is drawn only when it is the family's limit", {
  # The standard-error check's boundary rule for family parameters calls
  # two kinds of parameter "limit", held at variance 0 with sdreport()'s
  # own row kept as cov_fixed_raw: a family parameter at the end of its
  # range (the family's limiting case: se_limit_quiet() and, for the
  # parameter at the end, se_limit_draw_held()), and a mean or covariance
  # parameter that a fit with no maximum leaves undetermined (raw
  # variance 2.4e137 on an all-zero poisson). The first is drawn by the
  # side law as without the check; the second is never drawn from that
  # row. The check's verdict is put on the 3269 fit by hand, and its two
  # predicates are mocked, so this runs whether or not the build has the
  # rule
  real <- fit_draw_space
  with_limit <- function(fit) {
    ds <- real(fit)
    j <- match("zi_(Intercept)", ds$map$names)
    # a build whose check already holds it hands the raw row over itself
    ds$raw <- ds$raw %||% ds$V
    ds$V[j, ] <- 0
    ds$V[, j] <- 0
    ds$lost <- c(`zi_(Intercept)` = "limit")
    ds$jc <- list(null = diag(1, nrow(ds$V))[, j, drop = FALSE])
    ds$boundary_only <- TRUE
    ds
  }
  local_mocked_bindings(fit_draw_space = with_limit)
  # kind (a): the family's limit, restored and read side by side
  local_mocked_bindings(
    draw_limit_quiet = function(object, lost) TRUE,
    draw_limit_held = function(object, lost) {
      stats::setNames(lost == "limit", names(lost))
    })
  fa <- frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
            data = pde_zib_data(0))
  set.seed(1)
  ra <- pde_draw(predict_par_drawer(fa, TRUE, 400L))
  expect_match(ra$said, "flat toward one end of the range of `zi_(Intercept)`",
               fixed = TRUE)
  j <- match("zi_(Intercept)", outer_par_names(fa))
  dz <- vapply(seq_len(400L), function(s) attr(ra$p, "delta")(s)[j], 0)
  # the flat side is not drawn into, the identified one is
  expect_true(all(dz >= 0))
  expect_gt(sum(dz > 0), 100)
  # a family parameter the limit leaves undefined (quiet, not at its
  # end): held at its estimate, as the check holds it, and the rest drawn
  local_mocked_bindings(
    draw_limit_quiet = function(object, lost) TRUE,
    draw_limit_held = function(object, lost) {
      stats::setNames(logical(length(lost)), names(lost))
    })
  fu <- frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
            data = pde_zib_data(0))
  set.seed(1)
  ru <- pde_draw(predict_par_drawer(fu, TRUE, 400L))
  du <- vapply(seq_len(400L), function(s) attr(ru$p, "delta")(s), numeric(3))
  expect_true(all(du[j, ] == 0))
  expect_true(all(du[-j, ] != 0))
  # kind (b): not the family's limit, so no law: the existing fallback
  local_mocked_bindings(
    draw_limit_quiet = function(object, lost) FALSE,
    draw_limit_held = function(object, lost) {
      stats::setNames(logical(length(lost)), names(lost))
    })
  fb <- frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
            data = pde_zib_data(0))
  set.seed(1)
  rb <- allow_warnings(pde_draw(predict_par_drawer(fb, TRUE, 400L)),
                       "could not be recovered from the Hessian",
                       require = "could not be recovered from the Hessian")
  expect_length(rb$said, 0L)
  expect_null(attr(rb$p, "delta"))
})

test_that("an ordinal mixture's collapsed threshold gap draws no NA", {
  skip_on_cran()
  # 0.81.0 drew the second component's log gap (-13.1, standard error
  # 3840) from its normal law, so every threshold after it was infinite
  # in half the replicates: 42.75% of the draws were NA, with no word
  set.seed(7)
  d <- data.frame(x = stats::rnorm(400))
  k <- stats::rbinom(400, 1, 0.3)
  d$yo <- cut(d$x + ifelse(k == 1, 2, 0) + stats::rlogis(400),
              c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
  # the jittered starts of an ordinal mixture can find a step function
  # above this fit and warn; the fit they report is the default start's
  fit <- allow_warnings(
    frm(bf(yo ~ x), family = mixture(cumulative(), cumulative()), data = d),
    "jittered starts reached a higher log-likelihood")
  set.seed(13)
  r <- pde_draw(predict(fit, newdata = data.frame(x = c(-1, 0, 1)),
                        summary = FALSE, ndraws = 4000L))
  expect_length(r$said, 1L)
  expect_match(r$said, "log(mu2_Intercept[2] - mu2_Intercept[1])",
               fixed = TRUE)
  expect_false(anyNA(r$p))
  expect_true(all(r$p %in% 1:4))
  # the shares are NOT compared with fitted(): this mixture is weakly
  # identified, and the draw law's own spread moves them (the findings'
  # sweep, filed in the backlog)
})

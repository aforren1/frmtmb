# v0.35: conditional_effects() - likelihood-profile and
# parametric-bootstrap confidence bands, the default effect search
# (nonlinear, mo() and mi() fits), and the brms-compatibility refusals.
# Every test here either refits the model or runs a root search, so the
# whole file is CRAN-skipped.
skip_on_cran()

sim_ce_gauss <- function(n = 200, seed = 21) {
  set.seed(seed)
  dd <- data.frame(x = rnorm(n), f = factor(rep(c("a", "b"), n / 2)))
  dd$y <- rnorm(n, 1 + 0.8 * dd$x + (dd$f == "b") * 0.5, 1)
  dd
}

test_that("boot bands are percentiles of one shared bootstrap", {
  dd <- sim_ce_gauss(n = 120, seed = 31)
  fit <- frm(bf(y ~ x + f), family = gaussian(), data = dd)

  res <- 5
  ce <- conditional_effects(fit, effects = "x", resolution = res,
                            band = "boot", boot = 25, seed = 11)
  dx <- ce$x
  expect_equal(nrow(dx), res)
  # the estimate is still the fit's, not the bootstrap mean
  cw <- conditional_effects(fit, effects = "x", resolution = res)$x
  expect_equal(dx$estimate__, cw$estimate__)

  # exact identity with a bootstrap run by hand over the same grid and
  # the same seed: same draws, same refits, same percentiles
  grid <- data.frame(x = seq(min(dd$x), max(dd$x), length.out = res),
                     f = factor("a", levels = levels(dd$f)),
                     y = mean(dd$y))
  bsm <- frm_bootstrap(
    fit,
    FUN = function(f) {
      as.vector(predict(f, newdata = grid, type = "response", dpar = "mu",
                        resp = "y", re.form = NA))
    },
    nsim = 25, seed = 11
  )
  expect_true(all(bsm$converged))
  expect_equal(dx$lower__,
               unname(apply(bsm$t, 2, stats::quantile, 0.025)))
  expect_equal(dx$upper__,
               unname(apply(bsm$t, 2, stats::quantile, 0.975)))
  # se__ is the draw sd under band = "boot", not the Wald se
  expect_equal(dx$se__, unname(apply(bsm$t, 2, stats::sd)))
  expect_false(isTRUE(all.equal(dx$se__, cw$se__)))

  # ONE bootstrap covers every effect of the call: 2 grids, one t matrix
  ce2 <- conditional_effects(fit, effects = c("x", "f"), resolution = res,
                             band = "boot", boot = 25, seed = 11)
  bs <- attr(ce2, "boot")
  expect_s3_class(bs, "frmtmb_boot")
  expect_equal(bs$nsim, 25L)
  expect_equal(ncol(bs$t), res + nlevels(dd$f))
  expect_equal(ce2$x$lower__, dx$lower__)

  # and it can be handed back, which refits nothing
  ce3 <- conditional_effects(fit, effects = c("x", "f"), resolution = res,
                             band = "boot", boot = bs)
  expect_equal(ce3$x, ce2$x)
  expect_equal(ce3$f, ce2$f)
})

test_that("boot bands on a mixed model match a hand-run bootstrap", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  fs <- frm(bf(Reaction ~ Days + (1 | Subject)), family = gaussian(),
            data = sleepstudy)

  res <- 4
  ce <- conditional_effects(fs, effects = "Days", resolution = res,
                            band = "boot", boot = 20, seed = 3)$Days
  grid <- data.frame(Days = seq(min(sleepstudy$Days),
                                max(sleepstudy$Days), length.out = res),
                     Subject = sleepstudy$Subject[1],
                     Reaction = mean(sleepstudy$Reaction))
  bsm <- frm_bootstrap(
    fs,
    FUN = function(f) {
      as.vector(predict(f, newdata = grid, type = "response", dpar = "mu",
                        resp = "Reaction", re.form = NA))
    },
    nsim = 20, seed = 3
  )
  ok <- apply(bsm$t, 2, function(v) all(is.finite(v)))
  expect_true(all(ok))
  expect_equal(ce$lower__, unname(apply(bsm$t, 2, stats::quantile, 0.025)))
  expect_equal(ce$upper__, unname(apply(bsm$t, 2, stats::quantile, 0.975)))
  # the band is a population-level band: the group effects are excluded
  # from the predictions but present in every refit
  expect_true(all(ce$lower__ < ce$estimate__ & ce$estimate__ < ce$upper__))
})

test_that("boot bands announce the default run and refuse a foreign boot", {
  dd <- sim_ce_gauss(n = 60, seed = 32)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  # an explicit draw count is the user's decision, so it stays quiet
  expect_no_message(
    conditional_effects(fit, effects = "x", resolution = 3, band = "boot",
                        boot = 5, seed = 1)
  )
  # the default runs 200 refits and says so before it does
  expect_message(
    conditional_effects(fit, effects = "x", resolution = 3, band = "boot",
                        seed = 1),
    "refitting the model 200 times"
  )
  # a bootstrap of the coefficients is not a bootstrap of this grid
  bad <- frm_bootstrap(fit, nsim = 5, seed = 1)
  expect_error(
    conditional_effects(fit, effects = "x", resolution = 3, band = "boot",
                        boot = bad),
    "was not produced by a conditional_effects"
  )
  expect_error(
    conditional_effects(fit, effects = "x", band = "boot", boot = "lots"),
    "single number of bootstrap draws"
  )
})

test_that("a bootstrap is reused only for the grid it was run over", {
  dd <- sim_ce_gauss(n = 80, seed = 36)
  # a large level effect, so a band taken under the wrong condition is
  # visibly a band for some other curve
  dd$y <- dd$y + 4 * (dd$f == "b")
  fit <- frm(bf(y ~ x + f), family = gaussian(), data = dd)
  ce <- function(...) {
    conditional_effects(fit, effects = "x", resolution = 4,
                        band = "boot", seed = 2, ...)
  }

  ca <- ce(boot = 10, conditions = list(f = "a"))
  bs <- attr(ca, "boot")
  # the same call again reuses the draws and lands on the same numbers
  expect_equal(ce(boot = bs, conditions = list(f = "a"))$x, ca$x)
  # a different coverage is still the same draws: prob is read off them
  expect_no_error(ce(boot = bs, conditions = list(f = "a"), prob = 0.5))

  # the grid VALUES decide, not the grid's shape: the f = "b" call has
  # the same effect, row count and columns, and a different curve
  expect_error(ce(boot = bs, conditions = list(f = "b")),
               "predictions over a different grid")
  expect_error(
    conditional_effects(fit, effects = "x", resolution = 5,
                        band = "boot", seed = 2, boot = bs,
                        conditions = list(f = "a")),
    "predictions over a different grid"
  )
  expect_error(
    conditional_effects(fit, effects = "f", band = "boot", seed = 2,
                        boot = bs, conditions = list(f = "a")),
    "predictions over a different grid"
  )
  expect_error(
    conditional_effects(fit, effects = "x", resolution = 4,
                        band = "boot", seed = 2, boot = bs,
                        dpar = "sigma", conditions = list(f = "a")),
    "predictions over a different grid"
  )

  # what the refusal protects: run properly, the f = "b" band brackets
  # its own estimate, which the wrongly reused draws would not have
  cb <- ce(boot = 10, conditions = list(f = "b"))
  expect_gt(min(cb$x$estimate__ - ca$x$estimate__), 3)
  expect_true(all(cb$x$lower__ <= cb$x$estimate__ &
                    cb$x$estimate__ <= cb$x$upper__))
})

test_that("profile and Wald bands agree on a large-n gaussian", {
  set.seed(4)
  n <- 600
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + 0.8 * dd$x, 1)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)

  res <- 7
  cw <- conditional_effects(fit, effects = "x", resolution = res)$x
  cp <- conditional_effects(fit, effects = "x", resolution = res,
                            band = "profile", profile_points = res)$x
  expect_equal(cp$estimate__, cw$estimate__)
  # se__ stays the Wald standard error: the profile moves the endpoints
  expect_equal(cp$se__, cw$se__)
  width <- mean(cw$upper__ - cw$lower__)
  expect_lt(max(abs(cp$lower__ - cw$lower__)), 0.02 * width)
  expect_lt(max(abs(cp$upper__ - cw$upper__)), 0.02 * width)
  expect_true(all(cp$lower__ <= cp$estimate__ &
                    cp$estimate__ <= cp$upper__))

  # the same three bands nearly coincide, bootstrap included
  cb <- conditional_effects(fit, effects = "x", resolution = res,
                            band = "boot", boot = 60, seed = 7)$x
  expect_lt(max(abs(cb$lower__ - cw$lower__)), 0.5 * width)
  expect_lt(max(abs(cb$upper__ - cw$upper__)), 0.5 * width)
})

test_that("a profile band inverts the same likelihood as confint()", {
  # at the reference level of a treatment-coded factor the grid row IS
  # the intercept, so the band must equal the intercept's own interval
  set.seed(9)
  db <- data.frame(f = factor(rep(c("a", "b"), each = 20)))
  db$y <- rbinom(40, 1, c(0.7, 0.35)[as.integer(db$f)])
  fb <- frm(bf(y ~ f), family = bernoulli(), data = db)

  cp <- conditional_effects(fb, effects = "f", band = "profile")$f
  cw <- conditional_effects(fb, effects = "f")$f
  ci <- confint(fb, parm = "(Intercept)", method = "uniroot")
  expect_equal(stats::qlogis(cp$lower__[1]), unname(ci[1, "lwr"]))
  expect_equal(stats::qlogis(cp$upper__[1]), unname(ci[1, "upr"]))

  # Wald is symmetric on the link scale by construction; the profile is
  # not, and it leans the way the likelihood does
  wl <- stats::qlogis(cw$estimate__[1]) - stats::qlogis(cw$lower__[1])
  wu <- stats::qlogis(cw$upper__[1]) - stats::qlogis(cw$estimate__[1])
  expect_equal(wl, wu)
  pl <- stats::qlogis(cp$estimate__[1]) - stats::qlogis(cp$lower__[1])
  pu <- stats::qlogis(cp$upper__[1]) - stats::qlogis(cp$estimate__[1])
  expect_gt(pu, pl)
  expect_gt(pu, wu)
  expect_lt(pl, wl)
})

test_that("a profile band caps and interpolates a numeric grid", {
  dd <- sim_ce_gauss(n = 80, seed = 33)
  fit <- frm(bf(y ~ x + f), family = gaussian(), data = dd)
  res <- 9
  cp <- conditional_effects(fit, effects = "x", resolution = res,
                            band = "profile", profile_points = 3)$x
  expect_equal(nrow(cp), res)
  expect_true(all(is.finite(cp$lower__)))
  # profiled at 3 points (1, 5, 9), interpolated linearly between them
  # on the link scale, so an interior point sits on the chord
  expect_equal(cp$lower__[3], (cp$lower__[1] + cp$lower__[5]) / 2)
  expect_equal(cp$upper__[7], (cp$upper__[5] + cp$upper__[9]) / 2)
  # the cap does not touch a discrete grid: every level is profiled
  cf <- conditional_effects(fit, effects = "f", band = "profile",
                            profile_points = 1)$f
  expect_equal(nrow(cf), 2L)
  expect_true(all(is.finite(cf$lower__)))
})

test_that("a profile band covers a distributional parameter", {
  # the sigma predictor lives in the betad block, whose coefficients are
  # ranked past the fixed entries before they reach the outer vector
  set.seed(12)
  n <- 150
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + dd$x, exp(-0.5 + 0.3 * dd$x))
  fit <- frm(bf(y ~ x, sigma ~ x), family = gaussian(), data = dd)

  cw <- conditional_effects(fit, effects = "x", dpar = "sigma",
                            resolution = 3)$x
  cp <- conditional_effects(fit, effects = "x", dpar = "sigma",
                            resolution = 3, band = "profile",
                            profile_points = 3)$x
  expect_equal(cp$estimate__, cw$estimate__)
  expect_true(all(cp$lower__ < cp$estimate__ &
                    cp$estimate__ < cp$upper__))
  # a log-sd profile leans up relative to its symmetric Wald band
  expect_true(all(cp$upper__ > cw$upper__))
  expect_true(all(cp$lower__ > cw$lower__))
})

test_that("profile points that do not converge degrade to NA", {
  dd <- sim_ce_gauss(n = 60, seed = 34)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  testthat::local_mocked_bindings(
    tmbroot = function(...) c(NA_real_, NA_real_), .package = "TMB"
  )
  expect_warning(
    cp <- conditional_effects(fit, effects = "x", resolution = 5,
                              band = "profile", profile_points = 5)$x,
    "did not converge at 5 of 5"
  )
  expect_true(all(is.na(cp$lower__)))
  expect_true(all(is.na(cp$upper__)))
  # the estimate and the Wald se survive: only the band is missing
  expect_true(all(is.finite(cp$estimate__)))
  expect_true(all(is.finite(cp$se__)))
})

test_that("boot bands cover the ordinal per-category display", {
  set.seed(6)
  do <- data.frame(x = rnorm(100))
  do$y <- factor(cut(do$x + rnorm(100), 3), ordered = TRUE,
                 labels = c("l", "m", "h"))
  fo <- frm(bf(y ~ x), family = cumulative(), data = do)

  ob <- conditional_effects(fo, effects = "x", resolution = 4,
                            band = "boot", boot = 15, seed = 2)$x
  expect_equal(nrow(ob), 12L)
  expect_equal(levels(ob$cats__), c("l", "m", "h"))
  expect_true(all(is.finite(ob$lower__)))
  expect_true(all(ob$lower__ <= ob$estimate__ &
                    ob$estimate__ <= ob$upper__))
  expect_true(all(ob$lower__ >= 0 & ob$upper__ <= 1))
  # the estimates are still the fitted category probabilities
  ow <- conditional_effects(fo, effects = "x", resolution = 4)$x
  expect_equal(ob$estimate__, ow$estimate__)
})

test_that("bands other than wald refuse what they cannot cover", {
  set.seed(5)
  dz <- data.frame(x = rnorm(120))
  dz$y <- ifelse(runif(120) < 0.3, 0, rpois(120, exp(0.6 + 0.4 * dz$x)))
  fz <- frm(bf(y ~ x), family = zero_inflated_poisson(), data = dz)

  # method = "predict" bands are already simulation quantiles
  expect_error(
    conditional_effects(fz, effects = "x", method = "predict",
                        band = "profile"),
    "does not apply to method"
  )
  expect_error(
    conditional_effects(fz, effects = "x", method = "predict",
                        band = "boot"),
    "does not apply to method"
  )
  # a mean that runs through several dpars is not one lincomb
  expect_error(
    conditional_effects(fz, effects = "x", band = "profile"),
    "more than one distributional parameter"
  )
  # naming one predictor opts back in
  expect_s3_class(
    conditional_effects(fz, effects = "x", band = "profile", dpar = "mu",
                        resolution = 3, profile_points = 3),
    "frmtmb_conditional_effects"
  )

  set.seed(6)
  do <- data.frame(x = rnorm(80))
  do$y <- factor(cut(do$x + rnorm(80), 3), ordered = TRUE,
                 labels = c("l", "m", "h"))
  fo <- frm(bf(y ~ x), family = cumulative(), data = do)
  expect_error(
    conditional_effects(fo, effects = "x", band = "profile"),
    "ordinal category probability"
  )

  set.seed(7)
  dr <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dr$y <- rnorm(60, 1 + dr$x + rnorm(6, 0, 0.5)[dr$g], 1)
  fr <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dr,
            REML = TRUE)
  expect_error(
    conditional_effects(fr, effects = "x", band = "profile"),
    "requires an ML fit"
  )

  set.seed(8)
  ds <- data.frame(x = runif(80, -3, 3))
  ds$y <- rnorm(80, sin(ds$x), 0.3)
  fs <- frm(bf(y ~ s(x)), family = gaussian(), data = ds)
  expect_error(
    conditional_effects(fs, effects = "x", band = "profile"),
    "smooth, gp\\(\\) or hsgp\\(\\) term"
  )
  # the bootstrap has no such trouble: it refits
  sb <- conditional_effects(fs, effects = "x", resolution = 4,
                            band = "boot", boot = 8, seed = 3)$x
  expect_true(all(is.finite(sb$lower__)))
})

test_that("every band plots, prints and keeps cond__ working", {
  dd <- sim_ce_gauss(n = 80, seed = 35)
  fit <- frm(bf(y ~ x + f), family = gaussian(), data = dd)
  cnd <- data.frame(f = factor(c("a", "b"), levels = levels(dd$f)),
                    row.names = c("f = a", "f = b"))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (bnd in c("wald", "profile", "boot")) {
    ce <- conditional_effects(fit, effects = c("x", "x:f"),
                              resolution = 5, band = bnd, boot = 8,
                              profile_points = 3, seed = 5)
    expect_named(ce, c("x", "x:f"))
    expect_true(all(ce$x$lower__ <= ce$x$upper__))
    expect_equal(attr(ce$x, "band"), bnd)
    # the two-variable grid keeps its second predictor
    expect_equal(nrow(ce$`x:f`), 10L)
    expect_no_error(plot(ce, ask = FALSE))
    expect_no_error(plot(ce, ask = FALSE, points = TRUE))
    expect_no_error(print(ce))

    cc <- conditional_effects(fit, effects = "x", resolution = 4,
                              band = bnd, boot = 8, profile_points = 3,
                              seed = 5, conditions = cnd)
    # cond__ is a factor of the condition labels, as brms's is
    expect_equal(levels(cc$x$cond__), c("f = a", "f = b"))
    expect_equal(sort(unique(as.character(cc$x$cond__))),
                 c("f = a", "f = b"))
    expect_equal(nrow(cc$x), 8L)
    expect_no_error(plot(cc, ask = FALSE))
  }
})

# --- the default effect search (brms vignette port FN-3, FN-13) -------

test_that("a nonlinear fit finds the covariates of its nl body", {
  set.seed(1)
  d <- data.frame(x = rnorm(60))
  d$y <- 1 + 2 * d$x + rnorm(60)
  fnl <- frm(bf(y ~ a * x + b, a ~ 1, b ~ 1, nl = TRUE), data = d,
             family = gaussian(), start = list(beta = c(1, 1)))

  # the mu predictor has no terms of its own: x is read by the nl body
  ce <- conditional_effects(fnl, resolution = 5, band = "boot", boot = 8,
                            seed = 1)
  expect_named(ce, "x")
  expect_equal(nrow(ce$x), 5L)
  expect_true(all(is.finite(ce$x$estimate__)))
  expect_true(all(ce$x$lower__ <= ce$x$upper__))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(ce, ask = FALSE, points = TRUE))
  # a nonlinear predictor has no delta-method se, so the other bands say
  # so instead of failing somewhere inside predict()
  expect_error(conditional_effects(fnl, resolution = 5),
               "cannot put a wald band on a nonlinear predictor")
  expect_error(
    conditional_effects(fnl, resolution = 5, band = "profile"),
    "not available for a nonlinear predictor"
  )
})

test_that("an mo()-only fit finds its monotonic variable", {
  set.seed(1)
  d <- data.frame(o = factor(sample(1:4, 60, TRUE), ordered = TRUE))
  d$y <- 1 + as.integer(d$o) + rnorm(60)
  fmo <- frm(y ~ mo(o), data = d, family = gaussian())

  # the mo() design column is a placeholder, so the variable never
  # reaches terms(): it is read off the stored term instead
  ce <- conditional_effects(fmo)
  expect_named(ce, "o")
  expect_equal(nrow(ce$o), 4L)
  expect_true(all(is.finite(ce$o$se__)))
  expect_equal(as.character(ce$o$o), as.character(1:4))
  # naming it explicitly agrees with the search
  expect_equal(conditional_effects(fmo, effects = "o")$o, ce$o)
})

test_that("an mi() fit builds a grid from the complete cases", {
  set.seed(1)
  d <- data.frame(x = rnorm(60))
  d$y <- 1 + 2 * d$x + rnorm(60)
  d$x[1:5] <- NA
  fmi <- frm(bf(y | mi() ~ mi(x)) + bf(x | mi() ~ 1) + set_rescor(FALSE),
             data = d, family = gaussian())

  # the imputed predictor has gaps; the grid spans what is observed
  ce <- conditional_effects(fmi, resp = "y", resolution = 5)
  expect_named(ce, "x")
  expect_equal(range(ce$x$x), range(d$x, na.rm = TRUE))
  expect_true(all(is.finite(ce$x$estimate__)))
  expect_true(all(is.finite(ce$x$se__)))
  expect_true(all(ce$x$lower__ <= ce$x$upper__))
  # the raw-observation overlay carries the gaps and still draws
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(ce, ask = FALSE, points = TRUE))
})

# --- brms compatibility (FN-4) ---------------------------------------

test_that("surface = TRUE is refused and names the two-variable effect", {
  set.seed(1)
  d <- data.frame(x = rnorm(60), z = rnorm(60))
  d$y <- 1 + d$x * d$z + rnorm(60)
  fs <- frm(y ~ s(x, z), data = d, family = gaussian())

  expect_error(conditional_effects(fs, surface = TRUE),
               "surface = TRUE\\) is not implemented")
  # the alternative the message names actually works
  ce <- conditional_effects(fs, effects = "x:z", resolution = 5)
  expect_equal(nrow(ce$`x:z`), 15L)
  # and the brms default spelling is accepted, not refused
  expect_no_error(conditional_effects(fs, effects = "x", resolution = 4,
                                      surface = FALSE))
})

# --- per-category display contract (nominal families at merge) -------

test_that("the per-category display keys on the family contract", {
  stub <- function(type) list(family = list(type = type, family = type))
  # ordinal today, nominal categorical() at merge: both predict an
  # n x K probability matrix, and both take the cats__ display
  expect_true(frmtmb:::ce_cats_display(stub("ordinal"), NULL))
  expect_true(frmtmb:::ce_cats_display(stub("categorical"), NULL))
  expect_false(frmtmb:::ce_cats_display(stub("continuous"), NULL))
  expect_false(frmtmb:::ce_cats_display(stub("discrete"), NULL))
  # naming a dpar always opts back into the linear-predictor display
  expect_false(frmtmb:::ce_cats_display(stub("ordinal"), "mu"))
  expect_false(frmtmb:::ce_cats_display(stub("categorical"), "mu2"))
  expect_false(frmtmb:::ce_cats_display(list(family = list()), NULL))
})

test_that("a nominal family gets bootstrap bands and honest refusals", {
  skip_if_not(exists("categorical", envir = asNamespace("frmtmb")),
              "categorical() arrives with the families lane")
  set.seed(1)
  d <- data.frame(x = rnorm(120))
  # draw from a real multinomial logit: a deterministic threshold rule
  # separates the classes perfectly and the coefficients diverge
  eta2 <- -0.3 + 1.1 * d$x
  eta3 <- 0.4 - 0.8 * d$x
  pr <- cbind(1, exp(eta2), exp(eta3))
  pr <- pr / rowSums(pr)
  d$y <- factor(c("a", "b", "c")[
    apply(pr, 1, function(p) sample.int(3, 1, prob = p))])
  fc <- frm(y ~ x, data = d, family = get("categorical",
                                          envir = asNamespace("frmtmb"))())

  cb <- conditional_effects(fc, effects = "x", resolution = 4,
                            band = "boot", boot = 10, seed = 1)$x
  expect_equal(nrow(cb), 12L)
  expect_equal(levels(cb$cats__), levels(d$y))
  expect_true(all(cb$estimate__ >= 0 & cb$estimate__ <= 1))
  expect_true(all(cb$lower__ <= cb$estimate__ &
                    cb$estimate__ <= cb$upper__))
  # the delta method for a category probability is written for ordinal
  # thresholds, which a nominal family has none of
  expect_error(conditional_effects(fc, effects = "x"),
               "no analytic standard error")
  expect_error(conditional_effects(fc, effects = "x", band = "profile"),
               "no analytic standard error")
})

test_that("dharma_residuals() refuses a nominal response", {
  skip_if_not_installed("DHARMa")
  skip_if_not(exists("categorical", envir = asNamespace("frmtmb")),
              "categorical() arrives with the families lane")
  set.seed(1)
  d <- data.frame(x = rnorm(120))
  # draw from a real multinomial logit: a deterministic threshold rule
  # separates the classes perfectly and the coefficients diverge
  eta2 <- -0.3 + 1.1 * d$x
  eta3 <- 0.4 - 0.8 * d$x
  pr <- cbind(1, exp(eta2), exp(eta3))
  pr <- pr / rowSums(pr)
  d$y <- factor(c("a", "b", "c")[
    apply(pr, 1, function(p) sample.int(3, 1, prob = p))])
  fc <- frm(y ~ x, data = d, family = get("categorical",
                                          envir = asNamespace("frmtmb"))())
  # a quantile residual is a CDF evaluation, and a nominal category has
  # no CDF: the 1..K codes are an arbitrary labeling
  expect_error(dharma_residuals(fc, nsim = 5),
               "no meaning for a nominal response")
})

test_that("band = 'boot' carries a new group's variance, ordinal included", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  # re_formula = NULL conditions on a NEW group: its modes are zero, so
  # the CURVE is the population curve and the group variance belongs in
  # the BAND. The wald band gets it from lp_extra_var(); the bootstrap
  # has no variance to add, so each replicate DRAWS that group's
  # effects. Without the draw every refit predicted the new level's
  # zero modes and the interval was the population interval under
  # another name, while the frame's NA grouping column claimed
  # otherwise.
  ss <- lme4::sleepstudy
  fs <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(), data = ss)
  wid <- function(d) d$upper__ - d$lower__

  wp <- conditional_effects(fs, effects = "Days", resolution = 4)$Days
  wn <- conditional_effects(fs, effects = "Days", resolution = 4,
                            re_formula = NULL)$Days
  bp <- conditional_effects(fs, effects = "Days", resolution = 4,
                            band = "boot", boot = 100, seed = 11)$Days
  bn <- conditional_effects(fs, effects = "Days", resolution = 4,
                            re_formula = NULL, band = "boot", boot = 100,
                            seed = 11)$Days

  # the estimate is the population curve on every one of the four
  expect_equal(bn$estimate__, bp$estimate__, tolerance = 1e-12)
  expect_true(all(is.na(bn$Subject)))
  # and the band is NOT the population band any more: it was
  # bit-identical to it, ratio 1.00
  expect_true(all(wid(bn) > 2 * wid(bp)))
  # and it lands on the WALD band for the same quantity. The bounds
  # are wide because a percentile band at 100 refits is noisy: measured
  # 0.91-1.03 at this seed, 0.70-1.03 across seeds at 60 refits. What
  # is not noisy is the ratio to the population band above.
  expect_true(all(wid(bn) / wid(wn) > 0.6))
  expect_true(all(wid(bn) / wid(wn) < 1.5))
  # the population bootstrap band is untouched by the change
  expect_true(all(wid(bp) / wid(wp) > 0.8))
  expect_true(all(wid(bp) / wid(wp) < 1.2))

  # THE ORDINAL CASE, which is why this is not cosmetic: the
  # per-category delta method is refused with re_formula, by name, and
  # sends the user to band = "boot" - so boot is the ONLY band
  # available there and it was the one that dropped the variance.
  set.seed(8)
  n <- 300
  g <- factor(rep(seq_len(20), length.out = n))
  u <- stats::rnorm(20, 0, 1.2)
  x <- stats::rnorm(n)
  do <- data.frame(x = x, g = g)
  do$y <- ordered(cut(0.9 * x + u[g] + stats::rlogis(n),
                      c(-Inf, -1, 0.6, Inf), labels = 1:3))
  fo <- frm(bf(y ~ x + (1 | g)) + cumulative(), data = do)

  expect_error(conditional_effects(fo, re_formula = NULL),
               "written for the population-level curve")
  op <- conditional_effects(fo, effects = "x", resolution = 4,
                            band = "boot", boot = 50,
                            seed = 3)[["x:cats__"]]
  on <- conditional_effects(fo, effects = "x", resolution = 4,
                            re_formula = NULL, band = "boot", boot = 50,
                            seed = 3)[["x:cats__"]]
  expect_equal(on$estimate__, op$estimate__, tolerance = 1e-12)
  expect_true(all(is.na(on$g)))
  expect_true(all(wid(on) > wid(op)))
  expect_gt(median(wid(on) / wid(op)), 1.5)
})

test_that("effect2__ keeps distinct moderator values distinct", {
  skip_on_cran()

  # plot() groups on effect2__, so a level per CURVE is the contract.
  # Building the levels from round(v, 2) gave two curves one level
  # whenever two values rounded together, and dropped the user's own
  # names with them (the name-count guard stopped matching). Newly
  # reachable because int_conditions is newly implemented.
  set.seed(1)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- 1 + 0.8 * d$x - 0.4 * d$z + stats::rnorm(n)
  f <- frm(bf(y ~ x + z) + gaussian(), data = d)

  ce <- conditional_effects(
    f, effects = "x:z", resolution = 4,
    int_conditions = list(z = c(lo = 0.001, hi = 0.002)))[["x:z"]]
  expect_identical(sort(unique(ce$z)), c(0.001, 0.002))
  expect_identical(nlevels(ce$effect2__), 2L)
  expect_identical(levels(ce$effect2__), c("hi", "lo"))
  expect_identical(nrow(ce), 8L)

  # unnamed, the labels widen only as far as they must to stay distinct
  cu <- conditional_effects(
    f, effects = "x:z", resolution = 4,
    int_conditions = list(z = c(0.001, 0.002)))[["x:z"]]
  expect_identical(nlevels(cu$effect2__), 2L)
  expect_identical(levels(cu$effect2__), c("0.002", "0.001"))

  # and the ordinary case still rounds to two decimals, descending,
  # which is what brms puts there
  cd <- conditional_effects(f, effects = "x:z", resolution = 4)[["x:z"]]
  zz <- sort(unique(cd$z), decreasing = TRUE)
  expect_identical(levels(cd$effect2__), as.character(round(zz, 2)))
  expect_identical(nlevels(cd$effect2__), 3L)
})

test_that("a bootstrap cannot be reused across the new-group boundary", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  # The grid alone cannot tell a population call from a
  # re_formula = NULL one: ce_ref_value() holds an unvaried factor at
  # levels(col)[1] and the new-level placeholder is bk$levels[1], so on
  # sleepstudy both are "308" and the two grids are BYTE IDENTICAL.
  # Before the key carried the new-level spec, `boot =` reuse accepted a
  # population bootstrap for a new-group call and handed back the
  # population band - the very defect the per-replicate draw exists to
  # fix, through the reuse path the help page recommends - and the
  # converse gave a population call a band four times too wide.
  ss <- lme4::sleepstudy
  fs <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(), data = ss)
  wid <- function(d) d$upper__ - d$lower__
  base <- list(fs, effects = "Days", resolution = 4, band = "boot")

  cp <- do.call(conditional_effects, c(base, list(boot = 40, seed = 5)))
  cn <- do.call(conditional_effects,
                c(base, list(re_formula = NULL, boot = 40, seed = 5)))
  op <- attr(cp, "boot")
  on <- attr(cn, "boot")
  # the two bands are worth telling apart: about 4x
  expect_true(all(wid(cn$Days) > 2 * wid(cp$Days)))

  # BOTH directions are refused, by name and with the reason
  expect_error(
    do.call(conditional_effects,
            c(base, list(re_formula = NULL, boot = op))),
    "different group")
  expect_error(
    do.call(conditional_effects,
            c(base, list(re_formula = NULL, boot = op))),
    "conditions on a NEW group")
  expect_error(do.call(conditional_effects, c(base, list(boot = on))),
               "different group")
  expect_error(do.call(conditional_effects, c(base, list(boot = on))),
               "too wide")

  # and the MATCHING reuse still works, reproducing its own band exactly
  m1 <- do.call(conditional_effects, c(base, list(boot = op)))$Days
  m2 <- do.call(conditional_effects,
                c(base, list(re_formula = NULL, boot = on)))$Days
  expect_equal(wid(m1), wid(cp$Days), tolerance = 1e-12)
  expect_equal(wid(m2), wid(cn$Days), tolerance = 1e-12)
  # the grid-mismatch reason is still reachable, and is a different one
  expect_error(
    do.call(conditional_effects,
            c(list(fs, effects = "Days", resolution = 6, band = "boot"),
              list(boot = op))),
    "different grid")
})

# The draws method returns core's frame. Every assertion here is a
# comparison against the FIT method on the same model rather than
# against a written-out column list: what is being pinned is that one
# engine serves both surfaces, and a list would go stale the next time
# core adds a column.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

band_cols <- c("estimate__", "se__", "lower__", "upper__")

grid_of <- function(df) df[setdiff(names(df), band_cols)]

ce_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(9)
      dd <- data.frame(x = stats::rnorm(60), z = stats::rnorm(60),
                       g = factor(rep(1:6, 10)))
      dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x + 0.3 * dd$z +
                             stats::rnorm(6, 0, 0.5)[dd$g], 1)
      fit <- frm(bf(y ~ x + z + (1 | g)), family = gaussian(), data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 500, refresh = 0, seed = 1)))
      cache <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache
  }
})

ce_ord_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      dd <- v29_ordinal_data(46, n = 140)
      fit <- frm(bf(y ~ x) + cumulative(), data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 500, refresh = 0, seed = 2)))
      cache <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache
  }
})

test_that("the draws frame has the fit frame's columns and grid", {
  cs <- ce_case()
  cd <- conditional_effects(cs$ds, effects = "x", resolution = 8)
  cf <- conditional_effects(cs$fit, effects = "x", resolution = 8)
  expect_identical(names(cd), names(cf))
  expect_identical(names(cd$x), names(cf$x))
  expect_equal(grid_of(cd$x), grid_of(cf$x), ignore_attr = TRUE)
  expect_identical(attr(cd$x, "effects"), attr(cf$x, "effects"))
  # the held values and brms's plot() columns, which the old frame
  # dropped: without them a ported faceting call has nothing to facet
  expect_true(all(c("y", "z", "g", "cond__", "effect1__") %in%
                    names(cd$x)))
  # cond__ is present with ONE level when there is one condition set,
  # which is brms's rule and was not this method's
  expect_s3_class(cd$x$cond__, "factor")
  expect_length(levels(cd$x$cond__), 1L)
  expect_identical(attr(cd$x, "band"), "posterior")
})

test_that("the default effect list and its grids match the fit method", {
  cs <- ce_case()
  cd <- conditional_effects(cs$ds, resolution = 6)
  cf <- conditional_effects(cs$fit, resolution = 6)
  expect_identical(names(cd), names(cf))
  for (k in names(cf)) {
    expect_identical(names(cd[[k]]), names(cf[[k]]))
    expect_equal(grid_of(cd[[k]]), grid_of(cf[[k]]), ignore_attr = TRUE)
  }
})

test_that("conditions and int_conditions reach the draws grid", {
  cs <- ce_case()
  cnd <- data.frame(z = c(-1, 1))
  cd <- conditional_effects(cs$ds, effects = "x", resolution = 5,
                            conditions = cnd)
  cf <- conditional_effects(cs$fit, effects = "x", resolution = 5,
                            conditions = cnd)
  expect_identical(names(cd$x), names(cf$x))
  expect_equal(grid_of(cd$x), grid_of(cf$x), ignore_attr = TRUE)
  expect_length(levels(cd$x$cond__), 2L)

  # int_conditions was not a formal here at all, so a ported call
  # naming it hit "unused argument"
  ic <- list(z = c(lo = -1, hi = 1))
  cd2 <- conditional_effects(cs$ds, effects = "x:z", resolution = 5,
                             int_conditions = ic)
  cf2 <- conditional_effects(cs$fit, effects = "x:z", resolution = 5,
                             int_conditions = ic)
  expect_identical(names(cd2[["x:z"]]), names(cf2[["x:z"]]))
  expect_equal(grid_of(cd2[["x:z"]]), grid_of(cf2[["x:z"]]),
               ignore_attr = TRUE)
  # the user's own labels survive onto the moderator display
  expect_identical(levels(cd2[["x:z"]]$effect2__),
                   levels(cf2[["x:z"]]$effect2__))
})

test_that("re_formula = NULL conditions on a NEW group", {
  cs <- ce_case()
  pop <- conditional_effects(cs$ds, effects = "x", resolution = 6)
  new <- conditional_effects(cs$ds, effects = "x", resolution = 6,
                             re_formula = NULL, seed = 3)
  # the grouping column says the level is unobserved, as the fit
  # method's does; before this the call silently took group 1
  expect_true(all(is.na(new$x$g)))
  expect_identical(names(new$x), names(pop$x))
  # and the band carries that group's variance, so it is wider at every
  # grid point rather than the population band under another name
  expect_true(all(new$x$upper__ - new$x$lower__ >
                    pop$x$upper__ - pop$x$lower__))
  # seed = makes the drawn group reproducible
  again <- conditional_effects(cs$ds, effects = "x", resolution = 6,
                               re_formula = NULL, seed = 3)
  expect_equal(new$x$upper__, again$x$upper__)
})

test_that("an ordinal draws display is keyed and laid out like core's", {
  cs <- ce_ord_case()
  cd <- conditional_effects(cs$ds, effects = "x", resolution = 5)
  cf <- conditional_effects(cs$fit, effects = "x", resolution = 5,
                            band = "boot", boot = 20, seed = 5)
  # brms keys the per-category layout "x:cats__"; this method keyed it
  # "x", so brms's own plot() could not read the pair back
  expect_identical(names(cd), "x:cats__")
  expect_identical(names(cd), names(cf))
  expect_identical(names(cd[["x:cats__"]]), names(cf[["x:cats__"]]))
  expect_equal(grid_of(cd[["x:cats__"]]), grid_of(cf[["x:cats__"]]),
               ignore_attr = TRUE)
  expect_identical(attr(cd[["x:cats__"]], "effects"), c("x", "cats__"))
  # the probabilities of one grid row sum to one
  d <- cd[["x:cats__"]]
  s <- tapply(d$estimate__, d$x, sum)
  expect_true(all(abs(s - 1) < 1e-8))
})

test_that("categorical = FALSE draws the expected category number", {
  cs <- ce_ord_case()
  cd <- conditional_effects(cs$ds, effects = "x", resolution = 5,
                            categorical = FALSE)
  cf <- conditional_effects(cs$fit, effects = "x", resolution = 5,
                            categorical = FALSE)
  expect_identical(names(cd), "x")
  expect_identical(names(cd$x), names(cf$x))
  expect_equal(grid_of(cd$x), grid_of(cf$x), ignore_attr = TRUE)
  # sum_k k p_k, checked against the per-category display of the SAME
  # draws rather than against the range (1, 3), which holds for any
  # probability vector and so could not fail. The map is linear, so the
  # posterior mean of the weighted sum IS the weighted sum of the
  # posterior means, exactly.
  per_cat <- conditional_effects(cs$ds, effects = "x", resolution = 5)
  pc <- per_cat[["x:cats__"]]
  ncat <- nlevels(pc$cats__)
  ngrid <- nrow(pc) / ncat
  by_hand <- rowSums(vapply(seq_len(ncat), function(k) {
    k * pc$estimate__[(k - 1L) * ngrid + seq_len(ngrid)]
  }, numeric(ngrid)))
  expect_equal(cd$x$estimate__, by_hand)
  if (sampler_gates_on()) {
    expect_lt(max(abs(cd$x$estimate__ - cf$x$estimate__)), 0.05)
  }
  # categorical = TRUE is the default here, and refuses a dpar
  expect_error(conditional_effects(cs$ds, dpar = "mu", categorical = TRUE),
               "needs an ordinal or")
})

test_that("the draws method reports the arguments it cannot use", {
  cs <- ce_case()
  expect_warning(conditional_effects(cs$ds, effects = "x", resolution = 4,
                                     nonesuch = 1),
                 "ignoring unknown argument")
  expect_error(conditional_effects(cs$ds, effects = "x", re.form = NA),
               "spells this argument")
})


## ---- the plotted quantity, on families whose mean is not mu ---------
##
## Column and grid identity cannot see this: the NAMES match whether the
## curve is the expected response or the mu predictor. What separates
## them is the VALUE, so every block below compares estimate__ against
## the fit method, and each also carries one assertion that needs no
## chain agreement at all - a structural property the mu curve does not
## have.

mean_case <- local({
  cache <- list()
  function(which) {
    skip_sampler()
    if (is.null(cache[[which]])) {
      n <- 200L
      if (which == "zi") {
        set.seed(11)
        dd <- data.frame(x = stats::rnorm(n))
        mu <- exp(0.4 + 0.65 * dd$x)
        dd$y <- ifelse(stats::rbinom(n, 1, 0.265) == 1L, 0L,
                       stats::rpois(n, mu))
        fit <- frm(bf(y ~ x) + zero_inflated_poisson(), data = dd)
      } else if (which == "trunc") {
        set.seed(13)
        # drawn and REJECTED below the bound, so the density the model
        # fits is the one the data came from
        dd <- data.frame(x = stats::rnorm(3L * n))
        dd$y <- stats::rnorm(nrow(dd), 1.4 + 1.0 * dd$x, 1)
        dd <- dd[dd$y > 0, ][seq_len(n), ]
        rownames(dd) <- NULL
        fit <- frm(bf(y | trunc(lb = 0) ~ x) + gaussian(), data = dd)
      } else {
        set.seed(14)
        dd <- data.frame(x = stats::rnorm(n))
        z <- stats::rbinom(n, 1, 0.4)
        dd$y <- ifelse(z == 1L, stats::rnorm(n, -2 + 0.9 * dd$x, 0.7),
                       stats::rnorm(n, 2 + 0.9 * dd$x, 0.7))
        fit <- frm(bf(y ~ x) + mixture(gaussian(), gaussian()), data = dd,
                   start = list(beta = c(-2, 0.9, 2, 0.9)))
      }
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 800, refresh = 0, seed = 21)))
      cache[[which]] <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache[[which]]
  }
})

# max deviation from the fit curve, relative to the fit curve. The
# defect this guards was 36% on the zero-inflated fit and a whole
# different curve on the other two, so the gate is nowhere near it.
rel_gap <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-8))

test_that("ce_pred_dpar() is the seam that decides what is predicted", {
  # the export itself, before any sampling: NULL means the expected
  # response, and only a family whose mean IS the inverse link of mu
  # (with no truncation) keeps its dpar
  zi <- mean_case("zi")
  tr <- mean_case("trunc")
  expect_null(ce_pred_dpar(zi$fit$spec$responses[[1L]], "mu"))
  expect_null(ce_pred_dpar(tr$fit$spec$responses[[1L]], "mu"))
  # a named dpar is honored, and so is a category display
  expect_identical(ce_pred_dpar(zi$fit$spec$responses[[1L]], "zi",
                                dpar_given = TRUE), "zi")
  expect_identical(ce_pred_dpar(zi$fit$spec$responses[[1L]], "mu",
                                categorical = TRUE), "mu")
  # and a plain gaussian is unchanged
  cs <- ce_case()
  expect_identical(ce_pred_dpar(cs$fit$spec$responses[[1L]], "mu"), "mu")
})

test_that("a zero-inflated draws curve is the mean, not mu", {
  cs <- mean_case("zi")
  cd <- conditional_effects(cs$ds, effects = "x", resolution = 5)
  cf <- conditional_effects(cs$fit, effects = "x", resolution = 5)
  expect_identical(names(cd$x), names(cf$x))
  # STRUCTURAL, no chain agreement needed: the mean is (1 - zi) times
  # mu with zi > 0, so the default display must sit strictly below the
  # mu display. Before ce_pred_dpar() the two were the same curve.
  mu_only <- conditional_effects(cs$ds, effects = "x", resolution = 5,
                                 dpar = attr(cd$x, "dpar"))
  expect_true(all(cd$x$estimate__ < mu_only$x$estimate__))
  if (sampler_gates_on()) {
    expect_lt(rel_gap(cd$x$estimate__, cf$x$estimate__), 0.12)
  }
})

test_that("a truncated draws curve stays inside the support", {
  cs <- mean_case("trunc")
  cd <- conditional_effects(cs$ds, effects = "x", resolution = 5)
  cf <- conditional_effects(cs$fit, effects = "x", resolution = 5)
  expect_identical(names(cd$x), names(cf$x))
  # STRUCTURAL: a response truncated below at zero has a positive mean
  # everywhere. The mu predictor does not, and plot() drew the negative
  # value it reported.
  expect_true(all(cd$x$estimate__ > 0))
  expect_true(all(cd$x$lower__ > 0))
  mu_only <- conditional_effects(cs$ds, effects = "x", resolution = 5,
                                 dpar = attr(cd$x, "dpar"))
  expect_true(all(cd$x$estimate__ > mu_only$x$estimate__))
  if (sampler_gates_on()) {
    expect_lt(rel_gap(cd$x$estimate__, cf$x$estimate__), 0.12)
  }
})

test_that("a mixture draws curve is the mixture mean", {
  cs <- mean_case("mix")
  cd <- conditional_effects(cs$ds, effects = "x", resolution = 5)
  cf <- conditional_effects(cs$fit, effects = "x", resolution = 5)
  expect_identical(names(cd$x), names(cf$x))
  # STRUCTURAL: the components are separated by about 4, so the
  # theta-weighted mean cannot coincide with either component's mu
  mu_only <- conditional_effects(cs$ds, effects = "x", resolution = 5,
                                 dpar = attr(cd$x, "dpar"))
  expect_gt(min(abs(cd$x$estimate__ - mu_only$x$estimate__)), 0.5)
  if (sampler_gates_on()) {
    expect_lt(rel_gap(cd$x$estimate__, cf$x$estimate__), 0.12)
    # the explicit non-mu dpar path was already right and must stay so:
    # dpar_report_hook() softmax reporting lives inside
    # predict(type = "response", dpar = )
    tf <- conditional_effects(cs$fit, effects = "x", resolution = 5,
                              dpar = "theta1")
    td <- conditional_effects(cs$ds, effects = "x", resolution = 5,
                              dpar = "theta1")
    expect_lt(max(abs(td$x$estimate__ - tf$x$estimate__)), 0.05)
  }
})

test_that("allow_new_levels is accepted, not warned about", {
  cs <- ce_case()
  # core pulls both spellings out of the dots; the draws method used to
  # warn "ignoring unknown argument" for the same call
  expect_no_warning(
    conditional_effects(cs$ds, effects = "x", resolution = 4,
                        allow_new_levels = TRUE))
  expect_no_warning(
    conditional_effects(cs$ds, effects = "x", resolution = 4,
                        allow.new.levels = TRUE))
  expect_warning(conditional_effects(cs$ds, effects = "x", resolution = 4,
                                     nonesuch = 1),
                 "ignoring unknown argument")
})

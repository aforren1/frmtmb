# `route` on get_prior(): the table a caller gets must be a property of
# the call, not of the search path.
#
# The defect this file exists for: get_prior() consulted the
# prior-defaults registry unconditionally, so the same call on the same
# fit printed "(flat)" everywhere with frmtmb alone and brms's
# weakly-informative defaults after library(frmtmb.sample).
# dev/getprior-load-order-probe.R is that demonstration; the first test
# below is the same comparison with the attach replaced by the registry
# swap, because the two states have to exist inside one process.

# A stand-in for what a sampling package registers. Deliberately NOT
# frmtmb.sample's: this file must fail for the right reason whether or
# not that package is installed, so the provider is local and its
# densities are recognizable.
route_test_provider <- function(spec, frame) {
  set_prior("normal(0, 3)", class = "Intercept") +
    set_prior("exponential(2)", class = "sd")
}

# The registry is append-only for its real users; `swap_prior_defaults()`
# is the test seam that lets one process visit both states.
local_providers <- function(providers, env = parent.frame()) {
  old <- frmtmb:::swap_prior_defaults(providers)
  withr::defer(frmtmb:::swap_prior_defaults(old), envir = env)
  invisible(old)
}

route_test_data <- function() {
  set.seed(11)
  dd <- data.frame(y = 0, x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  # the prior table is a property of the call, not of the numbers, so
  # the response is drawn from the model the tests then fit
  # The draw takes its own seed, away from the fixture's: reusing
  # the fixture seed restarts the same random stream that made the
  # covariates, and the residuals come out equal to x.
  dd$y <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                       newparams = list(Intercept = 0, x = 0.5, sigma = 1,
                                        sd_g__Intercept = 0.5),
                       nsim = 1, seed = 1011)[[1]]
  dd
}

test_that("route = 'fit' answers the same with and without a provider", {
  dd <- route_test_data()
  form <- bf(y ~ x + (1 | g)) + gaussian()
  fit <- frm(form, data = dd)

  local_providers(list())
  bare_form <- get_prior(form, data = dd)
  bare_fit <- get_prior(fit)

  local_providers(list(route_test_provider))
  # the provider IS live: the sample route reads it
  expect_true(any(get_prior(form, data = dd, route = "sample")$prior ==
                    "normal(0, 3)"))
  # ... and the fit route is untouched by it, which is the whole point
  expect_identical(get_prior(form, data = dd), bare_form)
  expect_identical(get_prior(fit), bare_fit)
  expect_true(all(bare_form$prior == "(flat)"))
  expect_true(all(bare_fit$prior == "(flat)"))
})

test_that("route = 'sample' refuses when no package states the defaults", {
  dd <- route_test_data()
  form <- bf(y ~ x + (1 | g)) + gaussian()

  local_providers(list())
  expect_error(get_prior(form, data = dd, route = "sample"),
               "frmtmb.sample", fixed = TRUE)
  expect_error(get_prior(form, data = dd, route = "sample"),
               "library(frmtmb.sample)", fixed = TRUE)
  # and it says what to ask for instead
  expect_error(get_prior(form, data = dd, route = "sample"),
               "route = \"fit\"", fixed = TRUE)
  # a fit is refused on the same grounds, before any frame work
  fit <- frm(form, data = dd)
  expect_error(get_prior(fit, route = "sample"), "no loaded package")
})

test_that("route = 'sample' reports the registered defaults", {
  dd <- route_test_data()
  form <- bf(y ~ x + (1 | g)) + gaussian()
  local_providers(list(route_test_provider))

  gp <- get_prior(form, data = dd, route = "sample")
  # dpar narrows to mu's intercept. This model gives sigma no
  # predictor, so sigma has a row under its OWN class, which this
  # provider does not address and which therefore stays flat: a slot
  # key is a class plus its qualifiers
  expect_identical(gp$prior[gp$class == "Intercept" & gp$coef == "" &
                              gp$dpar == ""],
                   "normal(0, 3)")
  expect_identical(gp$prior[gp$class == "sigma"], "(flat)")
  expect_false(any(gp$class == "Intercept" & gp$dpar == "sigma"))
  expect_true(all(gp$prior[gp$class == "sd" & gp$coef == ""] ==
                    "exponential(2)"))
  # a default speaks for a class, not for one coefficient of it
  expect_true(all(gp$prior[nzchar(gp$coef)] == "(flat)"))
  # the slots the provider says nothing about stay flat rather than
  # borrowing a neighbor's density
  expect_true(all(gp$prior[gp$class == "b"] == "(flat)"))
})

test_that("an unknown route is refused rather than guessed", {
  dd <- route_test_data()
  expect_error(get_prior(bf(y ~ x) + gaussian(), data = dd,
                         route = "sampling"))
})

test_that("the printed table names its route", {
  dd <- route_test_data()
  form <- bf(y ~ x + (1 | g)) + gaussian()

  local_providers(list())
  gp <- get_prior(form, data = dd)
  expect_s3_class(gp, "frmtmb_prior_rows")
  expect_s3_class(gp, "data.frame")
  expect_identical(attr(gp, "route"), "fit")
  out <- utils::capture.output(print(gp))
  expect_match(out[1], "route = \"fit\"", fixed = TRUE)
  expect_match(out[1], "frm()", fixed = TRUE)

  local_providers(list(route_test_provider))
  gs <- get_prior(form, data = dd, route = "sample")
  expect_identical(attr(gs, "route"), "sample")
  outs <- utils::capture.output(print(gs))
  expect_match(outs[1], "route = \"sample\"", fixed = TRUE)
  expect_match(outs[1], "frm_sample()", fixed = TRUE)
  # the rows still print, under the label
  expect_true(any(grepl("Intercept", outs, fixed = TRUE)))
})

test_that("the row table is still an ordinary data frame to work with", {
  dd <- route_test_data()
  gp <- get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_true(is.data.frame(gp))
  expect_identical(names(gp),
                   c("prior", "class", "coef", "group", "dpar", "nlpar",
                     "resp", "lb", "ub"))
  sub <- gp[gp$class == "sd", ]
  expect_gt(nrow(sub), 0L)
  # fewer ROWS of a fit-route table are still a fit-route table, so the
  # label survives the subset rather than silently going missing
  expect_identical(attr(sub, "route"), "fit")
  expect_match(utils::capture.output(print(sub))[1], "route = \"fit\"",
               fixed = TRUE)

  # a COLUMN subset keeps the class but loses the attribute, because
  # `[.data.frame` drops it. That is wanted rather than a `[` method
  # that carries it through: a selection that has dropped the `prior`
  # column has no route left to describe. This is also the only
  # ordinary way to reach the print method's unlabeled branch.
  cols <- gp[, c("class", "coef")]
  expect_s3_class(cols, "frmtmb_prior_rows")
  expect_null(attr(cols, "route"))
  expect_false(grepl("route =", utils::capture.output(print(cols))[1],
                     fixed = TRUE))
})

test_that("coercion drops the class and the label together", {
  dd <- route_test_data()
  gp <- get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  plain <- as.data.frame(gp)
  # a coerced table that kept a stray `route` attribute would not be
  # identical() to the data frame a reader would write by hand, which
  # is the one thing coercion is asked for
  hand <- gp
  attr(hand, "route") <- NULL
  class(hand) <- "data.frame"
  expect_identical(class(plain), "data.frame")
  expect_null(attr(plain, "route"))
  expect_identical(plain, hand)
  expect_identical(plain[["prior"]], gp[["prior"]])
})

test_that("every row get_prior() lists is one set_prior() accepts", {
  # F1 of the punch round, and the test that stops it recurring. The
  # table is a promise: each row names a slot the caller can address.
  # A predictor-free dpar was routed to its own class without asking
  # whether that class is one the package REFUSES, so a mixture listed
  # class = "theta1", set_prior() refused it and named the Intercept
  # spelling, and the shape gate refused that and named the class back.
  # Advertised by the table, reachable by nothing.
  usable <- function(fit, gp) {
    for (i in seq_len(nrow(gp))) {
      r <- gp[i, ]
      lbl <- paste0("row ", i, ": ", r$class,
                    if (nzchar(r$coef)) paste0("/", r$coef),
                    if (nzchar(r$dpar)) paste0(" dpar=", r$dpar),
                    if (nzchar(r$nlpar)) paste0(" nlpar=", r$nlpar),
                    if (nzchar(r$group)) paste0(" group=", r$group))
      # class "cor" and its relatives take lkj() and nothing else, so
      # the density follows the class rather than the other way round
      dist <- if (r$class %in% c("cor", "cortime", "rescor")) {
        "lkj(2)"
      } else {
        "normal(0, 1)"
      }
      pl <- tryCatch(
        set_prior(dist, class = r$class, coef = r$coef, group = r$group,
                  dpar = r$dpar, nlpar = r$nlpar, resp = r$resp),
        error = function(e) {
          fail(paste0(lbl, " is refused by set_prior(): ",
                      conditionMessage(e)))
          NULL
        })
      if (is.null(pl)) next
      # and it must reach a parameter of the model it was listed for,
      # which is where the shape gate lives
      ok <- tryCatch({
        frmtmb:::resolve_prior_input(fit, pl)
        TRUE
      }, error = function(e) {
        fail(paste0(lbl, " does not resolve: ", conditionMessage(e)))
        FALSE
      })
      expect_true(ok)
    }
  }

  set.seed(19)
  dd <- data.frame(x = stats::rnorm(200), g = factor(rep(1:20, 10)))
  dd$y <- stats::rnorm(200, 1 + 0.5 * dd$x, 1.4)
  dd$cnt <- stats::rpois(200, 3)
  dd$p <- stats::rbeta(200, 2, 3)

  shapes <- list(
    ordinary = bf(y ~ x + (1 | g)) + gaussian(),
    distributional = bf(y ~ x, sigma ~ x) + gaussian(),
    intercept_only_dpar = bf(y ~ x, sigma ~ 1) + gaussian(),
    negbinomial = bf(cnt ~ x) + negbinomial(),
    zero_inflated = bf(cnt ~ x) + zero_inflated_poisson(),
    beta = bf(p ~ x) + Beta(),
    # the shape that produced the defect: a mixture proportion is a
    # predictor-free dpar whose class frmtmb refuses by name
    mixture = bf(y ~ x) + mixture(gaussian(), gaussian()))

  for (nm in names(shapes)) {
    form <- shapes[[nm]]
    gp <- get_prior(form, data = dd)
    fit <- suppressWarnings(frm(form, data = dd, dry_run = "objective"))
    usable(fit, gp)
  }

  # and the mixture's unreachable proportion is left OUT rather than
  # advertised: brms honors class = "theta" there and frmtmb has no
  # spelling at all, which the refusal says
  gpm <- get_prior(bf(y ~ x) + mixture(gaussian(), gaussian()), data = dd)
  expect_false(any(startsWith(gpm$class, "theta")))
  expect_false(any(startsWith(gpm$dpar, "theta")))
  # the component dispersions ARE listed, under their own classes
  expect_true(all(c("sigma1", "sigma2") %in% gpm$class))
})

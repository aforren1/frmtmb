# The brms prior interface, as a port of the brms nonlinear vignette
# hits it: the argument spelled `prior =` as brms spells it, brms's
# quoting `prior()` constructor, `nlpar =` addressing, and a brms-built
# prior object arriving because brms was attached and masked ours.
#
# The acceptance case is the vignette's own cumulative-loss model,
# `cum ~ ult * (1 - exp(-(dev/theta)^omega))`. Its data is not shipped
# with brms (the vignette reads a csv from GitHub), so the analogue
# below is simulated from the vignette's own parameter values: ten
# accident years, ten development lags, an ultimate loss around 5000
# with a between-year sd, omega 1.3 and theta 45.

loss_data <- function(seed = 903) {
  set.seed(seed)
  AY <- factor(rep(1988:1997, each = 10))
  dev <- rep(seq(6, 114, by = 12), 10)
  ult_g <- 5000 + stats::rnorm(10, 0, 400)
  cum <- ult_g[as.integer(AY)] * (1 - exp(-(dev / 45)^1.3)) +
    stats::rnorm(100, 0, 120)
  data.frame(cum = cum, dev = dev, AY = AY)
}

loss_form <- function() {
  bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
     ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE) + gaussian()
}

# the vignette's own starting region; a nonlinear body this shaped has
# no useful default start
loss_start <- list(beta = c(5000, 1, 45))

# ---- prior(), prior_(), prior_string() -------------------------------

test_that("prior() quotes its first argument, as brms's does", {
  expect_equal(unclass(prior(normal(5000, 1000), nlpar = "ult")),
               unclass(set_prior("normal(5000, 1000)", nlpar = "ult")))
  # every argument is deparsed, so the unquoted and quoted spellings of
  # a class agree
  expect_equal(unclass(prior(normal(0, 1), class = b)),
               unclass(prior(normal(0, 1), class = "b")))
  # combining works through the same c() / `+` methods
  pl <- c(prior(normal(5000, 1000), nlpar = "ult"),
          prior(normal(1, 2), nlpar = "omega"),
          prior(normal(45, 10), nlpar = "theta"))
  expect_s3_class(pl, "frmtmb_priorlist")
  expect_length(unclass(pl), 3L)
  expect_length(unclass(prior(normal(0, 1), class = "b") +
                          prior(exponential(1), class = "sd")), 2L)
  # a deparsed bound is a number by the time it is stored
  expect_identical(unclass(prior(normal(0, 1), lb = 0))[[1L]]$lb, 0)
})

test_that("prior_() and prior_string() are the programmatic spellings", {
  expect_equal(unclass(prior_(~normal(0, 10), class = ~b)),
               unclass(set_prior("normal(0, 10)", class = "b")))
  expect_equal(unclass(prior_(~student_t(3, 0, 2), class = "Intercept")),
               unclass(set_prior("student_t(3, 0, 2)",
                                 class = "Intercept")))
  # prior_string() takes a string computed at run time, which is
  # exactly what prior() cannot do
  s <- paste0("normal(0, ", 2 * 5, ")")
  expect_equal(unclass(prior_string(s, class = "b")),
               unclass(set_prior("normal(0, 10)", class = "b")))
  expect_error(prior_(list(1), class = "b"), "one-sided formulas")
})

test_that("prior() is frmtmb's only when frmtmb's is the one in scope", {
  # the masking claim this file rests on, asserted rather than assumed:
  # nothing but brms exports prior() into the search path, so the only
  # way frmtmb's is shadowed is brms being attached after it
  expect_identical(environmentName(environment(frmtmb::prior)),
                   "frmtmb")
  expect_s3_class(frmtmb::prior(normal(0, 1)), "frmtmb_priorlist")
})

# ---- set_prior(nlpar =) ----------------------------------------------

test_that("class b with nlpar covers the parameter's whole coefficient vector", {
  dd <- loss_data()
  fit <- frm(loss_form(), data = dd, start = loss_start)

  # the vignette's exact spelling lands on ult_(Intercept), because a
  # nonlinear parameter's sub-formula is not centered and brms holds
  # its intercept in the same coefficient vector as its slopes
  ri <- frmtmb:::resolve_prior_input(
    fit, set_prior("normal(5000, 1000)", nlpar = "ult"))
  expect_length(ri$entries, 1L)
  e <- ri$entries[[1L]]
  expect_identical(e$comp, "beta")
  expect_identical(names(fit$frame$par_template$beta)[e$idx],
                   "ult_(Intercept)")
  expect_equal(e$dist$location, 5000)

  # the quoting constructor resolves to the same entry
  ri_q <- frmtmb:::resolve_prior_input(
    fit, prior(normal(5000, 1000), nlpar = "ult"))
  expect_equal(ri_q$entries, ri$entries)

  # ... and it is one nonlinear parameter, not all three
  for (np in c("omega", "theta")) {
    ri_np <- frmtmb:::resolve_prior_input(
      fit, set_prior("normal(0, 1)", nlpar = np))
    expect_identical(names(fit$frame$par_template$beta)[
      ri_np$entries[[1L]]$idx], paste0(np, "_(Intercept)"))
  }
})

test_that("nlpar takes coef, class Intercept, and the sd/cor classes", {
  dd <- loss_data()
  fit <- frm(loss_form(), data = dd, start = loss_start)
  target <- function(pl) {
    names(fit$frame$par_template$beta)[
      frmtmb:::resolve_prior_input(fit, pl)$entries[[1L]]$idx]
  }
  # brms writes an intercept as "Intercept"; the design matrix spells
  # it "(Intercept)", and both name the same column
  expect_identical(target(set_prior("normal(5000, 1000)", nlpar = "ult",
                                    coef = "Intercept")),
                   "ult_(Intercept)")
  expect_identical(target(set_prior("normal(5000, 1000)", nlpar = "ult",
                                    coef = "(Intercept)")),
                   "ult_(Intercept)")
  expect_identical(target(set_prior("normal(5000, 1000)",
                                    class = "Intercept", nlpar = "ult")),
                   "ult_(Intercept)")

  # class "sd" narrows to the blocks of one nonlinear parameter
  ri <- frmtmb:::resolve_prior_input(
    fit, set_prior("exponential(0.01)", class = "sd", nlpar = "ult"))
  expect_length(ri$entries, 1L)
  expect_identical(ri$entries[[1L]]$scale, "sd")
  expect_error(frmtmb:::resolve_prior_input(
    fit, set_prior("exponential(1)", class = "sd", nlpar = "omega")),
    "No random-effect SDs")
})

test_that("nlpar separates two blocks on the SAME grouping factor", {
  # the case group= alone cannot address: two nonlinear parameters each
  # varying by g, so `group = "g"` names both and only nlpar tells them
  # apart. brms addresses it the same way
  set.seed(909)
  n <- 200L
  g <- factor(rep(1:20, each = 10))
  x <- stats::runif(n, 0, 5)
  u <- matrix(stats::rnorm(40), 2)
  dd <- data.frame(
    y = (2.5 + 0.4 * u[1, g]) * exp(-(0.7 + 0.1 * u[2, g]) * x) +
      stats::rnorm(n, 0, 0.1),
    x = x, g = g)
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1 + (1 | g),
                nl = TRUE) + gaussian(),
             data = dd, start = list(beta = c(2, 0.5)))
  expect_length(fit$frame$re_blocks, 2L)

  th_of <- function(pl) {
    sort(vapply(frmtmb:::resolve_prior_input(fit, pl)$entries,
                function(e) as.numeric(e$idx), 0))
  }
  both <- th_of(set_prior("exponential(1)", class = "sd", group = "g"))
  expect_length(both, 2L)
  only_a <- th_of(set_prior("exponential(1)", class = "sd",
                            nlpar = "a"))
  only_b <- th_of(set_prior("exponential(1)", class = "sd",
                            nlpar = "b"))
  expect_length(only_a, 1L)
  expect_length(only_b, 1L)
  expect_false(identical(only_a, only_b))
  expect_setequal(c(only_a, only_b), both)

  # get_prior() lists the two blocks apart on the same grounds
  gp <- get_prior(fit)
  sd_rows <- gp[gp$class == "sd" & gp$group == "g", ]
  expect_setequal(sd_rows$nlpar, c("a", "b"))
})

test_that("nlpar refuses what it cannot address, and never guesses", {
  dd <- loss_data()
  fit <- frm(loss_form(), data = dd, start = loss_start)
  # a typo names no nonlinear parameter, and the refusal lists the ones
  # the model has rather than resolving to nothing
  expect_error(frmtmb:::resolve_prior_input(
    fit, set_prior("normal(0, 1)", nlpar = "ULT")),
    "names no nonlinear parameter")
  # a DISTRIBUTIONAL parameter is not a nonlinear one, and the message
  # says which argument names it
  expect_error(frmtmb:::resolve_prior_input(
    fit, set_prior("normal(0, 1)", nlpar = "sigma")), "dpar =")
  # a linear model has no nonlinear parameters at all
  ld <- data.frame(y = stats::rnorm(40), x = stats::rnorm(40))
  lf <- frm(bf(y ~ x) + gaussian(), data = ld)
  expect_error(frmtmb:::resolve_prior_input(
    lf, set_prior("normal(0, 1)", nlpar = "a")), "nl = TRUE")
  # naming both dpar and nlpar is a question about intent
  expect_error(set_prior("normal(0, 1)", dpar = "sigma", nlpar = "ult"),
               "not both")
  # a coefficient that does not exist under this nlpar
  expect_error(frmtmb:::resolve_prior_input(
    fit, set_prior("normal(0, 1)", nlpar = "ult", coef = "zzz")),
    "not found")
})

test_that("get_prior lists a nonlinear parameter the way brms does", {
  dd <- loss_data()
  gp <- get_prior(loss_form(), data = dd)
  expect_true("nlpar" %in% names(gp))
  # class "b" with the parameter in the nlpar column, the intercept
  # among the coefficients rather than in its own class
  for (np in c("ult", "omega", "theta")) {
    expect_true(any(gp$class == "b" & gp$nlpar == np & gp$coef == ""))
    expect_true(any(gp$class == "b" & gp$nlpar == np &
                      gp$coef == "(Intercept)"))
    expect_false(any(gp$class == "Intercept" & gp$nlpar == np))
  }
  # the variance component is listed against the parameter that owns it
  expect_true(any(gp$class == "sd" & gp$group == "AY" &
                    gp$nlpar == "ult"))
  # and every listed nonlinear row round-trips into set_prior()
  fit <- frm(loss_form(), data = dd, start = loss_start)
  rows <- gp[gp$class == "b" & nzchar(gp$nlpar), ]
  for (i in seq_len(nrow(rows))) {
    pl <- set_prior("normal(0, 1000)", class = "b",
                    coef = rows$coef[i], nlpar = rows$nlpar[i])
    expect_gt(length(frmtmb:::resolve_prior_input(fit, pl)$entries), 0L)
  }
})

# ---- the argument is `prior`, and only `prior` ------------------------

test_that("every entry point spells the argument prior, brms's name", {
  # the contract, so an entry point added later cannot ship the retired
  # name: one spelling, brms's, with no alias behind it
  # frm_sample() is asserted the same way in frmtmb.sample's suite
  for (fn in c("frm", "frm_simulate")) {
    fo <- names(formals(getFromNamespace(fn, "frmtmb")))
    expect_true("prior" %in% fo, label = paste0(fn, " takes `prior`"))
    expect_false("priors" %in% fo,
                 label = paste0(fn, " does not take `priors`"))
  }
  # no dual-spelling machinery is wired to this pair anywhere
  expect_false(any(grepl("priors", deparse(frmtmb::frm), fixed = TRUE)))
})

test_that("the retired priors= spelling fails rather than passing", {
  dd <- loss_data()
  pl <- prior(normal(5000, 1000), nlpar = "ult")
  # R cannot partially match `priors` to `prior` (a longer name is not
  # a prefix), so a direct call fails on its own
  expect_error(frm(loss_form(), data = dd, priors = pl),
               "unused argument")
  expect_error(frm_simulate(bf(y ~ x) + gaussian(),
                            data.frame(x = 1:5, y = 0), priors = pl),
               "unused argument")
  # frm_sample()'s `...` WOULD have swallowed it, so that one is
  # refused by name rather than by R; frmtmb.sample asserts it
})

test_that("prior = takes the same specification at every entry point", {
  dd <- loss_data()
  pl <- prior(normal(5000, 1000), nlpar = "ult")
  a <- frm(loss_form(), data = dd, start = loss_start, prior = pl)
  b <- frm(loss_form(), data = dd, start = loss_start,
           prior = set_prior("normal(5000, 1000)", nlpar = "ult"))
  expect_equal(fixef(a)$ult, fixef(b)$ult)
  expect_equal(unclass(prior_summary(a)), unclass(prior_summary(b)))

  # no argument means no priors, plain ML, and the field the fit
  # carries is `prior`, as a brmsfit spells it
  ml <- frm(loss_form(), data = dd, start = loss_start)
  expect_null(ml$prior)
  expect_false("priors" %in% names(ml))
  expect_output(prior_summary(ml), "No priors were set")
  expect_s3_class(a$prior, "frmtmb_priorlist")

  set.seed(77)
  sdd <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)),
                    y = 0)
  spl <- prior(normal(0, 1), class = "b") +
    prior(normal(0, 2), class = "Intercept") +
    prior(exponential(1), class = "sd") +
    prior(exponential(1), class = "Intercept", dpar = "sigma")
  sims <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), sdd,
                       prior = spl, nsim = 2, seed = 1)
  expect_equal(dim(sims), c(60L, 2L))
  expect_false(is.null(attr(sims, "pars")))
})

# ---- the acceptance case: the vignette's loss model -------------------

test_that("the loss model takes the vignette's priors and they bind", {
  dd <- loss_data()
  form <- loss_form()
  vignette_priors <- c(prior(normal(5000, 1000), nlpar = "ult"),
                       prior(normal(1, 2), nlpar = "omega"),
                       prior(normal(45, 10), nlpar = "theta"))

  ml <- frm(form, data = dd, start = loss_start)
  map <- frm(form, data = dd, start = loss_start,
             prior = vignette_priors)
  expect_output(print(map), "MAP")

  # the vignette's own priors are weak next to 100 observations, so
  # what they must do is move the estimate toward their location and
  # leave the fit recognizable
  expect_lt(abs(fixef(map)$ult[[1L]] - 5000),
            abs(fixef(ml)$ult[[1L]] - 5000))
  expect_lt(abs(fixef(map)$ult[[1L]] - fixef(ml)$ult[[1L]]), 200)

  # a TIGHT prior in the same spelling proves the density is really on
  # ult and not merely accepted: the estimate follows it
  tight <- frm(form, data = dd, start = loss_start,
               prior = c(prior(normal(4000, 40), nlpar = "ult"),
                         prior(normal(1, 2), nlpar = "omega"),
                         prior(normal(45, 10), nlpar = "theta")))
  expect_lt(abs(fixef(tight)$ult[[1L]] - 4000), 250)
  expect_lt(fixef(tight)$ult[[1L]], fixef(ml)$ult[[1L]])

  # the penalized objective is the likelihood plus these three
  # densities, evaluated at the MAP solution
  nlp <- -sum(stats::dnorm(fixef(map)$ult[[1L]], 5000, 1000, log = TRUE),
              stats::dnorm(fixef(map)$omega[[1L]], 1, 2, log = TRUE),
              stats::dnorm(fixef(map)$theta[[1L]], 45, 10, log = TRUE))
  raw <- ml$obj$fn(map$opt$par)
  expect_lt(abs((-as.numeric(logLik(map))) - (raw + nlp)), 1e-6)

  # prior_summary() names the parameter each density landed on
  out <- utils::capture.output(prior_summary(map))
  expect_match(out[1L], "normal(5000, 1000) class=b nlpar=ult",
               fixed = TRUE)
  expect_match(out[2L], "nlpar=omega", fixed = TRUE)
  expect_match(out[3L], "nlpar=theta", fixed = TRUE)

  # set_prior()'s string spelling is the same call
  strs <- set_prior("normal(5000, 1000)", nlpar = "ult") +
    set_prior("normal(1, 2)", nlpar = "omega") +
    set_prior("normal(45, 10)", nlpar = "theta")
  expect_equal(fixef(frm(form, data = dd, start = loss_start,
                         prior = strs))$ult,
               fixef(map)$ult)
})


# ---- set_prior(resp =), the other addressing gap ----------------------

test_that("resp picks one response of a multivariate model", {
  set.seed(505)
  n <- 120L
  dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:12, 10)))
  dd$y1 <- stats::rnorm(n, 1 + 0.5 * dd$x, 1)
  dd$y2 <- stats::rnorm(n, -1 + 0.3 * dd$x, 1)
  fit <- frm(mvbf(bf(y1 ~ x + (1 | g)), bf(y2 ~ x)) + gaussian(),
             data = dd)

  ri <- frmtmb:::resolve_prior_input(
    fit, set_prior("normal(0, 0.1)", class = "b", resp = "y2"))
  expect_length(ri$entries, 1L)
  expect_identical(names(fit$frame$par_template$beta)[
    ri$entries[[1L]]$idx], "y2_x")
  # without it, the class covers both responses
  expect_length(frmtmb:::resolve_prior_input(
    fit, set_prior("normal(0, 0.1)", class = "b"))$entries, 2L)
  # class "sd" narrows to the response that owns the block
  expect_length(frmtmb:::resolve_prior_input(
    fit, set_prior("exponential(1)", class = "sd",
                   resp = "y1"))$entries, 1L)
  expect_error(frmtmb:::resolve_prior_input(
    fit, set_prior("exponential(1)", class = "sd", resp = "y2")),
    "No random-effect SDs")
  expect_error(frmtmb:::resolve_prior_input(
    fit, set_prior("normal(0, 1)", class = "b", resp = "zzz")),
    "not found")

  # get_prior() lists the same addressing it accepts
  gp <- get_prior(fit)
  expect_setequal(unique(gp$resp[gp$class == "b"]), c("y1", "y2"))
  expect_true(any(gp$class == "sd" & gp$group == "g" & gp$resp == "y1"))
})

# ---- a brms-built prior object ----------------------------------------

test_that("a brmsprior object is translated rather than refused", {
  skip_if_not_installed("brms")
  # exactly what a ported script produces once brms is attached and its
  # prior() masks frmtmb's
  bp <- c(brms::prior(normal(5000, 1000), nlpar = "ult"),
          brms::prior(normal(1, 2), nlpar = "omega"),
          brms::prior(normal(45, 10), nlpar = "theta"))
  expect_s3_class(bp, "brmsprior")

  pl <- frmtmb:::as_priorlist(bp)
  expect_s3_class(pl, "frmtmb_priorlist")
  # the rows arrive in brms's own order, which c() sorts by class
  spelled <- sort(vapply(unclass(pl), function(s) {
    paste0(s$dist$kind, s$dist$location, "|", s$class, "|", s$nlpar)
  }, ""))
  expect_identical(spelled,
                   sort(c("normal5000|b|ult", "normal1|b|omega",
                          "normal45|b|theta")))

  dd <- loss_data()
  map <- frm(loss_form(), data = dd, start = loss_start, prior = bp)
  own <- frm(loss_form(), data = dd, start = loss_start,
             prior = c(prior(normal(5000, 1000), nlpar = "ult"),
                       prior(normal(1, 2), nlpar = "omega"),
                       prior(normal(45, 10), nlpar = "theta")))
  expect_equal(fixef(map)$ult, fixef(own)$ult)
  # ... and through `priors =` as well, which is the same setting
  expect_equal(fixef(frm(loss_form(), data = dd, start = loss_start,
                         prior = bp))$ult, fixef(own)$ult)
})

test_that("brms prior rows frmtmb cannot mean are refused by name", {
  skip_if_not_installed("brms")
  # brms's bounds are strings in its frame; they arrive as numbers here
  pl <- frmtmb:::as_priorlist(brms::prior(normal(0, 1), class = "b",
                                          lb = 0))
  expect_identical(unclass(pl)[[1L]]$lb, 0)

  # a class whose word means different parameters in the two packages
  expect_error(frmtmb:::as_priorlist(brms::prior(normal(0, 1),
                                                 class = "theta")),
               "mixture proportion")
  # a distributional parameter's own class, whose density brms puts on
  # a different scale from frmtmb's nearest spelling
  expect_error(frmtmb:::as_priorlist(brms::prior(student_t(3, 0, 10),
                                                 class = "sigma")),
               "LINK scale")
  # a tag names a prior inside a Stan program
  expect_error(frmtmb:::as_priorlist(brms::prior(normal(0, 1),
                                                 class = "b",
                                                 tag = "mytag")),
               "Drop the tag")
  # a density frmtmb does not carry says which one it was
  expect_error(frmtmb:::as_priorlist(brms::prior(gamma(0.01, 0.01),
                                                 class = "sd")),
               "Unsupported prior distribution")

  # slot-listing rows carry no density and are dropped, not refused
  dd <- data.frame(y = stats::rnorm(40), x = stats::rnorm(40))
  gp <- brms::get_prior(y ~ x, data = dd, family = stats::gaussian())
  expect_message(pl2 <- frmtmb:::as_priorlist(gp), "as its own defaults")
  expect_true(is.null(pl2) || inherits(pl2, "frmtmb_priorlist"))
})

# ---- nothing that already worked changed ------------------------------

test_that("the shipped set_prior spellings still resolve as before", {
  set.seed(404)
  dd <- data.frame(x = stats::rnorm(120), z = stats::rnorm(120),
                   g = factor(rep(1:12, 10)))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x +
                         stats::rnorm(12, 0, 0.7)[dd$g], 1)
  fit <- frm(bf(y ~ x + z + (1 | g), sigma ~ x) + gaussian(), data = dd)

  # class b still EXCLUDES the intercept where there is no nlpar, which
  # is the centered parameterization frmtmb has always used
  ri <- frmtmb:::resolve_prior_input(fit,
    set_prior("normal(0, 5)", class = "b"))
  nms <- vapply(ri$entries, function(e) {
    names(fit$frame$par_template$beta)[e$idx]
  }, "")
  expect_setequal(nms, c("x", "z"))

  # dpar, group, coef and the bounds all address what they used to
  expect_identical(frmtmb:::resolve_prior_input(fit,
    set_prior("normal(0, 2)", class = "b",
              dpar = "sigma"))$entries[[1L]]$comp, "betad")
  expect_length(frmtmb:::resolve_prior_input(fit,
    set_prior("exponential(1)", class = "sd", group = "g"))$entries, 1L)
  ri_b <- frmtmb:::resolve_prior_input(fit,
    set_prior("normal(0, 5)", class = "b") +
      set_prior("", class = "b", coef = "x", lb = 0))
  expect_identical(unname(ri_b$lower["x"]), 0)

  # the new fields are present on every spec, so nothing downstream
  # has to guess whether an old object carries them
  s <- unclass(set_prior("normal(0, 1)", class = "b"))[[1L]]
  expect_true(all(c("resp", "nlpar") %in% names(s)))
  expect_identical(s$nlpar, "")
})

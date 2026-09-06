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
    # this model gives sigma no predictor, so sigma's own class is the
    # spelling that applies to it, here as in brms
    prior(exponential(1), class = "sigma")
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
  # and the structures frmtmb holds somewhere else entirely, each named
  # with the class that reaches the same parameters
  expect_error(frmtmb:::as_priorlist(brms::prior(student_t(3, 0, 1),
                                                 class = "sds")),
               "class = \"sd\" with group", fixed = TRUE)
  expect_error(frmtmb:::as_priorlist(brms::prior(student_t(3, 0, 1),
                                                 class = "sdgp")),
               "class = \"sd\" with group", fixed = TRUE)
  expect_error(frmtmb:::as_priorlist(brms::prior(normal(0, 1),
                                                 class = "lscale")),
               "class = \"theta\"", fixed = TRUE)
  expect_error(frmtmb:::as_priorlist(brms::prior(student_t(3, 0, 1),
                                                 class = "sdcar")),
               "class = \"sd\" with group", fixed = TRUE)
  expect_error(frmtmb:::as_priorlist(brms::prior(beta(1, 1),
                                                 class = "car")),
               "class = \"theta\"", fixed = TRUE)
  expect_error(frmtmb:::as_priorlist(brms::prior(dirichlet(1),
                                                 class = "simo")),
               "simo")
  # a tag names a prior inside a Stan program
  expect_error(frmtmb:::as_priorlist(brms::prior(normal(0, 1),
                                                 class = "b",
                                                 tag = "mytag")),
               "Drop the tag")
  # a density frmtmb does not carry says which one it was
  expect_error(frmtmb:::as_priorlist(brms::prior(uniform(0, 10),
                                                 class = "sd")),
               "Unsupported prior distribution")
  # brms's shrinkage priors reach that same message rather than the
  # generic parse failure. R2D2() needs both halves of the fix: an
  # upper-case name and an empty argument list
  for (p in c("R2D2()", "horseshoe(1)", "lasso(1)")) {
    expect_error(frmtmb:::parse_prior_dist(p),
                 "Unsupported prior distribution")
  }
})

test_that("a coef frmtmb cannot honor is refused, not applied wider", {
  skip_if_not_installed("brms")
  # brms narrows an sd row to one coefficient of a block and writes
  # exponential_lpdf(sd_1[2] | 1), keeping its default on sd_1[1].
  # frmtmb's class "sd" addresses a BLOCK and never reads `coef`, so
  # honoring the row without it would put the density on every standard
  # deviation of the block. Measured before this refusal: coef =
  # "Intercept", coef = "x" and no coef gave one bit-identical
  # objective. Silently widening a prior is the failure D1 exists to
  # remove, so the row is refused and the whole-block spelling named.
  expect_error(
    frmtmb:::as_priorlist(brms::prior(exponential(1), class = "sd",
                                      group = "g", coef = "x")),
    "addresses a whole random-effect BLOCK", fixed = TRUE)
  expect_error(
    frmtmb:::as_priorlist(brms::prior(lkj(2), class = "cor",
                                      group = "g", coef = "x")),
    "addresses a whole correlation matrix", fixed = TRUE)
  # the whole-block row, which is the spelling the message names, still
  # translates
  expect_s3_class(
    frmtmb:::as_priorlist(brms::prior(exponential(1), class = "sd",
                                      group = "g")),
    "frmtmb_priorlist")
  # and a coef frmtmb DOES honor is untouched
  expect_identical(
    unclass(frmtmb:::as_priorlist(brms::prior(normal(0, 1), class = "b",
                                              coef = "x")))[[1L]]$coef,
    "x")
})

test_that("coef and group narrow the classes that read them", {
  # the other half of the same question, and the reason the refusal
  # above is narrow: `coef` on class "b" and `group` on "sd"/"cor" DO
  # bite. A design with two slopes and two correlated blocks is what
  # makes that visible; with one of each every spelling picks the same
  # parameters and the test would pass while proving nothing.
  set.seed(77)
  n <- 300
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(1:20, 15)), h = factor(rep(1:15, 20)))
  d$y <- stats::rnorm(n, 1 + 0.5 * d$x - 0.3 * d$z +
                        stats::rnorm(20, 0, 0.7)[d$g], 1)
  fit <- frm(bf(y ~ x + z + (x | g) + (z | h)) + gaussian(), data = d)
  idx <- function(pl) {
    e <- frmtmb:::resolve_prior_input(fit, pl)$entries
    sort(unlist(lapply(e, function(z) paste0(z$comp, z$idx))))
  }
  b_all <- idx(set_prior("normal(0, 0.1)", class = "b"))
  expect_length(b_all, 2L)
  expect_identical(idx(set_prior("normal(0, 0.1)", class = "b",
                                 coef = "x")), b_all[1])
  expect_identical(idx(set_prior("normal(0, 0.1)", class = "b",
                                 coef = "z")), b_all[2])

  cor_all <- idx(set_prior("lkj(6)", class = "cor"))
  expect_length(cor_all, 2L)
  expect_length(idx(set_prior("lkj(6)", class = "cor", group = "g")), 1L)
  expect_length(idx(set_prior("lkj(6)", class = "cor", group = "h")), 1L)
  expect_false(identical(idx(set_prior("lkj(6)", class = "cor",
                                       group = "g")),
                         idx(set_prior("lkj(6)", class = "cor",
                                       group = "h"))))

  expect_length(idx(set_prior("exponential(1)", class = "sd")), 4L)
  expect_length(idx(set_prior("exponential(1)", class = "sd",
                              group = "g")), 2L)
})

test_that("every refused row of a table is named in one message", {
  skip_if_not_installed("brms")
  # A table is edited as a whole, so stopping at the first bad row costs
  # one round trip per bad row. y ~ gp(x) carries an `lscale` row AND an
  # `sdgp` row: naming only the first made it two edit rounds.
  set.seed(3)
  dd <- data.frame(x = stats::rnorm(60))
  dd$y <- stats::rnorm(60, dd$x)
  gp <- brms::get_prior(brms::bf(y ~ gp(x)), data = dd,
                        family = stats::gaussian())
  msg <- tryCatch(frmtmb:::as_priorlist(gp), error = conditionMessage)
  expect_match(msg, "2 rows", fixed = TRUE)
  expect_match(msg, "lscale", fixed = TRUE)
  expect_match(msg, "sdgp", fixed = TRUE)
  # and one bad row still reads as one row
  one <- tryCatch(frmtmb:::as_priorlist(
    brms::prior(student_t(3, 0, 1), class = "sds")),
    error = conditionMessage)
  expect_match(one, "1 row with", fixed = TRUE)
})

test_that("inv_gamma and beta are the densities brms's defaults need", {
  skip_if_not_installed("brms")
  # brms's default on a shape parameter is inv_gamma and on a zi/hu is
  # beta, so before these two arms every negbinomial, zero-inflated and
  # hurdle table stopped for want of a density rather than for want of
  # a parameter frmtmb holds.
  expect_identical(frmtmb:::parse_prior_dist("inv_gamma(0.4, 0.3)"),
                   prior_inv_gamma(0.4, 0.3))
  expect_identical(frmtmb:::parse_prior_dist("beta(1, 1)"),
                   prior_beta(1, 1))

  # the densities and their AD gradients, against R's own
  ad_grad <- function(dist, x) {
    tp <- RTMB::MakeTape(function(p) frmtmb:::prior_base_logdens(p[1],
                                                                 dist), x)
    c(as.numeric(tp(x)), as.numeric(tp$jacobian(x)))
  }
  for (x in c(0.15, 0.8, 3.2)) {
    d <- prior_inv_gamma(0.4, 0.3)
    got <- ad_grad(d, x)
    # dinvgamma equivalent: the gamma density of 1/x, times the
    # Jacobian of that reciprocal
    expect_equal(got[1],
                 stats::dgamma(1 / x, shape = 0.4, rate = 0.3,
                               log = TRUE) - 2 * log(x),
                 tolerance = 1e-15)
    expect_equal(got[2], -(0.4 + 1) / x + 0.3 / x^2, tolerance = 1e-15)
  }
  for (x in c(0.02, 0.5, 0.97)) {
    got <- ad_grad(prior_beta(2, 3), x)
    expect_equal(got[1], stats::dbeta(x, 2, 3, log = TRUE),
                 tolerance = 1e-15)
    expect_equal(got[2], 1 / x - 2 / (1 - x), tolerance = 1e-15)
  }
})

test_that("a zero-inflated table translates onto zi itself", {
  skip_if_not_installed("brms")
  # The row this removes from the stop set, end to end. brms declares
  # real<lower=0,upper=1> zi and puts beta(1, 1) on it; frmtmb holds zi
  # as a logit-scale intercept, so the density lands on zi through the
  # logit inverse link with that map's log-Jacobian, and brms's
  # lb = 0/ub = 1 become no constraint on the logit scale. That is the
  # non-log branch of the bound transform, exercised in both directions.
  set.seed(12)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  d$zi <- stats::rpois(n, exp(1 + 0.4 * d$x)) *
    stats::rbinom(n, 1, 0.75)
  gp <- brms::get_prior(brms::bf(zi ~ x), data = d,
                        family = brms::zero_inflated_poisson())
  pl <- frmtmb:::as_priorlist(gp)
  fit <- frm(bf(zi ~ x) + zero_inflated_poisson(), data = d, prior = pl)
  ri <- frmtmb:::resolve_prior_input(fit, pl)
  zi_e <- Filter(function(e) identical(e$comp, "betad"), ri$entries)
  expect_length(zi_e, 1L)
  expect_identical(zi_e[[1L]]$scale, "natural")
  expect_identical(zi_e[[1L]]$link$name, "logit")
  expect_identical(zi_e[[1L]]$dist$kind, "beta")
  expect_identical(unname(ri$lower[["zi_(Intercept)"]]), -Inf)
  expect_identical(unname(ri$upper[["zi_(Intercept)"]]), Inf)
})

test_that("a distributional class means the same in both spellings", {
  skip_if_not_installed("brms")
  # BEHAVIOR CHANGE. brms spells a prior on sigma itself class =
  # "sigma", and frmtmb's own set_prior() now takes that word with that
  # meaning: a density on the parameter, through its inverse link with
  # that map's log-Jacobian. The two routes are one code path, so this
  # compares them SPEC to SPEC rather than restating either.
  from_brms <- unclass(frmtmb:::as_priorlist(
    brms::prior(student_t(3, 0, 10), class = "sigma")))[[1L]]
  own <- unclass(set_prior("student_t(3, 0, 10)",
                           class = "sigma"))[[1L]]
  expect_identical(own, from_brms)
  # stored as the slot the resolver assigns to, and marked `natural`
  expect_identical(own$class, "Intercept")
  expect_identical(own$dpar, "sigma")
  expect_true(isTRUE(own$natural))
  # and printed as the word it was WRITTEN with, so what comes out can
  # be pasted back in
  expect_match(paste(utils::capture.output(print(set_prior(
    "student_t(3, 0, 10)", class = "sigma"))), collapse = ""),
    "class=sigma scale=natural", fixed = TRUE)

  # the LINK-scale spelling is still available and still means the log
  # scale; it carries no `natural` field at all, because a field
  # written as FALSE would itself be a change and frmtmb.sample reads
  # its absence
  lnk <- unclass(set_prior("student_t(3, 0, 10)", class = "Intercept",
                           dpar = "sigma"))[[1L]]
  expect_null(lnk$natural)

  # gamma is one of brms's dispersion defaults and now parses, so a
  # shape/phi/nu/kappa row translates rather than stopping at the parser
  gp <- frmtmb:::as_priorlist(brms::prior(gamma(0.01, 0.01),
                                          class = "phi"))
  expect_identical(unclass(gp)[[1L]]$dist$kind, "gamma")
  expect_identical(unclass(gp)[[1L]]$dpar, "phi")
  expect_identical(unclass(gp)[[1L]],
                   unclass(set_prior("gamma(0.01, 0.01)",
                                     class = "phi"))[[1L]])
})

test_that("each dpar spelling is refused on the model brms refuses it on", {
  # BEHAVIOR CHANGE, and the whole point of the flip. The two brms
  # spellings for a distributional parameter are mutually exclusive BY
  # MODEL SHAPE: measured off brms::make_stancode(), `class = "sigma"`
  # is accepted on `y ~ x` and refused on `bf(y ~ x, sigma ~ 1)`, and
  # `class = "Intercept", dpar = "sigma"` the other way round. frmtmb
  # now says the same, by name, in both directions.
  set.seed(91)
  dd <- data.frame(x = stats::rnorm(120))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x, exp(0.2 + 0.1 * dd$x))
  plain <- frm(bf(y ~ x) + gaussian(), data = dd)
  one <- frm(bf(y ~ x, sigma ~ 1) + gaussian(), data = dd)
  pred <- frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dd)
  nat <- set_prior("normal(0, 0.5)", class = "sigma")
  lnk <- set_prior("normal(0, 0.5)", class = "Intercept", dpar = "sigma")

  expect_length(frmtmb:::resolve_prior_input(plain, nat)$entries, 1L)
  expect_error(frmtmb:::resolve_prior_input(plain, lnk),
               "brms accepts only where sigma has one", fixed = TRUE)
  # `sigma ~ 1` is a PREDICTOR, which is brms's own reading of it: the
  # design is one intercept column either way, so the difference is who
  # wrote the formula, not what it contains
  expect_length(frmtmb:::resolve_prior_input(one, lnk)$entries, 1L)
  expect_error(frmtmb:::resolve_prior_input(one, nat),
               "brms accepts only where sigma has no predictor",
               fixed = TRUE)
  expect_length(frmtmb:::resolve_prior_input(pred, lnk)$entries, 1L)
  expect_error(frmtmb:::resolve_prior_input(pred, nat),
               "class = \"Intercept\", dpar = \"sigma\"", fixed = TRUE)

  # class "b" on a dpar with no slopes says which spelling has one
  expect_error(frmtmb:::resolve_prior_input(plain,
    set_prior("normal(0, 1)", class = "b", dpar = "sigma")),
    "addresses sigma's slopes", fixed = TRUE)
  # a class that names no dpar of this model lists the ones it has
  expect_error(frmtmb:::resolve_prior_input(plain,
    set_prior("normal(0, 1)", class = "shape")),
    "It has mu, sigma", fixed = TRUE)

  # and the bound travels with the placement it was written on: brms's
  # lb = 0 on a log-linked sigma is log(0) = -Inf, not a floor of 1
  ri <- frmtmb:::resolve_prior_input(plain,
    set_prior("", class = "sigma", lb = 0, ub = 3))
  expect_identical(unname(ri$lower[["sigma_(Intercept)"]]), -Inf)
  expect_equal(unname(ri$upper[["sigma_(Intercept)"]]), log(3))

  # the third spelling brms decides by shape, and the one the first
  # round left open. `sigma ~ 1` HAS a predictor, so it takes the
  # Intercept spelling, but it has no population-level slopes, so
  # class "b" there addressed an empty set: the penalty was
  # bit-identical to no prior at all, silently
  b_dp <- set_prior("normal(0, 1e-6)", class = "b", dpar = "sigma")
  expect_error(frmtmb:::resolve_prior_input(one, b_dp),
               "predictor is an intercept only, so it has none",
               fixed = TRUE)
  expect_error(frmtmb:::resolve_prior_input(one, b_dp),
               "class = \"Intercept\", dpar = \"sigma\"", fixed = TRUE)
  # and it is still honored where the predictor really has slopes
  expect_length(frmtmb:::resolve_prior_input(pred, b_dp)$entries, 1L)
})

test_that("brms refuses the same b/dpar row on an intercept-only dpar", {
  skip_if_not_installed("brms")
  # the refusal above is brms's rule, not a rule of frmtmb's own: brms
  # answers "do not correspond to any model parameter: b_sigma"
  set.seed(91)
  dd <- data.frame(x = stats::rnorm(60))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x, 1)
  msg <- tryCatch(brms::validate_prior(
    brms::prior(normal(0, 1), class = "b", dpar = "sigma"),
    brms::bf(y ~ x, sigma ~ 1), data = dd, family = gaussian()),
    error = conditionMessage)
  expect_true(is.character(msg))
  expect_match(msg, "do not correspond to any model parameter",
               fixed = TRUE)
  # and brms accepts it once sigma has slopes, as frmtmb does
  ok <- brms::validate_prior(
    brms::prior(normal(0, 1), class = "b", dpar = "sigma"),
    brms::bf(y ~ x, sigma ~ x), data = dd, family = gaussian())
  expect_s3_class(ok, "brmsprior")
})

test_that("a distributional class is refused by the same name on both routes", {
  # A class frmtmb keeps somewhere else used to be refused only when it
  # arrived on a brms frame; written by hand it fell through to a dpar
  # target that does not exist. Both paths now give the same sentence.
  expect_error(set_prior("student_t(3, 0, 1)", class = "sds"),
               "class = \"sd\" with group", fixed = TRUE)
  # the density gate runs first here, so the class is asked about with
  # a density set_prior() knows; brms's own dirichlet row is refused on
  # the translated route by the density instead, and either way the
  # call stops
  expect_error(set_prior("normal(0, 1)", class = "simo"),
               "Dirichlet on a mo() simplex", fixed = TRUE)
  expect_error(set_prior("normal(0, 1)", class = "theta2"),
               "mixture proportion", fixed = TRUE)
  expect_error(set_prior("normal(0, 1)", class = "car"),
               "spatial dependence parameter", fixed = TRUE)
  # a name that is not a class and cannot be a parameter says so at the
  # boundary, where there is still no model to check it against
  expect_error(set_prior("normal(0, 1)", class = "9x"),
               "nor the name of a distributional parameter",
               fixed = TRUE)
  # the two arguments that would name the parameter twice
  expect_error(set_prior("normal(0, 1)", class = "sigma",
                         dpar = "shape"), "names a second one",
               fixed = TRUE)
  expect_error(set_prior("normal(0, 1)", class = "sigma", nlpar = "a"),
               "nonlinear parameter's coefficients are class = \"b\"",
               fixed = TRUE)
})

test_that("set_prior() refuses an unhonored coef, as the brms route does", {
  # The review's residual R1. The translated route refused a `coef` on
  # class "sd" (frmtmb resolves that class per BLOCK, so the row would
  # have applied to every standard deviation of the block); the native
  # spelling accepted it silently and widened. Both refuse now.
  expect_error(set_prior("exponential(2)", class = "sd", group = "g",
                         coef = "x"),
               "addresses a whole random-effect BLOCK", fixed = TRUE)
  expect_error(set_prior("lkj(2)", class = "cor", coef = "x"),
               "names nothing it can narrow to", fixed = TRUE)
  # and `coef` on the classes that read it is untouched
  expect_silent(set_prior("normal(0, 1)", class = "b", coef = "x"))
  expect_silent(set_prior("normal(0, 1)", class = "Intercept",
                          coef = "Intercept"))

  # F7 of the punch round: the same rule, for the class family this
  # lane introduced. A distributional class is ONE parameter, so it
  # reads neither `coef` nor `group`; both used to reach the spec and
  # vanish while the prior applied anyway. brms refuses the same rows,
  # naming a parameter that does not exist (`sigma_x`)
  expect_error(set_prior("normal(0, 0.3)", class = "sigma", coef = "x"),
               "names nothing it can narrow to", fixed = TRUE)
  expect_error(set_prior("normal(0, 0.3)", class = "sigma",
                         group = "g"),
               "names nothing it can narrow to", fixed = TRUE)
  expect_error(set_prior("normal(0, 0.3)", class = "shape",
                         coef = "x"), "class = \"shape\"", fixed = TRUE)
  # `resp` is the one narrowing a distributional class does take
  expect_silent(set_prior("normal(0, 0.3)", class = "sigma",
                          resp = "y1"))
})

test_that("a brms table applies what its prior strings say", {
  skip_if_not_installed("brms")
  # BEHAVIOR CHANGE. Rows used to be dropped by the `source` column with
  # a message. `source` records who BUILT a row, not who wrote the
  # density in it, so the drop lost a prior the user had edited in place
  # and reported it as one brms had filled in. A row now applies
  # whatever its `prior` string says, which is brms's own rule.
  dd <- data.frame(y = stats::rnorm(40), x = stats::rnorm(40))
  gp <- brms::get_prior(y ~ x, data = dd, family = stats::gaussian())
  expect_silent(pl2 <- frmtmb:::as_priorlist(gp))
  expect_s3_class(pl2, "frmtmb_priorlist")
  # exactly the live rows, and nothing for a slot brms left flat
  live <- sum(nzchar(as.data.frame(gp)$prior))
  expect_identical(length(unclass(pl2)), as.integer(live))

  # a slot-listing row with an empty prior is not a prior to apply, so a
  # table with nothing live in it carries no density at all. What can
  # survive is brms echoing a parameter's own declared bound onto the
  # row, which frmtmb turns into the internal box that bound implies
  gp$prior <- ""
  blank <- frmtmb:::as_priorlist(gp)
  expect_true(is.null(blank) ||
                all(vapply(unclass(blank), function(s) is.null(s$dist),
                           TRUE)))
})

test_that("a flat row of a class frmtmb cannot name is not a refusal", {
  skip_if_not_installed("brms")
  # A smooth's wiggliness SD is class "sds" in brms and class "sd" with
  # a group here, so the LIVE row is refused by name with the spelling
  # that reaches the same parameters. Honoring the table exposes that
  # refusal where the row used to be dropped in silence, which is the
  # cost of reading the prior string instead of the `source` column.
  dd <- data.frame(y = stats::rnorm(60), x = stats::rnorm(60))
  gp <- brms::get_prior(y ~ s(x), data = dd, family = stats::gaussian())
  expect_error(frmtmb:::as_priorlist(gp), "class = \"sd\" with group",
               fixed = TRUE)

  g <- as.data.frame(gp)
  # brms puts a bound only on the row that carries the density, so its
  # own flat rows never reach the gate. A hand-written flat row with a
  # bound does, and it is a slot listing rather than a refusal: the
  # bound restates a declaration frmtmb's parameterization already makes
  g$prior[g$class == "sds"] <- ""
  g$lb[g$class == "sds"] <- "0"
  expect_silent(pl <- frmtmb:::as_priorlist(
    structure(g, class = c("brmsprior", "data.frame"))))
  expect_false(any(vapply(unclass(pl),
                          function(s) identical(s$class, "sds"), TRUE)))
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

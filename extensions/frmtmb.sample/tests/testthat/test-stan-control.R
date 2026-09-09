# frm_sample(control =): the SAMPLER's control list, spelled and meant
# as brms spells and means it, and frm_sample(fit_control =), which is
# what frm() calls `control`.
#
# Earlier releases had the two names the other way round, so a
# brms call that raised adapt_delta bound its list to the fit-time
# argument, which the fit route never reads: the call sampled with
# rstan's defaults and said nothing. The first block below is that
# silence turned into a refusal.

skip_on_cran()

sc_data <- function(seed = 9, n = 60L, ng = 6L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x +
                         stats::rnorm(ng, 0, 0.5)[dd$g], 1)
  dd
}

sc_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- frm(bf(y ~ x + (1 | g)) + gaussian(), data = sc_data())
    }
    cache
  }
})

## ---- the two lists are told apart by name ----------------------------

test_that("a frmtmb_control() list in control = is refused by name", {
  skip_sampler()
  fit <- sc_fit()
  expect_error(frm_sample(fit, control = frmtmb_control()),
               "fit_control")
  # the refusal names the fields it keyed on, so it reports what it saw
  expect_error(frm_sample(fit, control = frmtmb_control()),
               "grad_tol")
  # one field is enough: a partially written list is the same mistake
  expect_error(frm_sample(fit, control = list(profile = TRUE)),
               "profile")
  expect_error(frm_sample(bf(y ~ x) + gaussian(), data = sc_data(),
                          control = list(optimizer = "nlminb")),
               "fit_control")
})

test_that("the two vocabularies are disjoint, which is what makes the
          refusal decidable", {
  # the check tells a frmtmb_control() list from a sampler control list
  # by field name alone, so a shared name would make some lists
  # ambiguous. Asserted rather than assumed, because both vocabularies
  # belong to other packages and can grow.
  stan_names <- frmtmb.sample:::stan_control_names
  expect_length(intersect(stan_names, names(frmtmb_control())), 0L)
})

test_that("an option rstan does not have is refused by name", {
  skip_sampler()
  fit <- sc_fit()
  # rstan answers an unknown option by declining to sample and
  # returning an EMPTY fit, which would otherwise land on the
  # "returned no draws" refusal and be reported as a solver failure
  expect_error(frm_sample(fit, control = list(adapt_deltaa = 0.99)),
               "adapt_deltaa")
  expect_error(frm_sample(fit, control = list(adapt_delta = 0.9,
                                              nonsense = 1)),
               "nonsense")
  expect_error(frm_sample(fit, control = "adapt_delta = 0.99"),
               "NAMED list")
})

test_that("an empty control list is a no-op and an unnamed element is
           refused for being unnamed", {
  skip_sampler()
  withr::local_options(mc.cores = 1)
  # building the list programmatically and adding nothing to it is the
  # ordinary way to reach list(), and it means "rstan's defaults"
  ds <- suppressWarnings(suppressMessages(
    frm_sample(sc_fit(), chains = 1, iter = 40, refresh = 0, seed = 5,
               control = list())))
  expect_s3_class(ds, "frmtmb_draws")
  # a partially named list is a different mistake from a non-list, and
  # naming its class would not say which element is wrong
  expect_error(frm_sample(sc_fit(), control = list(0.99)), "no name")
  expect_error(
    frm_sample(sc_fit(), control = list(adapt_delta = 0.99, 12)),
    "1 of 2")
})

test_that("an abbreviation of control is refused, because it would
           reach rstan's own control through the dots", {
  skip_sampler()
  fit <- sc_fit()
  # `control` follows `...` in the formals so only the exact spelling
  # binds; the abbreviation lands in `...`, and rstan::sampling() takes
  # `control` BEFORE its own `...`, so partial matching binds it there.
  # Both of these used to reach the sampler unchecked, and the second
  # came back as the "no draws ... external solver" mis-report.
  expect_error(frm_sample(fit, contro = list(adapt_delta = 0.97)),
               "abbreviation of `control`")
  expect_error(frm_sample(fit, cont = frmtmb_control()),
               "abbreviation of `control`")
  # two abbreviations in one call: pmatch() reports only the first
  # unless duplicates.ok, and the second would then slip through
  err <- expect_error(frm_sample(fit, contro = list(), contr = list()))
  expect_match(conditionMessage(err), "contro")
  expect_match(conditionMessage(err), "contr`")
})

## ---- the list reaches rstan -------------------------------------------

test_that("control = travels to rstan and is recorded there", {
  skip_sampler()
  withr::local_options(mc.cores = 1)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(sc_fit(), chains = 1, iter = 60, refresh = 0, seed = 5,
               control = list(adapt_delta = 0.97, max_treedepth = 12))))
  rec <- ds$stanfit@stan_args[[1L]]$control
  expect_equal(rec$adapt_delta, 0.97)
  expect_equal(rec$max_treedepth, 12)
  # and the default leaves rstan's own defaults alone rather than
  # sending an empty list
  ds0 <- suppressWarnings(suppressMessages(
    frm_sample(sc_fit(), chains = 1, iter = 60, refresh = 0, seed = 5)))
  expect_null(ds0$stanfit@stan_args[[1L]]$control)
})

## ---- an empty stanfit is refused rather than returned ---------------

test_that("a run that produces no draws is refused rather than
           returned, and says it is the RUN and not the model", {
  skip_sampler()
  withr::local_options(mc.cores = 1)
  fit <- sc_fit()
  # An option rstan does not recognize makes it decline to sample and
  # hand back the shell of a stanfit with no error at all. It is the
  # cheapest deterministic way to reach that state and it involves no
  # solver, so it exercises the guard rather than the ODE story.
  # frm_sample() screens the option out before rstan sees it, which is
  # why this goes through as_tmbstan(), the escape hatch that does not.
  err <- expect_error(suppressWarnings(
    as_tmbstan(fit, chains = 1, iter = 40, refresh = 0,
               control = list(nonsense_option = 1))))
  expect_match(conditionMessage(err), "as_tmbstan[(][)]")
  expect_match(conditionMessage(err), "returned no draws")
  # the two ways to get nothing out of this package want different
  # fixes, so the message says which one this is
  expect_match(conditionMessage(err), "failure of THIS RUN")
  # and a run that DOES produce draws still gets its stanfit back
  sf <- suppressWarnings(suppressMessages(
    as_tmbstan(fit, chains = 1, iter = 40, refresh = 0)))
  expect_s4_class(sf, "stanfit")
  expect_gt(length(sf@sim[["samples"]]), 0L)
})

test_that("check_stan_draws() decides on the object, so an rstan that
           stops declining re-points one test instead of deleting the
           coverage", {
  skip_if_not_installed("rstan")
  chk <- frmtmb.sample:::check_stan_draws
  # the state the block above provokes through rstan is constructible
  # by hand, and this pins the guard to the SHAPE rather than to rstan
  # continuing to answer an unknown option by declining
  expect_error(chk(new("stanfit", sim = list()), "unit()"),
               "returned no draws")
  expect_error(chk(new("stanfit", sim = list(samples = list())),
                   "unit()"),
               "returned no draws")
  expect_silent(chk(new("stanfit", sim = list(samples = list(1))),
                    "unit()"))
  # the caller's name is what tells the two doors apart in a message
  expect_match(
    tryCatch(chk(new("stanfit", sim = list()), "as_tmbstan()"),
             error = conditionMessage),
    "as_tmbstan[(][)]")
})

## ---- fit_control still reaches frm() on the formula route ------------

test_that("fit_control = reaches frm() on the formula route", {
  skip_sampler()
  dd <- sc_data()
  dd$one <- factor(rep("a", nrow(dd)))
  form <- bf(y ~ x + (1 | one)) + gaussian()
  # check_nlev_1 is a frmtmb_control() field whose "stop" setting turns
  # a warning into an error during assembly, so the assembly's own
  # behavior reports whether the list arrived
  expect_error(
    frm_sample(form, data = dd, chains = 1, iter = 20, refresh = 0,
               fit_control = frmtmb_control(check_nlev_1 = "stop")),
    "single level")
  # the permissive settings travel too, so the argument is doing more
  # than turning one error on. A short chain warns about its own
  # mixing, which is not what is being read here, so the warnings are
  # collected and searched rather than matched one at a time.
  warned <- function(act) {
    w <- character(0)
    withCallingHandlers(
      suppressMessages(
        frm_sample(form, data = dd, chains = 1, iter = 20, refresh = 0,
                   fit_control = frmtmb_control(check_nlev_1 = act))),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      })
    any(grepl("single level", w))
  }
  expect_true(warned("warning"))
  expect_false(warned("ignore"))
})

test_that("an assembly argument on the FIT route is refused rather
           than accepted and discarded", {
  skip_sampler()
  fit <- sc_fit()
  dd <- sc_data()
  # fit_control is the one this release created: it is the name the
  # `control` refusal points a confused user at, so accepting a sampler
  # list under it and sampling with the defaults would put the silence
  # back under a new name
  expect_error(frm_sample(fit, fit_control = list(adapt_delta = 0.99)),
               "fit_control =")
  expect_error(frm_sample(fit, fit_control = frmtmb_control()),
               "formula interface")
  # and the message points at the argument that DOES tighten the
  # sampler, since that is what the caller wanted
  expect_error(frm_sample(fit, fit_control = list(adapt_delta = 0.99)),
               "adapt_delta")
  # the other four were accepted and discarded before this release too
  expect_error(frm_sample(fit, start = list(beta = c(999, 999))),
               "start =")
  expect_error(frm_sample(fit, data2 = list(z = 1)), "data2 =")
  expect_error(frm_sample(fit, na.action = stats::na.fail),
               "na.action =")
  expect_error(frm_sample(fit, REML = TRUE), "REML =")
  # several at once are named together rather than one per attempt
  err <- expect_error(frm_sample(fit, REML = TRUE, start = list()))
  expect_match(conditionMessage(err), "start =")
  expect_match(conditionMessage(err), "REML =")
  # passing the DEFAULT explicitly is still passing it: the argument
  # does nothing here either way, and missing() is what tells them apart
  expect_error(frm_sample(fit, REML = FALSE), "REML =")
  # none of this touches the formula route, where they do their work
  expect_s3_class(
    suppressWarnings(suppressMessages(
      frm_sample(bf(y ~ x) + gaussian(), data = dd, chains = 1,
                 iter = 40, refresh = 0, REML = FALSE,
                 fit_control = frmtmb_control()))),
    "frmtmb_draws")
})

## ---- what adapt_delta buys, measured by the run ----------------------

test_that("raising adapt_delta shortens the step and removes
           divergences on a centered funnel", {
  skip_sampler()
  skip_if_not(sampler_gates_on(), "chain-agreement gates are off")
  withr::local_options(mc.cores = 1)
  # Six groups of three observations and a group sd of 0.05 against a
  # residual sd of 1: the centered parameterization of this posterior is
  # Neal's funnel, which is the geometry adapt_delta exists for.
  # reparameterize = FALSE keeps it centered, and the default would
  # remove the funnel and with it the thing being measured. The
  # construction was chosen by measurement, not by eye: four candidates
  # were run over eight chain seeds each and every one moved the counts
  # the same way; this one has the widest margin (75 divergences at
  # 0.80 against 0 at 0.99). dev/stanctl-findings.md has the table.
  set.seed(2026)
  ng <- 6L
  dd <- data.frame(g = factor(rep(seq_len(ng), each = 3L)))
  dd$y <- stats::rnorm(nrow(dd), stats::rnorm(ng, 0, 0.05)[dd$g], 1)
  form <- bf(y ~ 1 + (1 | g)) + gaussian()
  fit <- frm(form, data = dd)

  run <- function(seed, ad) {
    ds <- suppressWarnings(suppressMessages(
      frm_sample(fit, chains = 1, iter = 1000, refresh = 0, seed = seed,
                 reparameterize = FALSE,
                 control = list(adapt_delta = ad))))
    np <- nuts_params(ds)
    list(div = sum(np$Value[np$Parameter == "divergent__"]),
         step = mean(np$Value[np$Parameter == "stepsize__"]))
  }
  seeds <- 11:18
  lo <- lapply(seeds, run, ad = 0.8)
  hi <- lapply(seeds, run, ad = 0.99)

  # Both assertions compare the two arms of the SAME runs, so neither
  # carries a threshold read off another machine. The step size is the
  # mechanism (adapt_delta is a target acceptance rate, and a higher
  # target is met with a shorter step) and holds per seed; the
  # divergence count is the consequence, is heavy-tailed across seeds,
  # and so is summed over the eight.
  for (i in seq_along(seeds)) {
    expect_lt(hi[[i]]$step, lo[[i]]$step)
  }
  n_lo <- sum(vapply(lo, function(z) z$div, 0))
  n_hi <- sum(vapply(hi, function(z) z$div, 0))
  # the low arm still diverging is a property of the CONSTRUCTION: if
  # this fails, the funnel has to be re-chosen before the comparison
  # below means anything
  expect_gt(n_lo, 0)
  expect_lt(n_hi, n_lo)
})

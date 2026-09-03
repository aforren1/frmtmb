# frm_sample() from a formula: assemble and tape the model, stop before
# the optimizer, and hand the objective to NUTS. Also the draws-side
# parameter-name convention, which this route is born with.

skip_on_cran()

sd_data <- function(seed = 9, n = 60L, ng = 6L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x +
                         stats::rnorm(ng, 0, 0.5)[dd$g], 1)
  dd
}

## ---- the unfitted object the formula route samples -------------------

test_that("dry_run = 'objective' stops before the optimizer", {
  dd <- sd_data()
  uf <- frm(bf(y ~ x + (1 | g)) + gaussian(), dd, dry_run = "objective")
  expect_s3_class(uf, "frmtmb_unfitted")
  expect_s3_class(uf, "frmtmb_fit")
  expect_null(uf$opt)
  # the tape is real and evaluable, which is all the sampler needs
  expect_true(is.finite(uf$obj$fn(uf$obj$par)))
  # and the frame is a normal frame, so the draws surface has its
  # structure
  expect_equal(stats::nobs(uf), nrow(dd))
  expect_length(uf$frame$re_blocks, 1L)

  # quadrature has no unfitted form: its tape is calibrated at an optimum
  dp <- dd
  dp$y <- stats::rpois(nrow(dp), exp(0.3 + 0.4 * dp$x))
  expect_error(frm(bf(y ~ x + (1 | g)) + poisson(), dp,
                   quadrature = TRUE, dry_run = "objective"),
               "no unfitted form")
})

test_that("methods needing an ML quantity refuse on an unfitted object", {
  dd <- sd_data()
  uf <- frm(bf(y ~ x + (1 | g)) + gaussian(), dd, dry_run = "objective")
  for (f in list(function() summary(uf), function() stats::vcov(uf),
                 function() stats::confint(uf), function() stats::logLik(uf),
                 function() stats::AIC(uf), function() fixef(uf),
                 function() ranef(uf), function() VarCorr(uf),
                 function() stats::predict(uf), function() stats::fitted(uf),
                 function() stats::residuals(uf),
                 function() stats::simulate(uf), function() print(uf))) {
    expect_error(f(), "needs a fitted model")
  }
  # what does NOT need one keeps working: the model description
  expect_equal(stats::nobs(uf), nrow(dd))
  expect_s3_class(stats::formula(uf), "formula")
  expect_equal(stats::family(uf)$family, "gaussian")
})

## ---- sampling ---------------------------------------------------------

test_that("the formula route samples the same posterior as the fit route", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  dd <- sd_data()
  form <- bf(y ~ x + (1 | g)) + gaussian()

  fit <- frm(form, data = dd)
  # both routes default to the same brms priors, so both sample the
  # SAME posterior; the fit route only starts the chains at the mode
  ds_fit <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 1000, refresh = 0, seed = 3)))
  ds_form <- suppressWarnings(suppressMessages(
    frm_sample(form, data = dd, chains = 1, iter = 1000,
               refresh = 0, seed = 3)))
  key <- function(pl) {
    vapply(unclass(pl), function(s) {
      paste0(s$class, ":", s$dpar, ":", s$dist$kind)
    }, "")
  }
  expect_equal(key(prior_summary(ds_fit)), key(prior_summary(ds_form)))

  a <- summary(ds_fit)[, "mean"]
  b <- summary(ds_form)[, "mean"]
  expect_setequal(names(a), names(b))
  # the sharp gate: both routes taped the same density, so the two
  # objectives agree pointwise; the draws comparison below then only
  # guards the sampler wiring
  p <- fit$obj$par
  expect_equal(ds_form$fit$obj$fn(p), fit$obj$fn(p), tolerance = 1e-8)
  expect_equal(ds_form$fit$obj$fn(p + 0.1), fit$obj$fn(p + 0.1),
               tolerance = 1e-8)
  # a seeded Stan run is not platform-deterministic, and the group-sd
  # components mix slowly with 6 groups (ESS in the tens on a 500-draw
  # chain): a 1.11-sd mean gap was observed on macOS CI, so the Monte
  # Carlo band is wide by design
  if (sampler_gates_on()) {
    expect_lt(max(abs(a - b[names(a)]) /
                    pmax(summary(ds_fit)[, "sd"], 1e-8)), 2)
  }

  # the whole draws surface runs off the formula-route object
  expect_s3_class(ds_form, "frmtmb_draws")
  expect_true(all(c("mean", "sd", "Rhat") %in% colnames(summary(ds_form))))
  expect_equal(nrow(fixef(ds_form)), 3L)
  expect_true(all(c("estimate", "lwr", "upr") %in% names(VarCorr(ds_form))))
  expect_named(ranef(ds_form), "1 | g")
  h <- hypothesis(ds_form, "x > 0")
  expect_s3_class(h, "frmtmb_hypothesis")
  expect_equal(dim(posterior_epred(ds_form, ndraws = 5)),
               c(5L, nrow(dd)))
  expect_equal(dim(posterior_predict(ds_form, ndraws = 5)),
               c(5L, nrow(dd)))
  expect_equal(dim(posterior_linpred(ds_form, ndraws = 5)),
               c(5L, nrow(dd)))

  # the embedded object still has no ML quantities, and says so
  expect_error(check_laplace(ds_form$fit), "needs a fitted model")
  expect_error(fixef(ds_form$fit), "needs a fitted model")
  expect_error(stats::logLik(ds_form$fit), "needs a fitted model")
})

## ---- default priors ---------------------------------------------------

# the priorlist the formula route would choose, as one string per spec,
# for comparison with brms's own default_prior() rows
def_strings <- function(form, data, family = NULL) {
  uf <- frm(form, data, family = family, dry_run = "objective")
  pl <- frmtmb:::default_priors_for(uf)
  vapply(unclass(pl), function(s) {
    d <- if (identical(s$dist$kind, "lkj")) {
      paste0("lkj(", s$dist$eta, ")")
    } else {
      paste0("student_t(", paste(c(s$dist$df, s$dist$location,
                                   s$dist$scale), collapse = ", "), ")")
    }
    paste0(d, "|", s$class, if (nzchar(s$dpar)) paste0(":", s$dpar))
  }, "")
}

sd_prior_data <- function() {
  set.seed(9)
  n <- 60L
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(1:6, length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x, 1)
  dd$y2 <- dd$y * 10 + 100
  dd$p <- stats::rpois(n, 3)
  dd$b <- stats::rbinom(n, 1L, 0.4)
  dd$w <- stats::rlnorm(n, 1.5, 0.7)
  # counts with MEDIAN 1 and zeros present: the case that separates
  # brms's per-element zero shift from shifting the whole vector, since
  # log(1) = 0 rounds to 0 and log(1.1) = 0.0953 rounds to 0.1
  repeat {
    dd$z <- stats::rpois(n, 1.2)
    if (stats::median(dd$z) == 1 && any(dd$z == 0)) break
  }
  dd$o <- factor(cut(0.8 * dd$x + stats::rnorm(n),
                     c(-Inf, -0.5, 0.5, Inf), labels = c("a", "b", "c")),
                 ordered = TRUE)
  dd
}

test_that("the default priors are brms 2.23's, read off brms itself", {
  skip_if_not_installed("brms")
  dd <- sd_prior_data()

  # brms's own answer, reduced to the class-level rows that carry a
  # distribution. Dropped: class b (a documented non-match: brms leaves
  # slopes flat and so do we), the shape/phi/nu dispersion classes
  # (gamma and inverse-gamma defaults set_prior() cannot carry), and the
  # Intercept row of an ordinal model (its thresholds). Each of those is
  # announced by the disclosure message and pinned by its own test
  # below. Class `cor` is NOT dropped any more: since 0.39 the lkj(1)
  # default matches brms's own, and the correlated case below is the row
  # that checks it. brms lists one sd row per dpar for a categorical
  # family where frmtmb has one shared block, so the rows are
  # deduplicated.
  brms_def <- function(f, family, drop = character(0)) {
    p <- brms::default_prior(f, data = dd, family = family)
    p <- p[nzchar(p$prior) & p$prior != "(flat)" & p$coef == "" &
             p$group == "" &
             !p$class %in% c("b", "shape", "phi", "nu", drop), ]
    cls <- ifelse(p$class == "sigma", "Intercept:sigma", p$class)
    sort(unique(paste0(p$prior, "|", cls)))
  }
  cases <- list(
    list(y ~ x + (1 | g), bf(y ~ x + (1 | g)) + gaussian(), gaussian()),
    # a CORRELATED block: brms adds lkj(1) on class cor, and so do we
    list(y ~ x + (1 + x | g), bf(y ~ x + (1 + x | g)) + gaussian(),
         gaussian()),
    list(y2 ~ x + (1 | g), bf(y2 ~ x + (1 | g)) + gaussian(), gaussian()),
    list(p ~ x + (1 | g), bf(p ~ x + (1 | g)) + poisson(), poisson()),
    list(b ~ x + (1 | g), bf(b ~ x + (1 | g)) + bernoulli(),
         brms::bernoulli()),
    # zeros present, median 1: the per-element shift case
    list(z ~ x + (1 | g), bf(z ~ x + (1 | g)) + poisson(), poisson()),
    list(z ~ x + (1 | g), bf(z ~ x + (1 | g)) + negbinomial(),
         brms::negbinomial()),
    # a log-SCALE family: mu's link is spelled identity, but the
    # response is logged inside the density and brms transforms it
    list(w ~ x + (1 | g), bf(w ~ x + (1 | g)) + lognormal(),
         brms::lognormal()),
    list(w ~ x + (1 | g), bf(w ~ x + (1 | g)) + Gamma(link = "log"),
         Gamma(link = "log")),
    # a factor response: no transform, so the integer codes never reach
    # median() or mad()
    list(o ~ x + (1 | g), bf(o ~ x + (1 | g)) + cumulative(),
         brms::cumulative())
  )
  for (cs in cases) {
    drop <- if (identical(deparse1(cs[[1L]][[2L]]), "o")) "Intercept"
    expect_equal(sort(unname(def_strings(cs[[2L]], dd))),
                 brms_def(cs[[1L]], cs[[3L]], drop %||% character(0)))
  }
})

test_that("the zero shift is per element and only under a log-like link", {
  dd <- sd_prior_data()
  # the reviewer's case made visible: median(z) is 1 and zeros are
  # present, so shifting the WHOLE vector would give log(1.1) = 0.1
  expect_equal(stats::median(dd$z), 1)
  expect_true(any(dd$z == 0))
  expect_equal(round(stats::median(log(dd$z + 0.1)), 1L), 0.1)

  scl <- function(form) {
    frmtmb:::default_prior_scale(frm(form, dd, dry_run = "objective"))
  }
  # log link: zeros only, so the location is log(median) = log(1) = 0
  expect_equal(scl(bf(z ~ x) + poisson())$location, 0)
  expect_equal(scl(bf(z ~ x) + negbinomial())$location, 0)
  # identity link never shifts, and never transforms
  s_id <- scl(bf(z ~ x) + gaussian())
  expect_equal(s_id$location, round(stats::median(dd$z), 1L))
  # a log-SCALE family transforms even though its mu link is identity
  expect_equal(scl(bf(w ~ x) + lognormal())$location,
               round(stats::median(log(dd$w)), 1L))
  expect_equal(scl(bf(w ~ x) + lognormal())$link, "log")
  # a factor response takes the untransformed fallback, so the category
  # codes never reach median() or mad()
  s_ord <- scl(bf(o ~ x) + cumulative())
  expect_equal(s_ord$location, 0)
  expect_equal(s_ord$scale, 2.5)
  expect_false(s_ord$centered)
})

test_that("brms's own conventions on slopes and on a modelled sigma hold", {
  dd <- sd_prior_data()
  # brms leaves slopes flat and so do we: no class "b" spec is produced
  pl <- frmtmb:::default_priors_for(
    frm(bf(y ~ x + (1 | g)) + gaussian(), dd, dry_run = "objective"))
  expect_false("b" %in% vapply(unclass(pl), `[[`, "", "class"))

  # a sigma with its own predictor takes the plain student_t(3, 0, 2.5)
  # on the LOG scale, which is brms's own distinction
  s2 <- unclass(frmtmb:::default_priors_for(
    frm(bf(y2 ~ x, sigma ~ x) + gaussian(), dd, dry_run = "objective")))
  sg <- Filter(function(s) identical(s$dpar, "sigma"), s2)[[1L]]
  expect_equal(sg$dist$scale, 2.5)
  expect_null(sg$natural)
})

test_that("every slot left flat is announced, never silently dropped", {
  dd <- sd_prior_data()
  notes <- function(form) {
    frmtmb:::default_prior_notes(frm(form, dd, dry_run = "objective"))
  }
  # an ordinal model's thresholds are not a design column, so
  # set_prior() cannot reach them; the priorlist really is empty here
  uf <- frm(bf(o ~ x) + cumulative(), dd, dry_run = "objective")
  expect_null(frmtmb:::default_priors_for(uf))
  expect_match(notes(bf(o ~ x) + cumulative()), "thresholds")
  # ... and it is still announced rather than passing in silence
  msg <- capture_messages(
    frmtmb:::sample_resolve_priors(uf, NULL))
  expect_match(paste(msg, collapse = ""), "thresholds")
  expect_match(paste(msg, collapse = ""), "see ?frm_sample", fixed = TRUE)
  # an ordinal model WITH random effects still gets its sd default
  uf2 <- frm(bf(o ~ x + (1 | g)) + cumulative(), dd,
             dry_run = "objective")
  cls <- vapply(unclass(frmtmb:::default_priors_for(uf2)), `[[`, "",
                "class")
  expect_equal(cls, "sd")
  expect_match(paste(capture_messages(
    frmtmb:::sample_resolve_priors(uf2, NULL)), collapse = ""),
    "thresholds")

  expect_match(notes(bf(z ~ x) + negbinomial()), "shape")
  # a correlated block is no longer a gap: it gets lkj(1), brms's own
  # default, so there is nothing left to announce
  expect_length(notes(bf(y ~ x + (1 + x | g)) + gaussian()), 0L)
  expect_length(notes(bf(y ~ x + (1 | g)) + gaussian()), 0L)
  # a structure with no LKJ density still is a gap, and is named by
  # name; test-lkj.R has that case, where a toep() block is built
})

test_that("the formula route discloses its defaults and reproduces them", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  dd <- sd_data()
  form <- bf(y ~ x + (1 | g)) + gaussian()

  msg <- capture_messages(suppressWarnings(
    ds <- frm_sample(form, data = dd, chains = 1, iter = 400,
                     refresh = 0, seed = 1)))
  m <- paste(msg, collapse = "")
  expect_match(m, "default priors")
  expect_match(m, "prior = \"flat\"", fixed = TRUE)
  # one compact line per class
  expect_match(m, "\n  Intercept  ")
  expect_match(m, "\n  sd  ")
  expect_match(m, "\n  b  ")

  # prior_summary() gives back exactly what was announced
  pl <- prior_summary(ds)
  expect_s3_class(pl, "frmtmb_priorlist")
  txt <- utils::capture.output(print(pl))
  expect_length(txt, 3L)
  expect_true(any(grepl("class=Intercept", txt, fixed = TRUE)))
  expect_true(any(grepl("class=sd", txt, fixed = TRUE)))
  expect_true(any(grepl("scale=natural", txt, fixed = TRUE)))
  for (s in unclass(pl)) {
    expect_equal(s$dist$kind, "t")
    expect_equal(s$dist$df, 3)
  }
  # and the message is suppressible
  expect_silent(suppressWarnings(suppressMessages(
    frm_sample(form, data = dd, chains = 1, iter = 200, refresh = 0,
               seed = 1))))
})

test_that("an ordinal formula-route call announces its threshold gap", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  dd <- sd_prior_data()

  # no random effects: the priorlist is empty, and the call must SAY so
  msg <- capture_messages(suppressWarnings(
    ds <- frm_sample(bf(o ~ x) + cumulative(), data = dd, chains = 1,
                     iter = 300, refresh = 0, seed = 6)))
  m <- paste(msg, collapse = "")
  expect_match(m, "default priors")
  expect_match(m, "thresholds")
  expect_null(prior_summary(ds))
  # the model still sampled, thresholds included
  expect_s3_class(ds, "frmtmb_draws")
  expect_true("x" %in% colnames(ds$draws))
  expect_gt(sum(grepl("^tau_raw", colnames(ds$draws))), 0L)
})

test_that("prior = 'flat' opts out and warns about propriety", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  dd <- sd_data()
  form <- bf(y ~ x + (1 | g)) + gaussian()

  w <- NULL
  msg <- capture_messages(suppressWarnings(withCallingHandlers(
    ds <- frm_sample(form, data = dd, prior = "flat", chains = 1,
                     iter = 400, refresh = 0, seed = 1),
    warning = function(x) w <<- c(w, conditionMessage(x)))))
  expect_length(grep("default priors", msg), 0L)
  expect_length(grep("flat", w), 1L)
  expect_null(prior_summary(ds))

  # it really is the old flat behavior, and the fit route reaches it
  # through the same opt-out: a seeded short run matches draw for draw
  fit <- frm(form, data = dd)
  ds_fit <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 1,
               init = "random", prior = "flat")))
  expect_equal(unname(ds$draws), unname(ds_fit$draws))

  # a model with no variance component has nothing improper to warn about
  w2 <- NULL
  suppressWarnings(suppressMessages(withCallingHandlers(
    frm_sample(bf(y ~ x) + gaussian(), data = dd, prior = "flat",
               chains = 1, iter = 200, refresh = 0, seed = 1),
    warning = function(x) w2 <<- c(w2, conditionMessage(x)))))
  expect_length(grep("flat", w2), 0L)

  expect_error(frm_sample(bf(y ~ x) + gaussian(), data = dd,
                          prior = "weak"),
               "the string .flat.")
  expect_error(frm_sample(fit, prior = "weak"), "the string .flat.")
})

test_that("the fit route defaults to the brms priors, and flat opts out", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  dd <- sd_data()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  # the disclosure banner is the formula route's, on the fit route
  msg <- capture_messages(suppressWarnings(
    ds <- frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 1)))
  expect_length(grep("default priors", msg), 1L)
  cls <- vapply(unclass(prior_summary(ds)), `[[`, "", "class")
  expect_true(all(c("Intercept", "sd") %in% cls))
  # and the sd default is what makes the block eligible to non-center
  expect_equal(ds$reparam$blocks, 1L)
  expect_length(grep("stays centered", msg), 0L)

  # prior = "flat" restores the bare objective: no banner, no defaults
  # in the summary, no reparameterization, and the propriety warning
  w <- NULL
  msg2 <- capture_messages(withCallingHandlers(
    ds2 <- frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 1,
                      prior = "flat"),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }))
  expect_length(grep("default priors", msg2), 0L)
  expect_length(grep("flat", w), 1L)
  expect_null(prior_summary(ds2))
  expect_null(ds2$reparam)
  expect_length(grep("stays centered", msg2), 1L)
})

test_that("the fit route unpins a chain from a boundary variance mode", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  # no group effect in the data at all, so maximum likelihood puts the
  # variance component on the boundary: the kidney pathology in
  # miniature (dev/brms-vignette-audit.md)
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(80),
                   g = factor(rep(seq_len(8), length.out = 80)))
  dd$y <- stats::rnorm(80, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_lt(fit$estimates$theta[[1L]], -4)

  flat <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 800, refresh = 0, seed = 2,
               prior = "flat")))
  def <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 800, refresh = 0, seed = 2)))
  # structural either way: the defaults are what open the gate
  expect_null(flat$reparam)
  expect_equal(def$reparam$blocks, 1L)

  if (sampler_gates_on()) {
    med <- function(ds) stats::median(exp(ds$draws[, "theta_1"]))
    # flat: the chain sits at the singular mode or walks the improper
    # tail below it, and reports a variance component of nothing
    expect_lt(med(flat), 0.01)
    # defaults: a variance component with a posterior on it
    expect_gt(med(def), 0.02)
    if (requireNamespace("posterior", quietly = TRUE)) {
      ess <- function(ds) posterior::ess_bulk(ds$draws[, "theta_1"])
      expect_gt(ess(def), 100)
      expect_gt(ess(def), 5 * ess(flat))
    }
  }
})

test_that("a MAP fit's prior stacks under a call prior and the defaults", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  dd <- sd_data()
  form <- bf(y ~ x + (1 | g)) + gaussian()
  mf <- frm(form, data = dd,
            prior = set_prior("normal(0, 0.3)", class = "b"))
  spec_of <- function(pl) {
    vapply(unclass(pl), function(s) {
      paste0(s$class, "=", s$dist$kind, "(",
             paste(unlist(s$dist[-1L]), collapse = ","), ")")
    }, "")
  }

  # the MAP prior survives, and the defaults fill the slots it left
  ds <- suppressWarnings(suppressMessages(
    frm_sample(mf, chains = 1, iter = 300, refresh = 0, seed = 1)))
  s <- spec_of(prior_summary(ds))
  expect_true("b=normal(0,0.3)" %in% s)
  expect_true(any(grepl("^Intercept=t", s)))
  expect_true(any(grepl("^sd=t", s)))

  # an explicit call-level prior takes the slot it addresses, and only
  # that one: the fit's own b prior is superseded, the defaults stay
  ds2 <- suppressWarnings(suppressMessages(
    frm_sample(mf, chains = 1, iter = 300, refresh = 0, seed = 1,
               prior = set_prior("normal(0, 0.05)", class = "b"))))
  s2 <- spec_of(prior_summary(ds2))
  expect_true("b=normal(0,0.05)" %in% s2)
  expect_false("b=normal(0,0.3)" %in% s2)
  expect_true(any(grepl("^Intercept=t", s2)))
  expect_true(any(grepl("^sd=t", s2)))
  # the tighter prior is the one the chain felt
  if (sampler_gates_on()) {
    expect_lt(stats::sd(ds2$draws[, "x"]), stats::sd(ds$draws[, "x"]))
  }

  # and the opt-out drops the MAP prior too, out loud: that prior is
  # taped INTO the fit's objective, so "flat" has to rebuild without it
  msg <- capture_messages(suppressWarnings(
    ds3 <- frm_sample(mf, chains = 1, iter = 300, refresh = 0, seed = 1,
                      prior = "flat")))
  expect_null(prior_summary(ds3))
  expect_length(grep("drops the prior this fit was made with", msg), 1L)
})

test_that("a prior is added to the sampled density exactly once", {
  # the trap: a MAP fit's penalty is taped INTO fit$obj at fit time, so
  # a sampling objective built by ADDING priors to fit$obj would count
  # that penalty twice. Every prior-carrying route rebuilds from the
  # bare likelihood instead, and this pins that down in objective
  # values, where a double count cannot hide
  set.seed(4)
  dd <- data.frame(x = stats::rnorm(120))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x, 1)
  pl <- set_prior("normal(0, 0.3)", class = "b")
  mf <- frm(bf(y ~ x) + gaussian(), data = dd, prior = pl)
  p <- mf$obj$par

  bare <- frmtmb:::prior_augmented_obj(mf, list())
  # the premise: the fit's own tape really does carry the penalty
  expect_gt(abs(mf$obj$fn(p) - bare$fn(p)), 1e-6)
  # rebuilding bare + the fit's own prior reproduces that tape exactly,
  # which is what makes the rebuild the safe way to add anything else
  own <- frmtmb:::resolve_prior_input(mf, pl)
  expect_equal(frmtmb:::prior_augmented_obj(mf, own$entries)$fn(p),
               mf$obj$fn(p), tolerance = 1e-8)

  # the sampling stack is defaults + the fit's prior, and the defaults
  # land on top of the penalized objective exactly once
  defs <- frmtmb:::default_priors_for(mf)
  all_e <- frmtmb:::resolve_prior_input(mf, defs + pl)$entries
  def_e <- frmtmb:::resolve_prior_input(mf, defs)$entries
  expect_equal(
    frmtmb:::prior_augmented_obj(mf, all_e)$fn(p) - mf$obj$fn(p),
    frmtmb:::prior_augmented_obj(mf, def_e)$fn(p) - bare$fn(p),
    tolerance = 1e-8)

  # and the resolved stack holds ONE entry per parameter position, so
  # there is nothing for a second density to attach to
  keys <- vapply(all_e, function(e) {
    paste0(e$comp, ".", paste(e$idx, collapse = ","))
  }, "")
  expect_equal(anyDuplicated(keys), 0L)
})

test_that("user priors override the defaults per class", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  dd <- sd_data()
  form <- bf(y ~ x + (1 | g)) + gaussian()

  ds <- suppressWarnings(suppressMessages(
    frm_sample(form, data = dd,
               prior = set_prior("exponential(1)", class = "sd"),
               chains = 1, iter = 400, refresh = 0, seed = 2)))
  pl <- unclass(prior_summary(ds))
  kinds <- vapply(pl, function(s) paste0(s$class, ":", s$dist$kind), "")
  # the sd default stepped aside; the Intercept ones stayed
  expect_true("sd:exponential" %in% kinds)
  expect_false("sd:t" %in% kinds)
  expect_true("Intercept:t" %in% kinds)

  # the legacy named-list spelling takes over the internal parameters it
  # names and leaves the remaining defaults in place
  ds2 <- suppressWarnings(suppressMessages(
    frm_sample(form, data = dd, prior = list(theta = prior_normal(0, 1)),
               chains = 1, iter = 400, refresh = 0, seed = 2)))
  pl2 <- prior_summary(ds2)
  expect_equal(names(attr(pl2, "overrides")), "theta")
  expect_true(any(grepl("class=Intercept",
                        utils::capture.output(print(pl2)))))
})

test_that("defaults tame a variance component flat priors cannot", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  set.seed(101)
  # three groups and no group-level signal: nothing stops a flat prior
  # on the log standard deviation running off toward minus infinity
  dd <- data.frame(g = factor(rep(1:3, each = 8L)))
  dd$y <- stats::rnorm(nrow(dd), 0, 1)
  form <- bf(y ~ 1 + (1 | g)) + gaussian()

  flat <- suppressWarnings(suppressMessages(
    frm_sample(form, data = dd, prior = "flat", chains = 1, iter = 1200,
               refresh = 0, seed = 4)))
  def <- suppressWarnings(suppressMessages(
    frm_sample(form, data = dd, chains = 1, iter = 1200, refresh = 0,
               seed = 4)))
  f <- flat$draws[, "theta_1"]
  d <- def$draws[, "theta_1"]
  # The flat chain visits log sd values the half-t simply does not, and
  # the gap is measured in the half-t chain's OWN spread rather than
  # against a fixed log sd. With three groups the posterior tail under
  # the half-t is genuinely long, and how far into it a chain gets is a
  # property of the sampler: the non-centered default reaches -3.7 where
  # the centered one stopped at -2.4, so a fixed floor would be testing
  # the parameterization instead of the prior.
  skip_if_not(sampler_gates_on(), "chain-agreement gates are off")
  expect_lt(stats::quantile(f, 0.025), -3.5)
  expect_gt((stats::quantile(d, 0.025) - stats::quantile(f, 0.025)) /
              stats::sd(d), 0.5)
  expect_gt(stats::sd(f) / stats::sd(d), 1.3)
  # and it mixes worse for it
  expect_gt(min(summary(def)[, "n_eff"]), min(summary(flat)[, "n_eff"]))
})

test_that("the formula route validates its own arguments", {
  dd <- sd_data()
  expect_error(frm_sample(bf(y ~ x) + gaussian()), "needs data =")
  expect_error(frm_sample(42, data = dd), "takes a frmtmb fit or a formula")
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  expect_error(frm_sample(fit, data = dd), "already fixed")
})

test_that("a mixture sampled from a formula uses random inits", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  set.seed(17)
  dd <- data.frame(y = c(stats::rnorm(60, -2, 0.6),
                         stats::rnorm(60, 2, 0.6)))
  form <- bf(y ~ 1) + mixture(gaussian(), gaussian())
  # the multimodality section of ?frm_sample recommends init = "random"
  # for a mixture; from a formula that is simply the default, since
  # there is no mode to anchor on
  ds <- suppressWarnings(suppressMessages(
    frm_sample(form, data = dd, chains = 1, iter = 400, refresh = 0,
               seed = 2)))
  expect_s3_class(ds, "frmtmb_draws")
  # the two component means come back separated, whichever label they
  # took
  fe <- fixef(ds)[, "Estimate"]
  # whether one short chain finds BOTH modes is chain luck on some
  # platforms, so the separation claim is a gated agreement assert
  if (sampler_gates_on()) {
    expect_gt(abs(fe[["mu1_Intercept"]] - fe[["mu2_Intercept"]]), 2)
  }
})

## ---- the draws-side name convention ----------------------------------

test_that("every draws accessor speaks the same parenthesis-free names", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  skip_if_not_installed("posterior")
  withr::local_options(mc.cores = 1)
  dd <- sd_data()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 400,
                                    refresh = 0, seed = 5))

  nm <- setdiff(colnames(ds$draws),
                c("lp__", grep("^b\\[", colnames(ds$draws), value = TRUE)))
  expect_equal(nm, c("Intercept", "x", "sigma_Intercept", "theta_1"))
  expect_false(any(grepl("[()]", colnames(ds$draws))))
  expect_equal(rownames(summary(ds)), nm)
  expect_equal(rownames(fixef(ds)), nm[1:3])
  expect_equal(setdiff(variables(ds), c("lp__")),
               setdiff(colnames(ds$draws), "lp__"))
  expect_equal(colnames(posterior::as_draws_matrix(as_draws(ds))),
               colnames(ds$draws))

  # the FIT side keeps its own canonical spelling
  expect_true("(Intercept)" %in% rownames(stats::vcov(fit)))

  # both spellings resolve in every draws-side lookup
  expect_equal(hypothesis(ds, "Intercept > 0")$estimate,
               hypothesis(ds, "`(Intercept)` > 0")$estimate)
  p1 <- frmtmb:::resolve_priors(fit, list(Intercept = prior_normal(0, 1)))
  p2 <- frmtmb:::resolve_priors(fit,
                                list(`(Intercept)` = prior_normal(0, 1)))
  expect_equal(p1[[1L]]$comp, p2[[1L]]$comp)
  expect_equal(p1[[1L]]$idx, p2[[1L]]$idx)

  # check_laplace() reports the same names
  cl <- suppressWarnings(check_laplace(fit, chains = 1, iter = 400,
                                       refresh = 0, seed = 5))
  expect_equal(cl$parameter, nm)
})

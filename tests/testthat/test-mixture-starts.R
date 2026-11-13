# The jittered starts of an ordinal mixture (mixture_starts()): the best
# maximum at which no component is a step function of its predictors, the
# degenerate boundary where maximum likelihood often has its supremum.
# The data are the two-class process of dev/optima-ordmix-starts.R
# (dev/ordstart-findings.md has the 40-fit tables).

oms_data <- function(seed, n = 300, k = 0) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  # a few ulps, so that each rounding path of the optimizer is tried
  d$x <- d$x * (1 + 2e-15 * k)
  d
}

# the fit and every warning it raised, as conditions
oms_fit <- function(d, fam, ..., f = bf(y ~ x + z)) {
  w <- list()
  fit <- withCallingHandlers(frm(f, family = fam, data = d, ...),
                             warning = function(cnd) {
                               w[[length(w) + 1L]] <<- cnd
                               invokeRestart("muffleWarning")
                             })
  list(fit = fit, warnings = w,
       msgs = vapply(w, conditionMessage, ""))
}

# mixture_starts_pick() with a count of its calls, one per run of the
# starts
oms_count_picks <- function(env = parent.frame()) {
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  pick <- frmtmb:::mixture_starts_pick
  local_mocked_bindings(mixture_starts_pick = function(...) {
    counter$n <- counter$n + 1L
    pick(...)
  }, .package = "frmtmb", .env = env)
  counter
}

test_that("the pick keeps a proper default and leaves a step function", {
  pick <- frmtmb:::mixture_starts_pick
  # proper default: a better proper point replaces it, a better flat
  # boundary or step function does not, and a gain within tol is none
  expect_identical(pick(c(10, 9), c(0L, 0L), 1e-3), 2L)
  expect_identical(pick(c(10, 9), c(0L, 1L), 1e-3), 1L)
  expect_identical(pick(c(10, 9), c(0L, 2L), 1e-3), 1L)
  expect_identical(pick(c(10, 10 - 5e-4), c(0L, 0L), 1e-3), 1L)
  # a flat boundary default stays, so that its warning does: the rule
  # flags one of two points of one likelihood by how far the optimizer
  # ran (cum_cum seed 10 under autoscale lost it)
  expect_identical(pick(c(10, 9, 11), c(1L, 0L, 0L), 1e-3), 1L)
  expect_identical(pick(c(10, 9), c(1L, 1L), 1e-3), 1L)
  expect_identical(pick(c(10, 10 - 5e-4), c(1L, 0L), 1e-3), 1L)
  expect_identical(pick(c(10, 9), c(1L, 2L), 1e-3), 1L)
  # a step default: the best proper point, whatever its likelihood, then
  # the best flat boundary, and itself when every start is a step
  expect_identical(pick(c(10, 9, 12, 13), c(2L, 2L, 0L, 0L), 1e-3), 3L)
  expect_identical(pick(c(10, 11, 12), c(2L, 1L, 0L), 1e-3), 3L)
  expect_identical(pick(c(10, 12, 13), c(2L, 1L, 1L), 1e-3), 2L)
  expect_identical(pick(c(10, 9), c(2L, 2L), 1e-3), 1L)
  expect_identical(pick(10, 2L, 1e-3), 1L)
  # a default the optimizer did not call converged: a step function gives
  # way to the best proper point only, and nothing else gives way
  expect_identical(pick(c(10, 12, 13), c(2L, 1L, 0L), 1e-3, FALSE), 3L)
  expect_identical(pick(c(10, 12, 13, 14), c(2L, 0L, 0L, 1L), 1e-3, FALSE),
                   2L)
  expect_identical(pick(c(10, 9), c(2L, 1L), 1e-3, FALSE), 1L)
  expect_identical(pick(c(10, 9), c(0L, 0L), 1e-3, FALSE), 1L)
  expect_identical(pick(c(10, 9), c(1L, 0L), 1e-3, FALSE), 1L)
})

test_that("an unconverged step-function default gives way to a proper one", {
  skip_on_cran()
  # cum_sratio seed 4: 0.81.0 ends at a step function with singular
  # convergence (code 1, logLik -383.63); the first build moved it to a
  # proper maximum at -388.08, and the fix round's code-0 rule kept the
  # step function, while the autoscaled fit of the same data moved.
  # Valid endings: a converged fit with no step function, or the default
  # itself when no converged proper start exists (the starts ran, and the
  # fit-end check says none escaped).
  fam <- mixture(cumulative(), sratio())
  for (k in 0:1) {
    d <- oms_data(4, k = k)
    off <- oms_fit(d, fam, control = frmtmb_control(starts = 0))$fit
    r <- oms_fit(d, fam)
    on <- r$fit
    off_step <- length(frmtmb:::mixture_ord_start_class(off, "y")$step) > 0L
    if (isTRUE(off$opt$convergence == 0) || !off_step) next
    expect_false(is.null(on$opt$mixture_starts), info = paste("k", k))
    if (identical(on$opt$par, off$opt$par)) {
      expect_true(any(grepl("None of 10 jittered starts", r$msgs,
                            fixed = TRUE)), info = paste("k", k))
    } else {
      expect_true(isTRUE(on$opt$convergence == 0), info = paste("k", k))
      expect_length(frmtmb:::mixture_ord_start_class(on, "y")$step, 0L)
      expect_length(frmtmb:::mixture_ord_start_class(on, "y")$boundary, 0L)
    }
  }
})

test_that("a step-function default gives way to a proper maximum", {
  skip_on_cran()
  # cum_sratio seed 9: 0.81.0 ended at a step function (logLik -374.88,
  # code 0); 10 of 40 jittered starts reach a point without one, the
  # best at -378.61 (dev/ordstart-classes.R). Valid endings: no step
  # function, with the warning that a higher step function was passed
  # over when one was; a step function that the fit-end check says no
  # start escaped; or a default the optimizer did not call converged,
  # which runs no starts and stays.
  fam <- mixture(cumulative(), sratio())
  for (k in 0:2) {
    r <- oms_fit(oms_data(9, k = k), fam)
    dg <- frmtmb:::mixture_ord_degeneracy(r$fit, "y")
    step <- dg$component[(dg$sharpen > -0.1) %in% TRUE]
    ms <- r$fit$opt$mixture_starts
    if (is.null(ms)) {
      expect_false(isTRUE(r$fit$opt$convergence == 0), info = paste("k", k))
      next
    }
    if (length(step)) {
      expect_true(any(grepl("None of 10 jittered starts", r$msgs,
                            fixed = TRUE)), info = paste("k", k))
    } else {
      expect_false(any(grepl("a step function of its predictors: doubling",
                             r$msgs, fixed = TRUE)), info = paste("k", k))
    }
    sw <- Filter(function(w) inherits(w, "frmtmb_mixture_step_start"),
                 r$warnings)
    # the warning names the step function's likelihood, above the fit's
    for (w in sw) {
      expect_match(conditionMessage(w), "is a step function of its predictors")
      expect_false(length(step) > 0L)
    }
  }
})

test_that("a proper default never loses likelihood or gains the warning", {
  skip_on_cran()
  # cum_cum seed 4: the default start ends at a proper maximum (-388.10)
  # and 13 of 40 jittered starts at step functions up to 6.1 higher
  # (dev/ordstart-classes.R). On every rounding the fit keeps a proper
  # maximum at least as high as the default start's.
  fam <- mixture(cumulative(), cumulative())
  for (k in 0:2) {
    d <- oms_data(4, k = k)
    off <- oms_fit(d, fam, control = frmtmb_control(starts = 0))
    on <- oms_fit(d, fam)
    expect_null(off$fit$opt$mixture_starts)
    tol <- 1e-3
    off_proper <- !any(grepl("ended at a degenerate boundary", off$msgs,
                             fixed = TRUE))
    if (off_proper) {
      expect_false(any(grepl("ended at a degenerate boundary", on$msgs,
                             fixed = TRUE)), info = paste("k", k))
      expect_gte(as.numeric(logLik(on$fit)),
                 as.numeric(logLik(off$fit)) - tol)
    }
    for (w in Filter(function(w) {
      grepl("jittered starts reached a higher", conditionMessage(w))
    }, on$warnings)) {
      expect_s3_class(w, "frmtmb_mixture_step_start")
    }
  }
})

test_that("the starts draw under a named generator, whatever the session's", {
  draws <- frmtmb:::mixture_starts_draws
  old <- RNGkind()
  on.exit(RNGkind(old[1L], old[2L], old[3L]), add = TRUE)
  ref <- draws(3, 2)
  for (k in list(c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
                 c("Mersenne-Twister", "Box-Muller", "Rejection"),
                 c("Mersenne-Twister", "Inversion", "Rounding"))) {
    suppressWarnings(RNGkind(k[1L], k[2L], k[3L]))
    set.seed(3)
    before <- .Random.seed
    expect_identical(draws(3, 2), ref, label = paste(k, collapse = " "))
    # the caller's kinds and seed are put back
    expect_identical(RNGkind(), k)
    expect_identical(.Random.seed, before)
  }
  RNGkind(old[1L], old[2L], old[3L])
  # a seed that did not exist does not exist after
  ge <- globalenv()
  saved <- get0(".Random.seed", envir = ge, inherits = FALSE)
  on.exit(if (!is.null(saved)) assign(".Random.seed", saved, envir = ge),
          add = TRUE)
  if (exists(".Random.seed", envir = ge, inherits = FALSE)) {
    rm(".Random.seed", envir = ge)
  }
  draws(2, 1)
  expect_false(exists(".Random.seed", envir = ge, inherits = FALSE))
})

test_that("the same call gives the same fit under any RNG kind", {
  skip_on_cran()
  # cum_cum seed 14: under the session's kind the starts gave -389.68 with
  # the default kind and -379.31 under L'Ecuyer-CMRG; 11 of the 40 fits
  # differed so (dev/reviews/2026-11-11-ordstart.md, finding 1)
  old <- RNGkind()
  on.exit(RNGkind(old[1L], old[2L], old[3L]), add = TRUE)
  d <- oms_data(14)
  fam <- mixture(cumulative(), cumulative())
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  a <- oms_fit(d, fam)$fit
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)
  st <- .Random.seed
  b <- oms_fit(d, fam)$fit
  expect_identical(RNGkind()[1L], "L'Ecuyer-CMRG")
  expect_identical(.Random.seed, st)
  expect_identical(a$opt$par, b$opt$par)
})

test_that("the starts neither read nor move the session's random state", {
  d <- oms_data(1)
  fam <- mixture(cumulative(), cumulative())
  set.seed(7)
  a <- stats::runif(3)
  set.seed(7)
  f1 <- oms_fit(d, fam)$fit
  b <- stats::runif(3)
  expect_identical(a, b)
  # and another state before the fit gives the same fit
  set.seed(99)
  f2 <- oms_fit(d, fam)$fit
  expect_identical(f1$opt$par, f2$opt$par)
  expect_identical(f1$opt$mixture_starts, f2$opt$mixture_starts)
})

test_that("the starts run only where they are asked for", {
  d <- oms_data(1)
  mix <- mixture(cumulative(), cumulative())
  # absent: a family that declares no degenerate boundary
  expect_null(oms_fit(d, cumulative())$fit$opt$mixture_starts)
  dg <- data.frame(x = d$x, z = d$z, y = d$x + stats::rnorm(nrow(d)))
  fg <- allow_warnings(frm(bf(y ~ x + z), family = mixture(gaussian(),
                                                           gaussian()),
                           data = dg),
                       c("converge", "Standard errors are not available"))
  expect_null(fg$opt$mixture_starts)
  # off, and a start of one's own, turn the default starts off; a count
  # runs them anyway
  expect_null(oms_fit(d, mix, control = frmtmb_control(starts = 0))$fit$opt$
                mixture_starts)
  f0 <- oms_fit(d, mix)$fit
  expect_identical(f0$opt$mixture_starts[["starts"]], 10)
  st <- f0$estimates[c("beta", "tau_raw1", "tau_raw2")]
  expect_null(oms_fit(d, mix, start = st)$fit$opt$mixture_starts)
  ms <- oms_fit(d, mix, start = st,
                control = frmtmb_control(starts = 2))$fit$opt$mixture_starts
  expect_identical(ms[["starts"]], 2)
})

test_that("an unmeasured sharpening counts as a step function", {
  # three cumulative components: the doubled likelihood is not finite at
  # points with coefficients of 740 to 9.5e4, and an NA sharpening ranked
  # them as sound (dev/reviews/2026-11-11-ordstart.md, finding 3)
  local_mocked_bindings(
    mixture_ord_degeneracy = function(fit, resp) {
      data.frame(component = 1:3, sharpen = c(NA, -5, 0.2),
                 near = c(1, 1, 1), collapsed = c(FALSE, FALSE, FALSE))
    }, .package = "frmtmb")
  cl <- frmtmb:::mixture_ord_start_class(NULL, "y")
  expect_identical(cl$step, c(1L, 3L))
  expect_identical(cl$boundary, integer(0))
})

test_that("only converged runs are ranked, and an unconverged default stays", {
  skip_on_cran()
  # the review's three-component process (dev/reviews/ordstart-review-
  # sweep.R, K = 5): seed 1 traded a code-0 default for a run stopped at
  # the evaluation limit, seed 7 for a singular-convergence point, and
  # seed 4's code-1 default for a point with coefficients of 1e4
  gen <- function(seed, n = 300, K = 5) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n), z = rnorm(n))
    cls <- stats::rbinom(n, 1, 0.4)
    lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
    cuts <- stats::quantile(lat, seq_len(K - 1) / K)
    d$y <- 1L + rowSums(outer(lat, cuts, ">"))
    d
  }
  fam <- mixture(cumulative(), cumulative(), cumulative())
  for (s in c(1, 4, 7)) {
    d <- gen(s)
    on <- oms_fit(d, fam)$fit
    off <- oms_fit(d, fam, control = frmtmb_control(starts = 0))$fit
    # no fit gains a convergence failure
    expect_false(isTRUE(off$opt$convergence == 0) &&
                   !isTRUE(on$opt$convergence == 0), info = paste("seed", s))
    # an unconverged default stays, unless a converged proper start
    # replaced its step function
    if (!isTRUE(off$opt$convergence == 0) &&
          !identical(on$opt$par, off$opt$par)) {
      expect_true(isTRUE(on$opt$convergence == 0), info = paste("seed", s))
    }
    # a point the starts chose is measured and sound
    if (isTRUE(on$opt$mixture_starts[["gain"]] != 0)) {
      expect_lt(max(abs(on$opt$par)), 740)
      expect_length(frmtmb:::mixture_ord_start_class(on, "y")$step, 0L)
    }
  }
})

test_that("autoscale runs the starts once, in the fit that is reported", {
  skip_on_cran()
  # x on a 1e-4 scale engages autoscale's pre-fit. The starts ran there
  # and again in the fit, and the pre-fit muffled the step-start warning
  # on 11 of 11 fits (dev/reviews/2026-11-11-ordstart.md, finding 2)
  picks <- oms_count_picks()
  fam <- mixture(cumulative(), cumulative())
  wfit <- function(d, f) {
    picks$n <- 0L
    r <- oms_fit(d, fam, f = f)
    list(fit = r$fit, picks = picks$n,
         step = any(vapply(r$warnings, inherits, NA,
                           "frmtmb_mixture_step_start")),
         degen = any(grepl("ended at a degenerate boundary", r$msgs,
                           fixed = TRUE)))
  }
  d <- oms_data(4)
  d$xs <- d$x * 1e-4
  plain <- wfit(d, bf(y ~ x + z))
  scaled <- wfit(d, bf(y ~ xs + z))
  expect_true(isTRUE(scaled$fit$autoscaled))
  expect_identical(plain$picks, 1L)
  expect_identical(scaled$picks, 1L)
  expect_false(is.null(scaled$fit$opt$mixture_starts))
  expect_identical(scaled$step, plain$step)
  # cum_cum seed 10: a flat-boundary default, whose degenerate warning
  # 0.81.0 gives with and without the rescaling
  d <- oms_data(10)
  d$xs <- d$x * 1e-4
  expect_true(wfit(d, bf(y ~ x + z))$degen)
  expect_true(wfit(d, bf(y ~ xs + z))$degen)
})

test_that("frm_allfit() refits run no starts", {
  skip_on_cran()
  f <- oms_fit(oms_data(4), mixture(cumulative(), cumulative()))$fit
  picks <- oms_count_picks()
  af <- frm_allfit(f, optimizers = list(nlminb = "nlminb", optim = "optim"),
                   start_from_mle = FALSE)
  expect_s3_class(af, "frmtmb_allfit")
  expect_identical(picks$n, 0L)
})

test_that("random effects run no starts unless a count asks", {
  cnt <- frmtmb:::mixture_starts_count
  expect_identical(cnt(list(), cold = TRUE, random = FALSE), 10L)
  expect_identical(cnt(list(), cold = TRUE, random = TRUE), 0L)
  expect_identical(cnt(list(), cold = FALSE, random = FALSE), 0L)
  expect_identical(cnt(list(starts = 3L), cold = FALSE, random = TRUE), 3L)
  skip_on_cran()
  # each start is a Laplace fit, 15 to 23 times the fit in all
  # (dev/ordstart-cost.R), so the default leaves them off
  d <- oms_data(1, n = 200)
  d$g <- factor(rep_len(1:10, nrow(d)))
  w <- list()
  f <- withCallingHandlers(
    frm(bf(y ~ x + (1 | g)), family = mixture(cumulative(), cumulative()),
        data = d),
    warning = function(cnd) {
      w[[length(w) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    })
  expect_null(f$opt$mixture_starts)
})

test_that("frmtmb_control(starts = ) takes a count or NULL", {
  expect_null(frmtmb_control()$starts)
  expect_identical(frmtmb_control(starts = 3)$starts, 3L)
  for (bad in list(-1, 1.5, "a", c(1, 2), NA)) {
    expect_error(frmtmb_control(starts = bad), "`starts` must be")
  }
})

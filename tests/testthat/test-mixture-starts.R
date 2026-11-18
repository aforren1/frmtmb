# The jittered and swapped starts of an ordinal mixture
# (mixture_starts()): the best maximum at which no component is a step
# function of its predictors, the degenerate boundary where maximum
# likelihood often has its supremum. The data are the two-class process
# of dev/optima-ordmix-starts.R (dev/ordstart-findings.md has the 40-fit
# tables, dev/ordstarts-findings.md the swapped start's).

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

# no fit runs jittered starts by default (backlog 6136), so a test of
# them asks for the count the default ran before
oms_ten <- frmtmb_control(starts = 10)

# the fit and every warning it raised, as conditions, the step-start
# messages (`steps`) and the kept-swap messages (`swaps`)
oms_fit <- function(d, fam, ..., f = bf(y ~ x + z)) {
  w <- list()
  s <- list()
  sw <- list()
  step_cls <- "frmtmb_mixture_step_start"
  swap_cls <- "frmtmb_mixture_swap"
  fit <- withCallingHandlers(frm(f, family = fam, data = d, ...),
                             warning = function(cnd) {
                               w[[length(w) + 1L]] <<- cnd
                               invokeRestart("muffleWarning")
                             },
                             message = function(cnd) {
                               if (inherits(cnd, step_cls)) {
                                 s[[length(s) + 1L]] <<- cnd
                                 invokeRestart("muffleMessage")
                               }
                               if (inherits(cnd, swap_cls)) {
                                 sw[[length(sw) + 1L]] <<- cnd
                                 invokeRestart("muffleMessage")
                               }
                             })
  # the step-start finding and a kept swap are messages (the user's
  # decisions of 2026-11-13 and 2026-11-17), never warnings
  expect_false(any(vapply(w, inherits, NA, step_cls)))
  expect_false(any(vapply(w, inherits, NA, swap_cls)))
  list(fit = fit, warnings = w, steps = s, swaps = sw,
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
    r <- oms_fit(d, fam, control = oms_ten)
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
  # function, with the message that a higher step function was passed
  # over when one was; a step function that the fit-end check says no
  # start escaped; or a default the optimizer did not call converged,
  # which runs no starts and stays.
  fam <- mixture(cumulative(), sratio())
  for (k in 0:2) {
    r <- oms_fit(oms_data(9, k = k), fam, control = oms_ten)
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
    # the message names the step function's likelihood, above the fit's
    for (m in r$steps) {
      expect_match(conditionMessage(m), "is a step function of its predictors")
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
    on <- oms_fit(d, fam, control = oms_ten)
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
    # the finding is a message of its class and nothing else
    expect_false(any(grepl("jittered starts reached a higher", on$msgs,
                           fixed = TRUE)))
    for (m in on$steps) {
      expect_match(conditionMessage(m), "jittered starts reached a higher",
                   fixed = TRUE)
    }
  }
})

test_that("the step-start finding is a message, not a warning", {
  # the user's decision of 2026-11-13: the point returned is the better
  # answer, so its signal is a message; the class is kept for a handler
  sig <- frmtmb:::mixture_starts_signal
  got <- NULL
  w <- NULL
  withCallingHandlers(sig("a step function was higher"),
                      message = function(m) {
                        got <<- m
                        invokeRestart("muffleMessage")
                      },
                      warning = function(x) {
                        w <<- x
                        invokeRestart("muffleWarning")
                      })
  expect_null(w)
  expect_s3_class(got, "frmtmb_mixture_step_start")
  expect_s3_class(got, "message")
  expect_identical(conditionMessage(got), "a step function was higher\n")
})

test_that("the starts draw under a named generator, whatever the session's", {
  draws <- frmtmb:::mixture_starts_draws
  # every kind RNGkind() reports, set by position: R-devel has a fourth,
  # binom.kind, after the three of R 4.6
  set_kinds <- function(k) suppressWarnings(do.call(RNGkind, as.list(k)))
  old <- RNGkind()
  on.exit(set_kinds(old), add = TRUE)
  ref <- draws(3, 2)
  cases <- list(c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
                c("Mersenne-Twister", "Box-Muller", "Rejection"),
                c("Mersenne-Twister", "Inversion", "Rounding"))
  # the fourth kind away from its default, where the R version has one
  if (length(old) >= 4L) {
    cases <- c(cases, list(c("Mersenne-Twister", "Inversion", "Rejection",
                             "Buggy BTPE")))
  }
  for (k in cases) {
    set_kinds(k)
    set.seed(3)
    want <- RNGkind()
    before <- .Random.seed
    expect_identical(draws(3, 2), ref, label = paste(k, collapse = " "))
    # the caller's kinds, all of them, and seed are put back
    expect_identical(RNGkind(), want)
    expect_identical(.Random.seed, before)
  }
  set_kinds(old)
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
  on.exit(suppressWarnings(do.call(RNGkind, as.list(old))), add = TRUE)
  d <- oms_data(14)
  fam <- mixture(cumulative(), cumulative())
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  a <- oms_fit(d, fam, control = oms_ten)$fit
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)
  st <- .Random.seed
  b <- oms_fit(d, fam, control = oms_ten)$fit
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
  f1 <- oms_fit(d, fam, control = oms_ten)$fit
  b <- stats::runif(3)
  expect_identical(a, b)
  # and another state before the fit gives the same fit
  set.seed(99)
  f2 <- oms_fit(d, fam, control = oms_ten)$fit
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
  # the default runs no jittered starts, and a mixture of one family has
  # no swap to run (backlog 6136), so no record; a count runs them, from
  # a start of one's own too
  expect_null(oms_fit(d, mix)$fit$opt$mixture_starts)
  expect_null(oms_fit(d, mix, control = frmtmb_control(starts = 0))$fit$opt$
                mixture_starts)
  f0 <- oms_fit(d, mix, control = oms_ten)$fit
  expect_identical(f0$opt$mixture_starts[["starts"]], 10)
  expect_identical(f0$opt$mixture_starts[["swaps"]], 0)
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
    },
    mixture_ord_thres_se = function(fit, resp, force = FALSE,
                                    lost_inf = FALSE, ...) {
      structure(c(1, 1, 1), lost = rep(FALSE, 3))
    }, .package = "frmtmb")
  cl <- frmtmb:::mixture_ord_start_class(NULL, "y")
  expect_identical(cl$step, c(1L, 3L))
  expect_identical(cl$boundary, integer(0))
})

test_that("a threshold above the cut, or lost, is flat to the ranking", {
  # the ranking called proper what the degenerate check warned about:
  # today's 10 starts moved 3 of 120 proper fits without random effects to
  # such a point (review of 2026-11-17), and 7 upward moves with (1 | g)
  # ended at a threshold the standard-error check lost
  # (dev/ordstarts-findings.md)
  dg <- data.frame(component = 1:3, sharpen = c(-5, -5, -5),
                   near = c(1, 1, 1), collapsed = c(FALSE, FALSE, FALSE))
  for (se in list(c(1, 101, 1), c(1, Inf, 1))) {
    local_mocked_bindings(
      mixture_ord_degeneracy = function(fit, resp) dg,
      mixture_ord_thres_se = function(fit, resp, force = FALSE,
                                      lost_inf = FALSE, ...) {
        expect_true(force)
        expect_true(lost_inf)
        structure(se, lost = is.infinite(se))
      }, .package = "frmtmb")
    cl <- frmtmb:::mixture_ord_start_class(NULL, "y")
    expect_identical(cl$step, integer(0))
    expect_identical(cl$unplaced, 2L)
    expect_identical(cl$boundary, integer(0))
  }
  # absent: every threshold placed, nothing flat
  local_mocked_bindings(
    mixture_ord_degeneracy = function(fit, resp) dg,
    mixture_ord_thres_se = function(fit, resp, force = FALSE,
                                    lost_inf = FALSE, ...) {
      structure(c(1, 99, 1), lost = rep(FALSE, 3))
    }, .package = "frmtmb")
  cl <- frmtmb:::mixture_ord_start_class(NULL, "y")
  expect_identical(cl$boundary, integer(0))
  expect_identical(cl$unplaced, integer(0))
  # the pick: a default flat by travel stays (sticky), one flat by its
  # standard errors gives way to a better proper point, and only to a
  # better one (review of 2026-11-18, item 1)
  pick <- frmtmb:::mixture_starts_pick
  expect_identical(pick(c(10, 9), c(1L, 0L), 1e-3, sticky = TRUE), 1L)
  expect_identical(pick(c(10, 9), c(1L, 0L), 1e-3, sticky = FALSE), 2L)
  expect_identical(pick(c(10, 11), c(1L, 0L), 1e-3, sticky = FALSE), 1L)
  expect_identical(pick(c(10, 9), c(1L, 1L), 1e-3, sticky = FALSE), 1L)
})

test_that("a lost threshold has no standard error, not its entry's", {
  # the covariance entry of a lost threshold read as 1.8 (cum_sratio n =
  # 300 seed 24 of dev/ordstarts-run.R), below the cut
  d <- oms_data(13)
  f <- suppressMessages(frm(bf(y ~ x + z), family = mixture(cumulative(),
                                                            sratio()),
                            data = d, control = frmtmb_control(starts = 0)))
  thse <- frmtmb:::mixture_ord_thres_se
  base <- thse(f, "y")
  expect_true(all(is.finite(base)))
  expect_identical(attr(base, "lost"), c(FALSE, FALSE))
  V <- f$cache$sdr$cov.fixed %||% frmtmb:::se_fit_analysis(f)$V
  # the check lost the middle raw threshold of component 1
  f$cache <- new.env(parent = emptyenv())
  f$cache$sdr <- list(cov.fixed = V, se_lost = c(tau_raw1_2 = NA))
  got <- thse(f, "y", lost_inf = TRUE)
  expect_identical(got[1], Inf)
  expect_identical(got[2], base[2])
  expect_identical(attr(got, "lost"), c(TRUE, FALSE))
  # the fit-end check keeps the entry (its lost list moves with the
  # weights' scale): the number, and the lost flag beside it
  kept <- thse(f, "y")
  expect_identical(unclass(kept)[1:2], unclass(base)[1:2])
  expect_identical(attr(kept, "lost"), c(TRUE, FALSE))
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
    on <- oms_fit(d, fam, control = oms_ten)$fit
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
    r <- oms_fit(d, fam, f = f, control = oms_ten)
    list(fit = r$fit, picks = picks$n,
         step = length(r$steps) > 0L,
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
  # the step-start signal is that of the one run whose record the fit
  # carries. How many starts reach a higher step function is each run's
  # optimizer path, which rescaling and rounding change: 0 to 6 of 10
  # over 240 fits of these data, every one at the same maximum, and with
  # x scaled by 1 + 48e-15 the scaled fit's starts reached none while the
  # plain fit's reached 5 (dev/ci082-findings.md)
  for (a in list(plain, scaled)) {
    rec <- a$fit$opt$mixture_starts
    expect_identical(a$step, rec[["degenerate_above"]] > 0 &&
                       rec[["tier"]] < 2)
  }
  # and both report the same maximum: every one of those 240 fits read
  # -388.095489 to 6 decimals, so within 1e-6 of each other, and 1e-7 of
  # the log-likelihood is 3.9e-5
  llp <- as.numeric(logLik(plain$fit))
  expect_lt(abs(as.numeric(logLik(scaled$fit)) - llp), 1e-7 * abs(llp))
  # cum_cum seed 10: a flat-boundary default, whose degenerate warning
  # 0.81.0 gives with and without the rescaling
  d <- oms_data(10)
  d$xs <- d$x * 1e-4
  expect_true(wfit(d, bf(y ~ x + z))$degen)
  expect_true(wfit(d, bf(y ~ xs + z))$degen)
})

test_that("frm_allfit() refits run no starts", {
  skip_on_cran()
  f <- oms_fit(oms_data(4), mixture(cumulative(), cumulative()),
               control = oms_ten)$fit
  picks <- oms_count_picks()
  af <- frm_allfit(f, optimizers = list(nlminb = "nlminb", optim = "optim"),
                   start_from_mle = FALSE)
  expect_s3_class(af, "frmtmb_allfit")
  expect_identical(picks$n, 0L)
})

test_that("no fit runs jittered starts unless a count asks", {
  # backlog 6136, the user's decision of 2026-11-17: a median 18 times
  # the fit with random effects, 5.6 to 7.1 times without
  # (dev/ordstarts-findings.md), so no model form runs them by default
  cnt <- frmtmb:::mixture_starts_count
  expect_identical(cnt(list(), cold = TRUE), 0L)
  expect_identical(cnt(list(), cold = FALSE), 0L)
  expect_identical(cnt(list(starts = 3L), cold = FALSE), 3L)
  expect_identical(cnt(list(starts = 0L), cold = TRUE), 0L)
})

test_that("the ranking keeps today's three proper fits proper", {
  skip_on_cran()
  # the review of 2026-11-17 (dev/reviews/ordstarts-review-nore.R): 10
  # starts moved these proper fits without the degenerate warning to a
  # higher point with it (largest threshold standard error 2705, 2853 and
  # 5190), which frmtmb 0.83.0 did by default. Valid endings: the fit with
  # starts = 0, or a higher point with no degenerate component
  fam <- mixture(cumulative(), sratio())
  for (cfg in list(c(1000, 13), c(1000, 34), c(300, 28))) {
    d <- oms_data(cfg[2], n = cfg[1])
    off <- oms_fit(d, fam, control = frmtmb_control(starts = 0))
    on <- oms_fit(d, fam, control = oms_ten)
    lab <- paste("n", cfg[1], "seed", cfg[2])
    expect_false(any(grepl("degenerate boundary", off$msgs, fixed = TRUE)),
                 info = lab)
    expect_false(any(grepl("degenerate boundary", on$msgs, fixed = TRUE)),
                 info = lab)
    tol <- 1e-3
    expect_gte(as.numeric(logLik(on$fit)),
               as.numeric(logLik(off$fit)) - tol)
    cl <- frmtmb:::mixture_ord_start_class(on$fit, "y")
    expect_length(cl$step, 0L)
    expect_length(cl$boundary, 0L)
  }
})

test_that("a mixture of two families keeps a swapped fit that is better", {
  skip_on_cran()
  # cum_sratio n = 1000 seed 13 (dev/ordstarts-smoke.R): the swapped
  # start reached a proper maximum 0.835 above the default's. Valid
  # endings on each rounding: that point, said in a message, or the
  # default's own when the swap is not better or does not hold
  fam <- mixture(cumulative(), sratio())
  for (k in 0:1) {
    d <- oms_data(13, n = 1000, k = k)
    off <- oms_fit(d, fam, control = frmtmb_control(starts = 0))
    on <- oms_fit(d, fam)
    ms <- on$fit$opt$mixture_starts
    lab <- paste("k", k)
    expect_identical(ms[["swaps"]], 1, info = lab)
    expect_identical(ms[["starts"]], 0, info = lab)
    ll_off <- as.numeric(logLik(off$fit))
    ll_on <- as.numeric(logLik(on$fit))
    expect_gte(ll_on, ll_off - 1e-3)
    if (ms[["swapped"]] == 1) {
      expect_length(on$swaps, 1L)
      expect_match(conditionMessage(on$swaps[[1L]]),
                   "swapped (component 1 from 2, component 2 from 1)",
                   fixed = TRUE)
      expect_gt(ll_on, ll_off + 1e-3)
      expect_identical(on$fit$opt$convergence, 0L)
      cl <- frmtmb:::mixture_ord_start_class(on$fit, "y")
      expect_length(cl$step, 0L)
      expect_length(cl$boundary, 0L)
      expect_false(any(grepl("degenerate boundary", on$msgs, fixed = TRUE)))
    } else {
      expect_length(on$swaps, 0L)
      expect_identical(on$fit$opt$par, off$fit$opt$par)
    }
  }
})

test_that("a swap that is the same model is not run: the objective says", {
  # permuting components of one family, link and layout is label
  # switching: from the swapped start the objective differed by at most
  # 2.8e-9 on 78 cum_cum fits (dev/ordstarts-sum4.R). The objective
  # decides, not the families' names (review of 2026-11-18, item 3)
  d <- oms_data(1)
  same <- frmtmb:::mixture_swap_distinct
  cc <- oms_fit(d, mixture(cumulative(), cumulative()))$fit
  expect_null(cc$opt$mixture_starts)
  expect_length(same(cc, frmtmb:::mixture_ord_swap_starts(cc, "y")), 0L)
  # sratio and cratio under the logit are one model
  # (P(Y = j | Y >= j) = F(tau_j - eta) in both)
  sc <- oms_fit(d, mixture(sratio(), cratio()),
                control = frmtmb_control(starts = 0))$fit
  expect_length(same(sc, frmtmb:::mixture_ord_swap_starts(sc, "y")), 0L)
  expect_null(oms_fit(d, mixture(sratio(), cratio()))$fit$opt$
                mixture_starts)
  # one family, but a prior on one component only: not the same model,
  # so the swap runs
  pr <- set_prior("normal(-0.8, 0.2)", class = "b", coef = "x",
                  dpar = "mu1")
  cp <- oms_fit(d, mixture(cumulative(), cumulative()), prior = pr,
                control = frmtmb_control(starts = 0))$fit
  expect_length(same(cp, frmtmb:::mixture_ord_swap_starts(cp, "y")), 1L)
  # absent: two families, one assignment
  cs <- oms_fit(d, mixture(cumulative(), sratio()),
                control = frmtmb_control(starts = 0))$fit
  sw <- frmtmb:::mixture_ord_swap_starts(cs, "y")
  expect_length(sw, 1L)
  expect_length(same(cs, sw), 1L)
  # every parameter moved to the other component, the weight's logit
  # negated
  p <- cs$opt$par
  nm <- frmtmb:::outer_par_map(cs)$names
  q <- sw[[1L]]$par
  expect_identical(unname(q[nm == "mu1_x"]), unname(p[nm == "mu2_x"]))
  expect_identical(unname(q[nm == "mu2_z"]), unname(p[nm == "mu1_z"]))
  expect_identical(unname(q[names(p) == "tau_raw1"]),
                   unname(p[names(p) == "tau_raw2"]))
  expect_identical(unname(q[nm == "theta1_(Intercept)"]),
                   -unname(p[nm == "theta1_(Intercept)"]))
  # K = 3: every assignment but the identity
  expect_length(frmtmb:::mixture_swap_assignments(3), 5L)
  expect_length(frmtmb:::mixture_swap_assignments(2), 1L)
  # from K = 4 the transpositions only (dev/ordstarts-sum5.R)
  expect_length(frmtmb:::mixture_swap_assignments(4), 6L)
})

test_that("warm refits, starts = 0 and a start of one's own run no swap", {
  skip_on_cran()
  d <- oms_data(13)
  fam <- mixture(cumulative(), sratio())
  f <- oms_fit(d, fam)$fit
  expect_identical(f$opt$mixture_starts[["swaps"]], 1)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  cand <- frmtmb:::mixture_swap_candidates
  local_mocked_bindings(mixture_swap_candidates = function(...) {
    counter$n <- counter$n + 1L
    cand(...)
  }, .package = "frmtmb")
  r <- suppressMessages(refit(f, d$y))
  expect_null(r$opt$mixture_starts)
  af <- suppressMessages(frm_allfit(
    f, optimizers = list(nlminb = "nlminb", optim = "optim")))
  expect_s3_class(af, "frmtmb_allfit")
  off <- oms_fit(d, fam, control = frmtmb_control(starts = 0))$fit
  expect_null(off$opt$mixture_starts)
  st <- f$estimates[c("beta", "tau_raw1", "tau_raw2")]
  own <- oms_fit(d, fam, start = st)$fit
  expect_null(own$opt$mixture_starts)
  expect_identical(counter$n, 0L)
})

test_that("the swap's check keeps its best proper end, not the default", {
  # review of 2026-11-18, item 5: a jittered end that climbs to a proper
  # point is the better answer, and one that climbs to a step function
  # is the supremum the pick avoids
  settle <- frmtmb:::mixture_swap_settle
  o <- list(par = c(1, 2), objective = 10)
  yes <- function(r) TRUE
  expect_identical(settle(o, function(p) NULL, 2L, 1e-3, yes), o)
  near <- function(p) list(par = p, objective = 10 - 5e-4)
  expect_identical(settle(o, near, 2L, 1e-3, yes), o)
  up <- function(p) list(par = p, objective = 9)
  expect_identical(settle(o, up, 2L, 1e-3, yes)$objective, 9)
  expect_identical(settle(o, up, 2L, 1e-3, function(r) FALSE), o)
  # four runs, the last one counts too (it alone caught the one point
  # that is not a maximum, dev/ordstarts-findings.md); the best is kept
  n <- 0L
  runs <- function(p) {
    n <<- n + 1L
    list(par = p, objective = c(9.5, 10, 9, 8)[n])
  }
  expect_identical(settle(o, runs, 2L, 1e-3, yes)$objective, 8)
  expect_identical(n, 4L)
  # the jitter draws neither read nor move the session's random state
  set.seed(5)
  a <- .Random.seed
  settle(o, function(p) NULL, 2L, 1e-3, yes)
  expect_identical(.Random.seed, a)
})

test_that("the starts named by the warning run on a default flat by its SE", {
  skip_on_cran()
  # review of 2026-11-18, item 1: a default flat only by a threshold's
  # standard error (or a lost one) kept the starts off, so the call the
  # degenerate warning names did nothing on 30 of 120 fits. cum_cum seeds
  # 3 and 16 (0.83.0's starts reached -385.4217 and -379.8336 there)
  fam <- mixture(cumulative(), cumulative())
  for (s in c(3, 16)) {
    d <- oms_data(s)
    off <- oms_fit(d, fam)
    on <- oms_fit(d, fam, control = oms_ten)
    lab <- paste("seed", s)
    expect_identical(on$fit$opt$mixture_starts[["starts"]], 10, info = lab)
    expect_gte(as.numeric(logLik(on$fit)),
               as.numeric(logLik(off$fit)) - 1e-3)
    # the default's warning names the call, which then runs
    if (any(grepl("degenerate boundary", off$msgs, fixed = TRUE))) {
      expect_true(any(grepl("frmtmb_control(starts = 10)", off$msgs,
                            fixed = TRUE)), info = lab)
    }
  }
  # test-compat.R's cell 2: 0.83.0's starts found a clean point 1.77
  # above the default's threshold with standard error 21900
  set.seed(20261005)
  n <- 400
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$o <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  fam <- mixture(cumulative(threshold = "sum_to_zero"), sratio())
  off <- oms_fit(d, fam, f = bf(o ~ x))
  on <- oms_fit(d, fam, f = bf(o ~ x), control = oms_ten)
  expect_identical(on$fit$opt$mixture_starts[["starts"]], 10)
  expect_gte(as.numeric(logLik(on$fit)), as.numeric(logLik(off$fit)) - 1e-3)
})

test_that("a default flat by its thresholds' travel names no starts", {
  skip_on_cran()
  # cum_cum seed 10: a threshold run off or collapsed, the flag the
  # optimizer's travel decides (under autoscale a gap of 1e-10 against
  # 0), keeps the default and runs no starts, so the warning does not
  # name them
  d <- oms_data(10)
  fam <- mixture(cumulative(), cumulative())
  off <- oms_fit(d, fam)
  degen <- grep("degenerate boundary", off$msgs, value = TRUE, fixed = TRUE)
  expect_true(length(degen) > 0L)
  on <- oms_fit(d, fam, control = oms_ten)
  ran <- on$fit$opt$mixture_starts[["starts"]] %||% 0
  if (ran == 0) {
    expect_false(any(grepl("frmtmb_control(starts = 10)", degen,
                           fixed = TRUE)))
  } else {
    expect_true(any(grepl("frmtmb_control(starts = 10)", degen,
                          fixed = TRUE)))
  }
})

test_that("constant weights do not change the maximum returned", {
  skip_on_cran()
  # review of 2026-11-18, item 2: gains were compared on the weighted
  # objective, so weights 1/1000 dropped the swap on 5 of 6 seeds and
  # the estimates moved by up to 4.39, with no warning. cum_sratio
  # n = 1000 seed 13, whose swap is kept
  d <- oms_data(13, n = 1000)
  fam <- mixture(cumulative(), sratio())
  got <- lapply(c(1 / 1000, 1, 10), function(w) {
    d$w <- w
    f <- oms_fit(d, fam, f = bf(y | weights(w) ~ x + z))$fit
    list(ll = as.numeric(logLik(f)) / w,
         swapped = f$opt$mixture_starts[["swapped"]])
  })
  # what the swap buys at weights 1, against the fit without it: the
  # per-unit logLik of every weighting must agree far inside that gain.
  # The parameters are not compared: along a flat direction they move
  # by up to 30 at one likelihood (short check of 2026-11-18)
  d$w <- 1
  off <- oms_fit(d, fam, f = bf(y | weights(w) ~ x + z),
                 control = frmtmb_control(starts = 0))$fit
  gain <- got[[2L]]$ll - as.numeric(logLik(off))
  expect_gt(gain, 0)
  for (g in got[-2L]) {
    expect_identical(g$swapped, got[[2L]]$swapped)
    expect_lt(abs(g$ll - got[[2L]]$ll), 1e-3 * gain)
  }
})

test_that("a swap moves the weights against their reference", {
  skip_on_cran()
  # review of 2026-11-18, item 4: theta2 ~ 1 moves the reference to
  # component 1, and the swapped start left theta2 as it was
  d <- oms_data(13, n = 1000)
  f <- oms_fit(d, mixture(cumulative(), sratio()),
               f = bf(y ~ x + z, theta2 ~ 1),
               control = frmtmb_control(starts = 0))$fit
  sw <- frmtmb:::mixture_ord_swap_starts(f, "y")
  expect_length(sw, 1L)
  p <- f$opt$par
  nm <- frmtmb:::outer_par_map(f)$names
  expect_identical(unname(sw[[1L]]$par[nm == "theta2_(Intercept)"]),
                   -unname(p[nm == "theta2_(Intercept)"]))
  # K = 3 with the reference on component 1: new eta_k = eta_s(k) -
  # eta_s(1) for the assignment s
  set.seed(3)
  n <- 600
  d3 <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- sample.int(3L, n, replace = TRUE)
  lat <- 2 * c(1.5, -0.8, 0.3)[cls] * d3$x + 2 * c(1, -1, 0)[cls] +
    stats::rlogis(n)
  d3$y <- 1L + (lat > -1) + (lat > 0) + (lat > 1)
  f3 <- suppressMessages(allow_warnings(
    frm(bf(y ~ x, theta2 ~ 1, theta3 ~ 1),
        family = mixture(cumulative(), sratio(), acat()), data = d3,
        control = frmtmb_control(starts = 0)),
    c("degenerate boundary", "converge", "Standard errors")))
  sw3 <- frmtmb:::mixture_ord_swap_starts(f3, "y")
  expect_length(sw3, 5L)
  p3 <- f3$opt$par
  nm3 <- frmtmb:::outer_par_map(f3)$names
  eta <- c(0, p3[nm3 == "theta2_(Intercept)"], p3[nm3 == "theta3_(Intercept)"])
  for (s in sw3) {
    sg <- s$sigma
    want <- eta[sg] - eta[sg[1L]]
    expect_equal(unname(s$par[nm3 == "theta2_(Intercept)"]), unname(want[2L]))
    expect_equal(unname(s$par[nm3 == "theta3_(Intercept)"]), unname(want[3L]))
  }
})

test_that("a swap moves the group sd entries with their component", {
  skip_on_cran()
  # the one unit check of a swapped start with (1 | g)
  d <- oms_data(13)
  d$g <- factor(rep_len(1:20, nrow(d)))
  f <- suppressMessages(allow_warnings(
    frm(bf(y ~ x + z + (1 | g)), family = mixture(cumulative(), sratio()),
        data = d, control = frmtmb_control(starts = 0)),
    c("degenerate boundary", "Standard errors", "converge")))
  sw <- frmtmb:::mixture_ord_swap_starts(f, "y")
  expect_length(sw, 1L)
  p <- f$opt$par
  th <- unname(p[names(p) == "theta"])
  expect_length(th, 2L)
  expect_identical(unname(sw[[1L]]$par[names(p) == "theta"]), rev(th))
  # and the modes with them, the inner problem's start
  bl <- f$frame$re_blocks
  b0 <- as.numeric(f$estimates$b)
  b1 <- sw[[1L]]$b
  expect_identical(b1[bl[[1L]]$b_idx], b0[bl[[2L]]$b_idx])
  expect_identical(b1[bl[[2L]]$b_idx], b0[bl[[1L]]$b_idx])
})

test_that("a kept swap is a message, not a warning", {
  sig <- frmtmb:::mixture_swap_signal
  got <- NULL
  w <- NULL
  withCallingHandlers(sig("swapped"),
                      message = function(m) {
                        got <<- m
                        invokeRestart("muffleMessage")
                      },
                      warning = function(x) {
                        w <<- x
                        invokeRestart("muffleWarning")
                      })
  expect_null(w)
  expect_s3_class(got, "frmtmb_mixture_swap")
  expect_s3_class(got, "message")
})

test_that("frmtmb_control(starts = ) takes a count or NULL", {
  expect_null(frmtmb_control()$starts)
  expect_identical(frmtmb_control(starts = 3)$starts, 3L)
  for (bad in list(-1, 1.5, "a", c(1, 2), NA)) {
    expect_error(frmtmb_control(starts = bad), "`starts` must be")
  }
})

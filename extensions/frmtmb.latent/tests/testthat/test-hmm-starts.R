## hmm_starts(): refits from jittered starting values.
##
## The construction at the top is probe D4 of dev/hmm-feasibility.md,
## the measurement this function exists for. Its point is not that
## hmm_starts() works, it is that WITHOUT it the fit is wrong and
## nothing says so: the cold start converges to a local optimum with
## `convergence == 0`, a positive definite Hessian, a gradient far
## inside the convergence test, and diagnose() reporting no problem at
## all. The first test asserts that clean-looking wrong answer, so that
## the second test's recovery is a recovery from something.

# probe D4's data, exactly: K = 2, 25 sequences of 30, stationary
# initial distribution, an independent random intercept per state mean.
d4_data <- function() {
  set.seed(2026)
  K <- 2L
  N <- 25L
  Tg <- 30L
  G <- matrix(c(0.85, 0.15, 0.20, 0.80), 2, 2, byrow = TRUE)
  mu <- c(0, 3)
  sg <- c(0.6, 0.6)
  sd_b <- c(0.7, 0.5)
  stat <- local({
    A <- rbind(t(diag(K) - G), 1)
    drop(qr.solve(A, c(rep(0, K), 1)))
  })
  b <- cbind(stats::rnorm(N, 0, sd_b[1L]), stats::rnorm(N, 0, sd_b[2L]))
  d <- do.call(rbind, lapply(seq_len(N), function(g) {
    s <- integer(Tg)
    s[1L] <- sample.int(K, 1L, prob = stat)
    for (t in seq_len(Tg - 1L)) {
      s[t + 1L] <- sample.int(K, 1L, prob = G[s[t], ])
    }
    data.frame(ID = g, t = seq_len(Tg),
               y = stats::rnorm(Tg, mu[s] + b[g, s], sg[s]))
  }))
  d$gf <- factor(d$ID)
  d
}

d4_form <- function() bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
d4_fam <- function() {
  hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")
}

# a small well-behaved chain, for the argument and plumbing tests
small_data <- function(seed = 4501L, N = 10L, Tl = 20L) {
  set.seed(seed)
  G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(N), function(id) {
    s <- integer(Tl)
    s[1L] <- 1L
    for (t in seq_len(Tl)[-1L]) {
      s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
    }
    data.frame(id = id, t = seq_len(Tl),
               y = stats::rnorm(Tl, c(0, 3)[s], 0.6))
  }))
}

## ---- the failure the function exists for ----------------------------

test_that("the cold start reaches a local optimum, diagnostics clean", {
  skip_on_cran()
  dd <- d4_data()
  fit <- frm(d4_form(), family = d4_fam(), data = dd)
  ll <- as.numeric(logLik(fit))
  dg <- frmtmb::diagnose(fit, quiet = TRUE)

  # every convergence test the package has, passing
  expect_identical(dg$convergence, 0L)
  expect_true(isTRUE(dg$pdHess))
  expect_length(dg$bad_se, 0L)
  expect_length(dg$flat, 0L)
  expect_lt(dg$max_grad / abs(ll), 1e-5)

  # and the answer is wrong. The reference is the hand-rolled
  # MakeADFun(random =) arbiter of probe D4, which hmmTMB reproduces to
  # the last digit; expressed as a ratio so nothing here is an absolute
  # tolerance. The gap is 8.099 log-likelihood units on a fit of
  # -1096.096, and it is 3 orders of magnitude larger than the
  # convergence test above.
  gap <- (-1087.99646521 - ll) / abs(ll)
  expect_gt(gap, 1e-3)
  expect_gt(gap / (dg$max_grad / abs(ll)), 1e3)

  # and this IS probe D4's own cold start, not merely a bad one: the
  # recorded figure is -1096.09575602 and pinning it relative to its own
  # magnitude is what makes the 8.099 on ?hmm a number a later reader
  # does not have to re-derive
  expect_lt(abs(ll - (-1096.09575602)) / abs(ll), 1e-6)
})

test_that("hmm_starts() escapes it and reports the two optima", {
  skip_on_cran()
  dd <- d4_data()
  fit <- frm(d4_form(), family = d4_fam(), data = dd)
  ll0 <- as.numeric(logLik(fit))
  ms <- hmm_starts(fit, n = 6, jitter = 2, seed = 4101)

  expect_s3_class(ms, "frmtmb_hmm_starts")
  llb <- as.numeric(logLik(ms$best))
  # the recovered optimum beats the cold start by far more than the
  # tolerance BOTH were declared converged at
  expect_gt((llb - ll0) / abs(ll0), ms$grad_tol)
  # and the surface really has more than one optimum, found by the run
  expect_gte(nrow(ms$modes), 2L)
  # the reported spread is the spread of the OPTIMA, not of the starts
  keep <- ms$table$status %in% c("original", "converged")
  expect_equal(ms$spread, diff(range(ms$table$logLik[keep])))
  expect_equal(ms$spread, llb - ll0)
  expect_gte(ms$spread_all, ms$spread)
  # the best fit is a usable fit, at the parameters it claims
  expect_s3_class(ms$best, "frmtmb_fit")
  expect_gt(as.numeric(logLik(ms$best)), ll0)
})

## ---- what happens when a refit does not converge --------------------

test_that("a refit that does not converge is counted and never wins", {
  skip_on_cran()
  dd <- small_data()
  # two optimizer iterations: every refit stops far from any optimum,
  # so the not-converged branch is exercised by construction rather
  # than waited for
  ctl <- frmtmb_control(optCtrl = list(iter.max = 2L, eval.max = 4L))
  fit <- suppressWarnings(frm(
    bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
    data = dd, control = ctl))
  ms <- suppressWarnings(hmm_starts(fit, n = 4, jitter = 2, seed = 7L))

  expect_gt(ms$n_not_converged, 0L)
  expect_identical(ms$n_not_converged + ms$n_converged + ms$n_error,
                   as.integer(ms$n))
  # the converged spread is taken over the converged rows only, and a
  # spread over ONE value is NA rather than 0: a comparison that was
  # never made must not print the same as perfect agreement
  keep <- ms$table$status %in% c("original", "converged")
  if (sum(keep) > 1L) {
    expect_equal(ms$spread, diff(range(ms$table$logLik[keep])))
  } else {
    expect_identical(ms$spread, NA_real_)
  }
  # and every finished row, converged or not, is inside the other one
  fin <- ms$table$status != "error"
  expect_equal(ms$spread_all, diff(range(ms$table$logLik[fin])))
  if (sum(keep) > 1L) expect_gte(ms$spread_all, ms$spread)
  # a refit that did not converge cannot be the answer, even when its
  # log-likelihood is the highest one in the table
  nc <- ms$table$logLik[ms$table$status == "not converged"]
  if (length(nc) && max(nc) > max(ms$table$logLik[keep])) {
    expect_lt(as.numeric(logLik(ms$best)), max(nc))
  }
  expect_true(as.numeric(logLik(ms$best)) %in% ms$table$logLik[keep])
  # the incumbent is held to the SAME test the refits are, and the
  # summary says so rather than quietly exempting it
  expect_false(ms$original_converged)
  expect_match(paste(utils::capture.output(print(ms)), collapse = " "),
               "ORIGINAL fit does not meet that test")
})

test_that("a unimodal fit is not told it found a local optimum", {
  skip_on_cran()
  # THE FALSE-ALARM CASE. Six independent well-behaved chains, one
  # hmm_starts() each. Every refit reaches the same optimum, so the
  # summary must say the original was the best found. Before the one
  # shared tolerance went in, `gap > 0` and `lli > logLik(best)` carried
  # none, and this printed "the original fit found a local optimum" on
  # all six, on gaps of 6e-11 to 1.5e-09, while the modes table three
  # lines below correctly said ONE optimum.
  n_false <- 0L
  n_moved <- 0L
  for (k in 1:6) {
    dd <- small_data(seed = 4500L + k)
    fit <- frm(bf(y ~ 1),
               family = hmm(K = 2, gaussian(), time = t, group = id),
               data = dd)
    ms <- suppressWarnings(
      hmm_starts(fit, n = 8, jitter = 2, seed = 8800L + k))
    out <- paste(utils::capture.output(print(ms)), collapse = " ")
    expect_identical(nrow(ms$modes), 1L)
    n_false <- n_false + grepl("found a local optimum", out)
    n_moved <- n_moved +
      !identical(as.numeric(logLik(ms$best)), ms$original_logLik)
  }
  expect_identical(n_false, 0L)
  expect_identical(n_moved, 0L)
})

# The invariant the three sites exist to keep, written once so all the
# tests below assert the same thing.
#
# IT USES NEITHER SITE'S THRESHOLD AS ITS YARDSTICK, deliberately. A
# first version compared the modes against `$mode_tol`, which is the
# unified value, and so could not see the state it exists to catch: a
# build whose modes table separates two optima while its verdict treats
# them as one still looks consistent when the modes are re-judged by
# the verdict's threshold. What is asked instead is what the TABLE
# itself asserts: if the table names the original as a distinct optimum
# AND names a better one beside it, the verdict must say the original
# found a local optimum. The only tolerance here is float equality.
starts_consistent <- function(ms) {
  out <- paste(utils::capture.output(print(ms)), collapse = " ")
  says_local <- grepl("found a local optimum", out)
  eps <- 1e-9 * max(abs(ms$original_logLik), 1)
  mm <- ms$modes$logLik
  orig_is_a_mode <- any(abs(mm - ms$original_logLik) <= eps)
  better <- orig_is_a_mode && any(mm > ms$original_logLik + eps)
  list(ok = identical(says_local, better),
       says_local = says_local, better = better,
       modes = nrow(ms$modes), out = out)
}

test_that("the three tolerances are one value", {
  skip_on_cran()
  dd <- small_data()
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  ms <- hmm_starts(fit, n = 4, jitter = 1, seed = 31L)
  # the threshold is a field, not three recomputations
  expect_true(is.numeric(ms$mode_tol) && length(ms$mode_tol) == 1L)
  expect_equal(ms$mode_tol,
               ms$grad_tol^2 * max(abs(ms$original_logLik), 1))
  ck <- starts_consistent(ms)
  expect_true(ck$ok)
  # and `$best` moved exactly when the table says it should have
  expect_identical(ck$better,
                   as.numeric(logLik(ms$best)) > ms$original_logLik)
})

test_that("no grad_tol makes the printout contradict itself", {
  skip_on_cran()
  # THE BAND. When the three sites called one helper with three
  # different references, the merge threshold was always the smaller
  # (a better optimum has the smaller magnitude for a negative
  # log-likelihood), so a band of grad_tol existed where the modes
  # table said two optima and the verdict said the incumbent was best.
  # On this fit that band is around grad_tol = 0.086, where the
  # threshold has grown to the size of the 8.099-unit gap itself. The
  # values below straddle it.
  dd <- d4_data()
  fit <- frm(d4_form(), family = d4_fam(), data = dd)
  for (gt in c(1e-3, 3e-2, 0.08, 0.085, 0.086, 0.09, 0.1)) {
    ms <- suppressWarnings(
      hmm_starts(fit, n = 3, jitter = 2, seed = 4101, grad_tol = gt))
    ck <- starts_consistent(ms)
    expect_true(ck$ok,
                info = paste0("grad_tol ", gt, ": modes ", ck$modes,
                              ", says_local ", ck$says_local,
                              ", better ", ck$better))
  }
})

test_that("a refit that errors is counted and is in neither spread", {
  skip_on_cran()
  dd <- small_data()
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  # a start hundreds of standard errors away puts sigma at exp(huge)
  # and the objective is not finite there
  ms <- suppressWarnings(hmm_starts(fit, n = 4, jitter = 1e4, seed = 3L))
  skip_if(ms$n_error == 0L,
          "no refit errored at this jitter; nothing to assert here")
  expect_identical(nrow(ms$table), as.integer(ms$n) + 1L)
  expect_identical(sum(is.na(ms$table$logLik)), as.integer(ms$n_error))
  fin <- ms$table$status != "error"
  if (sum(fin) > 1L) {
    expect_equal(ms$spread_all, diff(range(ms$table$logLik[fin])))
  } else {
    expect_identical(ms$spread_all, NA_real_)
  }
})

test_that("a run where every refit errored reports NA, not zero", {
  skip_on_cran()
  dd <- small_data()
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  ms <- suppressWarnings(hmm_starts(fit, n = 3, jitter = 1e6, seed = 9L))
  skip_if(ms$n_error < as.integer(ms$n),
          "not every refit errored at this jitter")
  # nothing was compared, so there is no spread. A zero here reads as
  # perfect agreement to a caller who never prints the object.
  expect_identical(ms$spread, NA_real_)
  expect_identical(ms$spread_all, NA_real_)
  out <- paste(utils::capture.output(print(ms)), collapse = " ")
  expect_match(out, "not measurable")
  expect_identical(as.numeric(logLik(ms$best)), ms$original_logLik)
})

## ---- plumbing -------------------------------------------------------

test_that("hmm_starts() is reproducible and leaves the random stream alone", {
  skip_on_cran()
  dd <- small_data()
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  set.seed(99)
  before <- stats::runif(1)
  seed_before <- get(".Random.seed", envir = globalenv())
  a <- hmm_starts(fit, n = 2, jitter = 1, seed = 12L)
  expect_identical(get(".Random.seed", envir = globalenv()), seed_before)
  after <- stats::runif(1)
  b <- hmm_starts(fit, n = 2, jitter = 1, seed = 12L)
  expect_equal(a$table$logLik, b$table$logLik)
  # the restored stream is the one that was interrupted
  set.seed(99)
  expect_equal(c(before, after), stats::runif(2))
})

test_that("the winning start is written into the returned fit's own call", {
  skip_on_cran()
  dd <- d4_data()
  fit <- frm(d4_form(), family = d4_fam(), data = dd)
  ms <- hmm_starts(fit, n = 6, jitter = 2, seed = 4101)
  skip_if(identical(as.numeric(logLik(ms$best)),
                    as.numeric(logLik(fit))),
          "no refit beat the original, so no call was rewritten")
  st <- ms$best$call[["start"]]
  expect_type(st, "list")
  again <- eval(ms$best$call)
  expect_equal(as.numeric(logLik(again)), as.numeric(logLik(ms$best)),
               tolerance = 1e-8)
})

test_that("keep = TRUE returns the refits and the default does not", {
  skip_on_cran()
  dd <- small_data()
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  a <- hmm_starts(fit, n = 2, jitter = 1, seed = 5L)
  expect_null(a$fits)
  b <- hmm_starts(fit, n = 2, jitter = 1, seed = 5L, keep = TRUE)
  expect_length(b$fits, 2L)
  for (f in b$fits) expect_s3_class(f, "frmtmb_fit")
})

test_that("the table carries a load-independent cost per refit", {
  skip_on_cran()
  dd <- small_data()
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  a <- hmm_starts(fit, n = 2, jitter = 1, seed = 21L)
  expect_true(all(c("fn_evals", "gr_evals") %in% names(a$table)))
  expect_true(all(a$table$gr_evals[-1L] > 0))
  # the same seed takes the same optimizer path, so the count is
  # identical across runs where the seconds are not
  b <- hmm_starts(fit, n = 2, jitter = 1, seed = 21L)
  expect_identical(a$table$gr_evals, b$table$gr_evals)
  expect_identical(a$table$fn_evals, b$table$fn_evals)
})

test_that("the scale falls back when confint and vcov disagree", {
  # THE BRANCHES `hmm_starts_scale()` CANNOT REACH THROUGH A REAL FIT.
  # A first version of this function claimed to CHECK the row order of
  # confint() against its `est` column; confint.frmtmb_fit() builds that
  # column as `object$opt$par`, so the test compared a vector with
  # itself and could never fail. What replaced it is a cross-check of
  # the standard errors against vcov(), which reaches the same
  # covariance by its own path and names its own rows, and this is
  # where that check is seen firing. A stub is used because no fit this
  # family produces can make the two disagree, which is the point.
  mk <- function(se_ci, se_vc, npar = 2L) {
    nm <- c("mu1_(Intercept)", "mu2_(Intercept)")[seq_len(npar)]
    structure(list(opt = list(par = stats::setNames(rep(0, npar),
                                                    rep("beta", npar))),
                   ci_se = se_ci, vc_se = se_vc, nm = nm),
              class = "latent_scale_stub")
  }
  q <- stats::qnorm(0.975)
  registerS3method("confint", "latent_scale_stub", function(object, ...) {
    if (is.null(object$ci_se)) stop("no confint for this object")
    cbind(lwr = -q * object$ci_se, upr = q * object$ci_se,
          est = rep(0, length(object$ci_se)))[, , drop = FALSE] |>
      `rownames<-`(object$nm)
  })
  registerS3method("vcov", "latent_scale_stub", function(object, ...) {
    if (is.null(object$vc_se)) stop("no vcov for this object")
    v <- diag(object$vc_se^2, nrow = length(object$vc_se))
    dimnames(v) <- list(object$nm, object$nm)
    v
  })
  f <- getFromNamespace("hmm_starts_scale", "frmtmb.latent")

  # agreeing: the standard errors are used
  ok <- f(mk(c(0.5, 0.25), c(0.5, 0.25)))
  expect_identical(ok$how, "se")
  expect_equal(ok$s, c(0.5, 0.25))

  # DISAGREEING, which is the case a real fit cannot produce: the two
  # paths name the same coefficients and report different errors
  bad <- f(mk(c(0.5, 0.25), c(0.25, 0.5)))
  expect_identical(bad$how, "unit")
  expect_equal(bad$s, c(1, 1))

  # confint() unavailable at all, the only entrance a real fit has
  none <- f(mk(NULL, c(0.5, 0.25)))
  expect_identical(none$how, "unit")
  # and vcov() unavailable, the second new entrance
  novc <- f(mk(c(0.5, 0.25), NULL))
  expect_identical(novc$how, "unit")
})

test_that("the perturbation scale is the fit's own standard errors", {
  skip_on_cran()
  dd <- small_data()
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  ms <- hmm_starts(fit, n = 1, jitter = 1, seed = 5L)
  expect_identical(ms$scale_how, "se")
})

test_that("print() names the counts, both spreads and the modes", {
  skip_on_cran()
  dd <- small_data()
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  ms <- hmm_starts(fit, n = 2, jitter = 1, seed = 5L)
  out <- paste(utils::capture.output(print(ms)), collapse = "\n")
  expect_match(out, "converged")
  expect_match(out, "logLik spread, the original")
  expect_match(out, "every one of the")
  expect_match(out, "distinct optima")
  expect_match(out, "seconds")
})

## ---- refusals -------------------------------------------------------

test_that("hmm_starts() refuses a fit that is not an hmm", {
  set.seed(31)
  dd <- data.frame(x = stats::rnorm(60))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  expect_error(hmm_starts(fit), "does not use an hmm\\(\\) family")
})

test_that("hmm_starts() refuses arguments it cannot use", {
  skip_on_cran()
  dd <- small_data()
  fit <- frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd)
  expect_error(hmm_starts(fit, n = 0), "whole number of refits")
  expect_error(hmm_starts(fit, n = 2.5), "whole number of refits")
  expect_error(hmm_starts(fit, n = c(1, 2)), "whole number of refits")
  expect_error(hmm_starts(fit, jitter = 0), "one positive number")
  expect_error(hmm_starts(fit, jitter = -1), "one positive number")
  expect_error(hmm_starts(fit, grad_tol = 0), "one positive number")
  expect_error(hmm_starts(fit, keep = "yes"), "must be TRUE or FALSE")
})

test_that("a fit whose data the caller cannot see is refused by name", {
  skip_on_cran()
  # the call names its data by SYMBOL, so a fit assembled in a scope
  # that has since gone has to fail once, loudly, rather than n times
  # as n bad starts
  fit <- local({
    dd_gone <- small_data()
    frm(bf(y ~ 1),
        family = hmm(K = 2, gaussian(), time = t, group = id),
        data = dd_gone)
  })
  expect_error(hmm_starts(fit, n = 1),
               "could not be re-evaluated from here")
})

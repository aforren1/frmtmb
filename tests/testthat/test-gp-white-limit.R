# An approximate gp() whose length scale ran under the spacing its basis
# resolves (frmtmb_gp_basis). Where the likelihood stays level from there to
# a length scale of 0, the white limit, sdgp and the length scales lie on a
# flat ridge (only sdgp^2 times their product is determined): a boundary
# fit, told by a message. Where it does not (a prior on lscale, a stop on
# the way there), the warning stays.

# Rows 3 (white) and 6 (resolved) of the gpfield review's 30 gaussian fits
# (dev/reviews/gpfield-review-basis.R): f from an exact squared exponential
# GP of length scale 0.7 on [0, 10], noise sd 0.3, 60 rows.
gpw_data <- function(row) {
  set.seed(9600 + row)
  x <- sort(stats::runif(60, 0, 10))
  K <- exp(-outer(x, x, "-")^2 / (2 * 0.7^2)) + diag(1e-8, 60)
  f <- as.vector(t(chol(K)) %*% stats::rnorm(60))
  data.frame(x = x, y = f + 0.3 * stats::rnorm(60))
}

# The fit, with the frmtmb_gp_basis condition it gave (a message is
# muffled; a warning is not, so that one nobody expected escapes).
gpw_fit <- function(...) {
  got <- list()
  fit <- withCallingHandlers(frm(...), frmtmb_gp_basis = function(cnd) {
    got[[length(got) + 1L]] <<- cnd
    if (inherits(cnd, "message")) invokeRestart("muffleMessage")
  })
  list(fit = fit, basis = got)
}

gpw_summary_text <- function(fit) {
  paste(utils::capture.output(print(summary(fit))), collapse = " ")
}

test_that("a field proved at its basis's white limit is a boundary fit", {
  r <- gpw_fit(bf(y ~ gp(x, k = 10)), data = gpw_data(3),
               family = gaussian())
  fit <- r$fit
  bk <- fit$frame$re_blocks[[1L]]
  expect_true(frmtmb:::gp_basis_white(bk, fit$estimates$theta))
  # one condition, a message of the boundary class, not a warning
  expect_length(r$basis, 1L)
  cnd <- r$basis[[1L]]
  expect_s3_class(cnd, "message")
  expect_s3_class(cnd, "frmtmb_boundary_fit")
  expect_match(conditionMessage(cnd), "white limit of its basis",
               fixed = TRUE)
  # sdgp and lscale have no standard error, and the summary says why; the
  # intercept and sigma keep theirs
  se <- sqrt(diag(frmtmb:::sdr_of(fit)$cov.fixed))
  th <- names(se) == "theta"
  expect_true(all(is.nan(se[th])))
  expect_true(all(is.finite(se[!th])))
  expect_identical(unname(frmtmb:::sdr_of(fit)$se_lost),
                   c("boundary", "boundary"))
  expect_match(gpw_summary_text(fit), "at the white limit of its basis",
               fixed = TRUE)
  # the product the data determine, as the message names it
  h <- regmatches(conditionMessage(cnd),
                  regexpr("\"[^\"]*\"", conditionMessage(cnd)))
  expect_length(h, 1L)
  hy <- hypothesis(fit, gsub("\"", "", h), class = "")$hypothesis
  bv <- frmtmb:::gp_brms_values(fit, fit$estimates$theta)
  expect_true(is.finite(hy$Est.Error) && hy$Est.Error > 0)
  expect_equal(hy$Estimate, 2 * log(bv[["sdgp_gpx"]]) +
                 log(bv[["lscale_gpx"]]))
})

test_that("a field whose ridge is not level keeps the warning", {
  d <- gpw_data(3)
  r <- gpw_fit(bf(y ~ gp(x, k = 10)), data = d, family = gaussian())
  fit <- r$fit
  bk <- fit$frame$re_blocks[[1L]]
  sp <- bk$gp_spacing
  # a prior holds the length scale at 0.05 of the spacing, an interior
  # maximum under the gate: the ridge rises both ways
  pr <- set_prior(sprintf("lognormal(%.10g, 0.05)", log(0.05 * sp)),
                  class = "lscale")
  rp <- allow_warnings(
    gpw_fit(bf(y ~ gp(x, k = 10)), data = d, family = gaussian(),
            prior = pr),
    "spacing its basis resolves", require = "spacing its basis resolves")
  expect_true(all(vapply(rp$basis, inherits, NA, "warning")))
  expect_true(frmtmb:::gp_basis_white(rp$fit$frame$re_blocks[[1L]],
                                      rp$fit$estimates$theta))
  se <- sqrt(diag(frmtmb:::sdr_of(rp$fit)$cov.fixed))
  expect_true(all(is.finite(se)))
  # the white fit moved back along its own ridge to 0.01 of the spacing,
  # where the likelihood still rises toward the limit (by about 1e-4,
  # dev/gpselimit-findings.md): no proof there
  om <- frmtmb:::outer_par_map(fit)
  pos <- which(om$comp == "theta")[bk$theta_idx]
  t <- log(0.01 * sp) - fit$estimates$theta[bk$theta_idx[2L]]
  mv <- c(-t / 2, t)
  f2 <- fit
  f2$cache <- new.env()
  f2$opt$par[pos] <- f2$opt$par[pos] + mv
  f2$estimates$theta[bk$theta_idx] <- f2$estimates$theta[bk$theta_idx] + mv
  expect_true(frmtmb:::gp_basis_white(bk, f2$estimates$theta))
  expect_false(frmtmb:::gp_white_ridge(f2, bk))
  expect_true(frmtmb:::gp_white_ridge(fit, bk))
})

test_that("a resolved field says nothing, and check_se governs the message", {
  r <- gpw_fit(bf(y ~ gp(x, k = 30)), data = gpw_data(6),
               family = gaussian())
  expect_length(r$basis, 0L)
  expect_false(frmtmb:::gp_basis_white(r$fit$frame$re_blocks[[1L]],
                                       r$fit$estimates$theta))
  expect_true(all(is.finite(sqrt(diag(
    frmtmb:::sdr_of(r$fit)$cov.fixed)))))
  d <- gpw_data(3)
  ri <- gpw_fit(bf(y ~ gp(x, k = 10)), data = d, family = gaussian(),
                control = frmtmb_control(check_se = "ignore"))
  expect_length(ri$basis, 0L)
  expect_error(frm(bf(y ~ gp(x, k = 10)), data = d, family = gaussian(),
                   control = frmtmb_control(check_se = "stop")),
               "white limit of its basis", fixed = TRUE)
})

test_that("a stop in the white plateau short of a maximum keeps the warning", {
  # row 6 is resolved (length scale 1.06 of the spacing); started deep in
  # the white plateau the optimizer stops there, 5.18 short, where every
  # gradient is about 0 and the ridge toward 0 is level: only the way back
  # tells it from the limit (dev/reviews/2026-11-17-gpselimit.md, finding
  # 2). From 1e-3 of the spacing it stops where the ridge still rises
  # toward 0. Neither is a boundary, and neither loses its SEs as one
  d <- gpw_data(6)
  f0 <- gpw_fit(bf(y ~ gp(x, k = 30)), data = d, family = gaussian())$fit
  sp <- f0$frame$re_blocks[[1L]]$gp_spacing
  ll0 <- as.numeric(logLik(f0))
  for (st in c(3e-5, 1e-3)) {
    r <- allow_warnings(
      gpw_fit(bf(y ~ gp(x, k = 30)), data = d, family = gaussian(),
              start = list(theta = c(log(30), log(st * sp)))),
      c("spacing its basis resolves", "Standard errors are not available"))
    fit <- r$fit
    bk <- fit$frame$re_blocks[[1L]]
    if (!frmtmb:::gp_basis_white(bk, fit$estimates$theta)) {
      # a rounding that leaves the plateau reaches the maximum: nothing to
      # say then
      expect_length(r$basis, 0L)
      next
    }
    expect_lt(as.numeric(logLik(fit)), ll0)
    expect_length(r$basis, 1L)
    expect_s3_class(r$basis[[1L]], "warning")
    expect_match(conditionMessage(r$basis[[1L]]),
                 "stopped short of its maximum", fixed = TRUE)
    lost <- frmtmb:::sdr_of(fit)$se_lost
    expect_false(any(lost[grepl("^theta", names(lost))] %in% "boundary"))
  }
})

test_that("the white sub-GPs of gp(by = f) have no standard error", {
  # 0.83.0 kept finite SEs on two of the four white hyperparameters here
  # (and on all four on sweep row 21, sdgp 197.7), with no word: the
  # Hessian did not remove the ridges (review finding 3)
  d <- gpw_data(3)
  d$f <- factor(rep(c("a", "b"), length.out = nrow(d)))
  r <- gpw_fit(y ~ gp(x, by = f, k = 10), data = d)
  fit <- r$fit
  se <- sqrt(diag(frmtmb:::sdr_of(fit)$cov.fixed))
  nm <- frmtmb:::outer_par_names(fit)
  thp <- which(names(se) == "theta")
  wh <- integer(0)
  for (bk in fit$frame$re_blocks) {
    if (frmtmb:::gp_basis_white(bk, fit$estimates$theta[bk$theta_idx])) {
      wh <- c(wh, thp[bk$theta_idx])
    }
  }
  expect_gt(length(wh), 0L)
  expect_true(all(is.nan(se[wh])))
  expect_true(all(frmtmb:::sdr_of(fit)$se_lost[nm[wh]] %in% "boundary"))
  expect_true(all(is.finite(se[setdiff(seq_along(se), thp)])))
  expect_true(all(vapply(r$basis, inherits, NA, "message")))
})

test_that("a white field at a doubted stop keeps the warning", {
  # a field in a nonlinear parameter whose fit stops at code 1, false
  # convergence, beside a white field (review finding 4, seed 2): the
  # ridge is level with the rest held at the stop, but the stop is not
  # known to be a maximum, so it is no boundary fit
  set.seed(52302)
  n <- 120L
  x <- sort(round(stats::runif(n, 0, 10), 3))
  z <- round(stats::rnorm(n), 3)
  f <- 2 * sin(3 * x) + x
  z <- abs(z) * 2
  y <- (2 + f / max(abs(f))) * (1 - exp(0.8 * z)) + 0.1 * stats::rnorm(n)
  d <- data.frame(x, z, y)
  seen <- character(0)
  basis <- list()
  allowed <- c("frmtmb_conv_code", "frmtmb_conv_gradient",
               "frmtmb_conv_gradient_nonfinite", "frmtmb_gp_basis")
  fit <- withCallingHandlers(
    frm(bf(y ~ a * (1 - exp(-b * z)), a ~ gp(x, k = 20), b ~ 1, nl = TRUE),
        data = d, family = gaussian()),
    frmtmb_gp_basis = function(cnd) basis[[length(basis) + 1L]] <<- cnd,
    warning = function(w) {
      if (class(w)[1L] %in% allowed) {
        seen <<- c(seen, class(w)[1L])
        invokeRestart("muffleWarning")
      }
    },
    message = function(m) invokeRestart("muffleMessage"))
  expect_s3_class(fit, "frmtmb_fit")
  # on a rounding without a convergence warning either condition is right
  if (any(startsWith(seen, "frmtmb_conv_"))) {
    expect_gt(length(basis), 0L)
    expect_true(all(vapply(basis, inherits, NA, "warning")))
  }
})

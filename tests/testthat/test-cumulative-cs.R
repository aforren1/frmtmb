# cs() on cumulative(), which brms 2.23.0 fits with a warning and which
# frmtmb refused up to 0.67.0 (lifted at 0.68.0, the user's decision of
# 2026-10-06). Row i reads the thresholds tau_k - cs_ik, brms's
# `Intercept - transpose(mucs[n])`, and a row whose offsets cross two
# thresholds has a negative category probability: its density is NaN,
# as brms's is, its fitted() probabilities NaN and its draw NA. The
# Stan-side identity is dev/rel068-cs-lpcheck.R (brms's compiled
# program, six shapes, at most 0.9 ulp).

ccs_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  # threshold-specific slopes small against the gaps, so that no row's
  # thresholds cross at the fit
  u <- stats::rlogis(n) + 0.5 * d$x + 0.4 * d$z
  d$y <- 1L + (u > -1 + 0.15 * d$x) + (u > 0.3) + (u > 1.5 - 0.15 * d$x)
  d$yh <- ifelse(stats::runif(n) < 0.2, 0L, d$y)
  d
}

ccs_warn <- paste("Category specific effects for this family should be",
                  "considered experimental")

# the n x K category probabilities, written out from the thresholds
ccs_probs <- function(eta, tau, cs, disc, cdf) {
  n <- length(eta)
  M <- disc * (matrix(tau, n, length(tau), byrow = TRUE) - cs - eta)
  Fm <- cbind(0, cdf(M), 1)
  Fm[, -1L] - Fm[, -ncol(Fm)]
}

test_that("cs() on cumulative() fits, with brms's warning", {
  d <- ccs_data(20261006)
  for (lk in c("logit", "probit")) {
    fit <- allow_warnings(
      frm(bf(y ~ z + cs(x), disc ~ 0 + z), family = cumulative(lk),
          data = d),
      ccs_warn, require = ccs_warn)
    est <- fit$estimates
    lp <- Filter(function(l) identical(l$dpar, "mu"),
                 fit$frame$linpreds)[[1L]]
    tau <- c(est$tau_raw[1], est$tau_raw[1] + cumsum(exp(est$tau_raw[-1])))
    cs <- outer(d$x, est[[lp$cs[[1L]]$par]])
    eta <- est$beta[["z"]] * d$z
    disc <- exp(est$betad[["disc_z"]] * d$z)
    P <- ccs_probs(eta, tau, cs, disc,
                   if (lk == "logit") stats::plogis else stats::pnorm)
    expect_true(all(P > 0))
    ref <- sum(log(P[cbind(seq_len(nrow(d)), d$y)]))
    expect_lt(abs(-fit$obj$fn(fit$opt$par) - ref),
              1e3 * .Machine$double.eps * abs(ref), label = lk)
    expect_lt(max(abs(unname(fitted(fit)[, "Estimate", ]) - P)),
              1e3 * .Machine$double.eps, label = lk)
    # brms's names for the per-threshold coefficients
    expect_true(all(paste0("bcs_x[", 1:3, "]") %in% variables(fit)))
  }
})

test_that("the density is brms's own R-side density", {
  skip_unless_brms()
  d <- ccs_data(20261007)
  fit <- allow_warnings(frm(y ~ cs(x), family = cumulative(), data = d),
                        ccs_warn, require = ccs_warn)
  est <- fit$estimates
  lp <- Filter(function(l) identical(l$dpar, "mu"),
               fit$frame$linpreds)[[1L]]
  tau <- c(est$tau_raw[1], est$tau_raw[1] + cumsum(exp(est$tau_raw[-1])))
  thr <- matrix(tau, nrow(d), 3, byrow = TRUE) -
    outer(d$x, est[[lp$cs[[1L]]$par]])
  P <- brms:::dcumulative(1:4, eta = rep(0, nrow(d)), thres = thr,
                          disc = 1, link = "logit")
  ref <- sum(log(P[cbind(seq_len(nrow(d)), d$y)]))
  expect_lt(abs(as.numeric(logLik(fit)) - ref),
            1e3 * .Machine$double.eps * abs(ref))
})

test_that("crossing cs() thresholds give NaN and NA, as brms's density does", {
  fam <- cumulative()
  raw <- c(-1, log(0.5), log(1))
  # row 2's offset lifts threshold 2 above threshold 3
  cs <- matrix(c(0, 0, 0, 0, -2, 0), 2, 3, byrow = TRUE)
  dp <- list(mu = c(0, 0), disc = c(1, 1), .cs = cs)
  out <- fam$lpdf(c(3, 3), dp, list(), list(tau_raw = raw))
  expect_true(is.finite(out[1]))
  expect_true(is.nan(out[2]))
  # the absent case: the same rows with no offsets are both finite
  out0 <- fam$lpdf(c(3, 3), list(mu = c(0, 0), disc = c(1, 1)), list(),
                   list(tau_raw = raw))
  expect_true(all(is.finite(out0)))
  set.seed(1)
  s <- fam$sim(dp, list(), 2L, list(tau_raw = raw))
  expect_true(!is.na(s[1]) && is.na(s[2]))
  # and on a fit, a newdata row far out on x crosses them
  d <- ccs_data(20261008)
  fit <- allow_warnings(frm(y ~ cs(x), family = cumulative(), data = d),
                        ccs_warn, require = ccs_warn)
  b <- fit$estimates[[Filter(function(l) identical(l$dpar, "mu"),
                             fit$frame$linpreds)[[1L]]$cs[[1L]]$par]]
  tau <- frmtmb:::ord_threshold_values(family(fit), fit$estimates$tau_raw)
  # the x at which thresholds 1 and 2 meet, and well past it
  x0 <- (tau[2] - tau[1]) / (b[2] - b[1])
  nd <- data.frame(x = c(0, 3 * x0))
  P <- fitted(fit, newdata = nd)[, "Estimate", ]
  expect_true(all(is.finite(P[1, ])))
  expect_true(all(is.nan(P[2, ])))
})

test_that("the one-step density reads the cs() offsets", {
  d <- ccs_data(20261009)
  fit <- allow_warnings(frm(y ~ z + cs(x), family = cumulative(), data = d),
                        ccs_warn, require = ccs_warn)
  ff <- family(fit)
  dp <- frmtmb:::eval_dpars(fit)[[1L]]
  lp <- Filter(function(l) identical(l$dpar, "mu"),
               fit$frame$linpreds)[[1L]]
  dp[[".cs"]] <- outer(d$x, fit$estimates[[lp$cs[[1L]]$par]])
  ex <- list(tau_raw = fit$estimates$tau_raw)
  osa <- methods::getClass("osa", where = asNamespace("RTMB"))
  a <- ff$lpdf(d$y, dp, list(), ex)
  b <- ff$lpdf(methods::new(osa, x = as.numeric(d$y),
                            keep = matrix(1, nrow(d), 1)), dp, list(), ex)
  expect_lt(max(abs(as.numeric(b) - a)) / max(abs(a)),
            1e3 * .Machine$double.eps)
  r <- residuals(fit, type = "osa")[, "Estimate"]
  expect_true(all(is.finite(r)))
})

test_that("the warning is brms's: per cs() predictor on the ordered families", {
  d <- ccs_data(20261010)
  seen <- function(expr) {
    w <- character()
    withCallingHandlers(expr, warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
    sum(grepl(ccs_warn, w, fixed = TRUE))
  }
  expect_identical(seen(frm(y ~ cs(x), family = cumulative(), data = d)), 1L)
  expect_identical(seen(frm(yh ~ cs(x), family = hurdle_cumulative("probit"),
                            data = d)), 1L)
  # once per cumulative() component of a mixture, as brms warns (two for
  # mixture(cumulative, cumulative), one beside sratio;
  # dev/rel068-cs-brms.R). The frame gives the warning, so the frame is
  # all these build: a cs() mixture's FIT stops on a NaN gradient on
  # some BLAS builds and not others (OpenBLAS 0.3.26 here,
  # dev/ciharden-findings.md), which says nothing about the warning
  expect_identical(seen(frm(y ~ cs(x), family = mixture(cumulative(),
                                                        sratio()),
                            data = d, dry_run = "frame")), 1L)
  expect_identical(seen(frm(y ~ cs(x), family = mixture(cumulative(),
                                                        cumulative()),
                            data = d, dry_run = "frame")), 2L)
  # absent: cs() on sratio(), and cumulative() without cs()
  expect_identical(seen(frm(y ~ cs(x), family = sratio(), data = d)), 0L)
  expect_identical(seen(frm(y ~ x, family = cumulative(), data = d)), 0L)
  # a family that takes no cs() is refused in brms's words
  d$yc <- rnorm(nrow(d))
  expect_error(frm(yc ~ cs(x), family = gaussian(), data = d),
               "Category specific effects are not supported for this family",
               fixed = TRUE)
})

test_that("predict() is NA only where the estimate crosses", {
  # predict() draws the parameters from their asymptotic law, and a
  # replicate is NA where its draw crosses two of the row's thresholds.
  # The user's decision of 2026-10-06: a row is NA only where the
  # plug-in estimate itself crosses, as fitted() is NaN there; elsewhere
  # the proportions are over the replicates that do not cross, and the
  # call warns once with the number dropped per row
  ccs_cross_x <- function(fit) {
    b <- fit$estimates[[Filter(function(l) identical(l$dpar, "mu"),
                               fit$frame$linpreds)[[1L]]$cs[[1L]]$par]]
    tau <- frmtmb:::ord_threshold_values(family(fit), fit$estimates$tau_raw)
    (tau[2] - tau[1]) / (b[2] - b[1])
  }
  d <- ccs_data(20261011)
  for (fam in list(cumulative(), hurdle_cumulative("probit"))) {
    hurdle <- identical(fam$family, "hurdle_cumulative")
    fit <- allow_warnings(frm(if (hurdle) yh ~ cs(x) else y ~ cs(x),
                              family = fam, data = d), ccs_warn)
    x0 <- ccs_cross_x(fit)
    # no crossing; near the crossing point, where some parameter draws
    # cross but the estimate does not; past it, where the estimate does
    nd <- data.frame(x = c(0, 0.9 * x0, 3 * x0))
    fe <- fitted(fit, newdata = nd)[, "Estimate", ]
    expect_identical(unname(rowSums(is.nan(fe)) > 0), c(FALSE, FALSE, TRUE))
    set.seed(2)
    raw <- predict(fit, newdata = nd, ndraws = 400, summary = FALSE)
    nas <- unname(colSums(is.na(raw)))
    expect_identical(nas[1L], 0)
    expect_gt(nas[2L], 0)
    expect_lt(nas[2L], 400)
    set.seed(2)
    p <- allow_warnings(predict(fit, newdata = nd, ndraws = 400),
                        "cs() offsets make two thresholds cross",
                        require = c(paste0("row 2 ", nas[2L]),
                                    "1 row(s) cross at the estimates"))
    expect_true(all(is.finite(p[1:2, ])), label = fam$family)
    # the proportions are over the replicates that are there
    ok <- raw[!is.na(raw[, 2L]), 2L]
    ref <- tabulate(as.integer(ok) + 1L - (if (hurdle) 0L else 1L),
                    ncol(p)) / length(ok)
    expect_identical(unname(p[2L, ]), ref, label = fam$family)
    expect_true(all(is.na(p[3L, ])), label = fam$family)
    # the absent case: no row crosses anywhere, and nothing warns
    set.seed(2)
    p0 <- predict(fit, newdata = nd[1L, , drop = FALSE], ndraws = 400)
    expect_true(all(is.finite(p0)))
  }
})

test_that("an ordered row without cs() keeps drawing a category", {
  # pnorm() is not monotone to the last ulp near +-0.6745 (7601 and 3897
  # negative steps in 2e5 adjacent doubles; dev/relrev-mono.R), so two
  # thresholds a few ulp apart there give a probability of -1.1e-16.
  # The NA a crossing row draws belongs to cs() alone; 0.67.0 drew a
  # category here, and so does this build (the release review, m7)
  fam <- cumulative("probit")
  a <- 0.67448975
  x <- a + (seq_len(200000) - 100000) * a * .Machine$double.eps
  hit <- NULL
  for (i in which(diff(stats::pnorm(x)) < 0)) {
    for (g in 1:4) {
      raw <- c(x[i], log(g * .Machine$double.eps), 0)
      tau <- frmtmb:::ord_threshold_values(fam, raw)
      if (tau[2] > tau[1] && stats::pnorm(tau[2]) < stats::pnorm(tau[1])) {
        hit <- raw
        break
      }
    }
    if (!is.null(hit)) break
  }
  expect_false(is.null(hit))
  set.seed(1)
  s <- vapply(1:200, function(j) {
    as.numeric(fam$sim(list(mu = 0, disc = 1), list(), 1L,
                       list(tau_raw = hit)))
  }, 0)
  expect_identical(sum(is.na(s)), 0L)
  expect_true(all(s %in% 1:4))
})

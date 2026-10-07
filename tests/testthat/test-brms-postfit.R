# brms post-processing calls on a maximum-likelihood fit that a ported
# script meets (lane surface, 2026-10-06): stancode(), standata(),
# pp_mixture(), add_criterion(), summary(waic =), plot()'s brms
# arguments, and the categorical refusal of conditional_effects().
#
# Seen to fail on 0.68.1 (rellib-r6): stancode() and standata() had no
# method for a fit ("no applicable method", and with brms loaded
# brms's default, "Data must be specified using the 'data' argument");
# pp_mixture() had none either; add_criterion() did not exist; plot()
# said "Did you mean `x`?" for `N` and named neither `variable` nor
# `regex`'s reason; summary(waic =) got the generic refusal; and the
# categorical refusal said "ordinal family" (dev/surface-repros.R).

fit_gauss <- function() {
  set.seed(1)
  dd <- data.frame(x = stats::rnorm(40))
  dd$y <- dd$x + stats::rnorm(40)
  frm(bf(y ~ x), family = gaussian(), data = dd)
}

test_that("stancode() and standata() on a fit refuse by name", {
  f <- fit_gauss()
  expect_error(stancode(f), "no meaning for a frmtmb fit: there is no Stan")
  expect_error(standata(f), "no meaning for a frmtmb fit: nothing is")
  # the closure the refusal names is the model
  nll <- build_objective(f$frame)
  expect_equal(nll(f$estimates), -as.numeric(logLik(f)),
               tolerance = 8 * .Machine$double.eps)
})

test_that("stancode() reaches the fit method with brms loaded", {
  skip_if_not_installed("brms")
  f <- fit_gauss()
  loadNamespace("brms")
  expect_error(brms::stancode(f), "no meaning for a frmtmb fit")
  expect_error(brms::standata(f), "no meaning for a frmtmb fit")
  expect_error(stancode(f), "no meaning for a frmtmb fit")
})

test_that("pp_mixture() on a fit is mixture_probs() in brms's array", {
  set.seed(4)
  dd <- data.frame(y = c(stats::rnorm(60, -2), stats::rnorm(60, 3)))
  fit <- frm(bf(y ~ 1), family = mixture(gaussian(), gaussian()),
             data = dd)
  pm <- pp_mixture(fit)
  expect_equal(dim(pm), c(120L, 4L, 2L))
  expect_equal(dimnames(pm)[[2L]], c("Estimate", "Est.Error", "Q2.5",
                                     "Q97.5"))
  expect_equal(dimnames(pm)[[3L]], c("P(K = 1 | Y)", "P(K = 2 | Y)"))
  expect_identical(unname(pm[, "Estimate", ]), unname(mixture_probs(fit)))
  # the band is on the logit, so it stays inside (0, 1) and holds the
  # estimate
  expect_true(all(pm[, "Q2.5", ] >= 0 & pm[, "Q97.5", ] <= 1))
  expect_true(all(pm[, "Q2.5", ] <= pm[, "Estimate", ] &
                    pm[, "Estimate", ] <= pm[, "Q97.5", ]))
  expect_true(all(is.finite(pm[, "Est.Error", ])))
  # the two components' errors are one error: P1 + P2 = 1
  expect_equal(pm[, "Est.Error", 1L], pm[, "Est.Error", 2L],
               tolerance = 1e-6)
  lg <- pp_mixture(fit, log = TRUE)
  expect_equal(lg[, "Estimate", ], log(pm[, "Estimate", ]))
  expect_equal(lg[, "Q97.5", ], log(pm[, "Q97.5", ]))
  pr <- pp_mixture(fit, probs = c(0.1, 0.9))
  expect_equal(dimnames(pr)[[2L]], c("Estimate", "Est.Error", "Q10", "Q90"))
  # brms's draws arguments, refused by name
  expect_error(pp_mixture(fit, summary = FALSE), "summary = FALSE")
  expect_error(pp_mixture(fit, robust = TRUE), "robust = TRUE")
  expect_error(pp_mixture(fit, ndraws = 10), "`ndraws`")
  expect_error(pp_mixture(fit, newdata = dd), "does not take newdata")
  expect_error(pp_mixture(fit, re_formula = NA), "does not take re_formula")
  expect_error(pp_mixture(fit, resp = "z"), "names no response")
})

test_that("pp_mixture()'s Est.Error is the SD of its interval's law", {
  # The review of 2026-10-07, B2: the delta-method error of the
  # probability was a median 0.36 of brms's posterior SD, 0.008 within
  # 1e-4 of an edge. Est.Error is now the SD of plogis(Z), Z normal
  # around the logit with the interval's standard error; checked here
  # against stats::integrate() of the same law, on the rows at the
  # middle, at the smallest probability and near 0.01 and 0.999.
  set.seed(4)
  dm <- data.frame(y = c(stats::rnorm(150, -1.5), stats::rnorm(150, 1.5)))
  fit <- frm(y ~ 1, family = mixture(gaussian(), gaussian()), data = dm)
  pm <- pp_mixture(fit)
  p <- pm[, "Estimate", 1L]
  lse <- frmtmb:::fit_fd_se(fit, mixture_probs)[, 1L] / (p * (1 - p))
  # the interval is built from that standard error on the logit
  expect_equal(unname(pm[, "Q97.5", 1L]),
               unname(stats::plogis(stats::qlogis(p) +
                                      stats::qnorm(0.975) * lse)))
  ref_sd <- function(i, log = FALSE) {
    # the complement above one half, where the values are small
    s <- if (p[i] > 0.5) -1 else 1
    g <- function(z) {
      if (log) stats::plogis(z, log.p = TRUE) else stats::plogis(s * z)
    }
    dn <- function(z) stats::dnorm(z, stats::qlogis(p[i]), lse[i])
    lo <- stats::qlogis(p[i]) - 12 * lse[i]
    hi <- stats::qlogis(p[i]) + 12 * lse[i]
    m1 <- stats::integrate(function(z) g(z) * dn(z), lo, hi,
                           rel.tol = 1e-12)$value
    stats::integrate(function(z) (g(z) - m1)^2 * dn(z), lo, hi,
                     rel.tol = 1e-12)$value^0.5
  }
  rows <- unique(c(which.min(abs(p - 0.5)), which.min(p),
                   which.min(abs(p - 0.01)), which.min(abs(p - 0.999))))
  for (i in rows) {
    expect_equal(pm[i, "Est.Error", 1L], ref_sd(i), tolerance = 1e-5,
                 info = paste("p =", p[i]))
  }
  lg <- pp_mixture(fit, log = TRUE)
  i <- which.min(abs(p - 0.01))
  expect_equal(lg[i, "Est.Error", 1L], ref_sd(i, log = TRUE),
               tolerance = 1e-5)
  # near an edge it is far above the first-order error, which is what
  # the law adds: at p near 0.01, more than twice
  se_delta <- frmtmb:::fit_fd_se(fit, mixture_probs)[, 1L]
  j <- which.min(abs(p - 0.01))
  expect_gt(pm[j, "Est.Error", 1L] / se_delta[j], 2)
})

test_that("pp_mixture() covers an ordinal mixture and refuses a plain fit", {
  set.seed(2)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x, -1.5 * d$x) + stats::rlogis(n)
  d$y <- cut(lat, c(-Inf, -1, 0, 1, Inf), labels = FALSE)
  fo <- allow_warnings(
    frm(y ~ x, family = mixture(cumulative(), cumulative()), data = d),
    "Standard errors are not available")
  po <- allow_warnings(pp_mixture(fo), "Standard errors are not available")
  expect_equal(dim(po), c(n, 4L, 2L))
  expect_equal(unname(po[, "Estimate", ]), unname(mixture_probs(fo)))
  expect_equal(unname(rowSums(po[, "Estimate", ])), rep(1, n))
  # brms's words on a fit that is not a mixture
  expect_error(pp_mixture(fit_gauss()),
               "can only be applied to mixture models", fixed = TRUE)
})

test_that("add_criterion() on a fit refuses each brms criterion by name", {
  f <- fit_gauss()
  for (k in c("loo", "waic", "kfold", "loo_subsample", "bayes_R2",
              "loo_R2", "marglik")) {
    expect_error(add_criterion(f, k),
                 paste0("criterion = \"", k, "\") has no meaning"),
                 fixed = TRUE)
  }
  expect_error(add_criterion(f, "loo"), "AIC() or BIC()", fixed = TRUE)
  # brms's own validation, in its words
  expect_error(add_criterion(f, "aic"), "should be a subset of")
  expect_error(add_criterion(f), "should be a subset of")
  allow_warnings(expect_error(add_criterion(f, "R2"), "bayes_R2"),
                 "Criterion 'R2' is deprecated",
                 require = "Criterion 'R2' is deprecated")
})

test_that("summary(waic =) is refused by name with the reason", {
  f <- fit_gauss()
  expect_error(summary(f, waic = TRUE), "cannot honor `waic`")
  expect_error(summary(f, waic = TRUE), "AIC and BIC")
  # the control: summary() itself still runs
  expect_s3_class(summary(f), "summary.frmtmb_fit")
})

test_that("plot() of a fit refuses brms's arguments by name", {
  f <- fit_gauss()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (a in c("N", "nvariables", "variable", "regex", "pars", "fixed",
              "combo", "bins", "theme", "plot", "newpage")) {
    args <- list(f, TRUE)
    names(args) <- c("x", a)
    expect_error(do.call(plot, args), paste0("cannot honor `", a, "`"),
                 fixed = TRUE)
  }
  expect_error(plot(f, N = 2, ask = FALSE), "frm_sample(fit)", fixed = TRUE)
  # an unknown name still gets the generic refusal, and the control
  # draws
  expect_error(plot(f, nosuch = 1), "has no argument `nosuch`")
  expect_identical(plot(f, ask = FALSE), f)
})

test_that("the categorical refusal of method = 'predict' names its family", {
  set.seed(1)
  d <- data.frame(x = stats::rnorm(200),
                  yc = factor(sample(c("a", "b", "c"), 200, TRUE)))
  fc <- frm(yc ~ x, family = categorical(), data = d)
  e <- tryCatch(conditional_effects(fc, "x", method = "predict"),
                error = conditionMessage)
  expect_match(e, "family 'categorical'", fixed = TRUE)
  expect_false(grepl("ordinal", e, fixed = TRUE))
  # the ordinal family keeps its wording
  d$yo <- sample(1:3, 200, TRUE)
  fo <- frm(yo ~ x, family = cumulative(), data = d)
  expect_error(conditional_effects(fo, "x", method = "predict"),
               "on an ordinal family")
})

test_that("a very wide logit law integrates without stopping", {
  # The re-check of 2026-10-07, n1: integrate() said "the integral is
  # probably divergent" at these three points, which stopped
  # pp_mixture() for the whole fit (dev/surface-rev2-lnsd3.R).
  lnsd <- frmtmb:::logitnormal_sd
  pts <- list(c(1e-3, 3000), c(1e-5, 3000), c(1e-30, 1e4))
  for (pt in pts) {
    v <- lnsd(pt[1L], 1 - pt[1L], pt[2L])
    expect_true(is.finite(v), info = paste(pt, collapse = " "))
    # a logit SD of thousands puts half the mass at each end: the SD of
    # the probability approaches its bound of one half
    expect_lte(v, 0.5)
    expect_gt(v / 0.5, 0.95)
  }
  # the log case at the same points is finite too
  expect_true(all(is.finite(vapply(pts, function(pt) {
    lnsd(pt[1L], 1 - pt[1L], pt[2L], log = TRUE)
  }, 0))))
})

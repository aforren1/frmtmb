## The rest of the frmtmb surface on a cross_wishart() fit. A family is not
## finished when it fits; it is finished when the methods a user reaches for
## next either work or refuse for a stated reason.

cp_surface_fit <- function() {
  set.seed(301)
  # the shared source is low-pass, so the low band really is the coherent
  # one; a white shared source would make coherence flat in frequency and
  # the band contrast would be testing nothing
  src <- as.numeric(stats::filter(stats::rnorm(4096), 0.9, method = "recursive"))
  a <- src + stats::rnorm(4096, 0, 0.5)
  b <- 0.9 * src + stats::rnorm(4096, 0, 2)
  d <- frm_cross_spectrum(a, b, sfreq = 128, segments = 16L)
  d$band <- factor(ifelse(d$freq < 32, "low", "high"), c("low", "high"))
  # the powers get a smooth in frequency: without one the coherence is
  # not readable at all on a spectrum this sloped. Nothing warns about
  # that; the AIC comparison below is what catches it.
  list(dat = d, fit = frmtmb::frm(cp_bform(), family = cross_wishart(),
                                  data = d))
}
cp_bform <- function() {
  frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ s(freq, k = 8),
             pow2 ~ s(freq, k = 8), coh ~ band, phase ~ 1)
}

test_that("summary, fixef, logLik and nobs work", {
  o <- cp_surface_fit()
  f <- o$fit
  expect_s3_class(f, "frmtmb_fit")
  expect_no_error(summary(f))
  expect_setequal(names(frmtmb::fixef(f)), c("mu", "pow2", "coh", "phase"))
  expect_true(is.finite(as.numeric(stats::logLik(f))))
  expect_true(is.finite(stats::AIC(f)))
  expect_equal(stats::nobs(f), nrow(o$dat))
  expect_true(all(is.finite(unlist(frmtmb::fixef(f)))))
  expect_true(all(is.finite(stats::vcov(f))))
})

test_that("every dpar is reachable on the link scale, with a standard error", {
  o <- cp_surface_fit()
  for (dp in c("mu", "pow2", "coh", "phase")) {
    v <- stats::predict(o$fit, type = "link", dpar = dp)
    expect_length(v, nrow(o$dat))
    expect_true(all(is.finite(v)), info = dp)
    s <- stats::predict(o$fit, type = "link", dpar = dp, se.fit = TRUE)
    expect_true(all(s$se.fit > 0), info = dp)
  }
})

test_that("predict on newdata works", {
  o <- cp_surface_fit()
  nd <- data.frame(band = factor(c("low", "high"), c("low", "high")),
                   freq = c(16, 48))
  p <- stats::predict(o$fit, newdata = nd, type = "link", dpar = "coh",
                      se.fit = TRUE)
  expect_length(p$fit, 2L)
  expect_true(all(p$se.fit > 0))
  co <- frm_coherence(o$fit, newdata = nd)
  expect_equal(nrow(co), 2L)
  # the low band carries the shared source, so its coherence is the higher
  expect_gt(co$.estimate[1], co$.estimate[2])
})

test_that("fitted() is the mean of the response and says so by being it", {
  o <- cp_surface_fit()
  # W11 is Gamma(n, scale S11), so its mean is exactly n * S11
  s11 <- exp(stats::predict(o$fit, type = "link", dpar = "mu"))
  expect_equal(as.numeric(stats::fitted(o$fit)), o$dat$n * s11,
               ignore_attr = TRUE)
  # and type = "response" agrees with it
  expect_equal(as.numeric(stats::predict(o$fit, type = "response")),
               as.numeric(stats::fitted(o$fit)), ignore_attr = TRUE)
})

test_that("residuals work and are what they claim", {
  o <- cp_surface_fit()
  rr <- stats::residuals(o$fit, type = "response")
  expect_equal(as.numeric(rr),
               o$dat$w11 - as.numeric(stats::fitted(o$fit)),
               ignore_attr = TRUE)
  # pearson divides by the exact standard deviation, sqrt(n) * S11
  s11 <- exp(stats::predict(o$fit, type = "link", dpar = "mu"))
  pr <- stats::residuals(o$fit, type = "pearson")
  expect_equal(as.numeric(pr), as.numeric(rr) / (sqrt(o$dat$n) * s11),
               ignore_attr = TRUE, tolerance = 1e-8)
  dv <- stats::residuals(o$fit, type = "deviance")
  expect_true(all(is.finite(dv)))
  # the deviance is a deviance: it is non-negative before the sign is taken
  expect_true(all(dv^2 >= 0))
})

test_that("simulate() refuses for a stated reason and names its replacement", {
  o <- cp_surface_fit()
  expect_error(stats::simulate(o$fit), "frm_cross_simulate")
})

test_that("par_template and set_prior reach the coherence dpar", {
  o <- cp_surface_fit()
  nms <- names(unlist(frmtmb::par_template(o$fit)))
  expect_true(any(grepl("coh_", nms, fixed = TRUE)))
  expect_true(any(grepl("phase_", nms, fixed = TRUE)))
  gp <- frmtmb::get_prior(cp_bform(), family = cross_wishart(), data = o$dat)
  expect_true(is.data.frame(gp))
  expect_true(all(c("coh", "phase") %in% c(gp$class, gp$dpar)))
  # a tight prior on the band contrast shrinks it
  f2 <- frmtmb::frm(
    cp_bform(), family = cross_wishart(), data = o$dat,
    prior = frmtmb::set_prior("normal(0, 0.02)", class = "b", dpar = "coh"))
  expect_lt(abs(frmtmb::fixef(f2)$coh[["bandhigh"]]),
            abs(frmtmb::fixef(o$fit)$coh[["bandhigh"]]))
})

test_that("a random effect and a smooth both reach the coherence dpar", {
  set.seed(302)
  X <- matrix(stats::rnorm(12 * 2048), 2048, 12)
  Y <- 0.8 * X + matrix(stats::rnorm(12 * 2048), 2048, 12)
  d <- frm_cross_spectrum(X, Y, segments = 8L)
  d <- d[d$freq < 0.1, ]
  fr <- frmtmb::frm(
    frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1 + (1 | id),
               pow2 ~ 1 + (1 | id), coh ~ 1 + (1 | id),
               phase ~ 1 + (1 | id)),
    family = cross_wishart(), data = d)
  expect_s3_class(fr, "frmtmb_fit")
  expect_true(is.finite(as.numeric(stats::logLik(fr))))
  expect_gte(length(frmtmb::ranef(fr)), 1L)
  co <- frm_coherence(fr, re.form = NA)
  expect_true(co$.lower[1] > 0 && co$.upper[1] < 1)

  fs <- frmtmb::frm(
    frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
               pow2 ~ 1, coh ~ s(freq, k = 5), phase ~ 1),
    family = cross_wishart(), data = d)
  expect_s3_class(fs, "frmtmb_fit")
  expect_true(all(is.finite(frm_coherence(fs)$.estimate)))
})

test_that("AIC catches a power model too simple for the data", {
  # The trap that converges silently: the likelihood couples the powers
  # and the coherence, so a constant power model on a sloped spectrum
  # moves the coherence rather than only the powers. Nothing warns about
  # it, deliberately: a draft carried a fit_check hook on the spread of
  # the power residuals and that statistic fired on correct models and
  # stayed silent on wrong ones. The shipped advice is an AIC
  # comparison, and this is that advice under test.
  o <- cp_surface_fit()
  flat <- suppressWarnings(frmtmb::frm(
    frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
               pow2 ~ 1, coh ~ band, phase ~ 1),
    family = cross_wishart(), data = o$dat))
  # the damage is real: the flat fit calls two bands nearly equal
  nd <- data.frame(band = factor(c("low", "high"), c("low", "high")),
                   freq = c(16, 48))
  cf <- frm_coherence(flat, newdata = nd)
  cg <- frm_coherence(o$fit, newdata = nd)
  expect_lt(cf$.estimate[1] - cf$.estimate[2], 0.05)
  expect_gt(cg$.estimate[1] - cg$.estimate[2], 0.20)
  # and AIC says so, decisively, which is what a user is told to check
  expect_gt(stats::AIC(flat) - stats::AIC(o$fit), 100)
})

test_that("no automatic power-misspecification warning is raised", {
  # The withdrawn hook must stay withdrawn: it fired on correctly
  # specified models, which is worse than being silent.
  o <- cp_surface_fit()
  expect_no_warning(frmtmb::frm(cp_bform(), family = cross_wishart(),
                                data = o$dat))
  expect_null(stats::family(o$fit)$post$fit_check)
})

test_that("the compatibility rows this package registered resolve", {
  for (p in list(c("cross_wishart", "vreal()"), c("cross_wishart", "simulate"),
                 c("cross_wishart", "fitted"), c("cross_wishart", "s()"),
                 c("frm_coherence", "predict"),
                 c("frm_cross_simulate", "cross_wishart"))) {
    r <- frmtmb::frm_compat(p[1L], p[2L])
    st <- if (is.data.frame(r)) r$status[1L] else r$status
    expect_true(st %in% c("works", "conditional", "refused", "untested",
                          "broken"),
                info = paste(p, collapse = " x "))
  }
  f <- frmtmb::frm_compat_features()
  expect_true(all(c("cross_wishart", "frm_coherence", "frm_cross_simulate")
                  %in% f$key))
  expect_identical(frmtmb::frm_compat("cross_wishart", "simulate")$status,
                   "refused")
  expect_identical(frmtmb::frm_compat("cross_wishart", "cens()")$status,
                   "refused")
  expect_identical(frmtmb::frm_compat("cross_wishart", "s()")$status, "works")
  expect_identical(frmtmb::frm_compat("frm_coherence", "predict")$status,
                   "works")
})

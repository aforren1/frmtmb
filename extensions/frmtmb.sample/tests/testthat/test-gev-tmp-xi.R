# A scalar gen_extreme_value() xi is sampled through brms's tmp_xi
# (frmtmb's gev_xi_chart_frame()). Sampled on its log1p scale instead,
# the support's edge was a wall: on this data set (seed 1402, n = 40,
# xi = -0.85) 3856 of 4000 transitions diverged under frm_sample(), Rhat
# 2.30, against 0 in brms (dev/reviews/2026-10-13-famuni.md, fix 5).

skip_on_cran()
skip_if_not_installed("tmbstan")

divergent <- function(dr) {
  sp <- rstan::get_sampler_params(dr$stanfit, inc_warmup = FALSE)
  sum(vapply(sp, function(s) sum(s[, "divergent__"]), 0))
}

test_that("a scalar xi samples without divergences, every draw in support", {
  rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
  set.seed(1402)
  d <- data.frame(y = rg(40, 0, 1, -0.85))
  fit <- allow_warnings(
    frm(y ~ 1, family = gen_extreme_value(), data = d),
    allowed = c("did not report convergence", "the edge at -1",
                "maximum absolute gradient", "Standard errors"))
  dr <- suppressMessages(frm_sample(fit, chains = 2, iter = 1000, seed = 1,
                                    refresh = 0))
  nd <- nrow(as.data.frame(dr))
  # the chart removes the wall: 0 of 4000 on 4 x 2000 at this seed
  # (dev/famuni-out/tmpxi-sample.txt); a stray divergence is not the
  # defect, a fraction of a percent is the bound
  expect_lte(divergent(dr), 0.005 * nd)
  dm <- as.data.frame(dr)
  xi <- dm[["xi"]]
  expect_true(all(is.finite(xi)))
  expect_true(all(xi > -1))
  # each draw's xi is the one its density used: every row inside the
  # support at that draw's mu and sigma, which is what scale_xi() is for
  tmin <- vapply(seq_along(xi), function(i) {
    min(1 + xi[i] * (d$y - dm[["b_Intercept"]][i]) / dm[["sigma"]][i])
  }, 0)
  expect_true(all(tmin > 0))
  # and the draws methods read it: log_lik is finite at every draw
  expect_true(all(is.finite(log_lik(dr))))
})

test_that("a prior on class xi is read on tmp_xi, and a cens() fit warns", {
  rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
  set.seed(1401)
  d <- data.frame(y = rg(40, 0, 1, -0.45))
  fit <- frm(y ~ 1, family = gen_extreme_value(), data = d)
  ent <- frmtmb::resolve_prior_input(
    fit, set_prior("normal(0, 2.5)", class = "xi"))$entries
  terms <- frmtmb::gev_xi_chart_terms(fit)
  expect_length(terms, 1L)
  expect_length(terms[[1L]]$blocked, 0L)
  moved <- frmtmb.sample:::mo_sample_entries(ent, fit)
  hit <- Filter(function(e) identical(e$scale, "tmp_xi"), moved)
  expect_length(hit, 1L)
  # the case where the chart must not be used: cens() needs a CDF the
  # chart does not carry, and the sampler says what that costs
  d$cc <- as.numeric(d$y > stats::quantile(d$y, 0.9))
  fc <- frm(y | cens(cc) ~ 1, family = gen_extreme_value(), data = d)
  expect_identical(frmtmb::gev_xi_chart_terms(fc)[[1L]]$blocked, "cens()")
  allow_warnings(
    suppressMessages(frm_sample(fc, chains = 1, iter = 300, seed = 1,
                                refresh = 0)),
    allowed = c("not through brms's tmp_xi", "divergent", "Bulk Effective",
                "Tail Effective", "R-hat", "treedepth", "Examine the pairs",
                "chains have not mixed"),
    require = "not through brms's tmp_xi")
})

test_that("a mixture of GEV components is named, not charted", {
  # the chart maps a plain GEV response's xi only; a two-GEV mixture at
  # xi -0.3 had 84 of 4000 transitions divergent and no word about why
  # (dev/reviews/2026-10-13-famuni.md, round 2, item 3)
  rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
  set.seed(35)
  d <- data.frame(y = c(rg(100, 0, 1, -0.3), rg(100, 6, 1, -0.3)))
  fm <- frm(y ~ 1, family = mixture(gen_extreme_value(),
                                    gen_extreme_value()), data = d)
  terms <- frmtmb::gev_xi_chart_terms(fm)
  expect_identical(vapply(terms, `[[`, "", "dpar"), c("xi1", "xi2"))
  expect_true(all(vapply(terms, function(tm) {
    identical(tm$blocked, "a mixture()")
  }, NA)))
  allow_warnings(
    suppressMessages(frm_sample(fm, chains = 1, iter = 300, seed = 1,
                                refresh = 0, init = "random")),
    allowed = c("not through brms's tmp_xi", "divergent", "Bulk Effective",
                "Tail Effective", "R-hat", "treedepth", "Examine the pairs",
                "chains have not mixed", "multimodal"),
    require = "xi1 of gen_extreme_value() is sampled on its log1p scale")
})

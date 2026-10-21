# loo_epred(), loo_linpred(), loo_predict(), loo_predictive_interval(),
# loo_R2() and loo_model_weights() on draws. Each is a function of the
# log_lik() matrix and one prediction matrix, so each is tested against
# loo::psis() and loo::E_loo() called by hand on those matrices.
#
# Seen to fail on frmtmb.sample 0.20.0 (rellib-r10): every function was
# missing ("could not find function"), add_criterion(x, "loo_R2")
# refused by name, and with brms attached loo_predict() on draws reached
# brms's generic with no method for frmtmb_draws. Log:
# dev/loodraws-log/tests-base-loo-predict.txt.

skip_on_cran()
skip_if_not_installed("loo")
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

lp_data <- function(seed = 21, n = 40L, ng = 8L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x +
                         stats::rnorm(ng, 0, 0.6)[dd$g])
  dd$y2 <- stats::rnorm(n, -1 + 0.3 * dd$x)
  dd$k <- stats::rpois(n, exp(0.2 + 0.4 * dd$x))
  dd$o <- factor(sample(1:3, n, TRUE), ordered = TRUE)
  dd
}

# draws jittered around the ML estimate on the unconstrained scale: no
# sampler is needed, because every quantity here is a function of
# whatever draws there are
lp_draws <- function(fit, nd = 120L, sd = 0.05, seed = 5L) {
  tpl <- fit$frame$par_template
  # a betad fixed by the family (an ordinal model's) is mapped, so the
  # sampler has no column for it, as brms_par_labels() has no name
  est <- unlist(lapply(names(tpl), function(cp) {
    v <- fit$estimates[[cp]]
    fx <- fit$frame[["betad_fixed_idx"]]
    if (identical(cp, "betad") && length(fx)) v <- v[-fx]
    v
  }))
  set.seed(seed)
  M <- matrix(rep(est, each = nd) + stats::rnorm(nd * length(est), 0, sd),
              nd, dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}

lp_case <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      dd <- lp_data()
      q <- function(e) suppressWarnings(suppressMessages(e))
      cache <<- list(
        dd = dd,
        ri = lp_draws(q(frm(bf(y ~ x + (1 | g)), family = gaussian(),
                            data = dd))),
        pois = lp_draws(q(frm(bf(k ~ x), family = poisson(), data = dd))),
        mv = lp_draws(q(frm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE),
                            data = dd))),
        ord = lp_draws(q(frm(bf(o ~ x), family = cumulative(), data = dd))))
    }
    cache
  }
})

# the weights by hand: brms's r_eff_log_lik() on the chains of the
# draws, which is one chain for draws with no stanfit
by_hand_psis <- function(ll) {
  suppressWarnings(loo::psis(-ll, r_eff = loo::relative_eff(
    exp(ll), chain_id = attr(ll, "chain_id"))))
}
by_hand_e <- function(x, ps, type, probs = 0.5, ll = NULL) {
  v <- loo::E_loo(x, ps, type = type, probs = probs, log_ratios = ll)$value
  if (is.matrix(v)) t(v) else matrix(v)
}
lp_ok <- c("Pareto k", "pareto_k", "Relative effective",
           "ESS has been capped", "Some Pareto k")
quiet <- function(e) allow_warnings(suppressMessages(e), lp_ok)

test_that("the weights are psis(), loo()'s and the hand-built ones", {
  cs <- lp_case()
  ll <- log_lik(cs$ri)
  ps <- by_hand_psis(ll)
  expect_identical(stats::weights(quiet(psis(cs$ri))), stats::weights(ps))
  expect_identical(
    stats::weights(quiet(loo(cs$ri, save_psis = TRUE))$psis_object),
    stats::weights(ps))
})

test_that("loo_epred() and loo_linpred() are E_loo() of the prediction", {
  cs <- lp_case()
  for (d in list(cs$ri, cs$pois)) {
    ll <- log_lik(d)
    ps <- by_hand_psis(ll)
    ep <- posterior_epred(d)
    lp <- posterior_linpred(d)
    for (ty in c("mean", "var", "quantile")) {
      pr <- if (ty == "quantile") c(0.1, 0.5, 0.9) else 0.5
      lty <- if (ty == "var") "variance" else ty
      lab <- if (ty == "quantile") c("q10", "q50", "q90") else ty
      got <- quiet(loo_epred(d, type = ty, probs = pr))
      ref <- by_hand_e(ep, ps, lty, pr)
      colnames(ref) <- lab
      expect_identical(got, ref)
      gotl <- quiet(loo_linpred(d, type = ty, probs = pr))
      refl <- by_hand_e(lp, ps, lty, pr)
      colnames(refl) <- lab
      expect_identical(gotl, refl)
    }
    # one quantile is one column, as brms's matrix(y)
    expect_identical(dim(quiet(loo_epred(d, type = "quantile"))),
                     c(ncol(ll), 1L))
  }
})

test_that("loo_predict() and the interval weight posterior_predict()", {
  cs <- lp_case()
  d <- cs$pois
  ps <- by_hand_psis(log_lik(d))
  set.seed(8)
  pp <- posterior_predict(d)
  set.seed(8)
  got <- quiet(loo_predict(d))
  ref <- by_hand_e(pp, ps, "mean")
  colnames(ref) <- "mean"
  expect_identical(got, ref)
  # brms's (1 - prob) / 2, which is not exactly 0.1
  a <- (1 - 0.8) / 2
  set.seed(8)
  iv <- quiet(loo_predictive_interval(d, prob = 0.8))
  set.seed(8)
  expect_identical(iv, quiet(loo_predict(d, type = "quantile",
                                         probs = c(a, 1 - a))))
  expect_identical(colnames(iv), c("q10", "q90"))
  expect_identical(dim(iv), c(nrow(cs$dd), 2L))
  expect_error(loo_predictive_interval(d, prob = c(0.5, 0.9)),
               "Argument 'prob' should be of length 1.", fixed = TRUE)
})

test_that("a given psis_object is used as it is, and checked", {
  cs <- lp_case()
  d <- cs$ri
  ps <- by_hand_psis(log_lik(d))
  # without a log_lik() of its own: a mock that fails proves none runs
  testthat::local_mocked_bindings(loo_draws_weights = function(...) {
    stop("computed the weights")
  })
  got <- loo_epred(d, psis_object = ps)
  ref <- by_hand_e(posterior_epred(d), ps, "mean")
  colnames(ref) <- "mean"
  expect_identical(got, ref)
  expect_error(loo_epred(d, psis_object = ps, ndraws = 50),
               "must be a psis object for the 50 draws")
  expect_error(loo_epred(d, psis_object = log_lik(d)),
               "must be a psis object")
})

test_that("ndraws and draw_ids pick the draws once, for both halves", {
  cs <- lp_case()
  d <- cs$ri
  # irregular rows, and enough of them for the Pareto fit of each column
  ids <- c(seq(2L, 60L, by = 2L), seq(61L, 119L, by = 3L))
  ll <- log_lik(d, draw_ids = ids)
  # a subset carries no chain structure: brms's one-chain r_eff
  expect_null(attr(ll, "chain_id"))
  ps <- suppressWarnings(loo::psis(-ll, r_eff = loo::relative_eff(
    exp(ll), chain_id = rep(1L, length(ids)))))
  ref <- by_hand_e(posterior_epred(d, draw_ids = ids), ps, "mean")
  colnames(ref) <- "mean"
  expect_identical(quiet(loo_epred(d, draw_ids = ids)), ref)
  rows <- frmtmb.sample:::draws_subsample(d, 60, NULL)
  expect_identical(quiet(loo_epred(d, ndraws = 60)),
                   quiet(loo_epred(d, draw_ids = rows)))
})

test_that("a subset of draws is weighted on one chain, everywhere", {
  # brms 2.23.0's get_chain_id(): every draw keeps the chains, a subset
  # is one chain, in loo(), psis(), the loo_*() functions and
  # pp_check(). frmtmb.sample 0.20.0 dropped r_eff (to 1) on a subset in
  # loo() and psis() while pp_check() took one chain, so the same draws
  # were weighted two ways; and psis() ignored draw_ids in silence.
  skip_if_not_installed("bayesplot")
  cs <- lp_case()
  d <- cs$ri
  # two chains, so the two rules give different efficiencies
  testthat::local_mocked_bindings(draws_nchains = function(x) 2L)
  ids <- c(seq(2L, 60L, by = 2L), seq(61L, 119L, by = 3L))
  ll <- log_lik(d, draw_ids = ids)
  one <- quiet(loo::relative_eff(exp(ll), chain_id = rep(1L, length(ids))))
  ps <- quiet(psis(d, draw_ids = ids))
  expect_identical(dim(ps), dim(ll))
  expect_identical(attr(ps, "r_eff"), one)
  lo <- quiet(loo(d, draw_ids = ids, save_psis = TRUE))
  expect_identical(stats::weights(lo$psis_object), stats::weights(ps))
  testthat::local_mocked_bindings(draws_bayesplot_fun = function(nm, what) {
    function(y, yrep, lw, ...) lw
  })
  lw <- quiet(pp_check(d, type = "loo_pit_overlay", draw_ids = ids))
  expect_identical(lw, stats::weights(ps, log = TRUE))
  # the absent case: every draw keeps both chains
  llall <- log_lik(d)
  expect_identical(attr(quiet(psis(d)), "r_eff"), quiet(loo::relative_eff(
    exp(llall), chain_id = rep(1:2, each = nrow(d$draws) / 2))))
  # and psis() refuses what it does not read
  expect_error(psis(d, junk = 1), "psis() has no argument `junk`",
               fixed = TRUE)
})

test_that("a multivariate model gives brms's layers, or one response", {
  cs <- lp_case()
  d <- cs$mv
  ps_joint <- by_hand_psis(log_lik(d))
  got <- quiet(loo_epred(d))
  expect_identical(dim(got), c(nrow(cs$dd), 1L, 2L))
  expect_identical(dimnames(got), list(NULL, "mean", c("y", "y2")))
  for (r in c("y", "y2")) {
    expect_identical(got[, , r],
                     by_hand_e(posterior_epred(d, resp = r), ps_joint,
                               "mean")[, 1])
  }
  # with resp, that response's own log-likelihood, as brms's loo(resp =)
  ps2 <- by_hand_psis(log_lik(d, resp = "y2"))
  ref <- by_hand_e(posterior_epred(d, resp = "y2"), ps2, "mean")
  colnames(ref) <- "mean"
  expect_identical(quiet(loo_epred(d, resp = "y2")), ref)
  expect_error(loo_epred(d, resp = "z"), "must name one response")
})

test_that("an ordinal loo_epred() keeps the category margin third", {
  cs <- lp_case()
  d <- cs$ord
  ps <- by_hand_psis(log_lik(d))
  ep <- posterior_epred(d)
  got <- quiet(loo_epred(d, type = "quantile", probs = c(0.25, 0.75)))
  expect_identical(dim(got), c(nrow(cs$dd), 2L, dim(ep)[3L]))
  expect_identical(dimnames(got)[[2L]], c("q25", "q75"))
  expect_identical(dimnames(got)[[3L]], dimnames(ep)[[3L]])
  ref <- by_hand_e(ep[, , 2L], ps, "quantile", c(0.25, 0.75))
  expect_identical(unname(got[, , 2L]), unname(ref))
})

test_that("loo_R2() is brms's bootstrap of the LOO residuals", {
  cs <- lp_case()
  d <- cs$ri
  y <- cs$dd$y
  ll <- log_lik(d)
  ps <- by_hand_psis(ll)
  ep <- posterior_epred(d)
  # the reference spells out brms's .loo_R2() with the weighted sums as
  # matrix products, which is the arithmetic the method uses
  ref <- local({
    set.seed(4)
    yl <- loo::E_loo(ep, ps, log_ratios = -ll)$value
    e <- yl - y
    S <- nrow(ep)
    N <- ncol(ep)
    w <- matrix(stats::rexp(S * N, rate = 1), nrow = S, ncol = N)
    w <- w / rowSums(w)
    s <- w %*% cbind(y, y^2, e, e^2)
    # brms's N / (N - 1) cancels in the ratio, but not in the rounding
    r <- 1 - ((N / (N - 1)) * (s[, 4] - s[, 3]^2)) /
      ((N / (N - 1)) * (s[, 2] - s[, 1]^2))
    pmin(pmax(r, -1), 1)
  })
  r2 <- quiet(loo_R2(d, summary = FALSE, seed = 4))
  expect_identical(dim(r2), c(nrow(d$draws), 1L))
  expect_identical(colnames(r2), "R2")
  expect_identical(r2[, 1], ref)
  # brms's seed handling: the caller's stream is put back
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  invisible(quiet(loo_R2(d, seed = 4)))
  expect_identical(stats::runif(1), before)
  sm <- quiet(loo_R2(d, seed = 4))
  expect_identical(dim(sm), c(1L, 4L))
  expect_identical(sm, posterior_summary(r2))
  expect_identical(colnames(quiet(loo_R2(cs$mv, summary = FALSE,
                                         seed = 1))),
                   c("R2y", "R2y2"))
  expect_error(loo_R2(cs$ord), "category distribution")
  expect_error(loo_R2(d, args_epred = list(ndraws = 5)),
               "give `ndraws` in `...`")
})

test_that("add_criterion() stores loo_R2, and loo_R2() reads it", {
  cs <- lp_case()
  d <- quiet(add_criterion(cs$ri, "loo_R2", seed = 4))
  expect_identical(d$criteria$loo_R2,
                   quiet(loo_R2(cs$ri, summary = FALSE, seed = 4)))
  # nothing asked beyond the summary: the stored draws, so no weights
  testthat::local_mocked_bindings(draws_psis_of = function(...) {
    stop("computed the weights")
  })
  expect_identical(loo_R2(d), posterior_summary(d$criteria$loo_R2))
  expect_identical(loo_R2(d, summary = FALSE), d$criteria$loo_R2)
  # a call that changes the draws behind it recomputes
  expect_error(loo_R2(d, ndraws = 10), "computed the weights")
})

test_that("loo_model_weights() is loo's, on each model's log_lik()", {
  cs <- lp_case()
  d1 <- cs$ri
  d2 <- lp_draws(suppressWarnings(frm(bf(y ~ 1), family = gaussian(),
                                      data = cs$dd)))
  lls <- list(log_lik(d1), log_lik(d2))
  reff <- lapply(lls, function(l) {
    loo::relative_eff(exp(l), chain_id = attr(l, "chain_id"))
  })
  w <- quiet(loo_model_weights(d1, d2))
  ref <- quiet(loo::loo_model_weights(lls, r_eff_list = reff))
  expect_identical(as.numeric(w), as.numeric(ref))
  expect_s3_class(w, "stacking_weights")
  expect_identical(names(w), c("d1", "d2"))
  expect_identical(names(quiet(loo_model_weights(d1, d1))), c("d1", "d1"))
  expect_identical(names(quiet(loo_model_weights(d1, b = d2))),
                   c("d1", "b"))
  set.seed(6)
  wp <- quiet(loo_model_weights(d1, d2, method = "pseudobma", BB_n = 200))
  set.seed(6)
  rp <- quiet(loo::loo_model_weights(lls, method = "pseudobma", BB_n = 200,
                                     r_eff_list = reff))
  expect_identical(as.numeric(wp), as.numeric(rp))
  # posterior_average()'s stacking is the same function
  pa <- quiet(posterior_average(d1, d2, variable = "sigma", seed = 1))
  expect_identical(unname(attr(pa, "weights")),
                   as.numeric(w) / sum(as.numeric(w)))
  expect_error(loo_model_weights(d1, d2, newdata = cs$dd),
               "cannot honor `newdata`")
  expect_error(loo_model_weights(d1, d2, optim = "x"),
               "has no argument `optim`")
  expect_error(loo_model_weights(d1, cs$dd), "takes models without a name")
})

test_that("pp_check()'s loo weights are psis()'s on every draw", {
  skip_if_not_installed("bayesplot")
  cs <- lp_case()
  d <- cs$ri
  testthat::local_mocked_bindings(draws_bayesplot_fun = function(nm, what) {
    function(y, yrep, lw, ...) lw
  })
  lw <- quiet(pp_check(d, type = "loo_pit_overlay"))
  expect_identical(lw, stats::weights(quiet(psis(d)), log = TRUE))
})

test_that("newdata and the refits are refused, as loo() refuses them", {
  cs <- lp_case()
  d <- cs$ri
  for (f in list(loo_epred, loo_linpred, loo_predict)) {
    expect_error(f(d, newdata = cs$dd[1:3, ]), "cannot honor `newdata`")
    expect_error(f(d, re_formula = NA), "cannot honor `re_formula`")
    expect_error(f(d, allow_new_levels = TRUE),
                 "cannot honor `allow_new_levels`")
    expect_error(f(d, moment_match = TRUE),
                 "cannot honor moment_match = TRUE")
    expect_error(f(d, save_psis = TRUE), "has no argument `save_psis`")
  }
  expect_error(loo_predictive_interval(d, newdata = cs$dd[1:3, ]),
               "loo_predictive_interval() cannot honor `newdata`",
               fixed = TRUE)
  expect_error(loo_R2(d, newdata = cs$dd), "cannot honor `newdata`")
  # brms 2.23.0 reads `scale` nowhere, so it is named rather than ignored
  expect_error(loo_linpred(d, scale = "response"),
               "has no argument `scale`")
  # the guard's absent case: brms's defaults pass
  expect_no_error(quiet(loo_epred(d, newdata = NULL, re_formula = NULL,
                                  moment_match = FALSE)))
})

test_that("on a maximum-likelihood fit each refuses and names frm_sample()", {
  fit <- suppressWarnings(frm(bf(y ~ x), family = gaussian(),
                              data = lp_data()))
  for (f in c("loo_epred", "loo_linpred", "loo_predict",
              "loo_predictive_interval", "loo_R2")) {
    expect_error(get(f)(fit), paste0(f, "(frmtmb.sample::frm_sample(fit))"),
                 fixed = TRUE)
  }
  expect_error(loo_model_weights(fit, fit), "frm_sample(fit1)",
               fixed = TRUE)
})

test_that("brms's own default values pass, in every function", {
  # brms's defaults are not all NULL: allow_new_levels = FALSE and
  # sample_new_levels = "uncertainty". A call that spells them out asks
  # for nothing log_lik() refuses, and frmtmb.sample 0.20.0's loo()
  # refused it, as the first spelling of these functions did.
  cs <- lp_case()
  d <- cs$ri
  dflt <- list(newdata = NULL, re_formula = NULL, allow_new_levels = FALSE,
               sample_new_levels = "uncertainty")
  for (f in list(loo_epred, loo_linpred, loo_predict,
                 loo_predictive_interval)) {
    set.seed(2)
    ref <- quiet(f(d))
    set.seed(2)
    expect_identical(quiet(do.call(f, c(list(d), dflt))), ref)
  }
  expect_identical(quiet(do.call(loo_R2, c(list(d, seed = 3), dflt))),
                   quiet(loo_R2(d, seed = 3)))
  d2 <- lp_draws(suppressWarnings(frm(bf(y ~ 1), family = gaussian(),
                                      data = cs$dd)))
  expect_identical(as.numeric(quiet(do.call(loo_model_weights,
                                            c(list(d, d2), dflt)))),
                   as.numeric(quiet(loo_model_weights(d, d2))))
  expect_identical(quiet(do.call(loo, c(list(d), dflt)))$estimates,
                   quiet(loo(d))$estimates)
  # draw_ids reaches log_lik() in loo_model_weights() as in the others
  ids <- seq(1L, 120L, by = 2L)
  lls <- list(log_lik(d, draw_ids = ids), log_lik(d2, draw_ids = ids))
  ref <- quiet(loo::loo_model_weights(lls, r_eff_list = lapply(lls,
    function(l) loo::relative_eff(exp(l), chain_id = rep(1L, nrow(l))))))
  expect_identical(as.numeric(quiet(loo_model_weights(d, d2,
                                                      draw_ids = ids))),
                   as.numeric(ref))
  # the inverse: a value that asks for something is still refused
  expect_error(loo_model_weights(d, d2, allow_new_levels = TRUE),
               "cannot honor `allow_new_levels`")
  expect_error(loo_epred(d, sample_new_levels = "gaussian"),
               "cannot honor `sample_new_levels`")
})

test_that("resp takes brms's spelling everywhere, and names as brms", {
  cs <- lp_case()
  dd <- cs$dd
  dd$y_a <- dd$y
  d <- lp_draws(suppressWarnings(frm(mvbf(bf(y_a ~ x), bf(y2 ~ x),
                                          rescor = FALSE), data = dd)))
  expect_identical(quiet(loo_epred(d, resp = "ya")),
                   quiet(loo_epred(d, resp = "y_a")))
  expect_identical(dimnames(quiet(loo_epred(d)))[[3L]], c("ya", "y2"))
  # loo_R2() follows brms's validate_resp() exactly, as bayes_R2() does
  # (test-brms-pins.R): brms's spelling only
  r2 <- quiet(loo_R2(d, resp = "ya", seed = 1))
  expect_identical(rownames(r2), "R2ya")
  expect_error(loo_R2(d, resp = "y_a"), "Valid response variables")
  expect_identical(
    quiet(loo_model_weights(d, d, resp = "ya")),
    quiet(loo_model_weights(d, d, resp = "y_a")))
  expect_error(loo_epred(d, resp = "yb"), "must name one response")
})

test_that("add_criterion() hands each criterion its own arguments", {
  # brms passes every argument to every criterion and each ignores what
  # it does not read; frmtmb's loo() refused loo_R2()'s `seed`
  cs <- lp_case()
  d <- quiet(add_criterion(cs$ri, c("loo", "loo_R2"), seed = 1))
  expect_setequal(names(d$criteria), c("loo", "loo_R2"))
  expect_identical(d$criteria$loo_R2,
                   quiet(loo_R2(cs$ri, summary = FALSE, seed = 1)))
  expect_identical(d$criteria$loo$estimates, quiet(loo(cs$ri))$estimates)
  expect_error(add_criterion(cs$ri, "loo", seed = 1),
               "no requested criterion (loo) takes the argument `seed`",
               fixed = TRUE)
})

test_that("posterior_average() weights do not read a stored loo", {
  # BREAKING at this release, as brms: brms's model_weights() computes
  # log_lik() for each model and never reads a stored criterion.
  # frmtmb.sample 0.20.0 used a loo() stored with other draws.
  cs <- lp_case()
  d1 <- cs$ri
  d2 <- lp_draws(suppressWarnings(frm(bf(y ~ 1), family = gaussian(),
                                      data = cs$dd)))
  pa <- quiet(posterior_average(d1, d2, variable = "sigma", seed = 1))
  d1s <- quiet(add_criterion(d1, "loo", ndraws = 40))
  d2s <- quiet(add_criterion(d2, "loo", ndraws = 40))
  pas <- quiet(posterior_average(d1s, d2s, variable = "sigma", seed = 1))
  expect_identical(unname(attr(pas, "weights")),
                   unname(attr(pa, "weights")))
  # the guard sees something: the stored criteria do differ
  expect_false(identical(d1s$criteria$loo$estimates,
                         quiet(loo(d1))$estimates))
})

test_that("on brms's own draws, every loo_*() is brms's to rounding", {
  # The same posterior draws on both sides: a brms 2.23.0 fit's draws,
  # in their chains, in a frmtmb_draws over the ML fit of the same
  # model. No Monte Carlo is left in the comparison: every value is
  # brms's exactly under reference BLAS, and to rounding otherwise.
  skip_unless_brms_fit()
  set.seed(3)
  n <- 80
  dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:8, 10)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + stats::rnorm(8, 0, 0.7)[dd$g])
  fitb <- function(f) {
    suppressWarnings(suppressMessages(brms::brm(
      f, data = dd, family = brms::brmsfamily("gaussian"), chains = 2,
      iter = 1000, seed = 1, refresh = 0, silent = 2)))
  }
  bg <- fitb(y ~ x)
  br <- fitb(y ~ x + (1 | g))
  mk <- function(bfit, f) {
    fit <- suppressWarnings(frm(f, family = gaussian(), data = dd))
    labs <- frmtmb::brms_par_labels(fit)
    M <- as.matrix(bfit)
    # the group SD: the conditional log_lik() and epred never read it
    if ("theta_1" %in% labs) M <- cbind(M, theta_1 = M[, "sd_g__Intercept"])
    structure(list(stanfit = bfit$fit,
                   draws = cbind(M[, labs, drop = FALSE], lp__ = 0),
                   fit = fit), class = "frmtmb_draws")
  }
  dg <- mk(bg, bf(y ~ x))
  dr <- mk(br, bf(y ~ x + (1 | g)))
  bq <- function(e) suppressWarnings(suppressMessages(e))
  # Every bound below is a stated multiple of what this run measures
  # between the two sides: dl, the largest difference of the
  # log-likelihoods (their rounding apart: brms forms the predictor with
  # X %*% b and frmtmb row by row, which under OpenBLAS differ in the
  # last bits; under reference BLAS dl is 0), and dx, the largest
  # relative difference of the prediction matrix. A leave-one-out
  # summary of a row is a self-normalized weighted sum over the S
  # draws: each weight is exp() of a log-likelihood, so its relative
  # error is at most 2 * dl once normalized, and each value's is at
  # most dx. Those errors do not add up over the draws: a weighted mean
  # errs by about (2 dl + dx) of the largest value. Only the rounding
  # of the sum grows with S, by eps per term. The bounds below multiply
  # every term by S, so the dl and dx parts carry a headroom of S (1000
  # here). The old r_eff = 1 rule on a subset moved
  # loo_epred(draw_ids =) by 1e-4 of the largest value, far outside
  # this bound (dev/loodraws-oldreff-check.R puts it back and shows
  # the subset assertion failing).
  eps <- .Machine$double.eps
  N <- nrow(dd)
  near <- function(a, b, tol) {
    expect_identical(dimnames(a), dimnames(b))
    expect_lte(max(abs(a - b)), tol)
  }
  for (k in list(list(bg, dg), list(br, dr))) {
    llf <- unname(unclass(quiet(log_lik(k[[2]])))[, ])
    llb <- unname(bq(brms::log_lik(k[[1]])))
    S <- nrow(llb)
    dl <- max(abs(llf - llb))
    # the inputs themselves agree to rounding: a guard on the guard. A
    # log_lik() value sums nothing over the draws, so S is headroom over
    # its few eps of rounding, not a count of terms
    expect_lte(dl, S * eps * max(abs(llb)))
    xb <- bq(brms::posterior_epred(k[[1]]))
    dx <- max(abs(posterior_epred(k[[2]]) - xb)) / max(abs(xb))
    per_term <- 2 * dl + dx + eps
    ps <- quiet(psis(k[[2]]))
    w <- stats::weights(ps, log = FALSE)
    for (f in c("epred", "linpred")) {
      fun_f <- if (f == "epred") loo_epred else loo_linpred
      fun_b <- if (f == "epred") brms::loo_epred else brms::loo_linpred
      mb <- bq(fun_b(k[[1]]))
      near(quiet(fun_f(k[[2]])), mb, S * per_term * max(abs(xb)))
      # variance: (sum w x^2 - (sum w x)^2) / (1 - sum w^2), two S-term
      # sums of squares, each term erring by 2 dl + 2 dx + eps
      vb <- bq(fun_b(k[[1]], type = "var"))
      near(quiet(fun_f(k[[2]], type = "var")), vb,
           2 * S * (2 * dl + 2 * dx + eps) * max(xb^2) /
             (1 - max(colSums(w^2))))
      # quantile: loo interpolates linearly in the cumulative weight,
      # whose error is S * (2 dl + eps); the slope between neighboring
      # draws is their gap over the upper one's weight, L at most
      pr <- c(0.025, 0.5, 0.975)
      L <- max(vapply(seq_len(ncol(xb)), function(i) {
        o <- order(xb[, i])
        max(diff(xb[o, i]) / w[o[-1L], i])
      }, 1))
      qb <- bq(fun_b(k[[1]], type = "quantile", probs = pr))
      near(quiet(fun_f(k[[2]], type = "quantile", probs = pr)), qb,
           S * (2 * dl + eps) * L + (dx + eps) * max(abs(xb)))
    }
    # a subset: brms's one-chain r_eff on the same rows, S draws fewer
    ids <- seq(1L, 1000L, by = 3L)
    near(quiet(loo_epred(k[[2]], draw_ids = ids)),
         bq(brms::loo_epred(k[[1]], draw_ids = ids)),
         length(ids) * per_term * max(abs(xb)))
    # loo_R2: 1 - var_err / var_y over N rows reweighted per bootstrap
    # draw (brms's rexp() stream at seed 7, replayed here). The LOO
    # residuals err by the mean's bound e1; var_err by two N-term sums
    # of terms erring by 2 max|err| e1 + N eps max(err^2) (brms sums
    # with rowSums(sweep()), this package with a matrix product); the
    # ratio by var_err's error over the smallest var_y, plus var_y's
    # own N-term rounding times var_err / var_y^2
    rf <- quiet(loo_R2(k[[2]], seed = 7, summary = FALSE))
    rb <- bq(brms::loo_R2(k[[1]], seed = 7, summary = FALSE))
    y <- dd$y
    ylo <- quiet(loo_epred(k[[2]]))[, 1]
    err <- ylo - y
    e1 <- S * per_term * max(abs(xb))
    set.seed(7)
    bw <- matrix(stats::rexp(S * N), S, N)
    bw <- bw / rowSums(bw)
    vy <- (N / (N - 1)) * (drop(bw %*% y^2) - drop(bw %*% y)^2)
    ve <- (N / (N - 1)) * (drop(bw %*% err^2) - drop(bw %*% err)^2)
    dve <- (N / (N - 1)) * 2 * (2 * max(abs(err)) * e1 +
                                  N * eps * max(err^2))
    dvy <- (N / (N - 1)) * 2 * N * eps * max(y^2)
    near(rf, rb, max(dve / vy + ve * dvy / vy^2))
  }
  # Model weights. Each is a function of the pointwise elpd_loo, which
  # differ between the two sides by dp at most.
  lf <- list(quiet(loo(dg)), quiet(loo(dr)))
  lb <- list(bq(brms::loo(bg)), bq(brms::loo(br)))
  pf <- sapply(lf, function(l) l$pointwise[, "elpd_loo"])
  pb <- sapply(lb, function(l) l$pointwise[, "elpd_loo"])
  dp <- max(abs(pf - pb))
  # stacking: optim() stops when the objective's relative change is
  # below its reltol, sqrt(eps), so the two optima are compared through
  # the objective, sum_i log(sum_k w_k exp(elpd_ik)), at frmtmb's elpd:
  # within reltol of its size, plus N * dp for brms's elpd being dp off
  obj <- function(wt) sum(log(exp(pf) %*% as.numeric(wt)))
  set.seed(5)
  wf <- quiet(loo_model_weights(dg, dr))
  set.seed(5)
  wb <- bq(brms::loo_model_weights(bg, br))
  expect_identical(class(wf), class(wb))
  expect_lte(abs(obj(wf) - obj(wb)),
             sqrt(eps) * (abs(obj(wb)) + sqrt(eps)) + 2 * N * dp)
  # not a degenerate pair, which any two implementations would match
  w2 <- as.numeric(wf)
  expect_true(all(w2 > 0.01 & w2 < 0.99))
  # pseudo-BMA: a softmax of sums over the N rows, so an exponent errs by
  # N dp and a weight by 2 N dp, plus eps per Bayesian-bootstrap draw
  set.seed(5)
  pbf <- quiet(loo_model_weights(dg, dr, method = "pseudobma"))
  set.seed(5)
  pbb <- bq(brms::loo_model_weights(bg, br, method = "pseudobma"))
  expect_identical(class(pbf), class(pbb))
  expect_lte(max(abs(as.numeric(pbf) - as.numeric(pbb))),
             2 * N * dp + 1000 * eps)
})

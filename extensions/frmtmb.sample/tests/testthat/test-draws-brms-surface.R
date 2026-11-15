# brms 2.23.0's surface on draws: a new level in conditional_effects()
# lends an observed level's draw (sample_new_levels = "uncertainty"), an
# ordinal posterior_epred() names its draw margin "1".."S", the loo_*()
# functions take a subset of the fitted rows, and a latent residual term
# carries brms's err[k]. brms's values are in dev/parity2-findings.md.
# The draws are built by hand around the ML estimates, so no sampler is
# needed and every comparison is exact.

srf_draws <- function(fit, n = 40L, seed = 1L, sd = 0.05) {
  tpl <- fit$frame[["par_template"]]
  # the sampler holds no column for a fixed dpar (an ordinal disc)
  fx <- fit$frame[["betad_fixed_idx"]]
  est <- unlist(lapply(names(tpl), function(cp) {
    v <- fit$estimates[[cp]]
    if (cp == "betad" && length(fx)) v <- v[-fx]
    v
  }))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, sd),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit, data = fit$data),
            class = "frmtmb_draws")
}

srf_ri <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(11)
      dg <- data.frame(x = stats::rnorm(48), g = factor(rep(1:6, each = 8)))
      dg$y <- 1 + 0.5 * dg$x + stats::rnorm(6, sd = 0.8)[dg$g] +
        stats::rnorm(48, sd = 0.5)
      fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dg)
      cache <<- list(dg = dg, fit = fit, ds = srf_draws(fit))
    }
    cache
  }
})

test_that("a new level lends an observed level's draw by default", {
  cs <- srf_ri()
  ds <- cs$ds
  ic <- list(x = c(-1, 0, 1))
  sp <- function(...) {
    ce <- conditional_effects(ds, "x", int_conditions = ic, spaghetti = TRUE,
                              seed = 4, ...)
    s <- attr(ce$x, "spaghetti")
    matrix(s$estimate__, ncol = 3L, byrow = TRUE)
  }
  pop <- sp(re_formula = NA)
  new_u <- sp(re_formula = NULL)
  new_g <- sp(re_formula = NULL, sample_new_levels = "gaussian")
  r <- ds$draws[, grep("^r_g[[]", colnames(ds$draws)), drop = FALSE]
  shift_u <- new_u - pop
  shift_g <- new_g - pop
  eps <- .Machine$double.eps
  scale <- max(abs(new_u))
  # per draw, the shift is the same at every grid point and is one of
  # that draw's six level effects
  expect_lt(max(abs(shift_u - shift_u[, 1L])), 64 * eps * scale)
  hit <- vapply(seq_len(nrow(r)), function(i) {
    min(abs(r[i, ] - shift_u[i, 1L])) <= 64 * eps * scale
  }, NA)
  expect_true(all(hit))
  # brms's "gaussian" draws a value that is no level's effect
  miss <- vapply(seq_len(nrow(r)), function(i) {
    min(abs(r[i, ] - shift_g[i, 1L])) > 64 * eps * scale
  }, NA)
  expect_true(all(miss))
})

test_that("an ordinal posterior_epred() names its draw margin 1..S", {
  set.seed(21)
  dd <- data.frame(x = stats::rnorm(80))
  dd$y <- cut(dd$x + stats::rnorm(80), c(-Inf, -0.5, 0.5, Inf),
              labels = c("lo", "mid", "hi"), ordered_result = TRUE)
  fit <- frm(bf(y ~ x), family = cumulative(), data = dd)
  ds <- srf_draws(fit, n = 12L)
  ep <- posterior_epred(ds, newdata = dd[1:4, ])
  # brms 2.23.0: list(c("1", ..., "S"), NULL, <categories>)
  expect_identical(dimnames(ep)[[1L]], as.character(1:12))
  expect_null(dimnames(ep)[[2L]])
  expect_identical(dimnames(ep)[[3L]], c("lo", "mid", "hi"))
  # brms's categorical draw margin stays unnamed
  fc <- frm(bf(y ~ x), family = categorical(),
            data = transform(dd, y = factor(y, ordered = FALSE)))
  epc <- posterior_epred(srf_draws(fc, n = 12L), newdata = dd[1:4, ])
  expect_null(dimnames(epc)[[1L]])
})

test_that("the loo_*() functions take a subset of the fitted rows", {
  skip_if_not_installed("loo")
  cs <- srf_ri()
  ds <- cs$ds
  rows <- c(5L, 3L, 20L, 3L)
  sub <- cs$dg[rows, ]
  # 40 draws jittered about the estimates are no posterior, and loo says
  # so; the comparisons here are of columns, which those words do not
  # touch
  quiet <- function(e) {
    suppressMessages(allow_warnings(e, c("The ESS has been capped",
                                         "Pareto k diagnostic values")))
  }
  # brms 2.23.0: loo_epred() and loo_linpred() of d[1:10, ] equal the
  # fitted rows 1:10, 0 difference
  expect_identical(quiet(loo_epred(ds, newdata = sub)),
                   quiet(loo_epred(ds))[rows, , drop = FALSE])
  expect_identical(quiet(loo_linpred(ds, newdata = sub)),
                   quiet(loo_linpred(ds))[rows, , drop = FALSE])
  set.seed(8)
  a <- quiet(loo_predict(ds, newdata = sub))
  set.seed(8)
  b <- quiet(loo_predict(ds))
  expect_identical(a, b[rows, , drop = FALSE])
  set.seed(8)
  a <- quiet(loo_predictive_interval(ds, newdata = sub, prob = 0.8))
  set.seed(8)
  b <- quiet(loo_predictive_interval(ds, prob = 0.8))
  expect_identical(a, b[rows, , drop = FALSE])
  # the fitted rows themselves are still the call without newdata
  expect_identical(quiet(loo_epred(ds, newdata = cs$dg)),
                   quiet(loo_epred(ds)))
  r2 <- quiet(loo_R2(ds, newdata = cs$dg[1:30, ], seed = 1))
  expect_identical(dim(r2), c(1L, 4L))
  # a row the model was not fitted on is still refused
  moved <- cs$dg[1:3, ]
  moved$x[2] <- moved$x[2] + 1
  expect_error(quiet(loo_epred(ds, newdata = moved)), "cannot honor")
})

test_that("a latent residual's draws carry brms's err[k]", {
  set.seed(4)
  dp <- data.frame(t = rep(1:8, 6), g = factor(rep(1:6, each = 8)),
                   x = stats::rnorm(48))
  e <- as.vector(replicate(6, stats::arima.sim(list(ar = 0.5), 8,
                                               sd = 0.6)))
  dp$y <- stats::rpois(48, exp(0.5 + 0.3 * dp$x + e))
  # what the fit says about sderr is not the point here, only the names
  # and values of the draws built around it
  fit <- allow_warnings(
    frm(bf(y ~ x + ar(time = t, gr = g, cov = FALSE)), family = poisson(),
        data = dp),
    c("The Laplace approximation of ar(", "Standard errors are not"))
  ds <- srf_draws(fit, n = 10L)
  v <- colnames(ds$draws)
  err <- paste0("err[", 1:48, "]")
  zerr <- paste0("zerr[", 1:48, "]")
  # brms 2.23.0 saves err[1..48] beside sderr on this model
  expect_true(all(err %in% v))
  vs <- variables(ds)
  expect_true(all(err %in% vs))
  # brms 2.23.0 sorts err, a class reorder_pars() does not list, last:
  # ... sderr Intercept lprior lp__ err[1] ... err[48]
  expect_true("lp__" %in% vs)
  expect_gt(min(match(err, vs)), match("lp__", vs))
  M <- ds$draws
  eps <- .Machine$double.eps
  expect_lt(max(abs(M[, err] - M[, "sderr"] * M[, zerr])),
            64 * eps * max(abs(M[, err])))
  # the model reads the internal columns, and err is not one of them
  mi <- frmtmb.sample:::draws_internal_matrix(ds)
  expect_false(any(err %in% colnames(mi)))
  expect_true(all(zerr %in% colnames(mi)))
})

# cox()'s baseline simplex in the chart of a mo() simplex
# (dev/moparity-findings.md, item 4). On 0.69.0 (rellib-r7) the
# baseline was a softmax: of the 200 data sets of
# dev/moparity-cox-study.R, 5 stopped at "singular convergence (7)"
# with weights near 1e-50 and a warning, and seed 162 stopped 1.86e-4
# below the maximum with a standard error lost.
#
# Moved from frmtmb with cox(). The test that a face the likelihood
# rises from is left, and the maximum is not, reads an internal of core
# (simplex_face_leave()), so it stays in core as
# tests/testthat/test-simplex-pars.R, on a test family that declares a
# simplex the same way cox() does.

cox_study_data <- function(seed) {
  set.seed(seed)
  n <- 300
  x <- rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
}

# The maximum by another route: the intercept folded into the scale of
# the weights, which then need no sum constraint, maximized by the
# multiplicative (EM) update, which rises at every step, alternated
# with Newton steps on the slope (concave in eta).
cox_polished_loglik <- function(fit, dd) {
  av <- fit$frame$aterm_values[[1L]]
  Zb <- av$Zbhaz
  Zc <- av$Zcbhaz
  ev <- dd$cens == 0
  x <- dd$x
  s <- pmax(as.numeric(cox_baseline(fit)), 1e-300) *
    exp(fit$estimates$beta[1])
  b <- fit$estimates$beta[2]
  Zbe <- Zb[ev, , drop = FALSE]
  for (outer in 1:300) {
    cvec <- colSums(Zc * exp(b * x))
    for (it in 1:2000) s <- s * colSums(Zbe / as.vector(Zbe %*% s)) / cvec
    b_old <- b
    for (it in 1:50) {
      m <- as.vector(Zc %*% s) * exp(b * x)
      step <- (sum(x[ev]) - sum(x * m)) / -sum(x^2 * m)
      b <- b - step
      if (abs(step) < 1e-14) break
    }
    if (abs(b - b_old) < 1e-13) break
  }
  eta <- log(sum(s)) + b * x
  s <- s / sum(s)
  sum((log(as.vector(Zb %*% s)) + eta)[ev]) -
    sum(as.vector(Zc %*% s) * exp(eta))
}

test_that("a cox() baseline at a face reaches the maximum, silently", {
  # seed 162 stopped 1.86e-4 short and seed 40 at singular convergence
  # with a weight of 6e-27 on 0.69.0
  for (seed in c(162, 40)) {
    dd <- cox_study_data(seed)
    w <- character()
    fit <- withCallingHandlers(
      frm(bf(time | cens(cens) ~ x), family = cox(), data = dd),
      warning = function(z) {
        w <<- c(w, conditionMessage(z))
        invokeRestart("muffleWarning")
      })
    expect_identical(w, character(), info = paste("seed", seed))
    expect_identical(fit$opt$convergence, 0L, info = paste("seed", seed))
    ref <- cox_polished_loglik(fit, dd)
    gap <- ref - as.numeric(logLik(fit))
    expect_lt(gap, 1e-8 * abs(ref), label = paste("gap, seed", seed))
    # a weight the maximum puts at 0 is at 0 here too
    expect_lt(min(cox_baseline(fit)), 1e-8 * max(cox_baseline(fit)))
  }
})

test_that("a cox() baseline is sampled as brms's sbhaz under dirichlet(1)", {
  dd <- cox_study_data(2)
  fit <- frm(bf(time | cens(cens) ~ x), family = cox(), data = dd)
  terms <- frmtmb::simplex_frame_terms(fit)
  expect_length(terms, 1L)
  expect_identical(terms[[1L]]$names, paste0("sbhaz[", 1:5, "]"))
  fr <- frmtmb::mo_chart_frame(fit$frame, "softmax")
  expect_identical(fr$aterm_values[[1L]]$sbhaz_chart, "softmax")
  # the same model in other coordinates: equal objectives at the same
  # weights
  obj2 <- RTMB::MakeADFun(frmtmb::build_objective(fr),
                          utils::modifyList(fit$estimates, list(
                            sbhaz_raw = frmtmb::mo_coords(
                              cox_baseline(fit) + 1e-3, "softmax"))),
                          silent = TRUE)
  p <- fit$opt$par
  w1 <- (cox_baseline(fit) + 1e-3) / sum(cox_baseline(fit) + 1e-3)
  p[names(p) == "sbhaz_raw"] <- frmtmb::mo_coords(w1)
  v1 <- fit$obj$fn(p)
  v2 <- obj2$fn(obj2$par)
  expect_lt(abs(v1 - v2), 1e-10 * abs(v1))
  # brms's class sbhaz, as a MAP prior and in the default table
  tab <- get_prior(bf(time | cens(cens) ~ x), data = dd, family = cox())
  expect_true("sbhaz" %in% tab$class)
  fm <- frm(bf(time | cens(cens) ~ x), family = cox(), data = dd,
            prior = set_prior("dirichlet(2)", class = "sbhaz"))
  expect_gt(min(cox_baseline(fm)), 1e-3 * max(cox_baseline(fm)))
  expect_error(set_prior("dirichlet(2)", class = "sbhaz", coef = "x"),
               "sbhaz")
})

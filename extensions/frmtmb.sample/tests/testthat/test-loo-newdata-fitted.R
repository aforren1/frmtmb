# brms's `loo(fit, newdata = fit$data)` (brms row brmsfit-methods:469).
# brms computes the leave-one-out functions on log_lik(newdata =); on
# the rows the model was fitted on that is its answer without newdata.
# Through frmtmb.sample 0.27.0 draws had no `$data`, so the row passed
# with newdata NULL, and with the data given every newdata was refused.
# Now the fitted rows are accepted and other rows are still refused.
skip_on_cran()

nf_data <- function(seed = 469, n = 48L, ng = 6L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)),
                   unused = stats::runif(n))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x +
                         stats::rnorm(ng, 0, 0.5)[dd$g], 1)
  dd
}

nf_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      dd <- nf_data()
      fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 300, refresh = 0, seed = 469)))
      cache <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache
  }
})

test_that("draws carry the fit's model frame as $data, as brms's fit does", {
  cs <- nf_case()
  expect_true("data" %in% names(cs$ds))
  expect_identical(cs$ds$data, cs$fit$data)
  expect_identical(cs$ds$data, model.frame(cs$ds))
})

test_that("the fitted rows as newdata give the answer without newdata", {
  cs <- nf_case()
  skip_if_not_installed("loo")
  ds <- cs$ds
  ll <- log_lik(ds)
  l0 <- suppressWarnings(loo(ds))
  w0 <- suppressWarnings(waic(ds))
  # brms's spelling, the model frame, and the data the model was fitted
  # to, which holds a column the model does not read
  for (nd in list(ds$data, model.frame(ds), cs$dd)) {
    expect_identical(log_lik(ds, newdata = nd), ll)
    expect_identical(suppressWarnings(loo(ds, newdata = nd))$estimates,
                     l0$estimates)
    expect_identical(suppressWarnings(waic(ds, newdata = nd))$estimates,
                     w0$estimates)
  }
  expect_identical(suppressWarnings(loo_epred(ds, newdata = cs$dd)),
                   suppressWarnings(loo_epred(ds)))
  expect_identical(stats::weights(suppressWarnings(psis(ds, cs$dd))),
                   stats::weights(suppressWarnings(psis(ds))))
})

test_that("rows other than the fitted ones are still refused", {
  cs <- nf_case()
  skip_if_not_installed("loo")
  ds <- cs$ds
  dd <- cs$dd
  y1 <- dd
  y1$y[3] <- y1$y[3] + 1
  ulp <- dd
  ulp$x[5] <- ulp$x[5] * (1 + 2e-15)
  others <- list(response = y1, ulp = ulp, subset = dd[1:20, ],
                 reordered = dd[c(2, 1, 3:nrow(dd)), ],
                 extra_row = rbind(dd, dd[1, ]),
                 no_response = dd[c("x", "g")])
  for (nm in names(others)) {
    nd <- others[[nm]]
    expect_error(loo(ds, newdata = nd), "cannot honor `newdata`",
                 info = nm)
    expect_error(log_lik(ds, newdata = nd), "does not take newdata",
                 info = nm)
    expect_error(loo_epred(ds, newdata = nd), "cannot honor `newdata`",
                 info = nm)
  }
  e <- tryCatch(loo(ds, newdata = y1), error = conditionMessage)
  expect_match(e, "holds exactly the rows the model was fitted on",
               fixed = TRUE)
  # a subset is fitted rows, so its refusal gives its own reason
  es <- tryCatch(loo(ds, newdata = dd[1:20, ]), error = conditionMessage)
  expect_match(es, "rows are matched by value", fixed = TRUE)
})

test_that("a transformed term is compared by the variables it reads", {
  # no sampler: the check reads the fit's frame, not the draws. A frame
  # rebuilt from the data differs from the fitted poly() basis by 2.8e-16,
  # so the check compares x, which the fit stores beside the frame
  set.seed(4691)
  dd <- data.frame(x = stats::rnorm(40), age = stats::runif(40, 1, 3),
                   g = factor(rep(1:5, 8)))
  dd$y <- stats::rpois(40, exp(0.2 + 0.3 * dd$x + log(dd$age)))
  fit <- frm(bf(y ~ poly(x, 2) + offset(log(age)) + (1 | g)),
             family = poisson(), data = dd)
  ds <- structure(list(stanfit = NULL, draws = matrix(0, 1, 1), fit = fit),
                  class = "frmtmb_draws")
  is_fitted <- frmtmb.sample:::draws_newdata_is_fitted
  # the frame holds poly(x, 2) and offset(log(age)), not x
  expect_false("x" %in% names(fit$data))
  expect_true(is_fitted(ds, dd))
  expect_true(is_fitted(ds, fit$data))
  shifted <- dd
  shifted$x[7] <- shifted$x[7] + 1e-9
  expect_false(is_fitted(ds, shifted))
  older <- dd
  older$age[2] <- older$age[2] + 0.5
  expect_false(is_fitted(ds, older))
  expect_false(is_fitted(ds, NULL))
  expect_false(is_fitted(ds, as.list(dd)))
  # the frame with a variable it does not hold beside it: that variable
  # is compared too. Through the first review the frame's own columns
  # were all that was compared, and a different x was taken for the
  # fitted rows
  with_x <- fit$data
  with_x$x <- dd$x
  expect_true(is_fitted(ds, with_x))
  with_x$x <- dd$x + 1
  expect_false(is_fitted(ds, with_x))
})

test_that("a variable read outside the frame is compared on every path", {
  # cens(cc, y2): the frame holds y, x and cc, and y2 is stored beside it
  # (frame$raw_vars). A different upper bound of an interval row is
  # another likelihood (review 2026-11-10, item 1: logLik -59.313127
  # against -58.130862), so it is not the fitted rows
  set.seed(2)
  dd <- data.frame(x = stats::rnorm(40), cc = rep(c(0, 0, 0, 2), 10))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(40)
  dd$y2 <- dd$y + 1
  fit <- frm(y | cens(cc, y2) ~ x, data = dd)
  ds <- structure(list(stanfit = NULL, draws = matrix(0, 1, 1), fit = fit),
                  class = "frmtmb_draws")
  is_fitted <- frmtmb.sample:::draws_newdata_is_fitted
  expect_false("y2" %in% names(fit$data))
  expect_true(is_fitted(ds, dd))
  expect_true(is_fitted(ds, fit$data))
  moved <- dd
  moved$y2[4] <- moved$y2[4] + 3
  expect_false(is_fitted(ds, moved))
  # and with the frame's columns beside the moved bound
  framed <- fit$data
  framed$y2 <- moved$y2
  expect_false(is_fitted(ds, framed))
  # the frame without the bound, but not the frame object itself, cannot
  # show the bound it was fitted with
  expect_false(is_fitted(ds, as.data.frame(as.list(fit$data))))
})

test_that("a column of another type is not the fitted column", {
  # a numeric k read as a factor or as text builds another design
  set.seed(3)
  di <- data.frame(y = stats::rnorm(40), k = rep(1:4, 10))
  fit <- frm(y ~ k, data = di)
  ds <- structure(list(stanfit = NULL, draws = matrix(0, 1, 1), fit = fit),
                  class = "frmtmb_draws")
  is_fitted <- frmtmb.sample:::draws_newdata_is_fitted
  as_double <- di
  as_double$k <- as.double(di$k)
  expect_true(is_fitted(ds, as_double))
  as_factor <- di
  as_factor$k <- factor(di$k)
  expect_false(is_fitted(ds, as_factor))
  as_text <- di
  as_text$k <- as.character(di$k)
  expect_false(is_fitted(ds, as_text))
  # a grouping factor given as text with the same labels is the same rows
  dg <- data.frame(y = stats::rnorm(40), g = factor(rep(letters[1:5], 8)))
  fg <- frm(y ~ g, data = dg)
  dsg <- structure(list(stanfit = NULL, draws = matrix(0, 1, 1), fit = fg),
                   class = "frmtmb_draws")
  as_chr <- dg
  as_chr$g <- as.character(dg$g)
  expect_true(is_fitted(dsg, as_chr))
})

# A GEV row outside a draw's support: log_lik() of draws reads the exact
# density through row_lpdf(exact = TRUE), -Inf, where the fitted
# objective keeps a finite barrier (-2.1e82 and below at these cells,
# dev/rel081-log/loo-exact-r20.txt). The cells come from draws moved by
# hand: sampled draws put no row there (0 of 44,000 in 22 fits that read
# the barrier, dev/reviews/2026-11-10-release.md, finding 2). The fitted
# rows as newdata reach the same reader. loo(), waic() and psis() stop
# with one classed error; through 0.80.0 loo() and psis() stopped in the
# loo package and waic() gave an elpd of NaN.
gev_off_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
      set.seed(349)
      d <- data.frame(x = stats::rnorm(40))
      d$y <- rg(40, 1 + 0.3 * d$x, 1, 0.3)
      fit <- suppressWarnings(frm(y ~ x, family = gen_extreme_value(),
                                  data = d))
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 400, seed = 349, refresh = 0)))
      # three draws with xi above 0 get the mean that puts the smallest y
      # below their lower end: t = 1 + xi (y - mu) / sigma = -xi / sigma
      i <- which.min(d$y)
      dr <- ds$draws
      k <- which(dr[, "xi"] > 0.05)[1:3]
      for (j in k) {
        mu <- d$y[i] + dr[j, "sigma"] / dr[j, "xi"] + 1
        dr[j, "b_Intercept"] <- mu - dr[j, "b_x"] * d$x[i]
      }
      moved <- ds
      moved$draws <- dr
      cache <<- list(d = d, ds = ds, moved = moved, i = i, k = k)
    }
    cache
  }
})

test_that("log_lik() of draws is -Inf off the GEV support, with newdata too", {
  cs <- gev_off_case()
  ll <- log_lik(cs$moved)
  expect_false(anyNA(cs$k))
  expect_true(all(ll[cs$k, cs$i] == -Inf))
  # the other draws are the sampled ones, and finite everywhere
  expect_true(all(is.finite(ll[-cs$k, ])))
  expect_identical(ll[-cs$k, ], log_lik(cs$ds)[-cs$k, ])
  # brms's newdata = fit$data and the data fitted reach the same reader
  expect_identical(log_lik(cs$moved, newdata = cs$moved$data), ll)
  expect_identical(log_lik(cs$moved, newdata = cs$d), ll)
})

test_that("loo(), waic() and psis() stop on a log_lik() cell that is -Inf", {
  cs <- gev_off_case()
  skip_if_not_installed("loo")
  ll <- log_lik(cs$moved)
  ncol_bad <- sum(colSums(!is.finite(ll)) > 0L)
  counts <- sprintf("in %d cells, at %d of %d observations and 3 of %d draws",
                    sum(!is.finite(ll)), ncol_bad, ncol(ll), nrow(ll))
  calls <- list(
    loo = function(nd) loo(cs$moved, newdata = nd),
    waic = function(nd) waic(cs$moved, newdata = nd),
    psis = function(nd) psis(cs$moved, nd))
  for (nm in names(calls)) {
    for (nd in list(NULL, cs$moved$data)) {
      e <- expect_error(suppressWarnings(calls[[nm]](nd)),
                        class = "frmtmb_loglik_not_finite", info = nm)
      expect_match(conditionMessage(e), counts, fixed = TRUE, info = nm)
      expect_match(conditionMessage(e), "outside the support",
                   fixed = TRUE, info = nm)
    }
  }
  # the same draws without the moved ones: the guard's condition absent
  keep <- setdiff(seq_len(nrow(ll)), cs$k)
  expect_s3_class(suppressWarnings(waic(cs$moved, draw_ids = keep)), "waic")
  expect_s3_class(suppressWarnings(loo(cs$moved, draw_ids = keep)), "loo")
  expect_s3_class(suppressWarnings(loo(cs$ds)), "loo")
  # pp_check()'s loo_* types build the matrix on their own path
  skip_if_not_installed("bayesplot")
  for (ty in c("loo_pit_overlay", "loo_intervals", "loo_pit_qq")) {
    expect_error(suppressWarnings(pp_check(cs$moved, type = ty)),
                 class = "frmtmb_loglik_not_finite", info = ty)
  }
})

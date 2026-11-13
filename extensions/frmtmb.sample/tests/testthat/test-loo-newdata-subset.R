# brms answers log_lik(), loo() and waic() of a subset of the fitted
# rows with those rows' pointwise values: on brms 2.23.0, 0 difference
# from its fitted columns on y ~ x + (1 | g), and up to 0.46 on an ar()
# model, whose rows read their neighbors in newdata
# (dev/sampfront-brms-probe.R, dev/sampfront-findings.md). Through
# frmtmb.sample 0.28.0 every such subset was refused. The draws here are
# built by hand around the ML estimates, so no sampler is needed and the
# columns are exact.

sub_data <- function(seed = 4717L) {
  set.seed(seed)
  dd <- data.frame(g = factor(rep(1:6, each = 8L)), t = rep(1:8, 6L),
                   x = stats::rnorm(48L), w = stats::runif(48L, 0.5, 2))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(6L, 0, 0.7)[dd$g] +
    stats::rnorm(48L)
  dd$z <- -1 + 0.3 * dd$x + stats::rnorm(48L)
  dd
}

# n draws around the estimates; the layout frm_sample() gives the matrix
sub_draws <- function(fit, n = 8L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit, data = fit$data),
            class = "frmtmb_draws")
}

sub_case <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      dd <- sub_data()
      fit <- frm(bf(y | weights(w) ~ x + (1 | g)), family = gaussian(),
                 data = dd)
      cache <<- list(dd = dd, fit = fit, ds = sub_draws(fit))
    }
    cache
  }
})

test_that("log_lik() of a subset of the fitted rows is their columns", {
  cs <- sub_case()
  ds <- cs$ds
  ll <- log_lik(ds)
  pick <- function(i) {
    out <- unclass(ll)[, i, drop = FALSE]
    attr(out, "chain_id") <- attr(ll, "chain_id")
    out
  }
  # the data's rows, the frame's rows, a reordering, and a repeated row,
  # which make.unique() names "3.1"
  for (i in list(1:10, c(20, 5, 3), c(3, 3, 7), 48:1)) {
    info <- paste(i, collapse = " ")
    expect_identical(log_lik(ds, newdata = cs$dd[i, ]), pick(i),
                     info = info)
    expect_identical(log_lik(ds, newdata = ds$data[i, ]), pick(i),
                     info = info)
  }
  # the guard's condition absent: the fitted rows themselves are still
  # the matrix without newdata
  expect_identical(log_lik(ds, newdata = cs$dd), ll)
  expect_identical(logLik(ds, newdata = cs$dd[2:4, ]), pick(2:4))
  expect_identical(log_lik(ds, newdata = cs$dd[2:4, ], ndraws = 4),
                   unclass(log_lik(ds, ndraws = 4))[, 2:4, drop = FALSE])
})

test_that("loo(), waic(), psis() of a subset are those rows' values", {
  skip_if_not_installed("loo")
  cs <- sub_case()
  ds <- cs$ds
  rows <- c(2, 9, 17, 30, 41)
  nd <- cs$dd[rows, ]
  q <- function(expr) suppressWarnings(expr)
  l0 <- q(loo(ds))
  l1 <- q(loo(ds, newdata = nd))
  expect_identical(nrow(l1$pointwise), length(rows))
  # PSIS smooths each column on its own, and the relative efficiency is
  # taken per column, so a row's value does not depend on the others
  expect_identical(unname(l1$pointwise[, "elpd_loo"]),
                   unname(l0$pointwise[rows, "elpd_loo"]))
  w0 <- q(waic(ds))
  w1 <- q(waic(ds, newdata = nd))
  expect_identical(unname(w1$pointwise[, "elpd_waic"]),
                   unname(w0$pointwise[rows, "elpd_waic"]))
  p0 <- q(psis(ds))
  p1 <- q(psis(ds, nd))
  expect_identical(unname(p1$log_weights),
                   unname(p0$log_weights[, rows]))
  # brms's yhash is newdata's response, so two models on one subset
  # compare without loo's warning about different responses
  expect_identical(attr(l1, "yhash"), list(cs$dd$y[rows]))
})

test_that("loo_model_weights() of a subset weights those rows", {
  skip_if_not_installed("loo")
  cs <- sub_case()
  d1 <- cs$ds
  f2 <- frm(bf(y | weights(w) ~ 1 + (1 | g)), family = gaussian(),
            data = cs$dd)
  d2 <- sub_draws(f2, seed = 2L)
  rows <- 1:20
  q <- function(expr) suppressWarnings(expr)
  lls <- list(log_lik(d1, newdata = cs$dd[rows, ]),
              log_lik(d2, newdata = cs$dd[rows, ]))
  reff <- lapply(lls, function(l) {
    loo::relative_eff(exp(l), chain_id = attr(l, "chain_id"))
  })
  w <- q(loo_model_weights(d1, d2, newdata = cs$dd[rows, ]))
  ref <- q(loo::loo_model_weights(lls, r_eff_list = reff))
  expect_identical(as.numeric(w), as.numeric(ref))
})

test_that("rows without the fitted row names are placed by their values", {
  # Review 2026-11-11 of lane sampfront, finding 2: a tibble has no row
  # names and dplyr's verbs number the rows anew, so tibble[11:20, ] and
  # dplyr::filter() were refused. Such rows are emulated here by
  # renumbering a data frame's rows, which is what those verbs leave
  cs <- sub_case()
  ds <- cs$ds
  dd <- cs$dd
  ll <- log_lik(ds)
  pick <- function(i) {
    out <- unclass(ll)[, i, drop = FALSE]
    attr(out, "chain_id") <- attr(ll, "chain_id")
    out
  }
  renumbered <- function(i) {
    out <- dd[i, ]
    rownames(out) <- NULL
    out
  }
  for (i in list(11:20, which(dd$g == "2"), c(30, 7, 7, 41))) {
    info <- paste(i, collapse = " ")
    expect_identical(log_lik(ds, newdata = renumbered(i)), pick(i),
                     info = info)
  }
  # names that point at other rows do not mislead the placement
  other <- dd[1:3, ]
  rownames(other) <- c("5", "6", "nope")
  expect_identical(log_lik(ds, newdata = other), pick(1:3))
  # a row equal to another in every variable the model reads has the
  # same column, which is why any of them may place it
  d2 <- dd
  d2[41:43, ] <- d2[1:3, ]
  f2 <- frm(bf(y | weights(w) ~ x + (1 | g)), family = gaussian(),
            data = d2)
  ds2 <- sub_draws(f2)
  l2 <- unclass(log_lik(ds2))
  expect_identical(unname(l2[, 41:43]), unname(l2[, 1:3]))
})

test_that("rows that are not a subset of the fitted ones are refused", {
  skip_if_not_installed("loo")
  cs <- sub_case()
  ds <- cs$ds
  dd <- cs$dd
  moved <- dd[1:10, ]
  moved$x[2] <- moved$x[2] + 1
  renamed <- moved
  rownames(renamed) <- NULL
  heavier <- dd[1:3, ]
  heavier$w[1] <- heavier$w[1] * 2
  ulp <- dd[1:3, ]
  ulp$x[3] <- ulp$x[3] * (1 + 2e-15)
  lacking <- dd[1:3, c("x", "g", "w")]
  why <- list(moved = "its row 2 equals no fitted row",
              renamed = "its row 2 equals no fitted row",
              heavier = "its row 1 equals no fitted row",
              ulp = "its row 3 equals no fitted row",
              lacking = "it does not hold the variable y")
  others <- list(moved = moved, renamed = renamed, heavier = heavier,
                 ulp = ulp, lacking = lacking)
  for (nm in names(others)) {
    nd <- others[[nm]]
    e <- tryCatch(log_lik(ds, newdata = nd), error = conditionMessage)
    expect_match(e, "does not take newdata", fixed = TRUE, info = nm)
    expect_match(e, why[[nm]], fixed = TRUE, info = nm)
    e <- tryCatch(loo(ds, newdata = nd), error = conditionMessage)
    expect_match(e, "cannot honor `newdata`", fixed = TRUE, info = nm)
    expect_match(e, why[[nm]], fixed = TRUE, info = nm)
    expect_error(psis(ds, nd), "does not take newdata", info = nm)
  }
  # the loo_*() functions take the fitted rows or none, and say so
  expect_error(loo_epred(ds, newdata = dd[1:5, ]),
               "columns; the loo_*() functions do not", fixed = TRUE)
})

test_that("a subset is refused where its columns are not brms's answer", {
  dd <- sub_data()
  q <- function(...) suppressWarnings(suppressMessages(frm(...)))
  # an ar() row reads its neighbors in newdata (brms: 0.46 from the
  # fitted columns on the rows with t >= 3)
  fa <- q(bf(y ~ x + ar(time = t, gr = g)), data = dd)
  e <- tryCatch(log_lik(sub_draws(fa), newdata = dd[dd$t >= 3, ]),
                error = conditionMessage)
  expect_match(e, "the autocorrelation term ar(time = t, gr = g)",
               fixed = TRUE)
  # a me() term's latent value per row
  dm <- dd
  dm$sx <- 0.2
  fm <- q(bf(y ~ me(x, sx)), data = dm)
  e <- tryCatch(log_lik(sub_draws(fm), newdata = dm[1:5, ]),
                error = conditionMessage)
  expect_match(e, "does not take newdata", fixed = TRUE)
  # named as written, not by its column label mexsx
  expect_match(e, "the term me(x, sx), whose latent value per row brms",
               fixed = TRUE)
  # responses fitted on different rows under subset()
  ds2 <- dd
  ds2$keep <- rep(c(TRUE, FALSE), 24L)
  fs <- q(mvbf(bf(y | subset(keep) ~ x), bf(z ~ x), rescor = FALSE),
          data = ds2)
  e <- tryCatch(log_lik(sub_draws(fs), newdata = ds2[1:5, ], resp = "z"),
                error = conditionMessage)
  expect_match(e, "fitted on different rows under subset()", fixed = TRUE)
  # the guard's condition absent: the same model without the term
  fp <- q(bf(y ~ x), data = dd)
  dp <- sub_draws(fp)
  expect_identical(log_lik(dp, newdata = dd[dd$t >= 3, ]),
                   structure(unclass(log_lik(dp))[, dd$t >= 3,
                                                  drop = FALSE],
                             chain_id = rep(1L, 8L)))
})

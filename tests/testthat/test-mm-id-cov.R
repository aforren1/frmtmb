# mm() with an |ID| key and with cov = A: brms 2.23.0 builds both. The
# identities against brms's log_prob are in test-mi-mm-brms-l1.R (gated);
# these check the block, the keys and the refusals.

mmid_data <- function(seed = 70, n = 160) {
  set.seed(seed)
  lv <- letters[1:8]
  d <- data.frame(x = rnorm(n), x1 = rnorm(n), x2 = rnorm(n),
                  g1 = factor(sample(lv, n, TRUE), levels = lv),
                  g2 = factor(sample(lv, n, TRUE), levels = lv),
                  w1 = runif(n, 0.2, 1), w2 = runif(n, 0.2, 2))
  u <- stats::setNames(rnorm(8, 0, 0.8), lv)
  d$y <- 1 + 0.5 * d$x + 0.5 * (u[as.character(d$g1)] + u[as.character(d$g2)]) +
    rnorm(n, 0, 0.5)
  d
}

mmid_A <- function() {
  lv <- letters[1:8]
  A <- 0.5^abs(outer(1:8, 1:8, "-"))
  dimnames(A) <- list(lv, lv)
  A
}

test_that("an |ID| key merges mm() terms across predictors", {
  dd <- mmid_data()
  fr <- frm(bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2))),
            data = dd, dry_run = "frame")
  expect_length(fr$re_blocks, 1L)
  bk <- fr$re_blocks[[1L]]
  expect_identical(bk$dim, 2L)
  expect_identical(bk$covstruct, "us")
  expect_identical(bk$group_name, "mmg1g2")
  # one correlation between mu's and sigma's effects over the pooled levels
  expect_length(bk$theta_idx, 3L)
})

test_that("id = is the |ID| key, as under gr()", {
  dd <- mmid_data()
  a <- frm(bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ 1 + (1 | p | mm(g1, g2))),
           data = dd)
  b <- frm(bf(y ~ x + (1 | mm(g1, g2, id = "p")),
              sigma ~ 1 + (1 | mm(g1, g2, id = "p"))), data = dd)
  expect_identical(as.numeric(logLik(a)), as.numeric(logLik(b)))
  expect_error(frm(bf(y ~ (1 | q | mm(g1, g2, id = "p"))), data = dd,
                   dry_run = "frame"), "names its correlation id twice")
})

test_that("a key takes one set of weights, and follows the members", {
  dd <- mmid_data()
  f <- function(form) frm(form, data = dd, dry_run = "frame")
  # brms reads the first term's for every term, so the two orders are two
  # models: refused in both
  expect_error(f(bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))),
                    sigma ~ 1 + (1 | p | mm(g1, g2)))),
               "reads the first term's weights", fixed = TRUE)
  expect_error(f(bf(y ~ x + (1 | p | mm(g1, g2)),
                    sigma ~ 1 + (1 | p | mm(g1, g2,
                                            weights = cbind(w1, w2))))),
               "reads the first term's weights", fixed = TRUE)
  expect_error(f(bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))),
                    sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2),
                                            scale = FALSE)))),
               "reads the first term's weights", fixed = TRUE)
  # one set of weights still fits
  fit <- frm(bf(y ~ x + (1 | p | mm(g1, g2, weights = cbind(w1, w2))),
                sigma ~ 1 + (1 | p | mm(g1, g2, weights = cbind(w1, w2)))),
             data = dd)
  expect_length(fit$frame$re_blocks, 1L)
  expect_true(is.finite(as.numeric(logLik(fit))))
  # different members are different groups, which brms refuses too
  expect_error(frm(bf(y ~ x + (1 | p | mm(g1, g2)),
                      sigma ~ 1 + (1 | p | mm(g2, g1))),
                   data = dd, dry_run = "frame"),
               "more than one grouping")
  expect_error(frm(bf(y ~ x + (1 | p | g1), sigma ~ 1 + (1 | p | mm(g1, g2))),
                   data = dd, dry_run = "frame"),
               "more than one grouping")
})

test_that("weights under a key are compared by value, not by their text", {
  dd <- mmid_data()
  dd$W <- cbind(dd$w1, dd$w2)
  a <- frm(bf(y ~ x + (1 | mm(g1, g2, id = "p", weights = W)),
              sigma ~ 1 + (1 | mm(g1, g2, id = "p", weights = cbind(w1, w2)))),
           data = dd)
  b <- frm(bf(y ~ x + (1 | p | mm(g1, g2, weights = W)),
              sigma ~ 1 + (1 | p | mm(g1, g2, weights = W))), data = dd)
  expect_length(a$frame$re_blocks, 1L)
  expect_identical(as.numeric(logLik(a)), as.numeric(logLik(b)))
  # the columns the other way round are other values, so brms's reading of the
  # first term's would drop what the second wrote
  expect_error(frm(bf(y ~ x + (1 | mm(g1, g2, id = "p", weights = W)),
                      sigma ~ 1 + (1 | mm(g1, g2, id = "p",
                                          weights = cbind(w2, w1)))),
                   data = dd, dry_run = "frame"),
               "reads the first term's weights", fixed = TRUE)
})

test_that("conflicting cor under a key over mm() is refused in both orders", {
  dd <- mmid_data()
  f <- function(form) frm(form, data = dd, dry_run = "frame")
  msg <- "ask for different correlation"
  expect_error(f(bf(y ~ x + (1 | mm(g1, g2, id = "p")),
                    sigma ~ 1 + (1 | mm(g1, g2, id = "p", cor = FALSE)))), msg)
  expect_error(f(bf(y ~ x + (1 | mm(g1, g2, id = "p", cor = FALSE)),
                    sigma ~ 1 + (1 | mm(g1, g2, id = "p")))), msg)
  # the |p| spelling drops its key under cor = FALSE, as brms's program does,
  # with two uncorrelated blocks in either order
  fp <- f(bf(y ~ x + (1 | p | mm(g1, g2)),
             sigma ~ 1 + (1 | p | mm(g1, g2, cor = FALSE))))
  expect_length(fp$re_blocks, 2L)
  fp <- f(bf(y ~ x + (1 | p | mm(g1, g2, cor = FALSE)),
             sigma ~ 1 + (1 | p | mm(g1, g2))))
  expect_length(fp$re_blocks, 2L)
  fg <- f(bf(y ~ x + (1 | p | gr(g1, cor = FALSE)),
             sigma ~ 1 + (1 | p | g1)))
  expect_length(fg$re_blocks, 2L)
  # every term saying so has nothing to correlate, as in brms
  fr <- f(bf(y ~ x + (1 | mm(g1, g2, id = "p", cor = FALSE)),
             sigma ~ 1 + (1 | mm(g1, g2, id = "p", cor = FALSE))))
  expect_length(fr$re_blocks, 2L)
})

test_that("a cor = FALSE key over another group is refused in brms's words", {
  dd <- mmid_data()
  f <- function(form) frm(form, data = dd, dry_run = "frame")
  msg <- "Can only combine group-level terms of the same grouping factors"
  expect_error(f(bf(y ~ x + (1 | gr(g1, id = "p")),
                    sigma ~ 1 + (1 | gr(g2, id = "p", cor = FALSE)))), msg)
  expect_error(f(bf(y ~ x + (1 | gr(g1, id = "p", cor = FALSE)),
                    sigma ~ 1 + (1 | gr(g2, id = "p")))), msg)
  expect_error(f(bf(y ~ x + (1 | gr(g1, id = "p", cor = FALSE)),
                    sigma ~ 1 + (1 | gr(g2, id = "p", cor = FALSE)))), msg)
  # brms takes the |p| spelling, where cor = FALSE drops the key, as this does
  fr <- f(bf(y ~ x + (1 | p | gr(g1, cor = FALSE)),
             sigma ~ 1 + (1 | p | gr(g2, cor = FALSE))))
  expect_length(fr$re_blocks, 2L)
})

test_that("a merged mm() block predicts as its terms do", {
  dd <- mmid_data()
  # the data have no slope and no sigma effect, so the fit can stop at
  # the edge of the correlation space: what is read here is the design
  fit <- allow_warnings(frm(bf(y ~ x + (1 + mmc(x1, x2) | p | mm(g1, g2)),
                               sigma ~ 1 + (1 | p | mm(g1, g2))),
                            data = dd),
                        "Optimizer did not report convergence")
  # on its own rows, newdata rebuilds the design the fit read
  a <- allow_warnings(fitted(fit), "covariance of the estimates")
  b <- allow_warnings(fitted(fit, newdata = dd),
                      "covariance of the estimates")
  expect_lt(max(abs(a[, "Estimate"] - b[, "Estimate"])) /
              max(abs(a[, "Estimate"])), 1e-10)
  # a member new in newdata drops to the population value, the other
  # member keeps its effect
  nd <- dd[1:4, ]
  nd$g1 <- factor("zz")
  expect_error(fitted(fit, newdata = nd), "New levels")
  n <- allow_warnings(fitted(fit, newdata = nd, allow_new_levels = TRUE),
                      "covariance of the estimates")
  expect_true(all(is.finite(n[, "Estimate"])))
  expect_identical(names(ranef(fit)), "mmg1g2")
})

test_that("mm(cov = ) with the identity matrix is the plain mm() model", {
  dd <- mmid_data()
  id <- diag(8)
  dimnames(id) <- list(letters[1:8], letters[1:8])
  a <- frm(bf(y ~ x + (1 | mm(g1, g2))), data = dd)
  b <- frm(bf(y ~ x + (1 | mm(g1, g2, cov = id))), data = dd,
           data2 = list(id = id))
  expect_identical(b$frame$re_blocks[[1L]]$covstruct, "gr_cov")
  expect_identical(b$frame$re_blocks[[1L]]$group_name, "mmg1g2")
  expect_lt(abs(as.numeric(logLik(a)) - as.numeric(logLik(b))) /
              abs(as.numeric(logLik(a))), 1e-7)
  # and a matrix that is not the identity is another model: the guard
  # above would pass on a cov = that nothing reads
  c <- frm(bf(y ~ x + (1 | mm(g1, g2, cov = A))), data = dd,
           data2 = list(A = mmid_A()))
  expect_gt(abs(as.numeric(logLik(a)) - as.numeric(logLik(c))),
            1e-4 * abs(as.numeric(logLik(a))))
})

test_that("mm(cov = ) takes several columns, which A %x% Sigma correlates", {
  dd <- mmid_data()
  fit <- frm(bf(y ~ x + (1 + mmc(x1, x2) | mm(g1, g2, cov = A))), data = dd,
             data2 = list(A = mmid_A()))
  bk <- fit$frame$re_blocks[[1L]]
  expect_identical(bk$dim, 2L)
  expect_length(bk$theta_idx, 3L)
  a <- allow_warnings(fitted(fit), "covariance of the estimates")
  b <- allow_warnings(fitted(fit, newdata = dd),
                      "covariance of the estimates")
  expect_lt(max(abs(a[, "Estimate"] - b[, "Estimate"])) /
              max(abs(a[, "Estimate"])), 1e-10)
})

test_that("mm(cov = ) composes with an |ID| key over one matrix", {
  dd <- mmid_data()
  fit <- frm(bf(y ~ x + (1 | p | mm(g1, g2, cov = A)),
                sigma ~ 1 + (1 | p | mm(g1, g2, cov = A))), data = dd,
             data2 = list(A = mmid_A()))
  expect_length(fit$frame$re_blocks, 1L)
  expect_identical(fit$frame$re_blocks[[1L]]$covstruct, "gr_cov")
  # one with a matrix and one without mix two structures
  expect_error(frm(bf(y ~ x + (1 | p | mm(g1, g2, cov = A)),
                      sigma ~ 1 + (1 | p | mm(g1, g2))), data = dd,
                   data2 = list(A = mmid_A()), dry_run = "frame"),
               "mix covariance")
})

test_that("mm(cov = ) under || is one block per column, as gr(cov = )", {
  dd <- mmid_data()
  A <- mmid_A()
  fit <- frm(bf(y ~ (1 + x || mm(g1, g2, cov = A))), data = dd,
             data2 = list(A = A))
  cs <- vapply(fit$frame$re_blocks, `[[`, "", "covstruct")
  expect_identical(cs, c("gr_cov", "gr_cov"))
  # the same term with cor = FALSE
  f2 <- frm(bf(y ~ (1 + x | mm(g1, g2, cov = A, cor = FALSE))), data = dd,
            data2 = list(A = A))
  expect_identical(as.numeric(logLik(fit)), as.numeric(logLik(f2)))
  # the two pieces by hand
  f3 <- frm(bf(y ~ (1 | mm(g1, g2, cov = A)) + (0 + x | mm(g1, g2, cov = A))),
            data = dd, data2 = list(A = A))
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(f3))) /
              abs(as.numeric(logLik(f3))), 1e-7)
  # several columns in one piece have no structure to keep: as for gr()
  dd$fc <- factor(rep(c("p", "q", "r"), length.out = nrow(dd)))
  expect_error(frm(bf(y ~ (1 + fc || mm(g1, g2, cov = A))), data = dd,
                   data2 = list(A = A), dry_run = "frame"),
               "has a term with", fixed = TRUE)
})

test_that("mm(cov = ) refuses what its matrix cannot cover", {
  dd <- mmid_data()
  bad <- mmid_A()[1:5, 1:5]
  expect_error(frm(bf(y ~ (1 | mm(g1, g2, cov = bad))), data = dd,
                   data2 = list(bad = bad), dry_run = "frame"),
               "covering all grouping levels")
  expect_error(frm(bf(y ~ (1 | mm(g1, g2, cov = bad))), data = dd,
                   data2 = list(bad = matrix(1, 3, 2)), dry_run = "frame"),
               "square matrix")
  # brms's model for a by-split with a known covariance has no block here
  expect_error(frm(bf(y ~ (1 | mm(g1, g2, by = cbind(w1, w2), cov = A))),
                   data = dd, data2 = list(A = mmid_A()),
                   dry_run = "frame"),
               "not supported")
})

test_that("the prior table and the names are brms's", {
  skip_if_not_installed("brms")
  dd <- mmid_data()
  A <- mmid_A()
  rows <- function(tab) {
    tab <- as.data.frame(tab)
    tab <- tab[tab$class %in% c("sd", "cor"), , drop = FALSE]
    out <- as.data.frame(lapply(tab[, c("class", "coef", "group", "dpar")],
                                as.character))
    rownames(out) <- NULL
    out
  }
  expect_equal(
    rows(default_prior(bf(y ~ x + (1 | p | mm(g1, g2)),
                          sigma ~ 1 + (1 | p | mm(g1, g2))), data = dd)),
    rows(brms::default_prior(brms::bf(y ~ x + (1 | p | mm(g1, g2)),
                                      sigma ~ 1 + (1 | p | mm(g1, g2))),
                             data = dd)))
  expect_equal(
    rows(default_prior(bf(y ~ x + (1 | mm(g1, g2, cov = A))), data = dd,
                       data2 = list(A = A))),
    rows(brms::default_prior(brms::bf(y ~ x + (1 | mm(g1, g2, cov = A))),
                             data = dd, data2 = list(A = A))))
})

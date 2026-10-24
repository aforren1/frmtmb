# car() blocks moved from frmtmb's tests/testthat/test-predict-re-uncertainty.R,
# test-re-formula-partial.R and test-data2.R, unchanged.

test_that("a block whose b is not positionwise refuses to batch", {
  # Punch round 2, B1. The batch assumes expand_b() carries b to the
  # coefficients POSITIONWISE. Two blocks break that: an rr block's
  # coefficient is that level's factors through the loadings, and
  # car(type = "esicar")'s is b with its connected component's mean
  # removed, so perturbing every level at once moves each row by its
  # own step MINUS the mean of all of them. The earlier test was
  # length(c_idx) != length(b_idx), which catches rr below full rank
  # and misses both esicar and rr at full rank. dev/reunc-log/fdcar.txt
  # has the wrong answer that produced: 0.583 relative in sample.
  skip_on_cran()
  lat <- expand.grid(r = 1:3, c = 1:3)
  W <- matrix(0L, nrow(lat), nrow(lat))
  for (i in seq_len(nrow(lat))) {
    for (j in seq_len(nrow(lat))) {
      if (abs(lat$r[i] - lat$r[j]) + abs(lat$c[i] - lat$c[j]) == 1L) {
        W[i, j] <- 1L
      }
    }
  }
  dimnames(W) <- list(as.character(seq_len(nrow(lat))),
                      as.character(seq_len(nrow(lat))))
  set.seed(31)
  m <- 20
  d <- data.frame(g = factor(rep(seq_len(nrow(lat)), each = m)),
                  x = stats::rnorm(nrow(lat) * m))
  u <- stats::rnorm(nrow(lat), 0, 0.8)
  u <- u - mean(u)
  eta <- 0.6 * d$x + u[as.integer(d$g)] + stats::rlogis(nrow(d))
  d$y <- factor(cut(eta, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
                ordered = TRUE)
  fit <- suppressWarnings(frm(bf(y ~ x + car(W, gr = g, type = "esicar")) +
                                cumulative(), data = d))
  nd <- d[c(1, 25, 45), c("x", "g")]
  gov <- frmtmb:::re_governed_b(fit)
  expect_length(gov, nrow(lat))
  # the preconditions, so this cannot pass vacuously: the block IS an
  # esicar block and the predicate DOES call it non-positionwise
  bk <- fit$frame[["re_blocks"]][[1L]]
  expect_true(frmtmb:::block_expands(bk))
  expect_false(frmtmb:::block_b_positionwise(bk))
  # and the two consequences, which are what the answer depends on
  expect_null(frmtmb:::re_b_batches(fit, NULL, "y", FALSE, gov))
  expect_identical(frmtmb:::re_used_b(fit, nd, "y", FALSE), gov)
  # the answer itself, against one coefficient at a time
  for (n in list(NULL, nd)) {
    ref <- frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x, n),
                              b_idx = gov, b_batch = NULL)
    got <- if (is.null(n)) fitted(fit) else fitted(fit, newdata = n)
    expect_identical(unname(as.vector(ref)),
                     unname(as.vector(got[, "Est.Error", ])))
  }
  # icar does NOT center, so it still batches and is still right: the
  # guard must be a property of the TYPE and not of car()
  fi <- suppressWarnings(frm(bf(y ~ x + car(W, gr = g, type = "icar")) +
                               cumulative(), data = d))
  gi <- frmtmb:::re_governed_b(fi)
  expect_true(frmtmb:::block_b_positionwise(fi$frame[["re_blocks"]][[1L]]))
  expect_gt(length(frmtmb:::re_b_batches(fi, NULL, "y", FALSE, gi)), 0)
  ref <- frmtmb:::fit_fd_se(fi, function(x) frmtmb:::fitted_point(x),
                            b_idx = gi, b_batch = NULL)
  got <- unname(as.vector(fitted(fi)[, "Est.Error", ]))
  expect_lt(max(abs(unname(as.vector(ref)) - got) /
                  (.Machine$double.eps * abs(unname(as.vector(ref))))), 8)
})

test_that("a car() field is addressed by its grouping factor", {
  # a car() field is a group-level term with the bar (1 | loc), so a
  # formula keeps or drops it by that name; it was never the bar-less
  # content the pre-0.65.0 refusal described
  set.seed(12)
  rc <- expand.grid(r = 1:4, c = 1:4)
  W <- 1 * (as.matrix(stats::dist(rc, method = "manhattan")) == 1)
  rownames(W) <- colnames(W) <- paste0("s", seq_len(nrow(W)))
  d <- data.frame(loc = factor(rep(rownames(W), each = 6),
                               levels = rownames(W)),
                  g = factor(rep(1:4, 24)), h = factor(rep(1:3, 32)),
                  x = stats::rnorm(96))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(16, 0, 0.8)[as.integer(d$loc)] +
    stats::rnorm(4, 0, 0.5)[d$g] + stats::rnorm(3, 0, 0.5)[d$h] +
    stats::rnorm(96, 0, 0.4)
  fit <- suppressWarnings(frm(bf(y ~ x + car(W, gr = loc, type = "icar") +
                                   (1 | g) + (1 | h)), data = d))
  full <- frm_linpred(fit)
  # naming all three terms is the full prediction
  expect_identical(frm_linpred(fit, re_formula = ~(1 | g) + (1 | h) +
                                 (1 | loc)), full)
  # leaving (1 | loc) out drops the field: what it adds back is one value
  # per location, and it is not zero
  add_loc <- full - frm_linpred(fit, re_formula = ~(1 | g) + (1 | h))
  tol <- 64 * .Machine$double.eps * max(abs(full))
  within_loc <- tapply(add_loc, d$loc, function(v) diff(range(v)))
  expect_true(all(within_loc <= tol))
  expect_gt(max(abs(add_loc)), 1e3 * tol)
})

test_that("car() resolves its adjacency matrix from data2", {
  # moved from test-data2.R's "every structural special resolves
  # from data2", which keeps the other structural specials
  d2_lattice_W <- function(r, c) {
    g <- expand.grid(r = seq_len(r), c = seq_len(c))
    n <- nrow(g)
    W <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) {
          W[i, j] <- 1
        }
      }
    }
    dimnames(W) <- list(paste0("L", seq_len(n)),
                        paste0("L", seq_len(n)))
    W
  }
  set.seed(12)
  W <- d2_lattice_W(2, 3)
  nl <- nrow(W)
  cd <- data.frame(loc = factor(rep(rownames(W), each = 6),
                                levels = rownames(W)),
                   x = stats::rnorm(nl * 6))
  cd$y <- 1 + 0.4 * cd$x + stats::rnorm(nl, 0, 0.8)[as.integer(cd$loc)] +
    stats::rnorm(nrow(cd), 0, 0.5)
  f_car <- frm(bf(y ~ x + car(W, gr = loc)) + gaussian(), data = cd)
  f_car2 <- frm(bf(y ~ x + car(Wd2, gr = loc)) + gaussian(), data = cd,
                data2 = list(Wd2 = W))
  expect_equal(as.numeric(logLik(f_car2)), as.numeric(logLik(f_car)),
               tolerance = 1e-10)

})


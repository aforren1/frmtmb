# predict() and fitted() carry the group effects' uncertainty at a
# grouping level the fit SAW, as brms's draws carry the posterior of r_g
# there (user decision, 2026-09-22). predict() draws each replicate's
# effects JOINTLY from their conditional law given the data, so two rows
# of one group share one draw; dev/reunc-findings.md has the
# construction and the coverage it claims.
#
# Every Monte Carlo comparison below is judged against its own Monte
# Carlo standard error, measured from the draws, never against a fixed
# number.

reunc_fixture <- function() {
  set.seed(4)
  G <- 6
  m <- 4
  d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                  x = stats::rnorm(G * m))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(G, 0, 1)[d$g] +
    stats::rnorm(G * m, 0, 0.8)
  d
}

reunc_fit <- function(mode) {
  d <- reunc_fixture()
  switch(mode,
    ML = frm(bf(y ~ x + (1 | g)), data = d),
    REML = frm(bf(y ~ x + (1 | g)), data = d, REML = TRUE),
    profile = frm(bf(y ~ x + (1 | g)), data = d,
                  control = frmtmb_control(profile = TRUE)))
}

# the sampling standard error of a covariance estimated from N draws of
# a bivariate normal pair
cov_se <- function(C, i, j, N) sqrt((C[i, i] * C[j, j] + C[i, j]^2) / N)

test_that("the group-effect draw is the analytic conditional law", {
  # The parameters have to be held at their estimates for the
  # conditional law of a one-way intercept to be analytic, and since
  # the 2026-09-23 rename `propagate_error = FALSE` holds the group
  # effects too, so there is no longer a public call that draws b
  # alone. The drawer is therefore driven directly, which is the thing
  # this claim is about: `predict()`'s use of it is covered by the
  # blocks below.
  N <- 8000
  for (mode in c("ML", "REML", "profile")) {
    fit <- reunc_fit(mode)
    # Var(b_g | y) = sigma^2 tau^2 / (sigma^2 + n_g tau^2), n_g = 4
    s2 <- sigma(fit)^2
    t2 <- exp(2 * fit$estimates[["theta"]][1])
    v <- s2 * t2 / (s2 + 4 * t2)
    set.seed(1)
    drawer <- frmtmb:::predict_b_drawer(fit, NULL, NULL, N, TRUE)
    expect_false(is.null(drawer))
    nb <- length(fit$estimates[["b"]])
    B <- t(vapply(seq_len(N), function(s) {
      drawer(fit, s)$estimates[["b"]]
    }, numeric(nb)))
    C <- stats::cov(B)
    # each level's own conditional variance, and levels independent
    expect_lt(abs(C[1, 1] - v) / (C[1, 1] * sqrt(2 / N)), 4, label = mode)
    expect_lt(abs(C[2, 2] - v) / (C[2, 2] * sqrt(2 / N)), 4, label = mode)
    expect_lt(abs(C[1, 2]) / cov_se(C, 1, 2, N), 4, label = mode)
  }
})

test_that("propagate_error = FALSE holds the group effects too", {
  # The 2026-09-23 decision. The argument names an axis, whether the
  # error in the estimates is propagated, and a group effect is
  # estimated, so FALSE holds it. Two rows of ONE group then carry no
  # shared random component at all and covary only through the noise,
  # which is zero. Under TRUE they share that group's draw, which the
  # block above and the one below measure.
  nd <- data.frame(x = c(0, 1, 0, -1),
                   g = factor(c(1, 1, 2, 3), levels = as.character(1:6)))
  N <- 8000
  for (mode in c("ML", "REML", "profile")) {
    fit <- reunc_fit(mode)
    s2 <- sigma(fit)^2
    set.seed(1)
    P <- predict(fit, newdata = nd, summary = FALSE, ndraws = N,
                 propagate_error = FALSE)
    C <- stats::cov(P)
    # nothing is shared between any pair of rows
    for (ij in list(c(1, 2), c(1, 3), c(3, 4))) {
      expect_lt(abs(C[ij[1], ij[2]]) / cov_se(C, ij[1], ij[2], N), 4,
                label = paste(mode, ij[1], ij[2]))
    }
    # and one row carries the observation noise ALONE, with no
    # conditional variance of the effect added on top
    expect_lt(abs(C[1, 1] - s2) / (C[1, 1] * sqrt(2 / N)), 4, label = mode)
    # The centre still carries the group's own effect, which is the
    # half re_formula cannot express. With everything held and the
    # same seed the two calls simulate from the same stream, so the
    # shift is not a mean over draws but an EXACT identity, draw by
    # draw: the only difference is whether b enters the predictor.
    set.seed(1)
    pop <- predict(fit, newdata = nd, summary = FALSE, ndraws = N,
                   propagate_error = FALSE, re_formula = NA)
    b1 <- fit$estimates[["b"]][1]   # row 1 is level 1 of the only block
    expect_lt(max(abs((P[, 1] - pop[, 1]) - b1)) / sqrt(s2), 1e-10,
              label = mode)
    # and that effect is not zero on this fixture, or the line above
    # would assert nothing
    expect_gt(abs(b1) / sqrt(s2), 0.05, label = mode)
  }
})

test_that("with the parameter draw, rows covary as the joint law says", {
  nd <- data.frame(x = c(0, 1, 0, -1),
                   g = factor(c(1, 1, 2, 3), levels = as.character(1:6)))
  N <- 8000
  for (mode in c("ML", "REML")) {
    fit <- reunc_fit(mode)
    set.seed(2)
    P <- predict(fit, newdata = nd, summary = FALSE, ndraws = N)
    C <- stats::cov(P)
    # a_i' V a_j over (beta, b): the fixed effects correlate rows of
    # DIFFERENT groups too, through the shared intercept
    lb <- frm_lp_basis(fit, newdata = nd)
    M <- lb$A %*% lb$V %*% t(lb$A)
    for (ij in list(c(1, 2), c(1, 3), c(3, 4))) {
      i <- ij[1]
      j <- ij[2]
      expect_lt(abs(C[i, j] - M[i, j]) / cov_se(C, i, j, N), 4,
                label = paste(mode, i, j))
    }
  }
})

test_that("fitted() at a known level carries the group effect", {
  fit <- reunc_fit("ML")
  nd <- data.frame(x = c(0, 1),
                   g = factor(c(1, 2), levels = as.character(1:6)))
  lb <- frm_lp_basis(fit, newdata = nd)
  # the delta method over (beta, b) jointly
  expect_equal(unname(fitted(fit, newdata = nd)[, "Est.Error"]),
               unname(sqrt(rowSums((lb$A %*% lb$V) * lb$A))))
  bcols <- lb$coef_pos %in% which(frm_joint_cov(fit)$names == "b")
  expect_true(any(bcols))
})

test_that("an ordinal fit's category probabilities carry it too", {
  set.seed(5)
  G <- 10
  m <- 12
  d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                  x = stats::rnorm(G * m))
  lat <- 0.8 * d$x + stats::rnorm(G, 0, 1)[d$g] + stats::rlogis(G * m)
  d$y <- factor(cut(lat, c(-Inf, -1, 0.5, Inf), labels = FALSE),
                ordered = TRUE)
  fit <- frm(bf(y ~ x + (1 | g)) + cumulative(), data = d)
  nd <- data.frame(x = 0, g = factor(1, levels = levels(d$g)))
  f <- fitted(fit, newdata = nd)[, "Est.Error", ]
  # brute force over the same joint law of (outer parameters, b),
  # pushed through fitted()'s own estimate
  ds <- frmtmb:::fit_draw_space(fit)
  jc <- frm_joint_cov(fit)
  pos <- c(unlist(lapply(unique(ds$map$comp),
                         function(cp) which(jc$names == cp))),
           which(jc$names == "b"))
  L <- chol(as.matrix(jc$V[pos, pos]))
  v0 <- c(frmtmb:::fit_outer_vector(fit, ds$map), fit$estimates[["b"]])
  p <- length(ds$map$names)
  N <- 1500
  set.seed(7)
  draws <- t(vapply(seq_len(N), function(s) {
    v <- v0 + drop(crossprod(L, stats::rnorm(length(v0))))
    fs <- frmtmb:::fit_set_outer(fit, v[seq_len(p)], ds$map)
    fs$estimates[["b"]][] <- v[-seq_len(p)]
    as.vector(frmtmb:::fitted_point(fs, nd))
  }, numeric(3)))
  mc <- apply(draws, 2, stats::sd)
  # the delta method is first order and a probability is curved, so the
  # two agree to a ratio, not to Monte Carlo error; before the b columns
  # joined the difference the ratio was about 0.6
  r <- unname(f / mc)
  expect_true(all(r > 0.85 & r < 1.18),
              info = paste(round(r, 3), collapse = " "))
})

test_that("a row at an unseen level is untouched by the known-level draw", {
  fit <- reunc_fit("ML")
  nd <- data.frame(x = c(0.5, 0, 1),
                   g = factor(c("new", "1", "2"), levels = c(1:6, "new")))
  set.seed(6)
  both <- predict(fit, newdata = nd, summary = FALSE, ndraws = 200,
                  allow_new_levels = TRUE)
  set.seed(6)
  alone <- predict(fit, newdata = nd[1, ], summary = FALSE, ndraws = 200,
                   allow_new_levels = TRUE)
  expect_identical(both[, 1], alone[, 1])
})

test_that("a quadrature fit says it cannot draw the group effects", {
  set.seed(10)
  d <- data.frame(g = factor(rep(1:10, each = 6)), x = stats::rnorm(60))
  d$y <- stats::rbinom(60, 1, stats::plogis(d$x + stats::rnorm(10)[d$g]))
  fit <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + bernoulli(), data = d,
                              quadrature = TRUE))
  expect_warning(predict(fit, ndraws = 20), class = "frmtmb_warning",
                 regexp = "quadrature")
  expect_no_warning(predict(fit, ndraws = 20, re_formula = NA))
})

test_that("the group-effect draw costs the caller no random numbers", {
  # Punch round 1, M2. The drawer takes one seed per replicate. Round 1
  # took them before predict() captured the .Random.seed it restores on
  # exit, so a conditional predict() returned identical values and left
  # the caller's stream ndraws uniforms further on, and the next
  # simulate() or sample() in a script gave different data than 0.61.0.
  # re_formula = NA is the reference because its drawer returns before
  # taking anything, so the two calls must leave the stream in the same
  # place. dev/reunc-log/stream.txt reconstructs the round-1 ordering
  # and shows this failing.
  fit <- reunc_fit("ML")
  nd <- reunc_fixture()[c(1, 5, 9), c("x", "g")]
  after <- function(expr) {
    set.seed(4242)
    force(expr)
    list(seed = get(".Random.seed", envir = globalenv()), u = stats::runif(1))
  }
  na <- after(predict(fit, newdata = nd, re_formula = NA, ndraws = 25))
  cond <- after(predict(fit, newdata = nd, ndraws = 25))
  expect_identical(cond$seed, na$seed)
  expect_identical(cond$u, na$u)
  # and at an unseen level, where the drawer runs but reaches nothing
  ndn <- data.frame(x = c(0.5, 0), g = factor(c("new", "1"),
                                              levels = c(1:6, "new")))
  new <- after(predict(fit, newdata = ndn, ndraws = 25,
                       allow_new_levels = TRUE))
  nan <- after(predict(fit, newdata = ndn, ndraws = 25, re_formula = NA,
                       allow_new_levels = TRUE))
  expect_identical(new$seed, nan$seed)
})

test_that("bounding and batching the b difference changes no answer", {
  # Punch round 1, M1. The finite-difference route used to perturb every
  # kept level, one at a time. It now perturbs only the levels the
  # predicted rows load (re_used_b()) and perturbs a block's levels
  # together (re_b_batches()). Both must give the SAME number, so this
  # compares the shipped call against the one-at-a-time full set.
  # dev/reunc-log/fdcost-lane.txt carries the cost this buys.
  skip_on_cran()
  set.seed(5)
  ng <- 20
  d <- data.frame(g = factor(rep(seq_len(ng), each = 6)),
                  x = stats::rnorm(ng * 6))
  lat <- 0.8 * d$x + stats::rnorm(ng, 0, 1)[d$g] + stats::rlogis(ng * 6)
  d$y <- factor(cut(lat, c(-Inf, -1, 0.5, Inf), labels = FALSE),
                ordered = TRUE)
  fit <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + cumulative(), data = d))
  nd <- d[c(1, 61), c("x", "g")]
  gov <- frmtmb:::re_governed_b(fit)
  used <- frmtmb:::re_used_b(fit, nd, "y", FALSE)
  # the test is worthless if the two sets are the same set
  expect_length(gov, ng)
  expect_length(used, 2L)
  # tolerated as a multiple of the values' OWN size: the batch adds the
  # same terms in a different order, so it agrees to a few ulps and not
  # bitwise
  ulp <- function(a, b) {
    a <- unname(as.vector(a))
    b <- unname(as.vector(b))
    max(abs(a - b) / (.Machine$double.eps * abs(a)))
  }
  full_nd <- frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x, nd),
                                b_idx = gov)
  expect_lt(ulp(full_nd, fitted(fit, newdata = nd)[, "Est.Error", ]), 8)
  full_in <- frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x),
                                b_idx = gov)
  expect_lt(ulp(full_in, fitted(fit)[, "Est.Error", ]), 8)
  # a multi-membership term puts two levels of one block on one row, so
  # the batch cannot attribute the difference and must refuse
  d$h <- factor(rep(1:5, length.out = nrow(d)))
  d$h2 <- factor(rep(c(2:5, 1), length.out = nrow(d)))
  fmm <- suppressWarnings(frm(bf(y ~ x + (1 | mm(h, h2))) + cumulative(),
                              data = d))
  bmm <- frmtmb:::re_governed_b(fmm)
  expect_null(frmtmb:::re_b_batches(fmm, NULL, "y", FALSE, bmm))
  expect_identical(
    unname(as.vector(frmtmb:::fit_fd_se(
      fmm, function(x) frmtmb:::fitted_point(x), b_idx = bmm))),
    unname(as.vector(fitted(fmm)[, "Est.Error", ])))
})

test_that("a factor smooth keeps its contribution on newdata", {
  # Punch round 1, M1, found in the punch round's own review. The batch
  # attributes a row's difference through re_row_support(), which on
  # newdata read only the ordinary group parts. A factor smooth's basis
  # is in the smooth parts, so its columns looked like coefficients no
  # row loads and the batch wrote them a derivative of zero. With an
  # ordinary group term beside it, to give the batch something to
  # attribute, Est.Error came out 0.676 of the way wrong and said
  # nothing. dev/reunc-log/fdsmooth-prefix.txt is that run.
  skip_on_cran()
  set.seed(21)
  ng <- 6
  d <- data.frame(g = factor(rep(seq_len(ng), each = 25)),
                  x = stats::runif(ng * 25, -2, 2))
  lat <- stats::rnorm(ng, 0, 1)[d$g] * sin(d$x) + 0.5 * d$x +
    stats::rlogis(nrow(d))
  d$y <- factor(cut(lat, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
                ordered = TRUE)
  d$h <- factor(rep(1:5, length.out = nrow(d)))
  fit <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | h)) +
                                cumulative(), data = d))
  nd <- d[c(3, 40, 90), c("x", "g", "h")]
  gov <- frmtmb:::re_governed_b(fit)
  used <- frmtmb:::re_used_b(fit, nd, "y", FALSE)
  # The preconditions. Without them a future change that made the batch
  # REFUSE on this fit would fall back to one coefficient at a time,
  # agree with the reference and pass, and the guard would stop
  # guarding anything. The defect lived in the ATTRIBUTION, which only
  # exists once a batch is formed, and it needed a smooth column to be
  # among the batched ones.
  bt <- frmtmb:::re_b_batches(fit, nd, "y", FALSE, used)
  expect_false(is.null(bt))
  expect_gt(length(bt), 0)
  sm <- fit$frame[["re_blocks"]][[1L]][["b_idx"]]
  expect_true(any(vapply(bt, function(b) any(b$idx %in% sm), logical(1))))
  ref <- frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x, nd),
                            b_idx = gov)
  got <- fitted(fit, newdata = nd)[, "Est.Error", ]
  rel <- max(abs(unname(as.vector(ref)) - unname(as.vector(got))) /
               abs(unname(as.vector(ref))))
  # relative to the values the run itself produced; the defect showed
  # here as 0.676, and on another design it OVERSTATED a cell instead
  expect_lt(rel, 8 * .Machine$double.eps)
})

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
  expect_true(frmtmb:::block_is_esicar(bk))
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

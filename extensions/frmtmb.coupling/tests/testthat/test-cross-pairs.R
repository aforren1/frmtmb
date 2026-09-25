## Item 3.2 of dev/extension-gaps-plan.md: frm_cross_pairs() stacks the
## channel pairs of one recording with a `pair` factor, so that one
## model with `(1 | pair)` can hold them all.
##
## As in test-epochs.R, calls are caught so that on a build without the
## function each assertion fails rather than the block stopping at the
## first error.

cls <- "frmtmb_coupling_error"

## A refusal, asserted so that on a build without this code it FAILS
## rather than stopping the block: expect_error(class = ) rethrows an
## error of another class, and every assertion after it would never run.
refuses <- function(expr, regexp) {
  cnd <- tryCatch(expr, error = identity)
  expect_s3_class(cnd, cls)
  expect_match(if (inherits(cnd, "condition")) conditionMessage(cnd) else
                 "", regexp)
}

xp_run <- function(expr) tryCatch(expr, error = identity)

## the largest distance on the logit scale, where frm_coherence()'s
## `.se` lives, in the reference's own standard errors; NA when the call
## did not return a frame, so the assertion fails
xp_dist <- function(a, ref) {
  if (!is.data.frame(a)) return(NA_real_)
  max(abs(stats::qlogis(a$.estimate) - stats::qlogis(ref$.estimate)) /
        ref$.se)
}

## the stacked fits stop at nlminb's relative tolerance with a gradient
## near core's 1e-3 warning line (measured 3.1e-03 and 2.1e-03 in
## dev/phase3a-log/phase3a-pairs-fit.txt), so that warning is let
## through and any OTHER warning still fails the test
xp_fit <- function(...) {
  w <- character(0)
  fit <- withCallingHandlers(
    xp_run(frmtmb::frm(...)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, other = w[!grepl("maximum absolute gradient", w)])
}

## Four channels driven by one source with different loadings, so the
## six pairs have six different coherences, flat in frequency.
xp_rec <- function(n = 4096L, load = c(1, 0.8, 0.5, 0.25)) {
  src <- stats::rnorm(n)
  X <- vapply(load, function(l) l * src + stats::rnorm(n), numeric(n))
  colnames(X) <- c("Fz", "Cz", "Pz", "Oz")
  X
}

test_that("each block is frm_cross_spectrum() on that pair", {
  set.seed(41)
  X <- xp_rec(2048L)
  xp <- xp_run(frm_cross_pairs(X, segments = 8L, sfreq = 256))
  expect_s3_class(xp, "data.frame")
  want <- c("Fz-Cz", "Fz-Pz", "Fz-Oz", "Cz-Pz", "Cz-Oz", "Pz-Oz")
  expect_identical(levels(xp$pair), want)
  expect_identical(names(xp)[1:3], c("pair", "ch1", "ch2"))
  for (p in want) {
    ab <- strsplit(p, "-", fixed = TRUE)[[1L]]
    ref <- frm_cross_spectrum(X[, ab[1L]], X[, ab[2L]], segments = 8L,
                              sfreq = 256)
    got <- if (is.data.frame(xp)) xp[xp$pair %in% p, -(1:3)]
    if (is.data.frame(got)) rownames(got) <- NULL
    expect_identical(got, ref, label = p)
  }
})

test_that("pairs can be chosen by name or by column number", {
  set.seed(42)
  X <- xp_rec(1024L)
  a <- xp_run(frm_cross_pairs(X, pairs = rbind(c("Pz", "Fz"),
                                                c("Cz", "Oz")),
                              segments = 4L))
  b <- xp_run(frm_cross_pairs(X, pairs = rbind(c(3, 1), c(2, 4)),
                              segments = 4L))
  expect_s3_class(a, "data.frame")
  expect_identical(a, b)
  expect_identical(levels(a$pair), c("Pz-Fz", "Cz-Oz"))
  # ch1 is the x signal, so the order in `pairs` is the phase's sign
  ref <- frm_cross_spectrum(X[, "Pz"], X[, "Fz"], segments = 4L)
  got <- if (is.data.frame(a)) a[a$pair %in% "Pz-Fz", -(1:3)]
  if (is.data.frame(got)) rownames(got) <- NULL
  expect_identical(got, ref)
  # a data frame of channels is the same recording
  d <- xp_run(frm_cross_pairs(as.data.frame(X), segments = 4L))
  m <- xp_run(frm_cross_pairs(X, segments = 4L))
  expect_s3_class(d, "data.frame")
  expect_identical(d, m)
})

test_that("a list of epochs goes through with its group", {
  set.seed(43)
  ep <- list(xp_rec(700L), xp_rec(900L), xp_rec(800L))
  g <- c("s1", "s1", "s2")
  xp <- xp_run(frm_cross_pairs(ep, pairs = rbind(c("Fz", "Cz")),
                               segments = 4L, group = g))
  ref <- xp_run(frm_cross_spectrum(lapply(ep, function(m) m[, "Fz"]),
                                   lapply(ep, function(m) m[, "Cz"]),
                                   segments = 4L, group = g))
  expect_s3_class(ref, "data.frame")
  got <- if (is.data.frame(xp)) xp[-(1:3)]
  expect_identical(got, ref)
})

test_that("frm_cross_pairs() refuses what it cannot pair", {
  set.seed(44)
  X <- xp_rec(512L)
  refuses(frm_cross_pairs(X, pairs = rbind(c("Fz", "Fz"))),
               "with itself")
  refuses(frm_cross_pairs(X, pairs = rbind(c("Fz", "Cz"),
                                                c("Cz", "Fz"))),
               "repeats the pair")
  refuses(frm_cross_pairs(X, pairs = rbind(c("Fz", "T7"))),
               "T7, which is not a channel")
  refuses(frm_cross_pairs(X, pairs = rbind(c(1, 5))),
               "not a column of `X`")
  refuses(frm_cross_pairs(X, pairs = c("Fz", "Cz")),
               "two-column matrix")
  refuses(frm_cross_pairs(X, segmnets = 4L),
               "unknown argument `segmnets`")
  refuses(frm_cross_pairs(X, NULL, 4L),
               "must be named")
  refuses(frm_cross_pairs(X, x = 1),
               "unknown argument `x`")
  refuses(frm_cross_pairs(X[, 1L, drop = FALSE]),
               "pair needs two")
  refuses(frm_cross_pairs(data.frame(a = 1:10, b = letters[1:10])),
               "non-numeric column b")
  refuses(frm_cross_pairs(list(X, X[, 4:1])),
               "same channels in the same order")
  Y <- X
  colnames(Y) <- c("a", "a", "b", "c")
  refuses(frm_cross_pairs(Y), "unique and non-empty")
})

test_that("four pairs fit with (1 | pair) and match their own fits", {
  skip_on_cran()
  set.seed(45)
  X <- xp_rec(8192L)
  pr <- rbind(c("Fz", "Cz"), c("Fz", "Pz"), c("Cz", "Pz"), c("Pz", "Oz"))
  xp <- xp_run(frm_cross_pairs(X, pairs = pr, segments = 16L))
  expect_s3_class(xp, "data.frame")
  form_sep <- frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
                         pow2 ~ 1, coh ~ 1, phase ~ 1)
  sep <- lapply(seq_len(nrow(pr)), function(k) {
    d <- frm_cross_spectrum(X[, pr[k, 1L]], X[, pr[k, 2L]], segments = 16L)
    fit <- frmtmb::frm(form_sep, family = cross_wishart(), data = d)
    frm_coherence(fit)[1L, ]
  })
  sep <- do.call(rbind, sep)
  lab <- paste(pr[, 1L], pr[, 2L], sep = "-")
  nd <- data.frame(pair = factor(lab, levels = lab))

  # every parameter by pair: the likelihood is a sum of one term per
  # pair, so each pair's coherence IS its own fit's, to the optimizer
  bp <- xp_fit(
    frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 0 + pair,
               pow2 ~ 0 + pair, coh ~ 0 + pair, phase ~ 0 + pair),
    family = cross_wishart(), data = xp)
  expect_identical(bp$other, character(0))
  by_pair <- bp$fit
  cb <- xp_run(frm_coherence(by_pair, newdata = nd))
  expect_s3_class(cb, "data.frame")
  # measured against each separate fit's own standard error, not an
  # absolute tolerance: 4.5e-04 of one at this seed, the optimizer's
  # tolerance and nothing else
  expect_lt(xp_dist(cb, sep), 1e-2)
  expect_lt(if (is.data.frame(cb)) max(abs(cb$.se / sep$.se - 1)) else NA,
            1e-3)

  # and the model the pair factor is for: one coherence per pair, drawn
  # toward the pairs' common level
  sh <- xp_fit(
    frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 0 + pair,
               pow2 ~ 0 + pair, coh ~ 1 + (1 | pair), phase ~ 0 + pair),
    family = cross_wishart(), data = xp)
  expect_identical(sh$other, character(0))
  shrunk <- sh$fit
  expect_s3_class(shrunk, "frmtmb_fit")
  expect_identical(if (inherits(shrunk, "frmtmb_fit"))
    shrunk$opt$convergence, 0L)
  cs <- xp_run(frm_coherence(shrunk, newdata = nd))
  expect_s3_class(cs, "data.frame")
  # the pairs differ by far more than their standard errors, so the
  # shrinkage is small: each pair stays within two of its own fit's
  # standard errors (0.29 at this seed), and the ORDER of the pairs is
  # kept
  expect_lt(xp_dist(cs, sep), 2)
  expect_identical(if (is.data.frame(cs)) order(cs$.estimate),
                   order(sep$.estimate))
})

# quadrature = TRUE when a row reads more than one random-effect block.
# One Gauss-Kronrod pass over every effect at once integrated the coarse
# block of (1 | ga/gb) and left the fine one to the Laplace
# approximation: the "quadrature" log-likelihood equaled the Laplace one
# to 1e-8, and summary() failed with "no 'dimnames' attribute for array"
# because sdreport() cannot mix the two methods. The gated fuzz plan
# found it on cumulative + weights (seed 20322452); the family and the
# weights play no part. Such fits are refused now, as lme4 refuses
# nAGQ > 1 for more than one term.

quad_nest_data <- function(seed = 20261020) {
  set.seed(seed)
  ng <- 12; nt <- 6; n <- ng * nt
  d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                  gb = factor(rep(c("i", "ii"), length.out = n)),
                  tt = factor(rep(seq_len(nt), ng)),
                  x = rnorm(n), w = runif(n, 0.5, 2))
  ab <- interaction(d$ga, d$gb, drop = TRUE)
  d$eta <- 0.2 + 0.5 * d$x + rnorm(ng, 0, 0.7)[d$ga] +
    rnorm(nlevels(ab), 0, 1)[ab]
  d$y <- rbinom(n, 1, plogis(d$eta))
  d$yo <- factor(cut(d$eta + stats::rlogis(n), c(-Inf, -0.8, 0.8, Inf),
                     labels = FALSE), ordered = TRUE)
  d$yg <- d$eta + rnorm(n, 0, 0.5)
  # a by-level that is constant within ga, so a by-split block reads the
  # rows of its own levels only
  d$gby <- factor(c("a", "b", "c")[(as.integer(d$ga) %% 3) + 1])
  d
}

quad_refusal <- "integrates one scalar random effect per row"

test_that("quadrature refuses a row that reads two blocks", {
  d <- quad_nest_data()
  # the fuzz spec's shape: weights and an ordinal family
  expect_error(
    frm(bf(yo | weights(w) ~ 1 + x + (1 | ga/gb)) + cumulative(),
        data = d, quadrature = TRUE),
    quad_refusal)
  expect_error(frm(bf(y ~ 1 + x + (1 | ga/gb)) + bernoulli(), data = d,
                   quadrature = TRUE),
               quad_refusal)
  # crossed
  expect_error(frm(bf(y ~ 1 + x + (1 | ga) + (1 | tt)) + bernoulli(),
                   data = d, quadrature = TRUE),
               quad_refusal)
  # one grouping factor, two dpars
  err_s <- tryCatch(frm(bf(yg ~ 1 + x + (1 | ga), sigma ~ (1 | ga)) +
                          gaussian(), data = d, quadrature = TRUE),
                    error = conditionMessage)
  expect_match(err_s, quad_refusal)
  # a dpar's block is named outside its term, as the formula writes it
  expect_match(err_s, "`sigma: (1 | ga)`", fixed = TRUE)
  # the message names the blocks as written
  err <- tryCatch(frm(bf(y ~ 1 + x + (1 | ga) + (1 | tt)) + bernoulli(),
                      data = d, quadrature = TRUE),
                  error = conditionMessage)
  expect_match(err, "1 | ga", fixed = TRUE)
  expect_match(err, "1 | tt", fixed = TRUE)
  # the same models fit under Laplace, which the refusal names
  expect_match(err, "quadrature = FALSE", fixed = TRUE)
  fl <- allow_boundary(frm(bf(yo | weights(w) ~ 1 + x + (1 | ga/gb)) +
                             cumulative(), data = d))
  sm <- allow_boundary(summary(fl))
  expect_s3_class(sm, "summary.frmtmb_fit")
})

test_that("one block per row is still integrated", {
  d <- quad_nest_data()
  # a by-split gives three blocks over ga, but each row reads one
  fb <- allow_boundary(frm(bf(y ~ x + (1 | gr(ga, by = gby))) +
                             bernoulli(), data = d, quadrature = TRUE))
  expect_length(fb$frame$re_blocks, 3L)
  expect_length(fb$obj$env$random, 0L)
  # a single block: every effect integrated, and summary() works
  fq <- allow_boundary(frm(bf(y | weights(w) ~ x + (1 | ga)) +
                             bernoulli(), data = d, quadrature = TRUE))
  expect_length(fq$obj$env$random, 0L)
  sm <- allow_boundary(summary(fq))
  expect_s3_class(sm, "summary.frmtmb_fit")
})

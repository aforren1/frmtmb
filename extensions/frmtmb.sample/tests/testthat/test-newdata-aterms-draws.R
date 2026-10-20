# Addition terms on newdata, on posterior draws. posterior_predict()
# read only trunc(), rate() and thres(gr = ) from newdata and handed the
# simulator an empty list otherwise, so a binomial with newdata N = 1000
# drew 0/1 counts, while posterior_epred() on the same newdata gave
# about 560 (dev/ndaterm-findings.md in the frmtmb repository).
# posterior_epred() reads every term through frmtmb, so it is the
# reference. Every tolerance is a multiple of a Monte Carlo standard
# error that the run measures from its own draws.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

nda_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(7)
      n <- 200
      d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE),
                      s = stats::runif(n, 0.2, 1))
      d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
      d$g <- 1 + d$x + stats::rnorm(n, 0, d$s)
      fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = d)
      fg <- frm(bf(g | se(s) ~ x), family = gaussian(), data = d)
      draw <- function(fit) {
        suppressWarnings(suppressMessages(
          frm_sample(fit, chains = 1, iter = 400, warmup = 200,
                     refresh = 0, seed = 7)))
      }
      cache <<- list(d = d, db = draw(fb), dg = draw(fg),
                     nd = data.frame(x = c(0, 0), N = c(1000, 2),
                                     s = c(0.05, 5)))
    }
    cache
  }
})

test_that("posterior_predict(newdata) draws trials() from newdata", {
  cs <- nda_case()
  set.seed(21)
  p <- posterior_predict(cs$db, newdata = cs$nd)
  e <- posterior_epred(cs$db, newdata = cs$nd)
  # 0.19.0: maximum 1 in both columns
  expect_true(all(p[, 1] <= 1000) && all(p[, 2] <= 2))
  expect_gt(max(p[, 1]), 2)
  # each predictive draw is a count around its own draw's expectation,
  # so the mean difference has the spread of the counts alone
  se <- apply(p - e, 2, stats::sd) / sqrt(nrow(p))
  expect_true(all(abs(colMeans(p) - colMeans(e)) <= 6 * se),
              info = paste(format(colMeans(p)), format(colMeans(e))))
})

test_that("predict(newdata) on draws summarizes those draws", {
  cs <- nda_case()
  set.seed(22)
  ps <- predict(cs$db, newdata = cs$nd)
  e <- colMeans(posterior_epred(cs$db, newdata = cs$nd))
  expect_true(all(abs(ps[, "Estimate"] - e) <=
                    6 * ps[, "Est.Error"] / sqrt(ndraws(cs$db))))
  expect_true(all(ps[, "Q97.5"] <= c(1000, 2)))
})

test_that("posterior_predict(newdata) draws at newdata's se()", {
  cs <- nda_case()
  set.seed(23)
  p <- posterior_predict(cs$dg, newdata = cs$nd)
  e <- posterior_epred(cs$dg, newdata = cs$nd)
  # the known sd beside the expected value's own posterior spread
  want <- sqrt(cs$nd$s^2 + apply(e, 2, stats::var))
  got <- apply(p, 2, stats::sd)
  rse <- 1 / sqrt(2 * (nrow(p) - 1))
  # 0.19.0: about 0.03 in both columns, the posterior spread of mu alone
  expect_true(all(abs(got / want - 1) <= 6 * rse),
              info = paste(format(got), format(want)))
})

test_that("draws of a model without a term need nothing from newdata", {
  cs <- nda_case()
  # the absent case: posterior_predict() of the fitted rows still reads
  # the frame's trials(), and nothing new is raised
  set.seed(24)
  expect_no_condition(p <- posterior_predict(cs$db, ndraws = 20))
  expect_true(all(sweep(p, 2, cs$d$N, "<=")))
  # newdata without the trials() variable is refused, naming it, as
  # brms refuses it
  expect_error(posterior_predict(cs$db, newdata = data.frame(x = 0)),
               "trials(N) could not be evaluated on newdata", fixed = TRUE)
})

test_that("posterior_predict() of a mixture reads theta and trials right", {
  skip_sampler()
  # a binomial mixture with N of 2 and 50 (the reviewer's construction,
  # dev/ndaterm-rev-probe2.R): each class must be drawn at its own
  # rows' trials(), and on newdata at the mixing weights the fit has
  set.seed(7)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), N = sample(c(2, 50), n, TRUE))
  cls <- stats::rbinom(n, 1, 0.5)
  d$y <- stats::rbinom(n, d$N, ifelse(cls == 1,
                                      stats::plogis(-1 + 0.5 * d$x),
                                      stats::plogis(1 + 0.5 * d$x)))
  fit <- frm(bf(y | trials(N) ~ 1, mu1 ~ x, mu2 ~ x),
             family = mixture(binomial(), binomial()), data = d)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 400, warmup = 200, refresh = 0,
               seed = 7)))
  set.seed(25)
  p <- posterior_predict(ds)
  # 0.19.0 on frmtmb 0.71.0: a quarter of the draws above N
  expect_identical(sum(sweep(p, 2, d$N, ">")), 0L)
  nd <- data.frame(x = c(0, 0), N = c(1000, 2))
  p <- posterior_predict(ds, newdata = nd)
  e <- posterior_epred(ds, newdata = nd)
  se <- apply(p - e, 2, stats::sd) / sqrt(nrow(p))
  expect_true(all(abs(colMeans(p) - colMeans(e)) <= 6 * se),
              info = paste(format(colMeans(p)), format(colMeans(e))))
})

test_that("posterior_predict(newdata) of a mixture reads theta unreported", {
  skip_sampler()
  # theta reports as the softmax probability; drawn from that as if it
  # were the predictor, the second component's weight was 0.31 where
  # the fit has 0.19 (0.19.0; dev/ndaterm-r1-theta-draws.R)
  set.seed(7)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  k <- stats::rbinom(n, 1, 0.2)
  d$y <- ifelse(k == 1, 6, 0) + d$x + stats::rnorm(n, 0, 0.5)
  fit <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = d)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 400, warmup = 200, refresh = 0,
               seed = 7)))
  nd <- data.frame(x = c(-1, 0, 1))
  set.seed(26)
  p <- posterior_predict(ds, newdata = nd)
  # the components sit 6 apart with sd 0.5, so a draw above 3 is a draw
  # of the second one; its share is the second weight, 1 - theta1
  w2 <- 1 - mean(posterior_epred(ds, newdata = nd, dpar = "theta1"))
  share <- mean(p > 3)
  expect_lt(abs(share - w2), 6 * sqrt(w2 * (1 - w2) / length(p)))
})

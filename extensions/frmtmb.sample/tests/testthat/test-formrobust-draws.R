# The draws side of two brms behaviors core now has
# (dev/formrobust-findings.md, sections 3 and 5):
# * posterior_predict(newdata = ) under ar(cov = FALSE) with the
#   response NA in some rows: brms fills each missing response with a
#   draw at that draw's parameters and runs the recursion over it.
#   Before, frmtmb.sample refused every newdata prediction of a
#   cov = FALSE model as a structured draw.
# * bernoulli() on a two-valued response that is not 0/1, which brms
#   codes by level order. The draws must agree with the 0/1 model's,
#   and a response read from newdata must be coded as the fit coded it.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

# the sampler's diagnostics on these short chains, and nothing else
fd_ess <- c("Effective Samples Size", "R-hat", "Rhat", "lp__")

fd_arma <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(31)
      G <- 30
      Tn <- 8
      d <- expand.grid(t = 1:Tn, g = factor(1:G))
      d$x <- stats::rnorm(nrow(d))
      e <- as.vector(apply(matrix(stats::rnorm(G * Tn), Tn, G), 2,
                           function(z) {
                             as.vector(stats::filter(z, 0.6, "recursive"))
                           }))
      d$y <- 1 + 0.5 * d$x + e
      ds <- NULL
      allow_warnings(ds <- suppressMessages(
        frm_sample(bf(y ~ x + ar(t, g, p = 1)), family = gaussian(),
                   data = d, chains = 1, iter = 600, refresh = 0,
                   seed = 3)), fd_ess)
      cache <<- list(d = d, ds = ds)
    }
    cache
  }
})

test_that("posterior_predict(newdata = ) fills a missing response", {
  cs <- fd_arma()
  nd <- cs$d[cs$d$g %in% c("1", "2"), ]
  nd$y[nd$g == "1" & nd$t >= 5] <- NA
  set.seed(2)
  pp <- posterior_predict(cs$ds, newdata = nd)
  expect_identical(dim(pp), c(300L, nrow(nd)))
  expect_false(anyNA(pp))
  # the rows before the first NA are drawn around the one-step mean,
  # which reads the observed residuals, as with every response observed
  full <- cs$d[cs$d$g %in% c("1", "2"), ]
  set.seed(2)
  pf <- posterior_predict(cs$ds, newdata = full)
  expect_false(anyNA(pf))
  ep <- posterior_epred(cs$ds, newdata = nd)
  expect_false(anyNA(ep))
  # a cov = FALSE term is drawn row by row, so re_formula is not refused
  # as it is for a structured draw
  set.seed(2)
  expect_false(anyNA(posterior_predict(cs$ds, re_formula = NA)))
  i4 <- which(nd$g == "1" & nd$t == 4)
  expect_identical(ep[, seq_len(i4)], posterior_epred(cs$ds,
                                                      newdata = full)[
                                                        , seq_len(i4)])
  # brms fills a missing response with a draw per posterior draw, so the
  # epred draws of the row after it carry the fill's spread, ar * sigma
  # (about 0.55 here), where a row with an observed past carries only
  # the parameters' (about 0.08): brms gives 0.58 against the
  # expected-value fill's 0.10 on these draws
  # (dev/formrobust-rev-log/sample-arma.txt), so the ratio is about 7
  # with the fill and about 1.2 without it
  sd5 <- stats::sd(ep[, i4 + 1L])
  sd6 <- stats::sd(ep[, i4 + 2L])
  expect_gt(sd6 / sd5, 3)
})

fd_bern <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(5)
      n <- 120
      d <- data.frame(x = stats::rnorm(n))
      yy <- stats::rbinom(n, 1, stats::plogis(-0.3 + 1.2 * d$x))
      d$y01 <- yy
      d$ym <- ifelse(yy == 1, -1, -2)
      s <- function(f) {
        out <- NULL
        allow_warnings(out <- suppressMessages(
          frm_sample(bf(f), family = bernoulli(), data = d, chains = 1,
                     iter = 300, refresh = 0, seed = 3)), fd_ess)
        out
      }
      cache <<- list(d = d, dm = s(ym ~ x), d0 = s(y01 ~ x))
    }
    cache
  }
})

test_that("bernoulli draws on -1/-2 are the 0/1 draws", {
  cs <- fd_bern()
  # the same codes reach the same tape, so the sampler takes the same
  # path: an identity
  expect_identical(as.matrix(cs$dm$draws), as.matrix(cs$d0$draws))
  expect_identical(log_lik(cs$dm), log_lik(cs$d0))
  set.seed(1)
  a <- posterior_predict(cs$dm)
  set.seed(1)
  b <- posterior_predict(cs$d0)
  expect_identical(a, b)
  # newdata holding only -2, which coded afresh would read as a 1
  nd <- cs$d[cs$d$ym == -2, ][1:6, ]
  set.seed(1)
  a <- predictive_error(cs$dm, newdata = nd)
  set.seed(1)
  b <- predictive_error(cs$d0, newdata = nd)
  expect_identical(a, b)
})

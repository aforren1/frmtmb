# Lane wt-predfix: on draws, an unseen level WITHOUT allow_new_levels
# reached core's refusal, whose remedy "Use allow_new_levels = TRUE" is
# the argument the draws methods then refuse. Seen to fail on the
# released 0.10.0 (dev/predfix-log/). Records: dev/predfix-findings.md.

pf_draws <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(20260922)
      n <- 60
      dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:6, each = 10)))
      dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x +
                             stats::rnorm(6, 0, 0.7)[dd$g], 1)
      fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
      cache <<- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 200, warmup = 100,
                   seed = 20260922, refresh = 0)))
    }
    cache
  }
})

pf_fns <- c("posterior_predict", "posterior_epred", "posterior_linpred",
            "predict", "fitted", "residuals")

pf_msg <- function(fn, ds, ...) {
  tryCatch({
    set.seed(1)
    do.call(fn, list(ds, ..., ndraws = 5))
    "answered"
  }, error = function(e) conditionMessage(e))
}

test_that("an unseen level without the flag gets no refused hint", {
  skip_on_cran()
  ds <- pf_draws()
  unseen <- data.frame(x = 0, g = factor("new"), y = 1)
  for (fn in pf_fns) {
    m <- pf_msg(fn, ds, newdata = unseen)
    # 0.10.0: "New levels in grouping factor `g`: new. Use
    # allow_new_levels = TRUE to predict them at the population level"
    expect_no_match(m, "Use allow_new_levels = TRUE", fixed = TRUE)
    expect_match(m, paste0("^", fn, "[(][)] on draws"), info = fn)
    expect_match(m, "refused here as well", fixed = TRUE, info = fn)
    expect_match(m, "`g`: new", fixed = TRUE, info = fn)
  }
})

test_that("sample_new_levels alone is the call without the flag", {
  skip_on_cran()
  ds <- pf_draws()
  unseen <- data.frame(x = 0, g = factor("new"), y = 1)
  for (fn in pf_fns) {
    m <- pf_msg(fn, ds, newdata = unseen, sample_new_levels = "gaussian")
    expect_no_match(m, "Use allow_new_levels = TRUE", fixed = TRUE)
    expect_match(m, "refused here as well", fixed = TRUE, info = fn)
  }
})

test_that("a newdata without the grouping column gets no refused hint", {
  skip_on_cran()
  ds <- pf_draws()
  nog <- data.frame(x = 0, y = 1)
  for (fn in pf_fns) {
    m <- pf_msg(fn, ds, newdata = nog)
    # 0.10.0: "... allow_new_levels = TRUE treats every row as an
    # unseen level ..."
    expect_no_match(m, "allow_new_levels = TRUE treats", fixed = TRUE)
    expect_match(m, "refused here as well", fixed = TRUE, info = fn)
  }
})

test_that("re_formula = NA answers an unseen level, flag or not", {
  skip_on_cran()
  ds <- pf_draws()
  unseen <- data.frame(x = 0, g = factor("new"), y = 1)
  # brms 2.23.0 answers both (dev/predfix-log/brms.txt); 0.10.0 refused
  # the flagged one because it asked about the level without re_formula
  for (fn in pf_fns) {
    expect_identical(pf_msg(fn, ds, newdata = unseen, re_formula = NA),
                     "answered", info = fn)
    expect_identical(pf_msg(fn, ds, newdata = unseen, re_formula = NA,
                            allow_new_levels = TRUE), "answered", info = fn)
  }
})

test_that("posterior_predict() on draws carries cs()", {
  skip_on_cran()
  skip_sampler()
  set.seed(11)
  n <- 300
  x <- stats::rnorm(n)
  p1 <- stats::plogis(-0.3 + 1.5 * x)
  p2 <- (1 - p1) * stats::plogis(0.8 - 1.2 * x)
  u <- stats::runif(n)
  d <- data.frame(x = x, yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
  fit <- frm(bf(yo ~ cs(x)), family = sratio(), data = d)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, warmup = 150, seed = 11,
               refresh = 0)))
  # at x = 3 the exact probability of category 1 is about 0.95; 0.10.0
  # drew it at about 0.42, the x = 0 value, because the cs() offsets
  # never reached the simulator
  # in sample, only the rows at large x: pooled over every row the
  # absent term and the present one average to nearly the same share
  hi <- which(d$x > 1.5)
  for (nd in list(data.frame(x = rep(3, 40)), NULL)) {
    set.seed(2)
    pp <- posterior_predict(ds, newdata = nd)
    ep <- posterior_epred(ds, newdata = nd)
    cols <- if (is.null(nd)) hi else seq_len(ncol(pp))
    obs <- mean(pp[, cols] == 1)
    ex <- mean(ep[, cols, 1L])
    expect_lt(abs(obs - ex) / sqrt(ex * (1 - ex) / length(pp[, cols])), 6,
              label = if (is.null(nd)) "in sample, x > 1.5" else "at x = 3")
  }
})

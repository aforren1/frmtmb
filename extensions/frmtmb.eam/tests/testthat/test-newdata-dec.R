# dec() on newdata. frmtmb 0.71.0's predict(newdata =) handed the
# simulator no addition terms but trunc(), rate() and thres(gr = ), so a
# wiener fit refused with "needs the boundary each row ended at, which
# is missing" although newdata carried it; fitted() and simulate() read
# it. conditional_effects() read dec() without this package's coercion,
# so a factor's levels reached the simulator swapped, and it held a
# numeric dec() at its mean, which every consumer read as the lower
# boundary (dev/ndaterm-findings.md in the frmtmb repository).
#
# bias = 0.3, so the two boundaries have different response-time
# distributions and a swapped or ignored boundary shows. Every tolerance
# is a multiple of a Monte Carlo standard error measured from the draws.

nd_dec <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(7)
      d <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.3)
      d$x <- stats::rnorm(300)
      d$side <- factor(ifelse(d$upper == 1, "up", "down"),
                       levels = c("up", "down"))
      cache <<- list(
        d = d,
        fit = frm(bf(rt | dec(upper) ~ x, bias = 0.3), family = wiener(),
                  data = d),
        fitf = frm(bf(rt | dec(side) ~ x, bias = 0.3), family = wiener(),
                   data = d),
        nd = data.frame(x = 0, upper = c(1, 0)))
    }
    cache
  }
})

test_that("predict(newdata, summary = FALSE) draws at newdata's dec()", {
  skip_if_not_installed("RWiener")
  cs <- nd_dec()
  set.seed(31)
  # 0.71.0: "wiener: simulating a response time needs the boundary each
  # row ended at, which is missing"
  p <- predict(cs$fit, newdata = cs$nd, summary = FALSE, ndraws = 1000)
  s <- as.matrix(simulate(cs$fit, nsim = 1000, newdata = cs$nd, seed = 31))
  se <- sqrt(apply(p, 2, stats::var) / nrow(p) +
               apply(s, 1, stats::var) / ncol(s))
  expect_true(all(abs(colMeans(p) - rowMeans(s)) <= 6 * se),
              info = paste(format(colMeans(p)), format(rowMeans(s))))
  # the boundaries differ by many standard errors, so a draw that
  # ignored dec() would not pass the line above
  expect_gt(abs(diff(rowMeans(s))), 10 * max(se))
})

test_that("posterior_predict(newdata) on draws reads dec() too", {
  skip_on_cran()
  skip_if_not_installed("RWiener")
  skip_if_not_installed("frmtmb.sample")
  skip_if_not_installed("tmbstan")
  skip_if(frmtmb.sample:::tmbstan_build_broken(),
          "this tmbstan samples a standard normal instead of the model")
  cs <- nd_dec()
  ds <- suppressMessages(allow_warnings(
    frmtmb.sample::frm_sample(cs$fit, chains = 1, iter = 400, warmup = 200,
                              seed = 3, refresh = 0),
    c("Effective Samples Size", "The largest R-hat is")))
  set.seed(32)
  p <- frmtmb.sample::posterior_predict(ds, newdata = cs$nd)
  e <- frmtmb.sample::posterior_epred(ds, newdata = cs$nd)
  # each draw's response time scatters around its own draw's
  # conditional mean, so the difference has the scatter alone
  se <- apply(p - e, 2, stats::sd) / sqrt(nrow(p))
  expect_true(all(abs(colMeans(p) - colMeans(e)) <= 6 * se),
              info = paste(format(colMeans(p)), format(colMeans(e))))
})

# Where `m` sits in the draws `ref`: the share of them below it, which
# is 1/2 at the median of the same distribution
share_below <- function(m, ref) mean(ref < m)

test_that("conditional_effects() codes a dec() factor as the fit did", {
  skip_if_not_installed("RWiener")
  cs <- nd_dec()
  for (lv in c("up", "down")) {
    set.seed(33)
    ce <- conditional_effects(
      cs$fitf, "x", method = "predict", ndraws = 1001,
      conditions = data.frame(side = factor(lv, c("up", "down"))))[[1]]
    i <- which.min(abs(ce$x))
    nd <- data.frame(x = ce$x[i], side = factor(lv, c("up", "down")))
    ref <- as.numeric(as.matrix(simulate(cs$fitf, nsim = 4000, newdata = nd,
                                         seed = 33)))
    # the estimate is the median of 1001 draws, so the share of the
    # reference below it is 1/2 within the two samples' binomial error;
    # 0.71.0 drew each level at the other boundary (shares off 1/2 by
    # 0.30 and 0.38, dev/ndaterm-out/t-eam-dec-r9.txt)
    tol <- 6 * sqrt(0.25 / 1001 + 0.25 / length(ref))
    expect_lt(abs(share_below(ce$estimate__[i], ref) - 0.5), tol)
  }
})

test_that("an explicit dec() value that is not a boundary is refused", {
  skip_if_not_installed("RWiener")
  cs <- nd_dec()
  expect_error(fitted(cs$fit, newdata = data.frame(x = 0, upper = 0.5)),
               "boundary indicator on newdata must be 0 (lower) or 1 (upper)",
               fixed = TRUE)
  # the absent case: a boundary of 0 or 1 is drawn, and an NA row stays
  # a missing row
  f <- fitted(cs$fit, newdata = data.frame(x = 0, upper = c(1, NA)))
  expect_true(is.finite(f[1, "Estimate"]))
  expect_no_condition(predict(cs$fit, summary = FALSE, ndraws = 3))
})

test_that("conditional_effects() holds a numeric dec() at a boundary", {
  skip_if_not_installed("RWiener")
  cs <- nd_dec()
  held <- "holding the addition-term variable(s) upper at 1"
  # the grid used to hold upper at its mean, 0.573, which every wiener
  # consumer read as the lower boundary; brms holds the variables of
  # its whole-valued terms at round(median()), 1 here
  expect_identical(round(stats::median(cs$d$upper)), 1)
  expect_message(ce <- conditional_effects(cs$fit, "x"), held,
                 fixed = TRUE)
  ref <- conditional_effects(cs$fit, "x", conditions = data.frame(upper = 1))
  expect_equal(ce[[1]]$estimate__, ref[[1]]$estimate__)
  expect_message(ce0 <- conditional_effects(cs$fit), held, fixed = TRUE)
  expect_true(all(is.finite(ce0[[1]]$estimate__)))
  for (eff in list(NULL, "x")) {
    set.seed(34)
    a <- c(list(cs$fit), eff, list(method = "predict", ndraws = 51))
    expect_message(cp <- do.call(conditional_effects, a), held,
                   fixed = TRUE)
    expect_true(all(is.finite(cp[[1]]$estimate__)))
  }
})

# fitted() on the factor fit at x = 0, for the given `side` column
dec_at <- function(cs, side) {
  unname(fitted(cs$fitf, newdata = data.frame(x = 0, side = side))[,
    "Estimate"])
}

dec_ref <- function(cs) {
  dec_at(cs, factor(c("down", "up"), c("up", "down")))
}

test_that("dec() on newdata is coded by the fit's levels, not its own", {
  skip_if_not_installed("RWiener")
  cs <- nd_dec()
  ref <- dec_ref(cs)
  # 0.71.0: relevel() and character swapped the boundaries in silence
  expect_equal(dec_at(cs, stats::relevel(factor(c("down", "up"),
                                                c("up", "down")), "down")),
               ref)
  expect_equal(dec_at(cs, c("down", "up")), ref)
})

test_that("one dec() value alone is coded, not refused for one level", {
  skip_if_not_installed("RWiener")
  cs <- nd_dec()
  ref <- dec_ref(cs)
  # 0.71.0: "could not be evaluated on newdata", from the coercion's
  # refusal of a one-level factor
  expect_equal(dec_at(cs, factor("down")), ref[1])
  expect_equal(dec_at(cs, "up"), ref[2])
})

test_that("a dec() value the fit did not see is refused by name", {
  skip_if_not_installed("RWiener")
  cs <- nd_dec()
  # 0.71.0 coded "left" and "right" as the two boundaries
  expect_error(dec_at(cs, factor(c("left", "right"))),
               "holds value(s) 'left', 'right' the fit did not see",
               fixed = TRUE)
})

test_that("a logical dec() is coded on one row", {
  skip_if_not_installed("RWiener")
  cs <- nd_dec()
  d <- cs$d
  d$lg <- d$upper == 1
  fl <- frm(bf(rt | dec(lg) ~ x, bias = 0.3), family = wiener(), data = d)
  one <- fitted(fl, newdata = data.frame(x = 0, lg = TRUE))[, "Estimate"]
  num <- fitted(cs$fit, newdata = data.frame(x = 0, upper = 1))[,
    "Estimate"]
  # 0.71.0: a one-row logical newdata was refused (one level)
  expect_equal(unname(one), unname(num))
})

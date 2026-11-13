# drop1() deletes every term of the main formula but the random-effect
# and autocorrelation terms. terms(fit) is the fixed design's terms,
# which the parser builds after mo(), me(), cs(), gp() and the smooths
# are taken out, so through 0.81.0 drop1(y ~ mo(k) + z) listed z alone
# and said nothing about mo(k). A me() term's measurements are part of
# the likelihood, so a row whose refit loses one has no statistics, as
# anova() refuses such a pair.
#
# Each row is checked against the refit written out by hand, as a ratio
# to the fit's own log-likelihood.

sp_data <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(11)
      n <- 200
      d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                      k = sample(1:4, n, TRUE),
                      g = factor(sample(letters[1:8], n, TRUE)),
                      t = rep(1:20, 10), id = factor(rep(1:10, each = 20)))
      d$xo <- d$x + stats::rnorm(n, 0, 0.3)
      d$sdx <- 0.3
      d$y <- 1 + 0.6 * d$k + 0.4 * d$z + sin(2 * d$x) +
        stats::rnorm(8)[d$g] + stats::rnorm(n, 0, 0.5)
      cache <<- d
    }
    cache
  }
})

ll <- function(f) as.numeric(logLik(f))
df_of <- function(f) attr(logLik(f), "df")
# the row of term `lab` against the hand refit `small`
expect_row <- function(tab, lab, full, small) {
  lrt <- 2 * (ll(full) - ll(small))
  expect_lt(abs(tab[lab, "LRT"] - lrt), 1e-8 * abs(ll(full)))
  expect_equal(tab[lab, "Df"], as.numeric(df_of(full) - df_of(small)))
}

test_that("drop1() drops a mo() term, which it used to leave out", {
  d <- sp_data()
  fit <- frm(y ~ mo(k) + z, data = d)
  tab <- drop1(fit, test = "Chisq")
  expect_identical(rownames(tab), c("<none>", "mo(k)", "z"))
  expect_row(tab, "mo(k)", fit, frm(y ~ z, data = d))
  expect_row(tab, "z", fit, frm(y ~ mo(k), data = d))
  # marginality protects a mo() main effect as it protects any other
  fi <- frm(y ~ mo(k) * z, data = d)
  expect_identical(rownames(drop1(fi)), c("<none>", "mo(k):z"))
})

test_that("drop1() drops a smooth with its unpenalized part", {
  d <- sp_data()
  fit <- frm(y ~ s(x) + z, data = d)
  tab <- drop1(fit, test = "Chisq")
  expect_identical(rownames(tab), c("<none>", "s(x)", "z"))
  expect_row(tab, "s(x)", fit, frm(y ~ z, data = d))
})

test_that("random-effect and autocorrelation terms stay out of scope", {
  d <- sp_data()
  fit <- frm(y ~ x + z + (1 | g), data = d)
  expect_identical(rownames(drop1(fit)), c("<none>", "x", "z"))
  # and add1() does not offer the fit's own random-effect term
  f1 <- frm(y ~ x + (1 | g), data = d)
  expect_identical(rownames(add1(f1, scope = ~ . + z + (1 | g))),
                   c("<none>", "z"))
  fa <- frm(y ~ x + ar(time = t, gr = id), data = d)
  expect_identical(rownames(drop1(fa)), c("<none>", "x"))
  # drop1() does not take a bar term at all
  expect_error(drop1(fit, scope = "1 | g"), "not a subset of the term")
})

test_that("add1() of a new random-effect term refits that term alone", {
  d <- sp_data()
  f0 <- frm(y ~ x, data = d)
  # through 0.81.0 the refit was update(f0, ~ . + 1 | g), which parses
  # as y ~ (x + 1 | g): a random slope in place of the fixed one
  hand <- frm(y ~ x + (1 | g), data = d)
  lrt <- 2 * (ll(hand) - ll(f0))
  for (sc in list(~ . + (1 | g), "1 | g")) {
    tab <- add1(f0, scope = sc, test = "Chisq")
    expect_identical(rownames(tab), c("<none>", "(1 | g)"))
    expect_lt(abs(tab["(1 | g)", "LRT"] - lrt), 1e-8 * abs(ll(hand)))
    expect_equal(tab["(1 | g)", "Df"], as.numeric(df_of(hand) - df_of(f0)))
  }
})

test_that("a refit that loses a me() call has no statistics, and says so", {
  d <- sp_data()
  fit <- frm(y ~ me(xo, sdx) + z, data = d)
  expect_message(tab <- drop1(fit, test = "Chisq"),
                 "no statistics for `me(xo, sdx)`", fixed = TRUE)
  expect_identical(rownames(tab), c("<none>", "me(xo, sdx)", "z"))
  expect_true(all(is.na(unlist(tab["me(xo, sdx)", ]))))
  # dropping z keeps the measurements, so that row is an ordinary test
  expect_row(tab, "z", fit, frm(y ~ me(xo, sdx), data = d))
})

test_that("add1() gives no statistics for a row that adds a me() call", {
  d <- sp_data()
  f0 <- frm(y ~ z, data = d)
  # through 0.81.0 the row read LRT -543.37: the likelihood of the
  # larger model holds the measurements of xo, and the smaller one not
  expect_message(tab <- add1(f0, scope = ~ . + me(xo, sdx) + mo(k),
                             test = "Chisq"),
                 "add1(): no statistics for `me(xo, sdx)`", fixed = TRUE)
  expect_true(all(is.na(unlist(tab["me(xo, sdx)", ]))))
  big <- frm(y ~ z + mo(k), data = d)
  lrt <- 2 * (ll(big) - ll(f0))
  expect_lt(abs(tab["mo(k)", "LRT"] - lrt), 1e-8 * abs(ll(big)))
  # the absent case: a me() call the fit already has adds no data
  fm <- frm(y ~ me(xo, sdx) + z, data = d)
  expect_no_message(tm <- add1(fm, scope = ~ . + me(xo, sdx):z,
                               test = "Chisq"))
  # through 0.81.0 add1() read the fixed design's terms, so it offered
  # me(xo, sdx) again with Df 0 and never the interaction
  expect_identical(rownames(tm), c("<none>", "me(xo, sdx):z"))
  big <- frm(y ~ me(xo, sdx) + z + me(xo, sdx):z, data = d)
  expect_lt(abs(tm["me(xo, sdx):z", "LRT"] - 2 * (ll(big) - ll(fm))),
            1e-8 * abs(ll(big)))
  fk <- frm(y ~ mo(k) + z, data = d)
  expect_identical(rownames(add1(fk, scope = ~ . + mo(k):z)),
                   c("<none>", "mo(k):z"))
})

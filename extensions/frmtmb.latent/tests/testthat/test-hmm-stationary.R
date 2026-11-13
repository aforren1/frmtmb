## A state's copy of a component dpar keeps the component's stationary
## point. skew_normal() declares alpha = 0, where the information is
## singular: a fit that stops there converges cleanly and reports a
## logLik tens of units below the maximum. hmm() starts every state's
## alpha from ONE pooled value, so a state skewed against the pooled
## residual starts on the wrong side of the point. Through 0.81.0 hmm()
## dropped the declaration and nothing sent the fit off it.

# two states skewed in opposite directions, so one of them is always
# skewed against the pooled residual
sn_hmm_data <- function(seed = 1L, nseq = 20L, Tl = 40L, alpha = 5) {
  set.seed(seed)
  rsn <- function(n, mu, sigma, a) {
    d <- a / sqrt(1 + a^2)
    mu + sigma * (d * abs(stats::rnorm(n)) +
                    sqrt(1 - d^2) * stats::rnorm(n))
  }
  G <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(nseq), function(id) {
    s <- integer(Tl)
    s[1L] <- sample.int(2L, 1L)
    for (t in seq_len(Tl)[-1L]) {
      s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
    }
    y <- ifelse(s == 1L, rsn(Tl, -3, 1, alpha), rsn(Tl, 3, 1, -alpha))
    data.frame(id = id, t = seq_len(Tl), y = y)
  }))
}

sn_hmm_fam <- function() hmm(K = 2, skew_normal(), time = t, group = id)

test_that("hmm() carries the component's stationary point to every state", {
  st <- sn_hmm_fam()$post$stationary
  expect_identical(names(st), c("alpha1", "alpha2"))
  own <- skew_normal()$post$stationary$alpha
  expect_identical(st$alpha1, own)
  expect_identical(st$alpha2, own)
  # a component that declares none gives a composition that declares none
  expect_null(hmm(K = 2, gaussian(), time = t, group = id)$post$stationary)
})

test_that("a state stopped at alpha = 0 is refit off it", {
  skip_on_cran()
  dd <- sn_hmm_data()
  fo <- bf(y ~ 1)
  u <- frm(fo, family = sn_hmm_fam(), data = dd, dry_run = "objective")
  bd <- u$estimates$betad
  a1 <- "alpha1_(Intercept)"
  a2 <- "alpha2_(Intercept)"
  # state 1 starts ON the stationary point; state 2 keeps its own start
  stuck <- bd
  stuck[[a1]] <- 0
  fit <- frm(fo, family = sn_hmm_fam(), data = dd,
             start = list(betad = stuck))
  esc <- fit$opt[["stationary_escape"]]
  expect_false(is.null(esc))
  # one stalled state, two sides
  expect_identical(unname(esc[["starts"]]), 2)
  ll <- as.numeric(logLik(fit))
  # the gain is far larger than the optimizer's own tolerance
  expect_gt(esc[["gain"]], abs(ll) * 1e-3)
  expect_gt(abs(fit$estimates$betad[[a1]]), 1)
  # the reference is this run's best fit from either side of the point
  ref <- max(vapply(c(2, -2), function(a) {
    s <- bd
    s[[a1]] <- a
    as.numeric(logLik(frm(fo, family = sn_hmm_fam(), data = dd,
                          start = list(betad = s))))
  }, 0))
  expect_lt(ref - ll, abs(ref) * 1e-8)
})

test_that("an hmm fit away from the point pays no refit", {
  skip_on_cran()
  dd <- sn_hmm_data()
  fo <- bf(y ~ 1)
  u <- frm(fo, family = sn_hmm_fam(), data = dd, dry_run = "objective")
  bd <- u$estimates$betad
  # each state starts on its own side, so neither stops near zero
  bd[["alpha1_(Intercept)"]] <- 2
  bd[["alpha2_(Intercept)"]] <- -2
  fit <- frm(fo, family = sn_hmm_fam(), data = dd,
             start = list(betad = bd))
  expect_null(fit$opt[["stationary_escape"]])
})

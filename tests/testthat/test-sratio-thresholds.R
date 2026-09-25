# sratio()'s thresholds are an unconstrained vector, as brms 2.23.0
# declares them (brms:::has_ordered_thres(sratio()) is FALSE). frmtmb
# used to hold them as (first threshold, log increments), like
# cumulative(), so where the optimum has crossing thresholds the fit sat
# on the ordering boundary instead of at brms's mode.
# dev/sratio-findings.md has the numbers; dev/sratio-brms-stan.R checks
# frmtmb's estimate against brms's compiled log density.

# Without predictors the stopping-ratio likelihood factorizes over the
# thresholds: threshold j sees n_j rows stop and m_j rows go on, so its
# maximum likelihood estimate is the link of the observed hazard
# n_j / (n_j + m_j), in closed form. The counts put the second hazard
# below the first, so the thresholds cross.
sr_counts <- c(50, 10, 40, 25)
sr_data <- function(counts = sr_counts) {
  data.frame(y = rep(seq_along(counts), counts))
}
sr_hazard <- function(counts) {
  counts[-length(counts)] / rev(cumsum(rev(counts)))[-length(counts)]
}
# the log-likelihood of threshold values `th` under distribution `Fd`
sr_loglik <- function(counts, th, Fd) {
  K1 <- length(th)
  go <- rev(cumsum(rev(counts)))[-1L]
  sum(counts[seq_len(K1)] * log(Fd(th)) + go * log(1 - Fd(th)))
}

# rows drawn from a stopping-ratio model with thresholds `th`: each row
# stops at category j with probability plogis(th_j - 0.8 x) once it has
# reached j
sr_sim <- function(th, seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  y <- rep(length(th) + 1L, n)
  for (j in rev(seq_along(th))) {
    y[stats::runif(n) < stats::plogis(th[j] - 0.8 * d$x)] <- j
  }
  d$y <- y
  d
}

test_that("crossing thresholds reach the closed-form maximum", {
  d <- sr_data()
  h <- sr_hazard(sr_counts)
  links <- list(logit = c(stats::qlogis, stats::plogis),
                probit = c(stats::qnorm, stats::pnorm),
                cloglog = c(function(p) log(-log(1 - p)),
                            function(x) 1 - exp(-exp(x))))
  for (lk in names(links)) {
    fit <- frm(y ~ 1, data = d, family = sratio(lk))
    fe <- fixef(fit)
    ref <- links[[lk]][[1L]](h)
    # the second threshold lies BELOW the first; an ordered storage
    # cannot reach it and stops where the two are equal
    expect_lt(fe["Intercept[2]", "Estimate"], fe["Intercept[1]", "Estimate"])
    expect_lt(max(abs(fe[, "Estimate"] - ref) / fe[, "Est.Error"]), 1e-3,
              label = lk)
    ll <- sr_loglik(sr_counts, ref, links[[lk]][[2L]])
    expect_lt(abs(as.numeric(logLik(fit)) - ll) / abs(ll),
              sqrt(.Machine$double.eps), label = lk)
    # the model is saturated, so every row's category probabilities are
    # the observed proportions
    P <- fitted(fit)[1L, "Estimate", ]
    p <- sr_counts / sum(sr_counts)
    expect_lt(max(abs(P - p) / p), 1e-4, label = lk)
  }
})

test_that("the thresholds are stored as they are reported", {
  fit <- frm(y ~ 1, data = sr_data(), family = sratio())
  raw <- fit$estimates[["tau_raw"]]
  expect_identical(unname(fixef(fit)[, "Estimate"]), unname(raw))
  map <- fit$spec$responses[[1L]]$family[["post"]][["ord_thresholds"]]
  expect_identical(map(raw), raw)
  # hypothesis() reads brms's names off the same values
  hy <- hypothesis(fit, "Intercept[2] < Intercept[1]")
  expect_equal(hy$hypothesis$Estimate, raw[2L] - raw[1L])
})

test_that("sratio and cratio are one model under the logit", {
  # F(t - eta) = 1 - F(eta - t) for a symmetric F, so the two sequential
  # families give every row the same hazard at the same thresholds.
  # Both now hold them unconstrained, so their maxima coincide even where
  # the thresholds cross, which the ordered storage broke
  d <- sr_sim(c(0.5, -1), seed = 3)
  fs <- frm(y ~ x, data = d, family = sratio())
  fc <- frm(y ~ x, data = d, family = cratio())
  th <- fixef(fs)[c("Intercept[1]", "Intercept[2]"), "Estimate"]
  expect_lt(th[[2L]], th[[1L]])
  ll <- as.numeric(logLik(fc))
  expect_lt(abs(as.numeric(logLik(fs)) - ll) / abs(ll),
            sqrt(.Machine$double.eps))
  se <- fixef(fc)[, "Est.Error"]
  expect_lt(max(abs(fixef(fs)[, "Estimate"] - fixef(fc)[, "Estimate"]) / se),
            1e-3)
})

test_that("with cs() the fit is K - 1 binomial regressions", {
  # the stopping-ratio model asks, at each category in turn, whether a
  # row that reached it stops there. With a slope per threshold those
  # questions share no parameter, so the maximum is that of K - 1
  # logistic regressions on shrinking subsets (Laara and Matthews 1985),
  # whose intercepts are unconstrained, as brms's thresholds are
  d <- sr_sim(c(0.5, -1), seed = 5)
  fit <- frm(y ~ cs(x), data = d, family = sratio())
  glms <- lapply(1:2, function(k) {
    sub <- d[d$y >= k, ]
    stats::glm(as.integer(sub$y == k) ~ x, data = sub,
               family = stats::binomial())
  })
  ll <- sum(vapply(glms, function(g) as.numeric(stats::logLik(g)), 0))
  expect_lt(abs(as.numeric(logLik(fit)) - ll) / abs(ll),
            sqrt(.Machine$double.eps))
  fe <- fixef(fit)
  a <- vapply(glms, function(g) stats::coef(g)[[1L]], 0)
  se <- vapply(glms, function(g) sqrt(stats::vcov(g)[1L, 1L]), 0)
  nm <- c("Intercept[1]", "Intercept[2]")
  expect_lt(a[2L], a[1L])
  expect_lt(max(abs(fe[nm, "Estimate"] - a) / se), 1e-3)
})

test_that("grouped thresholds cross within a level", {
  # thres(gr = g): one closed-form maximum per level
  ca <- sr_counts
  cb <- c(20, 30, 25)
  d <- data.frame(y = c(rep(seq_along(ca), ca), rep(seq_along(cb), cb)),
                  g = rep(c("a", "b"), c(sum(ca), sum(cb))))
  fit <- frm(y | thres(gr = g) ~ 1, data = d, family = sratio())
  fe <- fixef(fit)
  ref <- stats::qlogis(c(sr_hazard(ca), sr_hazard(cb)))
  nm <- c(paste0("Intercept[a,", 1:3, "]"), paste0("Intercept[b,", 1:2, "]"))
  expect_lt(max(abs(fe[nm, "Estimate"] - ref) / fe[nm, "Est.Error"]), 1e-3)
  ll <- sr_loglik(ca, ref[1:3], stats::plogis) +
    sr_loglik(cb, ref[4:5], stats::plogis)
  expect_lt(abs(as.numeric(logLik(fit)) - ll) / abs(ll),
            sqrt(.Machine$double.eps))
})

test_that("a class Intercept prior carries no Jacobian, as in brms", {
  # brms puts the prior on its unconstrained threshold vector, so the
  # posterior mode maximizes the log-likelihood plus the log prior
  # density and nothing else. Without predictors each threshold is its
  # own one-dimensional problem, solved here by optimize()
  d <- sr_data()
  go <- rev(cumsum(rev(sr_counts)))[-1L]
  ref <- vapply(1:3, function(j) {
    stats::optimize(function(t) {
      sr_counts[j] * stats::plogis(t, log.p = TRUE) +
        go[j] * stats::plogis(t, lower.tail = FALSE, log.p = TRUE) +
        stats::dnorm(t, 0, 2, log = TRUE)
    }, c(-10, 10), maximum = TRUE, tol = 1e-10)$maximum
  }, 0)
  pr <- set_prior("normal(0, 2)", class = "Intercept")
  for (fam in c("sratio", "cratio")) {
    ff <- get(fam, envir = asNamespace("frmtmb"))
    fit <- frm(y ~ 1, data = d, family = ff(), prior = pr)
    fe <- fixef(fit)
    expect_lt(max(abs(fe[, "Estimate"] - ref) / fe[, "Est.Error"]), 1e-3,
              label = fam)
    # the penalized objective at the reference point is minus the
    # log-likelihood minus the log prior, with no Jacobian term
    obj <- fit$obj
    par <- obj$par
    par[names(par) == "tau_raw"] <- ref
    want <- sr_loglik(sr_counts, ref, stats::plogis) +
      sum(stats::dnorm(ref, 0, 2, log = TRUE))
    expect_lt(abs(-obj$fn(par) - want), 64 * .Machine$double.eps *
                abs(want), label = fam)
  }
})

test_that("a multivariate model reaches each sratio response's own maximum", {
  set.seed(8)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y1 <- 1 + 0.5 * d$x + rnorm(n)
  d$o <- rep(seq_along(sr_counts), sr_counts)[sample(sum(sr_counts), n,
                                                     TRUE)]
  mv <- frm(bf(o ~ x) + sratio() + bf(y1 ~ x) + gaussian(), data = d)
  u <- frm(bf(o ~ x) + sratio(), data = d)
  fe <- fixef(mv)
  fu <- fixef(u)
  nm <- paste0("Intercept[", 1:3, "]")
  expect_lt(max(abs(fe[paste0("o_", nm), "Estimate"] - fu[nm, "Estimate"]) /
                  fu[nm, "Est.Error"]), 1e-3)
  expect_lt(fe["o_Intercept[2]", "Estimate"], fe["o_Intercept[1]", "Estimate"])
})

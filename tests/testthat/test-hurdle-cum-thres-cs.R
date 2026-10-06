# hurdle_cumulative() with thres(gr = ) and with cs(), which brms 2.23.0
# fits and 0.67.0 refused. The reference is brms's R side,
# brms:::posterior_epred_hurdle_cumulative(): the hurdle probability for
# category 0 and (1 - hu) times brms:::dcumulative() at the row's own
# thresholds for 1..K, each group's slice under thres(gr = ) and the
# thresholds less the cs() offsets under cs(). The Stan-side identity
# is test-ordinal-mixture-brms.R (gated) and dev/ordmix-lpcheck.R.

hcg_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  # the groups' cut points differ, and group b never reaches the top
  u <- stats::rlogis(n) + 0.8 * d$x
  top <- ifelse(d$g == "a", 1.5, Inf)
  y <- 1L + (u > -1) + (u > ifelse(d$g == "a", 0.3, 0.8)) + (u > top)
  d$y <- ifelse(stats::runif(n) < stats::plogis(-1 + 0.5 * d$z), 0L, y)
  d
}

hcg_epred <- function(eta, hu, disc, tau, link) {
  n <- length(eta)
  if (!is.matrix(tau)) tau <- matrix(tau, n, length(tau), byrow = TRUE)
  P <- brms:::dcumulative(seq_len(ncol(tau) + 1L), eta = eta, thres = tau,
                          disc = disc, link = link)
  cbind(hu, (1 - hu) * P)
}

hcg_ll <- function(P, y) sum(log(P[cbind(seq_along(y), y + 1L)]))

hcg_par <- function(fit, move = 0) {
  p <- fit$opt$par
  p + move * sin(seq_along(p) * 1.3)
}

test_that("thres(gr = ) keeps the hurdle in every group's density", {
  skip_unless_brms()
  d <- hcg_data(20261020)
  for (lk in c("logit", "probit")) {
    fit <- frm(bf(y | thres(gr = g) ~ x, hu ~ z),
               family = hurdle_cumulative(lk), data = d)
    th <- family(fit)[["thres"]]
    expect_identical(th[["nthres"]], c(3L, 2L))
    for (move in c(0, 0.1)) {
      est <- fit$obj$env$parList(hcg_par(fit, move))
      raw <- est$tau_raw
      # group a: three ordered thresholds; group b: two
      ta <- c(raw[1], raw[1] + cumsum(exp(raw[2:3])))
      tb <- c(raw[4], raw[4] + exp(raw[5]))
      eta <- est$beta[["x"]] * d$x
      hu <- stats::plogis(est$betad[["hu_(Intercept)"]] +
                            est$betad[["hu_z"]] * d$z)
      P <- matrix(0, nrow(d), 5)
      ia <- d$g == "a"
      P[ia, ] <- hcg_epred(eta[ia], hu[ia], 1, ta, lk)
      P[!ia, 1:4] <- hcg_epred(eta[!ia], hu[!ia], 1, tb, lk)
      ref <- hcg_ll(P, d$y)
      expect_lt(abs(-fit$obj$fn(hcg_par(fit, move)) - ref),
                1e3 * .Machine$double.eps * abs(ref), label = lk)
      if (move == 0) {
        Pf <- fitted(fit)[, "Estimate", ]
        expect_identical(colnames(Pf), paste0("P(Y = ", 0:4, ")"))
        # a column past a group's own categories is 0, as in brms
        expect_lt(max(abs(unname(Pf) - P)), 1e3 * .Machine$double.eps,
                  label = lk)
        nd <- d[c(1:3, which(!ia)[1:3]), ]
        Pn <- unname(fitted(fit, newdata = nd)[, "Estimate", ])
        expect_lt(max(abs(Pn - P[c(1:3, which(!ia)[1:3]), ])),
                  1e3 * .Machine$double.eps, label = lk)
      }
    }
  }
  expect_true(all(c("b_Intercept[a,3]", "b_Intercept[b,2]", "b_hu_z") %in%
                    variables(fit)))
})

test_that("thres(gr = ) with disc and with each structure", {
  skip_unless_brms()
  d <- hcg_data(20261021)
  fit <- frm(bf(y | thres(gr = g) ~ x, disc ~ 0 + z),
             family = hurdle_cumulative("probit"), data = d)
  est <- fit$obj$env$parList(hcg_par(fit, 0.1))
  raw <- est$tau_raw
  ta <- c(raw[1], raw[1] + cumsum(exp(raw[2:3])))
  tb <- c(raw[4], raw[4] + exp(raw[5]))
  eta <- est$beta[["x"]] * d$x
  disc <- exp(est$betad[["disc_z"]] * d$z)
  hu <- stats::plogis(est$betad[["hu_(Intercept)"]])
  ia <- d$g == "a"
  P <- matrix(0, nrow(d), 5)
  P[ia, ] <- hcg_epred(eta[ia], hu, disc[ia], ta, "probit")
  P[!ia, 1:4] <- hcg_epred(eta[!ia], hu, disc[!ia], tb, "probit")
  ref <- hcg_ll(P, d$y)
  expect_lt(abs(-fit$obj$fn(hcg_par(fit, 0.1)) - ref),
            1e3 * .Machine$double.eps * abs(ref))
  # equidistant: a first threshold and a delta per level, delta_<k> as
  # brms numbers them
  fe <- frm(bf(y | thres(gr = g) ~ x),
            family = hurdle_cumulative(threshold = "equidistant"), data = d)
  expect_true(all(c("delta_1", "delta_2") %in% variables(fe)))
  te <- fixef(fe)[c("Intercept[a,1]", "Intercept[a,2]", "Intercept[a,3]"),
                  "Estimate"]
  expect_lt(abs(diff(te)[2] - diff(te)[1]),
            1e3 * .Machine$double.eps * max(abs(te)))
  # sum_to_zero: each level's thresholds sum to zero
  fz <- frm(bf(y | thres(gr = g) ~ x),
            family = hurdle_cumulative(threshold = "sum_to_zero"), data = d)
  tz <- fixef(fz)[c("Intercept[b,1]", "Intercept[b,2]"), "Estimate"]
  expect_lt(abs(sum(tz)), 1e3 * .Machine$double.eps * max(abs(tz)))
})

test_that("simulate() draws each group's own categories and the hurdle", {
  d <- hcg_data(20261022, n = 300)
  fit <- frm(bf(y | thres(gr = g) ~ x), family = hurdle_cumulative(),
             data = d)
  sims <- simulate(fit, nsim = 300, seed = 3)
  codes <- vapply(sims, as.integer, integer(nrow(d)))
  # group b has categories 0..3 only
  expect_true(all(codes[d$g == "b", ] <= 3L))
  P <- fitted(fit)[, "Estimate", ]
  freq <- vapply(0:4, function(k) mean(codes == k), 0)
  se <- sqrt(colSums(P * (1 - P))) / (nrow(d) * sqrt(300))
  expect_true(all(abs(freq - colMeans(P)) < 5 * se))
})

test_that("cs() comes off each row's thresholds, as brms reads it", {
  skip_unless_brms()
  set.seed(20261023)
  n <- 400
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  # offsets small against the gaps between the thresholds, so that no
  # row's thresholds cross at the fit (the crossing case is below)
  u <- stats::rlogis(n) + 0.5 * d$x
  y <- 1L + (u > -1 + 0.15 * d$x) + (u > 0.3) + (u > 1.5 - 0.15 * d$x)
  d$y <- ifelse(stats::runif(n) < 0.2, 0L, y)
  for (f in list(bf(y ~ cs(x)), bf(y ~ cs(x), disc ~ 0 + z, hu ~ z))) {
    fit <- frm(f, family = hurdle_cumulative("probit"), data = d)
    lp <- Filter(function(l) identical(l$dpar, "mu"),
                 fit$frame$linpreds)[[1L]]
    for (move in c(0, 0.05)) {
      est <- fit$obj$env$parList(hcg_par(fit, move))
      tau <- c(est$tau_raw[1], est$tau_raw[1] + cumsum(exp(est$tau_raw[-1])))
      cs <- outer(d$x, est[[lp$cs[[1L]]$par]])
      bd <- est$betad
      disc <- if ("disc_z" %in% names(bd)) exp(bd[["disc_z"]] * d$z) else 1
      hu <- stats::plogis(bd[["hu_(Intercept)"]] +
                            if ("hu_z" %in% names(bd)) bd[["hu_z"]] * d$z
                          else 0)
      P <- hcg_epred(rep(0, n), hu, disc,
                     matrix(tau, n, 3, byrow = TRUE) - cs, "probit")
      ref <- hcg_ll(P, d$y)
      expect_lt(abs(-fit$obj$fn(hcg_par(fit, move)) - ref),
                1e3 * .Machine$double.eps * abs(ref),
                label = deparse1(f$formula))
      if (move == 0) {
        expect_lt(max(abs(unname(fitted(fit)[, "Estimate", ]) - P)),
                  1e3 * .Machine$double.eps)
      }
    }
  }
  expect_true(all(c("bcs_x[1]", "bcs_x[3]") %in% variables(fit)))
  # the simulator reads the same offsets
  sims <- simulate(fit, nsim = 300, seed = 4)
  codes <- vapply(sims, as.integer, integer(n))
  Pf <- fitted(fit)[, "Estimate", ]
  freq <- vapply(0:4, function(k) mean(codes == k), 0)
  se <- sqrt(colSums(Pf * (1 - Pf))) / (n * sqrt(300))
  expect_true(all(abs(freq - colMeans(Pf)) < 5 * se))
})

test_that("crossing cs() thresholds give NaN, as brms's density does", {
  # brms's log_inv_logit_diff() is undefined where the offsets reverse
  # two thresholds, and so is a probability below zero
  fam <- hurdle_cumulative()
  raw <- c(-1, log(0.5), log(1))
  cs <- matrix(c(0, 0, 0, 0, 2, 0), 2, 3, byrow = TRUE)
  out <- fam$lpdf(c(2, 2), list(mu = c(0, 0), hu = c(0.2, 0.2),
                                disc = c(1, 1), .cs = cs),
                  list(), list(tau_raw = raw))
  expect_true(is.finite(out[1]))
  expect_true(is.nan(out[2]))
  # and a negative probability is no distribution to draw from
  set.seed(1)
  s <- fam$sim(list(mu = c(0, 0), hu = c(0, 0), disc = c(1, 1), .cs = cs),
               list(), 2L, list(tau_raw = raw))
  expect_true(!is.na(s[1]) && is.na(s[2]))
})

test_that("cs() and thres(gr = ) are refused together, as in brms", {
  d <- hcg_data(20261024)
  expect_error(frm(bf(y | thres(gr = g) ~ cs(x)),
                   family = hurdle_cumulative(), data = d),
               "Cannot use category specific effects", fixed = TRUE)
  # absent: each alone fits
  expect_true(is.finite(as.numeric(logLik(
    frm(bf(y | thres(gr = g) ~ x), family = hurdle_cumulative(), data = d)))))
  expect_true(is.finite(as.numeric(logLik(
    frm(bf(y ~ cs(x)), family = hurdle_cumulative(), data = d)))))
  # cumulative() keeps its refusal of cs()
  d$y1 <- pmax(d$y, 1L)
  expect_error(frm(bf(y1 ~ cs(x)), family = cumulative(), data = d),
               "cs() needs an sratio, cratio, or acat family", fixed = TRUE)
})

test_that("conditional_effects() reads each group's thresholds and the hurdle", {
  d <- hcg_data(20261025)
  fit <- frm(bf(y | thres(gr = g) ~ x), family = hurdle_cumulative(),
             data = d)
  for (lv in c("a", "b")) {
    ce <- conditional_effects(fit, effects = "x",
                              conditions = data.frame(g = lv))[[1L]]
    expect_identical(levels(ce$cats__), as.character(0:4))
    nd <- ce[ce$cats__ == "0", c("x", "g")]
    P <- fitted(fit, newdata = nd)[, "Estimate", ]
    for (k in 0:4) {
      expect_lt(max(abs(ce$estimate__[ce$cats__ == as.character(k)] -
                          P[, k + 1L])), 1e3 * .Machine$double.eps)
    }
    # group b has no category 4
    if (lv == "b") expect_true(all(ce$estimate__[ce$cats__ == "4"] == 0))
  }
})

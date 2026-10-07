# Ordinal mixtures: brms's mixture() of cumulative(), sratio(), cratio(),
# acat() and hurdle_cumulative() components. The reference is brms
# 2.23.0's own R-side densities, brms:::dcumulative() and its three
# siblings, which are what brms's posterior_epred(), log_lik() and
# posterior_predict() read for each component, combined the way
# brms:::log_lik_mixture() and brms:::posterior_epred_mixture() combine
# them: the theta-weighted sum. The thresholds are rebuilt by hand from
# the internal vector, so the threshold maps are checked as well. The
# Stan-side identity is test-ordinal-mixture-brms.R (gated) and
# dev/ordmix-lpcheck.R.

omx_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) +
    stats::rlogis(n) / exp(0.3 * d$z * cls)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < stats::plogis(-1.2 + 0.5 * d$z), 0L, d$y)
  d
}

# the thresholds from one internal vector, written out per structure
omx_tau <- function(r, type, ordered, k1) {
  switch(type,
    flexible = if (ordered) c(r[1], r[1] + cumsum(exp(r[-1]))) else r,
    equidistant = r[1] + (seq_len(k1) - 1) * (if (ordered) exp(r[2]) else
      r[2]),
    sum_to_zero = if (ordered) {
      u <- c(0, cumsum(exp(r)))
      u - mean(u)
    } else c(r, -sum(r)))
}

# brms's n x K category probabilities of one component; `tau` is one
# vector or an n x (K - 1) matrix of row thresholds
omx_brms_probs <- function(family, eta, disc, tau, link) {
  dens <- get(paste0("d", family), asNamespace("brms"))
  n <- length(eta)
  if (!is.matrix(tau)) tau <- matrix(tau, n, length(tau), byrow = TRUE)
  dens(seq_len(ncol(tau) + 1L), eta = eta, thres = tau, disc = disc,
       link = link)
}

omx_ll <- function(P, y, code0 = 1L) {
  sum(log(P[cbind(seq_along(y), y + 1L - code0)]))
}

# the parameter vector a test evaluates the density at: the estimates,
# or the estimates moved by a fixed pattern, so that the check is not
# made only where the gradient vanishes
omx_par <- function(fit, move = 0) {
  p <- fit$opt$par
  p + move * sin(seq_along(p) * 1.7)
}

test_that("an ordinal mixture's density is the theta-weighted brms sum", {
  skip_unless_brms()
  d <- omx_data(20261005)
  shapes <- list(
    list(f = bf(y ~ x, disc1 ~ 0 + z, theta1 ~ z),
         fam = mixture(cumulative(), sratio()),
         comp = c("cumulative", "sratio"), type = c("flexible", "flexible"),
         link = c("logit", "logit"), shared = FALSE),
    list(f = bf(y ~ x), fam = mixture(cumulative("probit"), acat(),
                                      order = "mu"),
         comp = c("cumulative", "acat"), type = c("flexible", "flexible"),
         link = c("probit", "logit"), shared = TRUE),
    list(f = bf(y ~ x),
         fam = mixture(cratio(), cumulative(threshold = "sum_to_zero"),
                       order = "mu"),
         comp = c("cratio", "cumulative"),
         type = c("flexible", "sum_to_zero"), link = c("logit", "logit"),
         shared = TRUE),
    list(f = bf(y ~ x),
         fam = mixture(cumulative(threshold = "equidistant"),
                       sratio(threshold = "sum_to_zero")),
         comp = c("cumulative", "sratio"),
         type = c("equidistant", "sum_to_zero"),
         link = c("logit", "logit"), shared = FALSE)
  )
  for (s in shapes) {
    lab <- paste(s$comp, s$type, collapse = " + ")
    fit <- frm(s$f, family = s$fam, data = d)
    for (move in c(0, 0.15)) {
      p <- omx_par(fit, move)
      est <- fit$obj$env$parList(p)
      k1 <- 3L
      ordered <- s$comp %in% c("cumulative", "hurdle_cumulative")
      if (s$shared) {
        # one vector, ordered when any component is, and centered by
        # a sum-to-zero component
        tau0 <- omx_tau(est$tau_raw, "flexible", any(ordered), k1)
        taus <- lapply(s$type, function(ty) {
          if (ty == "sum_to_zero") tau0 - mean(tau0) else tau0
        })
      } else {
        taus <- lapply(1:2, function(k) {
          omx_tau(est[[paste0("tau_raw", k)]], s$type[k], ordered[k], k1)
        })
      }
      b <- est$beta
      bd <- est$betad
      disc1 <- if ("disc1_z" %in% names(bd)) exp(bd[["disc1_z"]] * d$z) else 1
      th <- bd[["theta1_(Intercept)"]] +
        if ("theta1_z" %in% names(bd)) bd[["theta1_z"]] * d$z else 0
      p1 <- stats::plogis(th)
      P <- p1 * omx_brms_probs(s$comp[1], b[["mu1_x"]] * d$x, disc1,
                               taus[[1]], s$link[1]) +
        (1 - p1) * omx_brms_probs(s$comp[2], b[["mu2_x"]] * d$x, 1,
                                  taus[[2]], s$link[2])
      ref <- omx_ll(P, d$y)
      got <- -fit$obj$fn(p)
      expect_lt(abs(got - ref), 1e3 * .Machine$double.eps * abs(ref),
                label = paste(lab, "move", move))
      if (move == 0) {
        Pf <- unname(fitted(fit)[, "Estimate", ])
        expect_lt(max(abs(Pf - P)), 1e3 * .Machine$double.eps, label = lab)
        nd <- d[1:7, ]
        Pn <- unname(fitted(fit, newdata = nd)[, "Estimate", ])
        expect_lt(max(abs(Pn - P[1:7, ])), 1e3 * .Machine$double.eps,
                  label = lab)
      }
    }
  }
})

test_that("cs() and thres(gr = ) reach each component, as brms reads them", {
  skip_unless_brms()
  d <- omx_data(20261006)
  # cs(): each component its own coefficients, off its own thresholds
  fit <- frm(bf(y ~ cs(x)), family = mixture(sratio(), acat()), data = d)
  expect_true(all(c("bcs_mu1_x[1]", "bcs_mu2_x[3]") %in% variables(fit)))
  for (move in c(0, 0.1)) {
    est <- fit$obj$env$parList(omx_par(fit, move))
    p1 <- stats::plogis(est$betad[["theta1_(Intercept)"]])
    cs <- function(k) {
      lp <- Filter(function(l) identical(l$dpar, paste0("mu", k)),
                   fit$frame$linpreds)[[1L]]
      outer(d$x, est[[lp$cs[[1L]]$par]])
    }
    # brms reads `Intercept - transpose(mucs[n])`: each row its own
    # thresholds, with the latent predictor at zero
    t1 <- matrix(est$tau_raw1, nrow(d), 3, byrow = TRUE) - cs(1)
    t2 <- matrix(est$tau_raw2, nrow(d), 3, byrow = TRUE) - cs(2)
    P <- p1 * omx_brms_probs("sratio", rep(0, nrow(d)), 1, t1, "logit") +
      (1 - p1) * omx_brms_probs("acat", rep(0, nrow(d)), 1, t2, "logit")
    ref <- omx_ll(P, d$y)
    expect_lt(abs(-fit$obj$fn(omx_par(fit, move)) - ref),
              1e3 * .Machine$double.eps * abs(ref))
  }
  expect_lt(max(abs(unname(fitted(fit)[, "Estimate", ]) - {
    est <- fit$estimates
    p1 <- stats::plogis(est$betad[["theta1_(Intercept)"]])
    lp1 <- Filter(function(l) identical(l$dpar, "mu1"),
                  fit$frame$linpreds)[[1L]]
    lp2 <- Filter(function(l) identical(l$dpar, "mu2"),
                  fit$frame$linpreds)[[1L]]
    t1 <- matrix(est$tau_raw1, nrow(d), 3, byrow = TRUE) -
      outer(d$x, est[[lp1$cs[[1L]]$par]])
    t2 <- matrix(est$tau_raw2, nrow(d), 3, byrow = TRUE) -
      outer(d$x, est[[lp2$cs[[1L]]$par]])
    p1 * omx_brms_probs("sratio", rep(0, nrow(d)), 1, t1, "logit") +
      (1 - p1) * omx_brms_probs("acat", rep(0, nrow(d)), 1, t2, "logit")
  })), 1e3 * .Machine$double.eps)

  # thres(gr = g): each component a threshold vector per level
  fg <- frm(bf(y | thres(gr = g) ~ x), family = mixture(cratio(), acat()),
            data = d)
  expect_true(all(c("b_mu1_Intercept[a,1]", "b_mu2_Intercept[b,3]") %in%
                    variables(fg)))
  for (move in c(0, 0.1)) {
    est <- fg$obj$env$parList(omx_par(fg, move))
    p1 <- stats::plogis(est$betad[["theta1_(Intercept)"]])
    gi <- as.integer(d$g)
    rowtau <- function(r) {
      M <- matrix(r, 2, 3, byrow = TRUE)
      M[gi, ]
    }
    P <- p1 * omx_brms_probs("cratio", est$beta[["mu1_x"]] * d$x, 1,
                             rowtau(est$tau_raw1), "logit") +
      (1 - p1) * omx_brms_probs("acat", est$beta[["mu2_x"]] * d$x, 1,
                                rowtau(est$tau_raw2), "logit")
    ref <- omx_ll(P, d$y)
    expect_lt(abs(-fg$obj$fn(omx_par(fg, move)) - ref),
              1e3 * .Machine$double.eps * abs(ref))
    if (move == 0) {
      expect_lt(max(abs(unname(fitted(fg)[, "Estimate", ]) - P)),
                1e3 * .Machine$double.eps)
    }
  }
})

test_that("hurdle components keep their hurdles inside the mixture", {
  skip_unless_brms()
  d <- omx_data(20261007)
  fit <- frm(bf(yh ~ x, hu1 ~ z),
             family = mixture(hurdle_cumulative(),
                              hurdle_cumulative("probit")),
             data = d)
  for (move in c(0, 0.1)) {
    est <- fit$obj$env$parList(omx_par(fit, move))
    p1 <- stats::plogis(est$betad[["theta1_(Intercept)"]])
    hu1 <- stats::plogis(est$betad[["hu1_(Intercept)"]] +
                           est$betad[["hu1_z"]] * d$z)
    hu2 <- stats::plogis(est$betad[["hu2_(Intercept)"]])
    # brms:::posterior_epred_hurdle_cumulative(): cbind(hu, (1 - hu) *
    # dcumulative()), for each component
    comp <- function(hu, eta, raw, link) {
      cbind(hu, (1 - hu) * omx_brms_probs("cumulative", eta, 1,
                                          omx_tau(raw, "flexible", TRUE, 3),
                                          link))
    }
    P <- p1 * comp(hu1, est$beta[["mu1_x"]] * d$x, est$tau_raw1, "logit") +
      (1 - p1) * comp(hu2, est$beta[["mu2_x"]] * d$x, est$tau_raw2,
                      "probit")
    ref <- omx_ll(P, d$yh, code0 = 0L)
    expect_lt(abs(-fit$obj$fn(omx_par(fit, move)) - ref),
              1e3 * .Machine$double.eps * abs(ref))
    if (move == 0) {
      Pf <- fitted(fit)[, "Estimate", ]
      expect_identical(colnames(Pf), paste0("P(Y = ", 0:4, ")"))
      expect_lt(max(abs(unname(Pf) - P)), 1e3 * .Machine$double.eps)
    }
  }
})

test_that("brms's names: per component, shared, equidistant, multivariate", {
  d <- omx_data(20261008)
  f1 <- frm(bf(y ~ x), family = mixture(cumulative(), cumulative()),
            data = d)
  v <- variables(f1)
  expect_true(all(c(paste0("b_mu1_Intercept[", 1:3, "]"),
                    paste0("b_mu2_Intercept[", 1:3, "]"),
                    "b_mu1_x", "b_mu2_x", "theta1", "theta2") %in% v))
  expect_true(all(c(paste0("mu1_Intercept[", 1:3, "]"), "mu2_x") %in%
                    rownames(fixef(f1))))
  # the reported thresholds are the ones the components read
  fx <- fixef(f1)
  expect_identical(unname(fx[paste0("mu1_Intercept[", 1:3, "]"),
                             "Estimate"]),
                   omx_tau(f1$estimates$tau_raw1, "flexible", TRUE, 3))
  h <- hypothesis(f1, "mu2_Intercept[2] - mu1_Intercept[2] = 0")
  expect_true(is.finite(h$hypothesis$Estimate))
  # order = "mu": both components report the one shared vector
  f2 <- frm(bf(y ~ x), family = mixture(cumulative(),
                                        sratio(threshold = "sum_to_zero"),
                                        order = "mu"), data = d)
  fx <- fixef(f2)
  t1 <- fx[paste0("mu1_Intercept[", 1:3, "]"), "Estimate"]
  t2 <- fx[paste0("mu2_Intercept[", 1:3, "]"), "Estimate"]
  expect_lt(max(abs(t2 - (t1 - mean(t1)))),
            1e3 * .Machine$double.eps * max(abs(t1)))
  # equidistant: one delta per component, as brms's transformed
  # parameters name it
  f3 <- frm(bf(y ~ x), family = mixture(cumulative(threshold = "equidistant"),
                                        acat(threshold = "equidistant")),
            data = d)
  expect_true(all(c("delta_mu1", "delta_mu2") %in% variables(f3)))
  tk <- fixef(f3)[paste0("mu2_Intercept[", 1:3, "]"), "Estimate"]
  # each threshold is tau_1 + (k - 1) delta, rounded on the way, so the
  # second difference is a few roundings of the largest threshold
  # (measured: 0, and 0.33 of that unit on component 1;
  # dev/ordmix-p1-m2.R)
  expect_lte(abs(diff(diff(tk))), 4 * .Machine$double.eps * max(abs(tk)))
  # the Links line names each component's distribution function and
  # hides a disc held at 1, as brms does. It lists theta1, which brms
  # leaves out when theta1 has no predictor: listing every mixing weight
  # is frmtmb's convention
  f4 <- frm(bf(y ~ x), family = mixture(cumulative("probit"), sratio()),
            data = d)
  expect_output(print(f4), "Links: mu1 = probit; mu2 = logit; theta1 =",
                fixed = TRUE)
  # multivariate: the response joins the prefix
  f5 <- frm(bf(y ~ x) + bf(yh ~ x),
            family = list(mixture(cumulative(), cumulative()),
                          hurdle_cumulative()), data = d)
  expect_true(all(c("b_mu1_y_Intercept[1]", "b_mu2_y_Intercept[3]",
                    "b_yh_Intercept[1]", "theta1_y") %in% variables(f5)))
})

test_that("default_prior() lists brms's threshold rows by component", {
  skip_unless_brms()
  d <- omx_data(20261009)
  rows <- function(p) {
    p <- as.data.frame(p)
    keep <- p$class %in% c("Intercept", "delta", "b") & !nzchar(p$coef) &
      !p$dpar %in% c("theta1", "theta2")
    sort(paste(p$class, p$group, p$dpar, sep = "|")[keep])
  }
  cases <- list(
    list(f = y ~ x, frm = mixture(cumulative(), sratio()),
         brm = brms::mixture(brms::cumulative(), brms::sratio())),
    list(f = y ~ x, frm = mixture(cumulative(), sratio(), order = "mu"),
         brm = brms::mixture(brms::cumulative(), brms::sratio(),
                             order = "mu")),
    list(f = y | thres(gr = g) ~ x, frm = mixture(cratio(), acat()),
         brm = brms::mixture(brms::cratio(), brms::acat())),
    list(f = y ~ x,
         frm = mixture(cumulative(threshold = "equidistant"), acat()),
         brm = brms::mixture(brms::cumulative(threshold = "equidistant"),
                             brms::acat()))
  )
  for (cs in cases) {
    ours <- default_prior(frmtmb::bf(cs$f), data = d, family = cs$frm)
    theirs <- suppressMessages(brms::default_prior(brms::bf(cs$f), data = d,
                                                   family = cs$brm))
    expect_identical(rows(ours), rows(theirs), label = deparse1(cs$f))
  }
})

test_that("a threshold prior lands where brms puts it", {
  d <- omx_data(20261010)
  lp_tot <- function(fit) {
    ent <- frmtmb:::resolve_prior_input(list(frame = fit$frame,
                                             spec = fit$spec),
                                        fit$prior)$entries
    -frmtmb:::neg_log_prior_fn(ent)(fit$estimates)
  }
  # order = "none", dpar = "mu2": brms's prior is on Intercept_mu2, the
  # thresholds of the CENTERED design, tau - mean(x) * b_mu2_x, with
  # the log-Jacobian of the ordered vector
  f1 <- frm(bf(y ~ x), family = mixture(cumulative(), cumulative()),
            data = d,
            prior = set_prior("normal(0, 1)", class = "Intercept",
                              dpar = "mu2"))
  est <- f1$estimates
  tau2 <- omx_tau(est$tau_raw2, "flexible", TRUE, 3)
  ref <- sum(stats::dnorm(tau2 - mean(d$x) * est$beta[["mu2_x"]], 0, 1,
                          log = TRUE)) + sum(est$tau_raw2[-1])
  expect_lt(abs(lp_tot(f1) - ref), 1e3 * .Machine$double.eps * abs(ref))
  # order = "mu": no dpar, on the shared vector, and brms centers no
  # design there; ordered because one component is cumulative
  f2 <- frm(bf(y ~ x), family = mixture(cumulative(), sratio(),
                                        order = "mu"), data = d,
            prior = set_prior("normal(0, 1)", class = "Intercept"))
  est <- f2$estimates
  tau <- omx_tau(est$tau_raw, "flexible", TRUE, 3)
  ref <- sum(stats::dnorm(tau, 0, 1, log = TRUE)) + sum(est$tau_raw[-1])
  expect_lt(abs(lp_tot(f2) - ref), 1e3 * .Machine$double.eps * abs(ref))
  # each spelling is refused where brms has no such row
  expect_error(frm(bf(y ~ x), family = mixture(cumulative(), cumulative()),
                   data = d, prior = set_prior("normal(0, 1)",
                                               class = "Intercept")),
               "names one of them with dpar =", fixed = TRUE)
  expect_error(frm(bf(y ~ x), family = mixture(cumulative(), cumulative(),
                                               order = "mu"), data = d,
                   prior = set_prior("normal(0, 1)", class = "Intercept",
                                     dpar = "mu1")),
               "Prior target not found", fixed = TRUE)
})

test_that("simulate() draws from the fitted category distribution", {
  d <- omx_data(20261011, n = 300)
  # this fit ends with component 1's thresholds collapsed, where the gap
  # floor binds (without it the fit stopped at a NaN gradient), and it
  # says so; the draws still come from the fitted distribution
  fit <- allow_warnings(
    frm(bf(y | thres(gr = g) ~ x), family = mixture(cumulative(), sratio()),
        data = d),
    "degenerate boundary", require = "degenerate boundary")
  # the collapsed thresholds have no standard error, so the rows whose
  # probabilities move with them get none either, said once (lane
  # setier: fitted() reads the covariance that holds them)
  P <- allow_warnings(fitted(fit), "move along a direction")[, "Estimate", ]
  sims <- simulate(fit, nsim = 400, seed = 7)
  codes <- vapply(sims, as.integer, integer(nrow(d)))
  freq <- vapply(1:4, function(k) mean(codes == k), 0)
  pbar <- colMeans(P)
  # the binomial standard error of each pooled frequency, from the
  # fitted probabilities themselves
  se <- sqrt(colSums(P * (1 - P))) / (nrow(d) * sqrt(400))
  expect_true(all(abs(freq - pbar) < 5 * se))
})

test_that("the refusals hold, and each absent case fits", {
  d <- omx_data(20261012)
  expect_error(mixture(cumulative(threshold = "equidistant"), cumulative(),
                       order = "mu"),
               "Cannot use equidistant and fixed thresholds", fixed = TRUE)
  f <- frm(bf(y ~ x), family = mixture(cumulative(threshold = "equidistant"),
                                       cumulative()), data = d)
  expect_true(is.finite(as.numeric(logLik(f))))
  expect_error(mixture(hurdle_cumulative(), cumulative()),
               "Cannot mix hurdle_cumulative() with an ordinal family",
               fixed = TRUE)
  expect_error(mixture(cumulative(), cumulative(), groups = ~g),
               "does not take ordinal components", fixed = TRUE)
  expect_error(mixture(cumulative(), gaussian()),
               "Cannot mix families with real and integer support",
               fixed = TRUE)
  # cs() on a cumulative() component fits since 0.68.0, with brms's
  # warning once per such component (test-cumulative-cs.R)
  # component 2 runs to a step function on these data, and says so
  f <- allow_warnings(frm(bf(y ~ x, mu2 ~ cs(x)),
                          family = mixture(cumulative(), sratio()),
                          data = d), "degenerate boundary")
  expect_true("bcs_mu2_x[1]" %in% variables(f))
  # order = "mu" stays refused for components that have a mu intercept
  expect_error(mixture(gaussian(), gaussian(), order = "mu"),
               "is not supported", fixed = TRUE)
})

test_that("cs() outside mu is refused: it used to move the thresholds", {
  # On 0.67.0 the offsets of every predictor went into the one slot the
  # densities read, so disc ~ cs(z) fitted y ~ x + cs(z) under disc's
  # name (dev/ordmix-base-behavior.R: the two logLiks agree to 1.1e-9,
  # and the fit reports the cs coefficients as disc_z[k]). brms refuses
  # it in these words.
  d <- omx_data(20261013)
  expect_error(frm(bf(y ~ x, disc ~ cs(z)), family = sratio(), data = d),
               "Category specific effects are only supported for the main",
               fixed = TRUE)
  expect_error(frm(bf(y ~ x, disc2 ~ cs(z)),
                   family = mixture(sratio(), sratio()), data = d),
               "Category specific effects are only supported for the main",
               fixed = TRUE)
  # brms's sentence on a family that is not ordinal as well, and the
  # mixture's predictors named only on a mixture
  d$yc <- d$x + d$z
  err <- tryCatch(frm(bf(yc ~ x, sigma ~ cs(z)), family = gaussian(),
                      data = d), error = conditionMessage)
  expect_match(err, "Category specific effects are only supported for the ",
               fixed = TRUE)
  expect_no_match(err, "mixture", fixed = TRUE)
  err <- tryCatch(frm(bf(y ~ x, disc ~ cs(z)), family = sratio(), data = d),
                  error = conditionMessage)
  expect_no_match(err, "mixture", fixed = TRUE)
})

test_that("shared thresholds start the components apart", {
  # two classes on one set of thresholds, apart in their slopes only:
  # the model order = "mu" describes
  set.seed(20261014)
  d <- data.frame(x = rnorm(400))
  cls <- stats::rbinom(400, 1, 0.4)
  lat <- ifelse(cls == 1, 2 * d$x, -1.5 * d$x) + stats::rlogis(400)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  fam <- mixture(cumulative(), cumulative(), order = "mu")
  fit <- frm(bf(y ~ x), family = fam, data = d)
  se <- sqrt(diag(vcov(fit)))
  b <- fixef(fit)[c("mu1_x", "mu2_x"), "Estimate"]
  expect_gt(abs(b[[1]] - b[[2]]), 3 * max(se[c("mu1_x", "mu2_x")]))
  # the reason for the spread: from equal starts the components stay
  # equal, a stationary point the optimizer does not leave
  stuck <- allow_warnings(frm(bf(y ~ x), family = fam, data = d,
                              start = list(beta = c(0, 0))),
                          c("NaN", "Hessian", "gradient", "converge",
                            "Standard errors are not available"))
  # the two components are the same function of the same numbers at
  # every step, so the slopes stay the same double (dev/ordmix-p1-m2.R)
  bs <- stuck$estimates$beta
  expect_identical(bs[[1]], bs[[2]])
  expect_lt(as.numeric(logLik(stuck)), as.numeric(logLik(fit)))
})

test_that("the flat directions of a hurdle mixture and of shared thresholds
          warn, and only those", {
  d <- omx_data(20261015)
  allow_warnings(
    frm(bf(yh ~ x), family = mixture(hurdle_cumulative(),
                                     hurdle_cumulative()), data = d),
    c("the likelihood sees only P(Y = 0)", "NaN", "Hessian", "gradient",
      "converge"),
    require = "the likelihood sees only P(Y = 0)")
  # absent: hu1 has a predictor, so P(Y = 0) moves with it
  expect_no_warning(frm(bf(yh ~ x, hu1 ~ z),
                        family = mixture(hurdle_cumulative(),
                                         hurdle_cumulative()), data = d))
  # absent: hu1 held at a value, 2 parameters for 2 numbers
  expect_no_warning(frm(bf(yh ~ x, hu1 = 0.1),
                        family = mixture(hurdle_cumulative(),
                                         hurdle_cumulative()), data = d))
  # shared thresholds, no mu predictor, a predictor on theta alone: the
  # components are one distribution
  allow_warnings(
    frm(bf(y ~ 1, theta1 ~ z), family = mixture(cumulative(), cumulative(),
                                                order = "mu"), data = d),
    c("the components are one", "NaN", "Hessian", "gradient", "converge"),
    require = "the components are one")
  # no predictor anywhere: one category distribution, either order
  for (o in c("none", "mu")) {
    allow_warnings(
      frm(bf(y ~ 1), family = mixture(cumulative(), cumulative(), order = o),
          data = d),
      c("every row has the same category", "NaN", "Hessian", "gradient",
        "converge"),
      require = "every row has the same category")
  }
  # absent: the same mixtures with a slope in mu
  expect_no_warning(frm(bf(y ~ x), family = mixture(cumulative(),
                                                    cumulative(),
                                                    order = "mu"), data = d))
  expect_no_warning(frm(bf(y ~ x), family = mixture(cumulative(),
                                                    cumulative()), data = d))
})

# The review's data (dev/ordmix-rev-falsealarm.R): theta1 moves with z,
# and at one threshold vector (-1, 0, 1) a cumulative and an sratio
# component give different category distributions; the hurdle
# components differ in hu.
omx_fa_data <- function(seed, n = 500) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  hu <- plogis(ifelse(cls == 1, -2, -0.4) + 0.8 * z)
  yh <- ifelse(runif(n) < hu, 0L, y)
  tau <- c(-1, 0, 1)
  pc <- diff(c(0, plogis(tau), 1))
  h <- plogis(tau)
  ps <- c(h[1], (1 - h[1]) * h[2], (1 - h[1]) * (1 - h[2]) * h[3],
          (1 - h[1]) * (1 - h[2]) * (1 - h[3]))
  c1 <- runif(n) < plogis(0.3 + 1.5 * z)
  ysh <- vapply(seq_len(n), function(i) {
    sample.int(4L, 1L, prob = if (c1[i]) pc else ps)
  }, 1L)
  data.frame(x, z, y, yh, ysh)
}

# every warning a fit raises, muffled and returned
omx_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  w
}

test_that("the flat-direction warnings stay off identified mixtures", {
  # Punch round 1, B1. On the first lane build each of these warned on
  # 20 of 20 seeds, though each is identified (dev/ordmix-rev-fa-log/).
  d <- omx_fa_data(3)
  # shared thresholds, no mu predictor, but the components are of
  # different families, or of one family with different links: two
  # distributions, so theta is in the likelihood
  for (fam in list(mixture(cumulative(), sratio(), order = "mu"),
                   mixture(cumulative(), cumulative("probit"),
                           order = "mu"))) {
    w <- omx_warnings(frm(bf(ysh ~ 1, theta1 ~ z), family = fam, data = d))
    expect_false(any(grepl("the components are one", w, fixed = TRUE)),
                 label = fam$family)
  }
  # hurdle components held apart by priors on hu1 and hu2, the warning's
  # own remedy
  w <- omx_warnings(
    frm(bf(yh ~ x), family = mixture(hurdle_cumulative(),
                                     hurdle_cumulative()), data = d,
        prior = set_prior("beta(4, 16)", class = "hu1") +
          set_prior("beta(8, 12)", class = "hu2")))
  expect_false(any(grepl("the likelihood sees only P(Y = 0)", w,
                         fixed = TRUE)))
  # and the no-predictor warning: priors on the thresholds place them
  w <- omx_warnings(
    frm(bf(y ~ 1), family = mixture(cumulative(), cumulative()), data = d,
        prior = set_prior("normal(0, 1)", class = "Intercept",
                          dpar = "mu1") +
          set_prior("normal(0, 1)", class = "Intercept", dpar = "mu2")))
  expect_false(any(grepl("every row has the same category", w,
                         fixed = TRUE)))
})

test_that("a component run off to a degenerate boundary warns", {
  # Punch round 1, B2: the review's data-generating process (seeds 1..40,
  # dev/ordmix-rev-degen.R); seed 19 at n = 300 is one of its degenerate
  # fits, a component a step function of x, latent distance 7817
  # (dev/ordmix-p1-degen-log/rev_cum_cum_300.txt)
  gen <- function(seed) {
    set.seed(seed)
    x <- rnorm(300)
    cls <- rbinom(300, 1, 0.4)
    lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(300)
    data.frame(x = x, y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5))
  }
  w <- omx_warnings(frm(bf(y ~ x), family = mixture(cumulative(),
                                                    cumulative()),
                        data = gen(19)))
  expect_true(any(grepl("degenerate boundary", w, fixed = TRUE)))
  # absent: an identified fit of the same process (seed 3, distance 6.3)
  w <- omx_warnings(frm(bf(y ~ x), family = mixture(cumulative(),
                                                    cumulative()),
                        data = gen(3)))
  expect_false(any(grepl("degenerate boundary", w, fixed = TRUE)))
})

# the review's strong-predictor process (dev/ordmix-rev2-margin.R): two
# classes, slopes 1.5 c and -0.8 c, latent noise from the link's own
# distribution
omx_margin_data <- function(seed, link, cc, n) {
  noise <- switch(link, logit = stats::rlogis, cauchit = stats::rcauchy)
  set.seed(seed)
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * cc * x + 1, -0.8 * cc * x - 1) + noise(n)
  cut_at <- c(-1.5, 0, 1.5) * max(1, cc / 2)
  data.frame(x = x, y = 1L + (lat > cut_at[1]) + (lat > cut_at[2]) +
               (lat > cut_at[3]))
}

test_that("strong predictors are not a degenerate boundary", {
  # Punch round 2, B3. The round-1 cut, a latent distance of 50, fired
  # on 13 of 16 sound fits at slopes of 15 and -8 (n = 500) and on 20 of
  # 20 at n = 2000. Seed 1 is sound: every standard error finite, the
  # slopes 12% from the truth, latent distance 58
  # (dev/ordmix-p2-crit-log/margin_logit_10_500.txt).
  w <- omx_warnings(frm(bf(y ~ x), family = mixture(cumulative(),
                                                    cumulative()),
                        data = omx_margin_data(1, "logit", 10, 500)))
  expect_false(any(grepl("degenerate boundary", w, fixed = TRUE)))
  # a cauchit step function: latent distance 2.4e5, doubling its
  # discrimination gains 0.0025 in log-likelihood
  w <- omx_warnings(frm(bf(y ~ x), family = mixture(cumulative("cauchit"),
                                                    cumulative("cauchit")),
                        data = omx_margin_data(4, "cauchit", 1, 500)))
  expect_true(any(grepl("a step function of its predictors", w,
                        fixed = TRUE)))
})

test_that("a threshold prior does not place the hu and theta direction", {
  # Punch round 2, B4: round 1 silenced the flat-direction warnings on
  # any prior at all (dev/ordmix-rev2-b1.R, seeds 1..10: 0 of 10 warned,
  # with standard errors of 59.5 to 296 or NaN in all 10)
  d <- omx_fa_data(1)
  hm <- mixture(hurdle_cumulative(), hurdle_cumulative())
  w <- omx_warnings(frm(bf(yh ~ x), family = hm, data = d,
                        prior = set_prior("normal(0, 2)",
                                          class = "Intercept",
                                          dpar = "mu1")))
  expect_true(any(grepl("the likelihood sees only P(Y = 0)", w,
                        fixed = TRUE)))
  # absent: a prior on hu1 places it (K = 2: 3 parameters, 2 numbers)
  w <- omx_warnings(frm(bf(yh ~ x), family = hm, data = d,
                        prior = set_prior("beta(4, 16)", class = "hu1")))
  expect_false(any(grepl("the likelihood sees only P(Y = 0)", w,
                         fixed = TRUE)))
})

test_that("with no predictor, priors are counted against the proportions", {
  # Punch round 2, B4: a prior on one component's thresholds holds 3 of
  # 7 parameters, and the data identify 3 proportions; round 1 went
  # silent, with fixef() standard errors NaN in 7 of 10 fits
  d <- omx_fa_data(1)
  fam <- mixture(cumulative(), cumulative())
  w <- omx_warnings(frm(bf(y ~ 1), family = fam, data = d,
                        prior = set_prior("normal(0, 2)",
                                          class = "Intercept",
                                          dpar = "mu1")))
  expect_true(any(grepl("no distributional parameter has a predictor", w,
                        fixed = TRUE)))
  # absent: both components' thresholds held leaves 1 free parameter
  w <- omx_warnings(frm(bf(y ~ 1), family = fam, data = d,
                        prior = set_prior("normal(0, 2)",
                                          class = "Intercept",
                                          dpar = "mu1") +
                          set_prior("normal(0, 2)", class = "Intercept",
                                    dpar = "mu2")))
  expect_false(any(grepl("no distributional parameter has a predictor", w,
                         fixed = TRUE)))
  # c1: a hurdle mixture with thres(gr = ) identifies one zero
  # proportion for all groups, as hu and theta are shared: 5
  # thresholds (3 + 2) and one zero, Hessian rank 6 of 13
  # (dev/ordmix-rev3-rank.R). Priors on mu1's 5 thresholds and on hu1
  # leave 7 free, one more than that; round 2 counted 7 identified and
  # went silent.
  set.seed(20261007)
  n <- 1500
  g <- factor(sample(c("a", "b"), n, TRUE))
  y4 <- sample(1:4, n, TRUE, prob = c(0.3, 0.2, 0.25, 0.25))
  y3 <- sample(1:3, n, TRUE, prob = c(0.4, 0.35, 0.25))
  dg <- data.frame(g = g, y = ifelse(g == "a", y4, y3))
  runif(n)
  dg$yh2 <- ifelse(runif(n) < ifelse(dg$g == "a", 0.1, 0.35), 0L, dg$y)
  w <- omx_warnings(frm(bf(yh2 | thres(gr = g) ~ 1),
                        family = mixture(hurdle_cumulative(),
                                         hurdle_cumulative()), data = dg,
                        prior = set_prior("normal(0, 2)",
                                          class = "Intercept",
                                          dpar = "mu1") +
                          set_prior("beta(2, 8)", class = "hu1")))
  expect_true(any(grepl("no distributional parameter has a predictor", w,
                        fixed = TRUE)))
})

# the review's degenerate-fit process (dev/ordmix-rev-degen.R), with
# case weights
omx_wdata <- function(seed, n, w) {
  set.seed(seed)
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  data.frame(x = x, y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5),
             w = w)
}

test_that("constant weights do not change the degenerate verdict", {
  # Punch round 2b, B5: the sharpening cut is an absolute change in
  # log-likelihood, and the check summed weights as given, so weights
  # 1/n warned on 20 of 20 fits, 11 of them sound
  # (dev/ordmix-rev3-b3.R, w_sum1). Seed 3 is sound: slopes 20% from
  # the truth, largest standard error 11.8.
  f <- function(seed, n, w) {
    frm(bf(y | weights(w) ~ x), family = mixture(cumulative(),
                                                 cumulative()),
        data = omx_wdata(seed, n, w))
  }
  w <- omx_warnings(f(3, 500, 1 / 500))
  expect_false(any(grepl("degenerate boundary", w, fixed = TRUE)))
  # weights 1 and 10 give the same verdict, on a degenerate fit (seed
  # 19, n = 300) and on a sound one. The two optimizer paths end within
  # 1e-5 of each other, not at the same doubles, so the sharpening is
  # compared at ONE fit's parameters, its weights scaled by 10, 1/300
  # and 1
  for (s in c(19, 3)) {
    w1 <- character(0)
    w10 <- character(0)
    f1 <- withCallingHandlers(f(s, 300, 1), warning = function(cnd) {
      w1 <<- c(w1, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
    f10 <- withCallingHandlers(f(s, 300, 10), warning = function(cnd) {
      w10 <<- c(w10, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
    expect_identical(any(grepl("degenerate boundary", w1, fixed = TRUE)),
                     any(grepl("degenerate boundary", w10, fixed = TRUE)))
    sh <- function(c) {
      g <- f1
      g$frame$aterm_values$y$weights <- c * f1$frame$aterm_values$y$weights
      frmtmb:::mixture_ord_degeneracy(g, "y")$sharpen
    }
    expect_equal(sh(10), sh(1))
    expect_equal(sh(1 / 300), sh(1))
  }
})

test_that("the degenerate check says only what it found", {
  # three components on two-class data (dev/ordmix-cum3two.R): one
  # component gives NaN for a category on every row, and the check's
  # max() over no number warned "no non-missing arguments to max"
  # beside its own warning. Any warning but the check's fails here, with
  # one exception: three components on two classes have a ridge, and
  # whether nlminb stops on it with "singular convergence (7)" is
  # rounding (it does with OpenBLAS 0.3.26, dev/ciharden-findings.md).
  # That is the optimizer's verdict, not the check's, so it may come.
  set.seed(20261005 + 11)
  n <- 400
  x <- rnorm(n)
  rnorm(n)
  sample(c("a", "b"), n, TRUE)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  d <- data.frame(y = as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf))),
                  x = x)
  allow_warnings(
    frm(bf(y ~ x), family = mixture(cumulative(), cumulative(),
                                    cumulative()), data = d),
    c("degenerate boundary", "Optimizer did not report convergence"),
    require = "degenerate boundary")
})

test_that("thres(x = ) above the data warns once per threshold block", {
  d <- omx_data(20261016)
  # allowed by its threshold, not by "bound": the degenerate warning
  # says "boundary", and a threshold above the data is not a degenerate
  # component (it runs off whatever the component does, and the check
  # leaves it out; its latent distance was 2.2e7 before)
  allow_warnings(
    frm(bf(y | thres(4) ~ x), family = mixture(cumulative(), sratio()),
        data = d),
    "Intercept[4]", require = "Intercept[4]")
  w <- character()
  withCallingHandlers(
    frm(bf(y | thres(4) ~ x), family = mixture(cumulative(), sratio(),
                                               order = "mu"), data = d),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  expect_identical(sum(grepl("Intercept[4]", w, fixed = TRUE)), 1L)
})

test_that("conditional_effects() shows the mixture's category probabilities", {
  d <- omx_data(20261017)
  fit <- frm(bf(y ~ x, theta1 ~ z), family = mixture(cumulative(), acat()),
             data = d)
  ce <- conditional_effects(fit, effects = "x")[[1L]]
  expect_identical(levels(ce$cats__), as.character(1:4))
  nd <- ce[ce$cats__ == "1", c("x", "z")]
  P <- fitted(fit, newdata = nd)[, "Estimate", ]
  for (k in 1:4) {
    expect_lt(max(abs(ce$estimate__[ce$cats__ == as.character(k)] - P[, k])),
              1e3 * .Machine$double.eps)
  }
  expect_true(all(is.finite(ce$se__)))
})

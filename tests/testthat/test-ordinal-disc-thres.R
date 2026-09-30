# brms's discrimination parameter `disc` and its threshold structures
# (`threshold = "equidistant"` and `"sum_to_zero"`) on the ordinal
# families. The reference is brms 2.23.0's own R-side densities,
# brms:::dcumulative() and its three siblings, which are what brms's
# posterior_epred() and log_lik() read. They take `disc` and the
# thresholds as arguments, so they check frmtmb's density, simulator and
# threshold maps without compiling anything. The Stan-side identity is
# row 12e of test-brms-likelihood.R (gated).

ordt_data <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
  d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
  d$yh <- ifelse(runif(n) < stats::plogis(-1 + 0.3 * d$x), 0L, d$y)
  d
}

# brms's category probabilities at frmtmb's estimates, one row per
# observation: brms reads a draws-by-thresholds matrix, and the rows
# here are the observations instead
ordt_brms_probs <- function(family, eta, disc, tau, link) {
  dens <- get(paste0("d", family), asNamespace("brms"))
  n <- length(eta)
  dens(seq_len(length(tau) + 1L), eta = eta,
       thres = matrix(tau, n, length(tau), byrow = TRUE), disc = disc,
       link = link)
}

# the log prior a MAP fit carries, at its estimates, entry by entry
ordt_log_prior <- function(fit) {
  ent <- frmtmb:::resolve_prior_input(list(frame = fit$frame,
                                           spec = fit$spec),
                                      fit$prior)$entries
  -frmtmb:::neg_log_prior_fn(ent)(fit$estimates)
}

ordt_parts <- function(fit) {
  fam <- family(fit)
  list(eta = as.numeric(frm_linpred(fit, dpar = "mu", type = "link")),
       disc = as.numeric(frm_linpred(fit, dpar = "disc",
                                     type = "response")),
       tau = frmtmb:::ord_threshold_values(fam, fit$estimates$tau_raw),
       link = fam[["ord_link"]][["name"]])
}

test_that("disc and every threshold structure are brms's density", {
  skip_unless_brms()
  d <- ordt_data(20260930)
  for (fam in c("cumulative", "sratio", "cratio", "acat")) {
    for (th in c("flexible", "equidistant", "sum_to_zero")) {
      lab <- paste(fam, th)
      fit <- frm(bf(y ~ x, disc ~ 0 + z),
                 family = get(fam)(threshold = th), data = d)
      p <- ordt_parts(fit)
      P <- ordt_brms_probs(fam, p$eta, p$disc, p$tau, p$link)
      Pf <- fitted(fit)[, "Estimate", ]
      expect_lt(max(abs(Pf - P)), 1e3 * .Machine$double.eps, label = lab)
      ll <- sum(log(P[cbind(seq_len(nrow(d)), d$y)]))
      ours <- as.numeric(logLik(fit))
      expect_lt(abs(ours - ll) / abs(ours), 1e-12, label = lab)
      # the simulator's probabilities are the density's
      Ps <- frmtmb:::ord_cat_probs(fam, p$eta, p$tau, NULL, p$link, p$disc)
      expect_lt(max(abs(Ps - P)), 1e3 * .Machine$double.eps, label = lab)
      # the structure itself holds
      if (th == "equidistant") {
        dt <- diff(p$tau)
        expect_lt(max(abs(dt - dt[1L])) / abs(dt[1L]), 1e-12, label = lab)
      }
      if (th == "sum_to_zero") {
        expect_lt(abs(sum(p$tau)) / max(abs(p$tau)), 1e-12, label = lab)
      }
    }
  }
})

test_that("acat off the logit is brms's second form", {
  skip_unless_brms()
  d <- ordt_data(20260939)
  osa <- methods::getClass("osa", where = asNamespace("RTMB"))
  for (lk in c("probit", "probit_approx", "cloglog", "cauchit", "softit")) {
    fit <- frm(bf(y ~ x, disc ~ 0 + z),
               family = acat(lk, threshold = "equidistant"), data = d)
    p <- ordt_parts(fit)
    P <- if (lk == "probit_approx") {
      # brms's R side reads probit_approx as the exact pnorm while its
      # Stan program, the likelihood, reads Phi_approx; frmtmb follows
      # the likelihood, so the reference is brms's acat formula
      # (brms:::inv_link_acat()) written out with Phi_approx
      X <- p$disc * (p$eta - matrix(p$tau, nrow(d), 4L, byrow = TRUE))
      Fx <- stats::plogis(0.07056 * X^3 + 1.5976 * X)
      U <- vapply(1:5, function(k) {
        a <- if (k > 1L) apply(Fx[, seq_len(k - 1L), drop = FALSE], 1L,
                               prod) else 1
        b <- if (k < 5L) apply(1 - Fx[, k:4, drop = FALSE], 1L, prod) else 1
        a * b
      }, numeric(nrow(d)))
      U / rowSums(U)
    } else {
      ordt_brms_probs("acat", p$eta, p$disc, p$tau, lk)
    }
    expect_lt(max(abs(fitted(fit)[, "Estimate", ] - P)),
              1e3 * .Machine$double.eps, label = lk)
    ours <- as.numeric(logLik(fit))
    ll <- sum(log(P[cbind(seq_len(nrow(d)), d$y)]))
    expect_lt(abs(ours - ll) / abs(ours), 1e-12, label = lk)
    Ps <- frmtmb:::ord_cat_probs("acat", p$eta, p$tau, NULL, lk, p$disc)
    expect_lt(max(abs(Ps - P)), 1e3 * .Machine$double.eps, label = lk)
    ff <- family(fit)
    dp <- frmtmb:::eval_dpars(fit)[[1L]]
    ex <- list(tau_raw = fit$estimates$tau_raw)
    a <- ff$lpdf(d$y, dp, list(), ex)
    b <- ff$lpdf(methods::new(osa, x = as.numeric(d$y),
                              keep = matrix(1, nrow(d), 1)), dp, list(), ex)
    expect_lt(max(abs(as.numeric(b) - a)) / max(abs(a)),
              1e3 * .Machine$double.eps, label = lk)
  }
  # grouped thresholds read the same form on each level's own slice
  fit <- frm(y | thres(gr = g) ~ x, family = acat("cloglog"), data = d)
  p <- ordt_parts(fit)
  nth <- family(fit)$thres$nthres
  Pf <- fitted(fit)[, "Estimate", ]
  gi <- as.integer(d$g)
  for (g in 1:2) {
    tg <- p$tau[sum(nth[seq_len(g - 1L)]) + seq_len(nth[g])]
    rows <- which(gi == g)
    P <- ordt_brms_probs("acat", p$eta[rows], p$disc[rows], tg, "cloglog")
    expect_lt(max(abs(Pf[rows, seq_len(nth[g] + 1L)] - P)),
              1e3 * .Machine$double.eps)
  }
})

test_that("disc reaches the grouped thresholds and the hurdle", {
  skip_unless_brms()
  d <- ordt_data(20260931)
  fit <- frm(bf(y | thres(gr = g) ~ x, disc ~ 0 + z),
             family = acat(threshold = "equidistant"), data = d)
  p <- ordt_parts(fit)
  nth <- family(fit)$thres$nthres
  Pf <- fitted(fit)[, "Estimate", ]
  gi <- as.integer(d$g)
  for (g in 1:2) {
    tg <- p$tau[sum(nth[seq_len(g - 1L)]) + seq_len(nth[g])]
    rows <- which(gi == g)
    P <- ordt_brms_probs("acat", p$eta[rows], p$disc[rows], tg, "logit")
    expect_lt(max(abs(Pf[rows, seq_len(nth[g] + 1L)] - P)),
              1e3 * .Machine$double.eps)
  }
  fh <- frm(bf(yh ~ x, disc ~ 0 + z),
            family = hurdle_cumulative(threshold = "sum_to_zero"), data = d)
  p <- ordt_parts(fh)
  hu <- as.numeric(frm_linpred(fh, dpar = "hu", type = "response"))
  P <- cbind(hu, (1 - hu) * ordt_brms_probs("cumulative", p$eta, p$disc,
                                            p$tau, "logit"))
  expect_lt(max(abs(fitted(fh)[, "Estimate", ] - P)),
            1e3 * .Machine$double.eps)
  expect_lt(abs(sum(p$tau)) / max(abs(p$tau)), 1e-12)
})

test_that("the structures are nested in the flexible model", {
  d <- ordt_data(20260932)
  ll <- vapply(c("flexible", "equidistant", "sum_to_zero"), function(th) {
    as.numeric(logLik(frm(y ~ x, family = sratio(threshold = th),
                          data = d)))
  }, 0)
  # an equidistant or zero-sum threshold vector is one the flexible
  # model can also take, so its maximum is never above the flexible one
  expect_gte(ll[["flexible"]], ll[["equidistant"]])
  expect_gte(ll[["flexible"]], ll[["sum_to_zero"]])
  # and they are fits of different models, not the same one relabelled
  expect_gt(ll[["flexible"]] - ll[["equidistant"]], 0)
})

test_that("a structure on too few thresholds is refused, and only then", {
  d <- ordt_data(20260933)
  d$y2 <- ifelse(d$y > 2, 2L, 1L)
  d$y3 <- pmin(d$y, 3L)
  # two categories: flexible fits, the other two have nothing to hold
  # brms's message suggesting bernoulli() is not the point here
  expect_s3_class(suppressMessages(frm(y2 ~ x, family = cumulative(),
                                       data = d)),
                  "frmtmb_fit")
  expect_error(frm(y2 ~ x, family = cumulative(threshold = "equidistant"),
                   data = d), "delta, the distance")
  # a single zero-sum threshold is 0: no threshold parameter at all,
  # and a model all the same (glm(y == 1 ~ 0 + I(-x)) is the same one)
  f0 <- suppressMessages(frm(y2 ~ x, family = cumulative(threshold =
                                                            "sum_to_zero"),
                             data = d))
  expect_length(f0$estimates$tau_raw, 0L)
  g0 <- stats::glm(I(y2 == 1) ~ 0 + I(-x), family = binomial, data = d)
  expect_lt(abs(as.numeric(logLik(f0)) - as.numeric(logLik(g0))) /
              abs(as.numeric(logLik(g0))), 1e-10)
  # three categories: both structures fit
  for (th in c("equidistant", "sum_to_zero")) {
    expect_s3_class(frm(y3 ~ x, family = cumulative(threshold = th),
                        data = d), "frmtmb_fit")
  }
  # a thres(gr = ) level of one threshold: sum_to_zero holds it at zero
  # and fits, as brms does; equidistant has no delta to place there
  d$y3g <- ifelse(d$g == "a", pmin(d$y, 2L), d$y3)
  fs <- frm(y3g | thres(gr = g) ~ x,
            family = cumulative(threshold = "sum_to_zero"), data = d)
  tau <- frmtmb:::ord_threshold_values(family(fs), fs$estimates$tau_raw)
  expect_identical(tau[1L], 0)
  expect_length(tau, 3L)
  expect_error(frm(y3g | thres(gr = g) ~ x,
                   family = cumulative(threshold = "equidistant"), data = d),
               "level(s) 'a' of thres(gr = )", fixed = TRUE)
  expect_error(cumulative(threshold = "equi"), "takes one of")
})

test_that("sum-to-zero with one threshold in every thres(gr = ) level", {
  skip_unless_brms()
  # brms declares a grouped count int<lower=1>, so this is a model it
  # runs: every threshold is 0, and b and disc are estimated
  d <- ordt_data(20261001)
  d$y2 <- 1L + (stats::rlogis(nrow(d)) + 0.8 * d$x > 0.3)
  for (fam in c("cumulative", "sratio", "acat")) {
    f <- suppressMessages(
      frm(bf(y2 | thres(gr = g) ~ x, disc ~ 0 + z),
          family = get(fam)(threshold = "sum_to_zero"), data = d))
    expect_length(f$estimates$tau_raw, 0L)
    expect_true(all(c("b_Intercept[a,1]", "b_Intercept[b,1]") %in%
                      variables(f)), label = fam)
    p <- ordt_parts(f)
    expect_identical(p$tau, c(0, 0))
    P <- ordt_brms_probs(fam, p$eta, p$disc, 0, "logit")
    expect_lt(max(abs(fitted(f)[, "Estimate", ] - P)),
              1e3 * .Machine$double.eps, label = fam)
    ll <- sum(log(P[cbind(seq_len(nrow(d)), d$y2)]))
    expect_lt(abs(as.numeric(logLik(f)) - ll) / abs(ll), 1e-12, label = fam)
    expect_s3_class(summary(f), "summary.frmtmb_fit")
  }
})

# The density at a CHOSEN parameter vector, disc away from 1, against
# brms's R-side density there. At an optimum a density that ignored disc
# would still agree with brms, because the disc coefficient would stay at
# its start of 0 (disc = 1): the review's mutants that dropped disc from
# acat's logit path, from cratio and from the grouped sequential density
# passed every test that compared at the fit's own optimum.
ordt_fixed_point <- function(fit, disc_coef) {
  par <- fit$opt$par
  j <- which(names(par) == "betad")
  par[j] <- disc_coef
  pl <- fit$obj$env$parList(par)
  list(par = par, pl = pl)
}

test_that("each density reads disc at a fixed point, not only at the optimum", {
  skip_unless_brms()
  d <- ordt_data(20261002)
  shapes <- list(
    list(fam = "cumulative", link = "logit"),
    list(fam = "cumulative", link = "probit"),
    list(fam = "sratio", link = "logit"),
    list(fam = "sratio", link = "cauchit"),
    list(fam = "cratio", link = "logit"),
    list(fam = "cratio", link = "probit"),
    list(fam = "acat", link = "logit"),
    list(fam = "acat", link = "cloglog"))
  for (s in shapes) {
    for (grouped in c(FALSE, TRUE)) {
      lab <- paste(s$fam, s$link, if (grouped) "thres(gr = g)")
      form <- if (grouped) bf(y | thres(gr = g) ~ x, disc ~ 0 + z) else
        bf(y ~ x, disc ~ 0 + z)
      fit <- frm(form, family = get(s$fam)(s$link), data = d)
      fp <- ordt_fixed_point(fit, 0.7)
      tau <- frmtmb:::ord_threshold_values(family(fit), fp$pl$tau_raw)
      eta <- d$x * fp$pl$beta[["x"]]
      disc <- exp(0.7 * d$z)
      ll <- 0
      if (grouped) {
        nth <- family(fit)$thres$nthres
        for (g in 1:2) {
          rows <- which(as.integer(d$g) == g)
          tg <- tau[sum(nth[seq_len(g - 1L)]) + seq_len(nth[g])]
          P <- ordt_brms_probs(s$fam, eta[rows], disc[rows], tg, s$link)
          ll <- ll + sum(log(P[cbind(seq_along(rows), d$y[rows])]))
        }
      } else {
        P <- ordt_brms_probs(s$fam, eta, disc, tau, s$link)
        ll <- sum(log(P[cbind(seq_len(nrow(d)), d$y)]))
      }
      ours <- -fit$obj$fn(fp$par)
      expect_lt(abs(ours - ll) / abs(ll), 1e-12, label = lab)
    }
  }
})

test_that("every threshold map inverts to the internal vector", {
  # the inverse is what frmtmb.sample hands a stored draw back through,
  # so a wrong one moves posterior_epred() while every fit agrees
  set.seed(20261003)
  for (type in c("flexible", "equidistant", "sum_to_zero")) {
    for (ordered in c(TRUE, FALSE)) {
      lay <- frmtmb:::thres_layout(c(4L, 2L, 3L), type)
      raw <- stats::rnorm(sum(lay$rlen))
      tau <- frmtmb:::thres_tau(raw, lay, ordered)
      back <- frmtmb:::thres_raw_from_tau(tau, lay, ordered)
      expect_lt(max(abs(back - raw)) / max(abs(raw)),
                1e3 * .Machine$double.eps,
                label = paste(type, if (ordered) "ordered" else "unordered"))
    }
  }
})

test_that("an intercept in disc warns on every ordinal family", {
  d <- ordt_data(20260934)
  for (fam in c("cumulative", "sratio", "cratio", "acat")) {
    allow_warnings(
      frm(bf(y ~ x, disc ~ z), family = get(fam)(), data = d),
      paste0(fam, ": disc has an intercept"),
      require = paste0(fam, ": disc has an intercept"))
  }
  # absent without the intercept, and absent when a prior holds it
  expect_no_warning(frm(bf(y ~ x, disc ~ 0 + z), family = acat(), data = d))
  expect_no_warning(frm(bf(y ~ x, disc ~ z), family = acat(), data = d,
                        prior = set_prior("normal(0, 1)", class = "Intercept",
                                          dpar = "disc")))
  # cell means span the intercept without a column of that name
  d$h <- factor(sample(c("p", "q"), nrow(d), TRUE))
  allow_warnings(
    frm(bf(y ~ x, disc ~ 0 + h), family = sratio(), data = d),
    "sratio: disc has an intercept (its columns, hp, hq, add up to one)",
    require = "sratio: disc has an intercept (its columns, hp, hq")
  # and a prior on them holds it
  expect_no_warning(frm(bf(y ~ x, disc ~ 0 + h), family = sratio(),
                        data = d, prior = set_prior("normal(0, 1)",
                                                    class = "b",
                                                    dpar = "disc")))
})

test_that("disc held at 1 is shown nowhere, as brms shows it nowhere", {
  d <- ordt_data(20261004)
  fit <- frm(y ~ x, family = cumulative(), data = d)
  s <- summary(fit)
  expect_identical(s$links, "cdf = logit")
  expect_length(s$fixed_dpars, 0L)
  expect_false("disc" %in% names(s$coefficients))
  expect_identical(names(fixef(fit, flatten = TRUE)), "x")
  out <- utils::capture.output(print(fit))
  expect_false(any(grepl("disc", out, fixed = TRUE)))
  # fixef(flatten = TRUE) lines up with confint() again
  expect_true(all(names(fixef(fit, flatten = TRUE)) %in%
                    rownames(confint(fit))))
  # the mapped coefficient is still there for the objective
  expect_identical(fit$frame$linpreds[["y.disc"]]$constant, 1)
  # modeled, it shows, and so does a value the user fixed
  fm <- frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(), data = d)
  expect_identical(summary(fm)$links, "cdf = logit; disc = log")
  expect_true("disc_z" %in% names(fixef(fm, flatten = TRUE)))
  f2 <- frm(bf(y ~ x, disc = 2), family = cumulative(), data = d)
  expect_identical(summary(f2)$fixed_dpars[["disc"]], 2)
  # a multivariate coef() lists no disc block
  d$y2 <- pmin(d$y, 3L)
  mv <- frm(bf(y ~ x) + bf(y2 ~ x), data = d,
            family = list(cumulative(), sratio()))
  expect_false(any(grepl("disc", names(coef(mv)), fixed = TRUE)))
})

test_that("default_prior() shows brms's delta bound and per-level rows", {
  d <- ordt_data(20261005)
  dp <- default_prior(y | thres(gr = g) ~ x, data = d,
                      family = cumulative(threshold = "equidistant"))
  dl <- dp[dp$class == "delta", ]
  expect_identical(sort(dl$group), c("", "a", "b"))
  expect_true(all(dl$lb == 0))
  expect_identical(sort(dp$group[dp$class == "Intercept"]), c("", "a", "b"))
  # sratio holds delta itself, unbounded
  ds <- default_prior(y ~ x, data = d,
                      family = sratio(threshold = "equidistant"))
  expect_true(is.na(ds$lb[ds$class == "delta"]))
  # each listed row is one the resolver takes
  expect_s3_class(frm(y | thres(gr = g) ~ x, data = d,
                      family = cumulative(),
                      prior = set_prior("normal(0, 3)", class = "Intercept",
                                        group = "a")),
                  "frmtmb_fit")
})

test_that("class delta and class Intercept follow brms's parameters", {
  d <- ordt_data(20260935)
  # brms's default_prior() for sratio(threshold = "equidistant") with a
  # cse() term lists b four times, delta and one Intercept (priors:14)
  dp <- default_prior(y ~ x + z + cse(g), data = d,
                      family = sratio(threshold = "equidistant"))
  expect_setequal(dp$class, c("b", "delta", "Intercept"))
  expect_identical(sum(dp$class == "delta"), 1L)
  expect_identical(sum(dp$class == "Intercept"), 1L)
  # under sum_to_zero brms's Intercept is a vector frmtmb does not have
  dz <- default_prior(y ~ x, data = d,
                      family = cumulative(threshold = "sum_to_zero"))
  expect_false("Intercept" %in% dz$class)
  expect_error(frm(y ~ x, family = cumulative(threshold = "sum_to_zero"),
                   data = d, prior = set_prior("normal(0, 1)",
                                               class = "Intercept")),
               "sum_to_zero")
  # and delta exists only under equidistant
  expect_error(frm(y ~ x, family = cumulative(), data = d,
                   prior = set_prior("normal(0, 1)", class = "delta")),
               "class \"delta\" is the distance")
  # the MAP penalty is brms's density: delta on its own scale, with the
  # log-Jacobian of the log that cumulative() holds it on, and the first
  # threshold at the mean of the predictors
  pl <- set_prior("normal(1.5, 0.2)", class = "delta") +
    set_prior("normal(-1, 0.5)", class = "Intercept")
  fit <- frm(y ~ x, family = cumulative(threshold = "equidistant"),
             data = d, prior = pl)
  tau <- frmtmb:::ord_threshold_values(family(fit), fit$estimates$tau_raw)
  delta <- tau[2L] - tau[1L]
  first <- tau[1L] - mean(d$x) * fixef(fit)["x", "Estimate"]
  lprior <- ordt_log_prior(fit)
  ref <- stats::dnorm(delta, 1.5, 0.2, log = TRUE) + log(delta) +
    stats::dnorm(first, -1, 0.5, log = TRUE)
  expect_lt(abs(lprior - ref) / abs(ref), 1e-10)
  # an unordered family holds delta itself, with no Jacobian
  fs <- frm(y ~ x, family = sratio(threshold = "equidistant"), data = d,
            prior = set_prior("normal(0.5, 0.1)", class = "delta"))
  tau <- frmtmb:::ord_threshold_values(family(fs), fs$estimates$tau_raw)
  lprior <- ordt_log_prior(fs)
  ref <- stats::dnorm(tau[2L] - tau[1L], 0.5, 0.1, log = TRUE)
  expect_lt(abs(lprior - ref) / abs(ref), 1e-10)
})

test_that("delta is named as brms names it and reported", {
  d <- ordt_data(20260936)
  fit <- frm(y ~ x, family = cumulative(threshold = "equidistant"), data = d)
  v <- variables(fit)
  expect_true(all(c(paste0("b_Intercept[", 1:4, "]"), "delta") %in% v))
  tau <- frmtmb:::ord_threshold_values(family(fit), fit$estimates$tau_raw)
  h <- hypothesis(fit, "delta = 0", class = NULL)$hypothesis
  expect_lt(abs(h$Estimate - (tau[2L] - tau[1L])) / (tau[2L] - tau[1L]),
            1e-12)
  sp <- summary(fit)$spec_pars
  expect_true("delta" %in% rownames(sp))
  expect_lt(abs(sp["delta", "Estimate"] - h$Estimate) / h$Estimate, 1e-12)
  # its interval is the log-scale Wald interval mapped back, so it cannot
  # reach zero
  expect_gt(sp["delta", 3L], 0)
  fg <- frm(y | thres(gr = g) ~ x,
            family = sratio(threshold = "equidistant"), data = d)
  expect_true(all(c("delta_1", "delta_2") %in% variables(fg)))
  # the threshold count, not the parameter count, sizes a cs() term and
  # pins a refit
  fc <- frm(y ~ x + cs(z), family = sratio(threshold = "equidistant"),
            data = d)
  expect_true(all(paste0("bcs_z[", 1:4, "]") %in% variables(fc)))
  expect_identical(frmtmb:::thres_pin_of_fit(fc)[[1L]]$nthres, 4L)
})

test_that("the fitted paths read the structure on new data", {
  d <- ordt_data(20260937)
  fit <- frm(bf(y ~ x, disc ~ 0 + z),
             family = cratio(threshold = "sum_to_zero"), data = d)
  nd <- d[1:5, ]
  expect_identical(dim(fitted(fit, newdata = nd)), c(5L, 4L, 5L))
  expect_lt(max(abs(fitted(fit, newdata = nd)[, "Estimate", ] -
                      fitted(fit)[1:5, "Estimate", ])),
            1e3 * .Machine$double.eps)
  sim <- simulate(fit, nsim = 2, seed = 1)
  expect_true(all(unlist(sim) %in% 1:5))
  r <- residuals(fit, type = "osa")[, "Estimate"]
  expect_true(all(is.finite(r)))
})

test_that("the one-step density reads disc and the structure", {
  # with every row kept, the one-step branch of a density (the response
  # on the tape, the category picked arithmetically) is the data branch
  d <- ordt_data(20260938)
  osa <- methods::getClass("osa", where = asNamespace("RTMB"))
  for (fam in c("cumulative", "sratio", "cratio", "acat")) {
    for (th in c("flexible", "equidistant", "sum_to_zero")) {
      f <- frm(bf(y ~ x, disc ~ 0 + z), family = get(fam)(threshold = th),
               data = d)
      ff <- family(f)
      dp <- frmtmb:::eval_dpars(f)[[1L]]
      ex <- list(tau_raw = f$estimates$tau_raw)
      a <- ff$lpdf(d$y, dp, list(), ex)
      b <- ff$lpdf(methods::new(osa, x = as.numeric(d$y),
                                keep = matrix(1, nrow(d), 1)),
                   dp, list(), ex)
      expect_lt(max(abs(as.numeric(b) - a)) / max(abs(a)),
                1e3 * .Machine$double.eps, label = paste(fam, th))
    }
  }
})

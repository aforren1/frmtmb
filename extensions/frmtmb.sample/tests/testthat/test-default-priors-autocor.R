# frm_sample()'s default prior on the residual-correlation classes is
# brms 2.23.0's (flat) on the NATURAL parameter: on (-1, 1) for a
# cov = TRUE ar or ma coefficient, on [0, 1] for cosy, and nothing for
# the cov = FALSE coefficients, which are the sampled parameters
# themselves (brms's default_prior() and stancode(),
# dev/acprior-brms-defaults.R). Through frmtmb.sample 0.22.0 the default
# was flat on frmtmb's internal parameter instead: (1 - rho^2)^(-3/2)
# on a first-order rho, and a logit-flat cosy over its whole window.
# On the data of dev/acprior-agree.R the ar posterior mean was 42.3
# Monte Carlo standard errors from brms's own, and 48 percent of the
# cosy draws were below brms's bound (dev/acprior-findings.md).

dpa_data <- function(seed = 4, ng = 12, k = 5) {
  set.seed(seed)
  d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                  x = stats::rnorm(ng * k))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(ng)[d$g] * 0.3 + stats::rnorm(ng * k)
  d
}

dpa_ri <- function(fit, prior = NULL) {
  suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, prior))$ri
}

# Flat on the natural coefficients means the sampled density on the
# internal ones is the map's Jacobian alone, so -log prior + log |J| is
# the same at every internal point. The spread is measured against the
# size of the two parts, which is what a cancellation is relative to;
# flat on the internal parameter leaves the spread of log |J| itself
dpa_flat_spread <- function(fit, ri, ths, logjac) {
  nlp <- neg_log_prior_fn(ri$entries)
  parts <- vapply(ths, function(th) {
    p <- fit$frame$par_template
    p[["thetaac"]][] <- th
    c(nlp(p), logjac(th))
  }, c(0, 0))
  v <- colSums(parts)
  diff(range(v)) / max(abs(parts))
}

# log |det| of the map from the internal parameters to the natural
# ones, by central differences of core's autocor_natural(): independent
# of the closed form the prior itself adds (ac_trans_logjac()), so a
# wrong Jacobian there cannot cancel itself here
dpa_fd_logjac <- function(ac) {
  function(th) {
    f <- function(t) as.numeric(frmtmb:::autocor_natural(t, ac))
    J <- vapply(seq_along(th), function(j) {
      h <- 1e-5 * max(1, abs(th[j]))
      e <- replace(numeric(length(th)), j, h)
      (f(th + e) - f(th - e)) / (2 * h)
    }, numeric(length(th)))
    log(abs(det(as.matrix(J))))
  }
}

dpa_messages <- function(expr) {
  out <- character(0)
  withCallingHandlers(expr, message = function(m) {
    out <<- c(out, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  out
}

test_that("the default is flat on a cov = TRUE coefficient itself", {
  d <- dpa_data()
  for (fo in list(y ~ x + ar(t, g, cov = TRUE), y ~ x + ma(t, g, cov = TRUE),
                  y ~ x + ar(t, g, p = 2, cov = TRUE),
                  y ~ x + arma(t, g, cov = TRUE))) {
    fit <- frm(fo, data = d, dry_run = "objective")
    ri <- dpa_ri(fit)
    # the density on thetaac alone: drop sigma's, which reads thetaac too
    # on these terms (the innovation sd, tested below)
    ri$entries <- Filter(function(e) identical(e$comp, "thetaac"),
                         ri$entries)
    ac <- fit$frame$autocor[[1L]]
    set.seed(1)
    ths <- replicate(5, stats::rnorm(ac$npar, 0, 1.5), simplify = FALSE)
    expect_lt(dpa_flat_spread(fit, ri, ths, dpa_fd_logjac(ac)), 1e-6,
              label = deparse(fo))
    # and no bound narrows it: the uniform spans the map's whole image
    expect_length(ri$lower, 0L)
  }
})

test_that("cov = FALSE coefficients get no density, as brms's flat", {
  d <- dpa_data()
  fit <- frm(y ~ x + arma(t, g), data = d, dry_run = "objective")
  ri <- dpa_ri(fit)
  hit <- Filter(function(e) identical(e$comp, "thetaac"), ri$entries)
  expect_length(hit, 0L)
})

test_that("cosy is flat on brms's [0, 1], and a user prior keeps lb = 0", {
  d <- dpa_data()
  fit <- frm(y ~ x + cosy(t, g), data = d, dry_run = "objective")
  a <- 1 / (5 - 1)
  ri <- dpa_ri(fit)
  # the image of 0 under the window's logistic map
  expect_equal(ri$lower[["thetaac_1"]], stats::qlogis(a / (1 + a)))
  lj <- dpa_fd_logjac(fit$frame$autocor[[1L]])
  ri$entries <- Filter(function(e) identical(e$comp, "thetaac"), ri$entries)
  expect_lt(dpa_flat_spread(fit, ri, list(-1, 0, 0.5, 2, 4), lj), 1e-6)
  # a density of the user's own keeps brms's declared bound
  ru <- dpa_ri(fit, set_prior("normal(0, 0.1)", class = "cosy"))
  expect_equal(ru$lower[["thetaac_1"]], stats::qlogis(a / (1 + a)))
  # unless it sets its own
  rl <- dpa_ri(fit, set_prior("normal(0, 0.1)", class = "cosy", lb = -0.1))
  expect_equal(rl$lower[["thetaac_1"]],
               stats::qlogis((-0.1 + a) / (1 + a)))
})

test_that("a bounds-only prior keeps the flat density on the coefficient", {
  # brms's set_prior("", class = "ar", lb = 0.5) is flat on [0.5, 1];
  # without the default's density the slot would be flat on frmtmb's
  # internal parameter instead
  d <- dpa_data()
  fit <- frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = "objective")
  ri <- dpa_ri(fit, set_prior("", class = "ar", lb = 0.5))
  expect_equal(ri$lower[["thetaac_1"]], 0.5 / sqrt(1 - 0.25))
  lj <- dpa_fd_logjac(fit$frame$autocor[[1L]])
  ri$entries <- Filter(function(e) identical(e$comp, "thetaac"), ri$entries)
  expect_lt(dpa_flat_spread(fit, ri, list(0.6, 1, 2, 5), lj), 1e-6)
  # the absent case: a density of the user's own is the one applied, on
  # the coefficient, with no bound added, as brms's normal_lpdf(ar | ...)
  ru <- dpa_ri(fit, set_prior("normal(0.5, 0.2)", class = "ar"))
  e <- Filter(function(e) identical(e$comp, "thetaac"), ru$entries)
  expect_length(e, 1L)
  expect_identical(e[[1L]]$dist$inner$kind, "normal")
  expect_length(ru$lower, 0L)
})

test_that("the default on sigma is brms's, on sigma itself", {
  # brms's cholesky_cor_ar1(), _ma1() and _arma1() are covariances, so
  # its sigma under a cov = TRUE term is the innovation sd, and so is
  # frmtmb's since sigma became the innovation sd. The default
  # student_t(3, 0, s) is brms's density at sigma = exp(eta), with the
  # log link's Jacobian and no factor of the correlation parameters.
  # frmtmb.sample 0.23.0 read it at sigma * c with the Jacobian of c,
  # because frmtmb's sigma was the marginal sd.
  d <- dpa_data()
  for (fo in list(y ~ x + ar(t, g, cov = TRUE), y ~ x + ma(t, g, cov = TRUE),
                  y ~ x + arma(t, g, cov = TRUE),
                  y ~ x + ar(t, g, p = 2, cov = TRUE), y ~ x + cosy(t, g))) {
    fit <- frm(fo, data = d, dry_run = "objective")
    ri <- dpa_ri(fit)
    ac <- fit$frame$autocor[[1L]]
    sg <- Filter(function(e) !identical(e$comp, "thetaac") &&
                   identical(e$dist$kind, "t") && identical(e$scale, "sd"),
                 ri$entries)
    expect_length(sg, 1L)
    expect_null(sg[[1L]]$innov)
    sc <- sg[[1L]]$dist$scale
    nlp <- neg_log_prior_fn(sg)
    set.seed(2)
    v <- vapply(1:5, function(k) {
      p <- fit$frame$par_template
      p[["thetaac"]][] <- stats::rnorm(ac$npar)
      eta <- stats::rnorm(1, -0.5, 0.5)
      p[[sg[[1L]]$comp]][sg[[1L]]$idx] <- eta
      # brms's density on sigma = exp(eta), carried to eta
      ref <- stats::dt(exp(eta) / sc, 3, log = TRUE) - log(sc) + eta
      c(nlp(p), ref)
    }, c(0, 0))
    expect_lt(diff(range(colSums(v))) / max(abs(v)), 1e-10,
              label = deparse(fo))
  }
  # under sigma ~ x the default is on brms's Intercept_sigma, the log
  # innovation sd at the predictor means, which no correlation
  # parameter moves
  fx <- frm(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x), data = d,
            dry_run = "objective")
  ex <- Filter(function(e) identical(e$comp, "betad"), dpa_ri(fx)$entries)
  expect_length(ex, 1L)
  nlx <- neg_log_prior_fn(ex)
  set.seed(4)
  p <- fx$frame$par_template
  p[["betad"]][] <- c(-0.4, 0.3)
  vx <- vapply(1:5, function(k) {
    p[["thetaac"]][] <- stats::rnorm(1, 0, 2)
    nlx(p)
  }, 0)
  expect_identical(diff(range(vx)), 0)
  msg <- dpa_messages(frmtmb.sample:::sample_resolve_priors(
    frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = "objective"),
    NULL))
  expect_false(any(grepl("innovation sd", msg, fixed = TRUE)))
})

test_that("a bounded parameter is sampled through one logistic", {
  # Stan's log bound on the internal cosy, stacked on the window's
  # logistic, gave 5 or 6 divergent transitions per 4 chains of 2000 at
  # cosy near 0, where brms and this chart have none
  # (dev/acprior-cosydiv.R)
  d <- dpa_data()
  fit <- frm(y ~ x + cosy(t, g), data = d, dry_run = "objective")
  ri <- dpa_ri(fit)
  ch <- autocor_bound_chart(fit, ri$lower, ri$upper)
  expect_length(ch, 1L)
  tm <- ch[[1L]]
  expect_identical(c(tm$L, tm$U), c(0, 1))
  # v -> internal -> natural is the logistic onto (0, 1), and back
  v <- c(-6, -1, 0, 2, 7)
  th <- tm$to_internal(v)
  expect_equal(vapply(th, tm$to_natural, 0), stats::plogis(v))
  expect_equal(vapply(th, tm$from_internal, 0), v)
  # its log Jacobian against central differences
  fd <- vapply(v, function(x) {
    h <- 1e-5 * max(1, abs(x))
    log(abs((tm$to_internal(x + h) - tm$to_internal(x - h)) / (2 * h)))
  }, 0)
  lj <- vapply(v, tm$logjac, 0)
  expect_lt(max(abs(lj - fd)) / max(abs(fd)), 1e-6)
  # no chart without a bound, or for a cov = FALSE coefficient, whose
  # bound is Stan's own transform onto it already
  fa <- frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = "objective")
  ra <- dpa_ri(fa)
  expect_length(autocor_bound_chart(fa, ra$lower, ra$upper), 0L)
  fb <- frm(y ~ x + ar(t, g), data = d, dry_run = "objective")
  rb <- dpa_ri(fb, set_prior("", class = "ar", lb = 0))
  expect_length(autocor_bound_chart(fb, rb$lower, rb$upper), 0L)
})

test_that("a mode below cosy's default bound says where the bound is from", {
  # the default lb = 0 is brms's declaration, which the user did not
  # write, so the warning names it and the spelling of frm()'s window
  set.seed(11)
  d <- data.frame(t = rep(1:4, 30), g = rep(1:30, each = 4),
                  x = stats::rnorm(120))
  # residuals less half their group mean: a negative cosy
  e <- stats::rnorm(120)
  d$y <- 1 + 0.5 * d$x + e - 0.5 * stats::ave(e, d$g)
  fit <- frm(y ~ x + cosy(t, g), data = d)
  est <- as.numeric(frmtmb:::autocor_natural(fit$estimates[["thetaac"]],
                                             fit$frame$autocor[[1L]]))
  expect_lt(est, 0)
  ri <- dpa_ri(fit)
  ch <- autocor_bound_chart(fit, ri$lower, ri$upper)
  w <- tryCatch(frmtmb.sample:::chart_mode_warning(fit, ch),
                warning = conditionMessage)
  expect_match(w, "lower bound brms declares for cosy", fixed = TRUE)
  expect_match(w, "set_prior(\"\", class = \"cosy\", lb = -1/3)",
               fixed = TRUE)
  # the absent case: a mode inside the interval warns of nothing
  dp <- d
  dp$y <- 1 + 0.5 * d$x + e + stats::rnorm(30)[d$g]
  fp <- frm(y ~ x + cosy(t, g), data = dp)
  rp <- dpa_ri(fp)
  expect_gt(as.numeric(frmtmb:::autocor_natural(fp$estimates[["thetaac"]],
                                                fp$frame$autocor[[1L]])), 0)
  expect_null(frmtmb.sample:::chart_mode_warning(
    fp, autocor_bound_chart(fp, rp$lower, rp$upper)))
  # and a bound the user wrote is called the requested one
  ru <- dpa_ri(fit, set_prior("", class = "cosy", lb = est + 0.05))
  wu <- tryCatch(frmtmb.sample:::chart_mode_warning(
    fit, autocor_bound_chart(fit, ru$lower, ru$upper)),
    warning = conditionMessage)
  expect_match(wu, "outside the requested interval", fixed = TRUE)
})

test_that("the sampling route's table and message name brms's interval", {
  d <- dpa_data()
  tab <- default_prior(y ~ x + cosy(t, g), data = d, route = "sample")
  expect_equal(unlist(tab[tab$class == "cosy", c("lb", "ub")],
                      use.names = FALSE), c(0, 1))
  ta <- default_prior(y ~ x + ar(t, g, cov = TRUE), data = d,
                      route = "sample")
  expect_equal(unlist(ta[ta$class == "ar", c("lb", "ub")],
                      use.names = FALSE), c(-1, 1))
  expect_identical(ta$prior[ta$class == "ar"], "uniform(-1, 1)")
  # cov = FALSE: brms's unbounded (flat), which needs no density
  tc <- default_prior(y ~ x + ar(t, g), data = d, route = "sample")
  expect_identical(tc$prior[tc$class == "ar"], "(flat)")
  expect_true(all(is.na(unlist(tc[tc$class == "ar", c("lb", "ub")]))))
  fit <- frm(y ~ x + cosy(t, g), data = d, dry_run = "objective")
  msg <- dpa_messages(frmtmb.sample:::sample_resolve_priors(fit, NULL))
  expect_true(any(grepl("cosy +uniform[(]-0.25, 1[)], lb = 0", msg)))
})

test_that("a registered term's bounded parameter follows the same rule", {
  # the rule is core's, read through autocor_flat_priors(), so a term an
  # extension registers (frmtmb.spatial's sar()) gets it with no list of
  # names here: a uniform on its interval, kept under a bounds-only prior
  frmtmb::frmtmb_local_registries()
  lo <- -0.4
  hi <- 0.9
  frmtmb::frmtmb_register_autocor(
    "acpterm", args = "V",
    build = function(spec, rows, n_all, lookup) list(),
    npar = 1L,
    natural = function(theta, aux) {
      c(acprho = lo + (hi - lo) / (1 + exp(-theta[1])))
    },
    bounds = function(aux) cbind(lower = lo, upper = hi),
    whiten = function(y, mu, sigma, theta, aux) {
      list(u = (y - mu) / sigma, logjac = -length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) mu + sigma * u,
    sigma_varies = FALSE)
  set.seed(3)
  d <- data.frame(x = stats::rnorm(30))
  d$y <- 1 + d$x + stats::rnorm(30)
  fit <- frm(y ~ x + acpterm(V), data = d, data2 = list(V = 1),
             dry_run = "objective")
  bnd <- dpa_fd_logjac(fit$frame$autocor[[1L]])
  only_ac <- function(ri) {
    ri$entries <- Filter(function(e) identical(e$comp, "thetaac"),
                         ri$entries)
    ri
  }
  ri <- only_ac(dpa_ri(fit))
  expect_lt(dpa_flat_spread(fit, ri, list(-2, 0, 1, 3), bnd), 1e-6)
  rb <- only_ac(dpa_ri(fit, set_prior("", class = "acprho", lb = 0)))
  expect_equal(rb$lower[["thetaac_1"]], stats::qlogis(-lo / (hi - lo)))
  expect_lt(dpa_flat_spread(fit, rb, list(-0.2, 0, 1, 3), bnd), 1e-6)
  # the bound is sampled through the chart too, onto (0, 0.9)
  ch <- autocor_bound_chart(fit, rb$lower, rb$upper)
  expect_identical(c(ch[[1L]]$L, ch[[1L]]$U), c(0, hi))
  # and its internal spelling is refused by the class's name
  expect_error(dpa_ri(fit, set_prior("normal(0, 1)", class = "theta",
                                     coef = "thetaac_1")),
               "internal parameter of acprho", class = "frmtmb_error")
})

test_that("the sampling table passed back is the defaults' density", {
  # default_prior(route = "sample") passed back as prior = (and the
  # default's text typed as a user prior) is the density the defaults
  # are: through fix round 1 of this change the sigma row read the
  # marginal sd as a user prior and the innovation sd as a default, up
  # to 0.72 nats apart (dev/review-acprior-r1-trip.R)
  d <- dpa_data()
  for (fo in list(y ~ x + ar(t, g, cov = TRUE), y ~ x + ma(t, g, cov = TRUE),
                  y ~ x + arma(t, g, cov = TRUE))) {
    fit <- frm(fo, data = d, dry_run = "objective")
    tab <- default_prior(fo, data = d, route = "sample")
    # one density per parameter, keyed by the parameter, so the
    # comparison is exact and does not depend on the order the entries
    # are summed in (which differs between the two stacks by an ulp)
    per <- function(ri, p) {
      v <- vapply(ri$entries, function(e) neg_log_prior_fn(list(e))(p), 0)
      names(v) <- vapply(ri$entries, function(e) {
        paste0(e$comp, ".", paste(e$idx, collapse = ","))
      }, "")
      v[order(names(v))]
    }
    r0 <- dpa_ri(fit)
    r1 <- dpa_ri(fit, tab)
    sg <- dpa_ri(fit, set_prior(tab$prior[tab$class == "sigma"],
                                class = "sigma"))
    set.seed(6)
    for (k in 1:4) {
      p <- fit$frame$par_template
      p[["thetaac"]][] <- stats::rnorm(length(p[["thetaac"]]))
      p[["betad"]][] <- stats::rnorm(length(p[["betad"]]), 0, 0.5)
      p[["beta"]][] <- stats::rnorm(length(p[["beta"]]))
      d0 <- per(r0, p)
      expect_identical(per(r1, p), d0, label = deparse(fo))
      # the default's text typed: its sigma density is the default's
      e_sg <- Filter(function(e) identical(e$dist$kind, "t") &&
                       identical(e$scale, "sd"), sg$entries)
      expect_length(e_sg, 1L)
      sig <- paste0(e_sg[[1L]]$comp, ".", e_sg[[1L]]$idx)
      expect_identical(per(sg, p)[[sig]], d0[[sig]], label = deparse(fo))
    }
  }
})

test_that("frm()'s MAP under the sampling table is the sampler's mode", {
  # the same prior on both routes: frm(prior = the sampling table) is
  # the mode of the density frm_sample() samples by default
  set.seed(7)
  d <- dpa_data(seed = 7, ng = 20)
  fo <- y ~ x + ar(t, g, cov = TRUE)
  tab <- default_prior(fo, data = d, route = "sample")
  fmap <- frm(fo, data = d, prior = tab)
  fit <- frm(fo, data = d)
  rp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))
  fit$prior <- rp$effective
  obj <- frmtmb.sample:::prior_augmented_obj(fit, rp$ri$entries)
  at_map <- fmap$opt$par
  g0 <- max(abs(obj$gr(at_map)))
  g1 <- max(abs(obj$gr(at_map + 0.05)))
  # stationary at the MAP, against the gradient a small step away
  expect_lt(g0 / g1, 1e-3)
  # and the two objectives are the same function there
  expect_equal(obj$fn(at_map), fmap$obj$fn(at_map))
})

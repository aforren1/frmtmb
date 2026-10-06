# The draws surface on core's ordinal mixtures and on hurdle_cumulative()
# with thres(gr = ) and cs(). A draw stores each component's thresholds
# under brms's names, b_mu1_Intercept[k] and b_mu2_Intercept[k], both
# from the one shared vector under order = "mu"; these check the names,
# the map both ways, and posterior_epred() and log_lik() at a draw
# against brms's own R-side densities, theta-weighted as
# brms:::posterior_epred_mixture() weights them.

omd_data <- function(seed, n = 200) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < 0.2, 0L, d$y)
  d
}

# Stan's own sampler diagnostics: these short single chains raise them
# (measured per fit in dev/ordmix-p1-m3.R), and the tests read names
# and densities at the draws, not the posterior summaries. Any other
# warning escapes and fails.
omd_stan <- c("divergent transitions", "maximum treedepth",
              "Examine the pairs() plot", "R-hat", "Effective Samples Size")

# brms's default student_t on the thresholds, written by dpar as brms
# keys its rows: flat thresholds leave the sampler an improper
# direction, and frm_sample() sets no default there
omd_draws <- local({
  cache <- list()
  function(which) {
    skip_sampler()
    if (is.null(cache[[which]])) {
      d <- omd_data(71)
      st <- "student_t(3, 0, 2.5)"
      nb <- set_prior("normal(0, 2)", class = "b", dpar = "mu1") +
        set_prior("normal(0, 2)", class = "b", dpar = "mu2")
      # the "none" fit's component 1 ends at a degenerate boundary
      # (disc1's slope has no prior) and says so; the fit is only the
      # sampler's start, and the tests read the draws
      fit <- allow_warnings(switch(which,
        none = frm(bf(y ~ x, disc1 ~ 0 + z),
                   family = mixture(cumulative(), sratio()), data = d,
                   prior = nb +
                     set_prior(st, class = "Intercept", dpar = "mu1") +
                     set_prior(st, class = "Intercept", dpar = "mu2")),
        mu = frm(bf(y ~ x),
                 family = mixture(cumulative(),
                                  acat(threshold = "sum_to_zero"),
                                  order = "mu"),
                 data = d, prior = nb + set_prior(st, class = "Intercept")),
        hgr = frm(bf(yh | thres(gr = g) ~ x), family = hurdle_cumulative(),
                  data = d, prior = set_prior(st, class = "Intercept")),
        hcs = frm(bf(yh ~ cs(x)), family = hurdle_cumulative("probit"),
                  data = d,
                  prior = set_prior(st, class = "Intercept") +
                    set_prior("normal(0, 1)", class = "b"))),
        "degenerate boundary")
      ds <- allow_warnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3)),
        omd_stan)
      cache[[which]] <<- list(d = d, fit = fit, ds = ds)
    }
    cache[[which]]
  }
})

omd_probs <- function(family, eta, tau, disc = 1, link = "logit") {
  dens <- get(paste0("d", family), asNamespace("brms"))
  n <- length(eta)
  if (!is.matrix(tau)) tau <- matrix(tau, n, length(tau), byrow = TRUE)
  dens(seq_len(ncol(tau) + 1L), eta = eta, thres = tau, disc = disc,
       link = link)
}

test_that("each component's thresholds are drawn under brms's names", {
  cs <- omd_draws("none")
  v <- variables(cs$ds)
  expect_true(all(c(paste0("b_mu1_Intercept[", 1:3, "]"),
                    paste0("b_mu2_Intercept[", 1:3, "]"), "b_mu1_x",
                    "b_disc1_z", "theta1", "theta2") %in% v))
  # the internal matrix maps forward onto the stored columns again
  int <- frmtmb.sample:::draws_internal_matrix(cs$ds)
  for (o in frmtmb.sample:::draws_natural_cols(cs$ds$fit)$ordinal) {
    fwd <- t(apply(int[, o$names[seq_along(o$internal)], drop = FALSE],
                   1L, o$map))
    expect_lt(max(abs(fwd - cs$ds$draws[, o$names])) /
                max(abs(cs$ds$draws[, o$names])),
              1e3 * .Machine$double.eps)
  }
  # disc1's intercept is held at 0 here, so frm_sample() has no default
  # to add; the fit's own priors are the only ones
  expect_identical(nrow(fixef(cs$ds)), 9L)
})

test_that("a draw's category probabilities are the theta-weighted brms sum", {
  skip_if_not_installed("brms")
  for (which in c("none", "mu")) {
    cs <- omd_draws(which)
    M <- posterior::as_draws_matrix(cs$ds)
    e <- posterior_epred(cs$ds)
    ll <- log_lik(cs$ds)
    fams <- if (which == "none") c("cumulative", "sratio") else
      c("cumulative", "acat")
    for (s in c(1L, nrow(M))) {
      P <- 0
      for (k in 1:2) {
        th <- as.numeric(M[s, paste0("b_mu", k, "_Intercept[", 1:3, "]")])
        eta <- as.numeric(M[s, paste0("b_mu", k, "_x")]) * cs$d$x
        disc <- if (which == "none" && k == 1L) {
          exp(as.numeric(M[s, "b_disc1_z"]) * cs$d$z)
        } else 1
        P <- P + as.numeric(M[s, paste0("theta", k)]) *
          omd_probs(fams[k], eta, th, disc)
      }
      expect_lt(max(abs(e[s, , ] - P)), 1e3 * .Machine$double.eps,
                label = paste(which, s))
      ref <- log(P[cbind(seq_len(nrow(cs$d)), cs$d$y)])
      expect_lt(max(abs(ll[s, ] - ref) / abs(ref)),
                1e3 * .Machine$double.eps, label = paste(which, s))
    }
  }
  # order = "mu": the sum-to-zero component reads the shared vector
  # centered, as brms's Intercept_mu2_stz
  cs <- omd_draws("mu")
  M <- posterior::as_draws_matrix(cs$ds)
  t1 <- M[, paste0("b_mu1_Intercept[", 1:3, "]")]
  t2 <- M[, paste0("b_mu2_Intercept[", 1:3, "]")]
  expect_lt(max(abs(t2 - (t1 - rowMeans(t1)))) / max(abs(t1)),
            1e3 * .Machine$double.eps)
})

test_that("hurdle_cumulative() draws under thres(gr = ) and cs()", {
  skip_if_not_installed("brms")
  cs <- omd_draws("hgr")
  v <- variables(cs$ds)
  expect_true(all(c("b_Intercept[a,1]", "b_Intercept[b,3]", "hu") %in% v))
  M <- posterior::as_draws_matrix(cs$ds)
  e <- posterior_epred(cs$ds)
  gi <- as.integer(cs$d$g)
  for (s in c(1L, nrow(M))) {
    eta <- as.numeric(M[s, "b_x"]) * cs$d$x
    hu <- as.numeric(M[s, "hu"])
    P <- matrix(0, nrow(cs$d), 5)
    for (g in 1:2) {
      lv <- c("a", "b")[g]
      th <- as.numeric(M[s, paste0("b_Intercept[", lv, ",", 1:3, "]")])
      r <- gi == g
      P[r, ] <- cbind(hu, (1 - hu) * omd_probs("cumulative", eta[r], th))
    }
    expect_lt(max(abs(e[s, , ] - P)), 1e3 * .Machine$double.eps)
  }
  cs <- omd_draws("hcs")
  expect_true(all(c("bcs_x[1]", "bcs_x[3]") %in% variables(cs$ds)))
  M <- posterior::as_draws_matrix(cs$ds)
  e <- posterior_epred(cs$ds)
  for (s in c(1L, nrow(M))) {
    th <- as.numeric(M[s, paste0("b_Intercept[", 1:3, "]")])
    bcs <- as.numeric(M[s, paste0("bcs_x[", 1:3, "]")])
    tau <- matrix(th, nrow(cs$d), 3, byrow = TRUE) - outer(cs$d$x, bcs)
    hu <- as.numeric(M[s, "hu"])
    P <- cbind(hu, (1 - hu) * omd_probs("cumulative", rep(0, nrow(cs$d)),
                                        tau, link = "probit"))
    # a draw whose offsets cross two thresholds on a row has a negative
    # "probability" there in brms, and NaN here
    ok <- rowSums(P < 0) == 0
    expect_gt(sum(ok), 0)
    expect_lt(max(abs(e[s, ok, ] - P[ok, ])), 1e3 * .Machine$double.eps)
    expect_true(all(is.nan(e[s, !ok, ])))
  }
})

test_that("a model with no location column keeps its other names", {
  # On 0.67.0 a location predictor with no column (y ~ 1 on an ordinal
  # family, or y ~ cs(x)) dropped the label of every betad entry, so the
  # draws were named one column off: disc's coefficient was stored as
  # b_Intercept[1] and the last threshold as tau_raw[3]
  # (dev/ordmix-emptybeta.R)
  skip_sampler()
  set.seed(20261005)
  d <- data.frame(z = rnorm(200))
  u <- stats::rlogis(200) / exp(0.4 * d$z)
  d$y <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
  fit <- frm(bf(y ~ 1, disc ~ 0 + z), family = cumulative(), data = d,
             prior = set_prior("student_t(3, 0, 2.5)", class = "Intercept"))
  ds <- allow_warnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3)),
    omd_stan)
  expect_setequal(setdiff(variables(ds), "lp__"), variables(fit))
  M <- posterior::as_draws_matrix(ds)
  th <- M[, paste0("b_Intercept[", 1:3, "]")]
  # thresholds, so ordered in every draw; the disc slope sits near its
  # ML estimate, in units of its own posterior spread
  expect_true(all(th[, 2] > th[, 1] & th[, 3] > th[, 2]))
  dz <- as.numeric(M[, "b_disc_z"])
  expect_lt(abs(mean(dz) - fixef(fit)["disc_z", "Estimate"]) / stats::sd(dz),
            1)
  # and hurdle_cumulative() with cs() alone keeps its hu
  expect_true("hu" %in% variables(omd_draws("hcs")$ds))
})

test_that("frm_sample() gives each modeled disc<k> brms's default", {
  d <- omd_data(72)
  uf <- frm(bf(y ~ x, disc2 ~ z), family = mixture(cumulative(),
                                                    cumulative()),
            data = d, dry_run = "objective")
  defs <- unclass(frmtmb.sample:::default_priors_for(uf))
  s <- Filter(function(s) identical(s$class, "Intercept") &&
                identical(s$dpar, "disc2"), defs)
  expect_length(s, 1L)
  expect_identical(c(s[[1L]]$dist$mean %||% s[[1L]]$dist$location,
                     s[[1L]]$dist$sd %||% s[[1L]]$dist$scale), c(0, 1))
  # absent: a disc held at 1 has no intercept to place
  expect_false(any(vapply(defs, function(s) identical(s$dpar, "disc1"), NA)))
})

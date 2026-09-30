# The draws surface on core's ordinal discrimination `disc` and its
# threshold structures. A draw stores the thresholds under brms's names
# and on brms's scale, `b_Intercept[k]`, plus brms's `delta` under
# equidistant thresholds, while the sampler holds fewer internal
# parameters than that; these check the names, the map both ways, and
# the density at a draw against brms's own R-side densities.

ordd_data <- function(seed, n = 150) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
  d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
  d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
  d
}

ordd_draws <- local({
  cache <- list()
  function(which) {
    skip_sampler()
    if (is.null(cache[[which]])) {
      d <- ordd_data(62)
      fit <- switch(which,
        equi = frm(bf(y ~ x, disc ~ 0 + z),
                   family = cumulative(threshold = "equidistant"), data = d),
        stz = frm(y ~ x, family = sratio(threshold = "sum_to_zero"),
                  data = d),
        hurdle = frm(yh ~ x, family = hurdle_cumulative(), data = d),
        # no prior of its own, so the sampler's default fills the slot;
        # the ML fit warns that nothing places the intercept
        hdisc = allow_warnings(
          frm(bf(yh ~ x, disc ~ z), family = hurdle_cumulative(),
              data = d),
          "disc has an intercept"))
      msg <- character()
      ds <- withCallingHandlers(
        suppressWarnings(frm_sample(fit, chains = 1, iter = 300,
                                    refresh = 0, seed = 3)),
        message = function(m) {
          msg <<- c(msg, conditionMessage(m))
          invokeRestart("muffleMessage")
        })
      cache[[which]] <<- list(d = d, fit = fit, ds = ds, msg = msg)
    }
    cache[[which]]
  }
})

test_that("equidistant draws carry every threshold and delta", {
  cs <- ordd_draws("equi")
  v <- variables(cs$ds)
  thr <- paste0("b_Intercept[", 1:4, "]")
  expect_true(all(c(thr, "delta") %in% v))
  m <- as.matrix(cs$ds, variable = c(thr, "delta"))
  # the thresholds of every draw are equally spaced by that draw's delta
  sp <- m[, thr] - m[, thr[1L]] -
    outer(m[, "delta"], 0:3)
  expect_lt(max(abs(sp)) / max(abs(m)), 1e3 * .Machine$double.eps)
  expect_true(all(m[, "delta"] > 0))
  # the internal matrix maps forward onto the stored columns again
  o <- frmtmb.sample:::draws_natural_cols(cs$ds$fit)$ordinal[[1L]]
  R <- length(o$internal)
  expect_identical(R, 2L)
  int <- frmtmb.sample:::draws_internal_matrix(cs$ds)
  expect_false(any(o$names[-seq_len(R)] %in% colnames(int)))
  fwd <- t(apply(int[, o$names[seq_len(R)], drop = FALSE], 1L, o$map))
  expect_lt(max(abs(fwd - cs$ds$draws[, o$names])) /
              max(abs(cs$ds$draws[, o$names])), 1e3 * .Machine$double.eps)
  expect_identical(nrow(fixef(cs$ds)), 6L)
})

test_that("a draw's density is brms's, disc included", {
  skip_if_not_installed("brms")
  cs <- ordd_draws("equi")
  ll <- log_lik(cs$ds)
  idx <- frmtmb.sample:::draws_par_index(cs$ds$fit)
  for (s in c(1L, nrow(ll))) {
    f <- frmtmb.sample:::draws_fit_at(cs$ds, s, idx)
    dp <- frmtmb:::eval_dpars(f)[[1L]]
    tau <- frmtmb:::ord_threshold_values(family(f), f$estimates$tau_raw)
    n <- nrow(cs$d)
    P <- brms:::dcumulative(1:5, eta = dp$mu,
                            thres = matrix(tau, n, 4L, byrow = TRUE),
                            disc = dp$disc, link = "logit")
    ref <- log(P[cbind(seq_len(n), cs$d$y)])
    expect_lt(max(abs(ll[s, ] - ref) / pmax(1, abs(ref))),
              1e3 * .Machine$double.eps)
  }
})

test_that("sum-to-zero draws sum to zero under brms's names", {
  cs <- ordd_draws("stz")
  thr <- paste0("b_Intercept[", 1:4, "]")
  expect_true(all(thr %in% variables(cs$ds)))
  m <- as.matrix(cs$ds, variable = thr)
  expect_lt(max(abs(rowSums(m))) / max(abs(m)), 1e3 * .Machine$double.eps)
  ep <- posterior_epred(cs$ds, ndraws = 5)
  expect_identical(dim(ep), c(5L, nrow(cs$d), 5L))
})

test_that("a hurdle ordinal's thresholds are stored under brms's names", {
  cs <- ordd_draws("hurdle")
  v <- variables(cs$ds)
  expect_true(all(paste0("b_Intercept[", 1:4, "]") %in% v))
  expect_false(any(grepl("^tau_raw", v)))
})

test_that("an intercept in disc gets brms's normal(0, 1) by default", {
  cs <- ordd_draws("hdisc")
  ps <- utils::capture.output(print(prior_summary(cs$ds)))
  expect_true(any(grepl("Intercept_disc ~ normal(0, 1)", ps, fixed = TRUE)))
})

# A draw's category probabilities against brms's R-side density at the
# STORED columns, so a wrong inverse map (the columns handed back to the
# model) or a wrong forward map (the columns written) both show.
ordd_brms_P <- function(fam, eta, thres, disc = 1, link = "logit") {
  n <- length(eta)
  get(paste0("d", fam), asNamespace("brms"))(
    seq_len(length(thres) + 1L), eta = eta,
    thres = matrix(thres, n, length(thres), byrow = TRUE), disc = disc,
    link = link)
}

ordd_mv_data <- function(seed, n = 150) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  cut_at <- function(br) {
    u <- stats::rlogis(n) + 0.8 * d$x
    1L + rowSums(outer(u, br, ">"))
  }
  d$y <- cut_at(c(-1.2, -0.2, 0.8, 1.8))   # 5 categories
  d$y2 <- cut_at(c(-0.5, 1))              # 3 categories
  d
}

ordd_mv_draws <- local({
  cache <- list()
  function(which) {
    skip_sampler()
    if (is.null(cache[[which]])) {
      d <- ordd_mv_data(63)
      fam <- switch(which,
        flex = list(cumulative(), sratio()),
        equi = list(cumulative(threshold = "equidistant"),
                    sratio(threshold = "equidistant")),
        mixed = list(cumulative(threshold = "equidistant"),
                     acat(threshold = "sum_to_zero")))
      fit <- frm(bf(y ~ x) + bf(y2 ~ x), data = d, family = fam)
      ds <- withCallingHandlers(
        suppressWarnings(frm_sample(fit, chains = 1, iter = 300,
                                    refresh = 0, seed = 3)),
        message = function(m) invokeRestart("muffleMessage"))
      cache[[which]] <<- list(d = d, fit = fit, ds = ds,
                              fams = vapply(fam, `[[`, "", "family"))
    }
    cache[[which]]
  }
})

test_that("two ordinal responses each keep their own threshold maps", {
  skip_if_not_installed("brms")
  for (which in c("flex", "equi", "mixed")) {
    cs <- ordd_mv_draws(which)
    M <- as.matrix(cs$ds)
    for (r in 1:2) {
      rn <- c("y", "y2")[r]
      K1 <- c(4L, 2L)[r]
      thr <- paste0("b_", rn, "_Intercept[", seq_len(K1), "]")
      expect_true(all(thr %in% colnames(M)), label = paste(which, rn))
      ep <- posterior_epred(cs$ds, resp = rn)
      expect_identical(dim(ep), c(nrow(M), nrow(cs$d), K1 + 1L))
      for (s in c(1L, nrow(M))) {
        P <- ordd_brms_P(cs$fams[r], cs$d$x * M[s, paste0("b_", rn, "_x")],
                         M[s, thr])
        expect_lt(max(abs(ep[s, , ] - P)), 1e3 * .Machine$double.eps,
                  label = paste(which, rn, s))
      }
      dl <- paste0("delta_", rn)
      if (dl %in% colnames(M)) {
        # the stored delta is the stored spacing, not its internal value
        sp <- M[, thr[2L]] - M[, thr[1L]]
        expect_lt(max(abs(M[, dl] - sp)) / max(abs(sp)),
                  1e3 * .Machine$double.eps, label = paste(which, dl))
      }
    }
  }
})

test_that("sum-to-zero draws hand the right thresholds back to the model", {
  skip_if_not_installed("brms")
  # sratio holds the first K - 2 thresholds and the last is minus their
  # sum; posterior_epred() reads the internal vector the inverse map
  # rebuilds from the stored columns, so a wrong inverse shows here
  cs <- ordd_draws("stz")
  M <- as.matrix(cs$ds)
  thr <- paste0("b_Intercept[", 1:4, "]")
  ep <- posterior_epred(cs$ds)
  for (s in c(1L, nrow(M))) {
    P <- ordd_brms_P("sratio", cs$d$x * M[s, "b_x"], M[s, thr])
    expect_lt(max(abs(ep[s, , ] - P)), 1e3 * .Machine$double.eps)
  }
})


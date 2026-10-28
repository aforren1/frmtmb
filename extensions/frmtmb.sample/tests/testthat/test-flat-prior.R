# brms's (flat) is flat on brms's own parameter, inside brms's bounds.
# brms 2.23.0 has no prior = "flat"; its nearest model is "" in every
# slot, and its program then carries no lprior term at all: Stan's
# implicit uniform on `vector<lower=0> sd_1`, `real<lower=0> sigma`,
# `cholesky_factor_corr L_1`, `vector<lower=-1,upper=1> ar`,
# `real<lower=0> delta`, `ordered Intercept` (brms 2.23.0's stancode()).
# frmtmb holds these through maps, so flat on
# brms's parameter is the map's Jacobian on the internal one. Through
# frmtmb.sample 0.23.0 prior = "flat" added nothing, which is flat on
# the log standard deviation and (1 - rho^2)^(-3/2) on a first-order
# correlation. Against brms on y ~ 1 + (1 | g), 6 groups of 5
# (data seed 11, 4 chains of 3000 draws): sd_g__Intercept
# mean 1.458 in brms, 1.455 now, 1.197 before.

fp_data <- function(seed = 5, G = 8, n = 6) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(G), each = n)),
                  t = rep(seq_len(n), G), x = stats::rnorm(G * n),
                  z = stats::rnorm(G * n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(G, 0, 0.5)[d$g] + stats::rnorm(G * n)
  d$ord <- cut(d$y, 4, labels = FALSE)
  d
}

fp_rp <- function(fit, prior) {
  suppressWarnings(suppressMessages(
    frmtmb.sample:::sample_resolve_priors(fit, prior)))
}

# -log prior of the resolved entries at the fit's template with `set`
# positions replaced
fp_nlp <- function(fit, entries, set) {
  p <- fit$frame$par_template
  for (cmp in names(set)) {
    v <- set[[cmp]]
    p[[cmp]][seq_along(v)] <- v
  }
  neg_log_prior_fn(entries)(p)
}

test_that("prior = \"flat\" is flat on the sd and on sigma themselves", {
  d <- fp_data()
  fit <- frm(y ~ x + (1 | g), data = d, dry_run = "objective")
  rp <- fp_rp(fit, "flat")
  # theta_1 = log sd and betad_1 = log sigma: flat on sd and sigma is
  # their Jacobian, so -log prior + theta_1 + betad_1 is 0 everywhere
  set.seed(1)
  pts <- replicate(6, stats::rnorm(2, 0, 2), simplify = FALSE)
  res <- vapply(pts, function(v) {
    fp_nlp(fit, rp$ri$entries, list(theta = v[1], betad = v[2])) +
      v[1] + v[2]
  }, 0)
  scale <- max(abs(unlist(pts)))
  expect_lt(max(abs(res)), 1e-12 * scale)
  # written as brms writes the rows: (flat)
  pe <- as.data.frame(rp$effective)
  expect_setequal(pe$class, c("Intercept", "sigma", "sd"))
  expect_true(all(pe$prior == ""))
})

test_that("prior = \"flat\" on a correlation matrix is flat on its factor", {
  # brms declares cholesky_factor_corr L_1, and "" leaves Stan's implicit
  # uniform over L's free entries; lkj(1) is that only at d = 2. The
  # reference here is independent of the prior code: the log
  # determinant of the map from the three correlation parameters to
  # L's free entries, by central differences through us_sigma()
  d <- fp_data()
  fit <- frm(y ~ x + (1 + x + z | g), data = d, dry_run = "objective")
  rp <- fp_rp(fit, "flat")
  cor_e <- Filter(function(e) identical(e$dist$kind, "lkj"),
                  rp$ri$entries)
  expect_length(cor_e, 1L)
  expect_identical(cor_e[[1L]]$idx, 4:6)
  freeL <- function(t) {
    L <- t(chol(stats::cov2cor(as.matrix(us_sigma(c(0, 0, 0, t), 3)))))
    L[lower.tri(L)]
  }
  logdet <- function(t) {
    h <- 1e-5
    J <- vapply(1:3, function(k) {
      e <- replace(numeric(3), k, h)
      (freeL(t + e) - freeL(t - e)) / (2 * h)
    }, numeric(3))
    log(abs(det(J)))
  }
  set.seed(2)
  pts <- replicate(6, stats::rnorm(3, 0, 0.8), simplify = FALSE)
  gap <- vapply(pts, function(t) {
    -fp_nlp(fit, cor_e, list(theta = c(0, 0, 0, t))) - logdet(t)
  }, 0)
  lj <- vapply(pts, logdet, 0)
  expect_lt(diff(range(gap)), 1e-6 * max(abs(lj)))
  # and it is not lkj(1), whose density differs from it by a factor
  # in L's second diagonal entry at d = 3
  lkj1 <- fp_rp(fit, NULL)$ri$entries
  lkj1 <- Filter(function(e) identical(e$dist$kind, "lkj"), lkj1)
  gap1 <- vapply(pts, function(t) {
    -fp_nlp(fit, lkj1, list(theta = c(0, 0, 0, t))) - logdet(t)
  }, 0)
  expect_gt(diff(range(gap1)), 1e3 * diff(range(gap)))
})

test_that("prior = \"flat\" keeps the flat defaults on an interval", {
  d <- fp_data()
  fit <- frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = "objective")
  flat <- fp_rp(fit, "flat")$ri$entries
  dflt <- fp_rp(fit, NULL)$ri$entries
  ac <- function(es) Filter(function(e) identical(e$comp, "thetaac"), es)
  expect_identical(ac(flat), ac(dflt))
  # brms's sigma under a first-order cov = TRUE ar is the innovation sd,
  # and (flat) on it is read there, as every density on class sigma is
  sg <- Filter(function(e) identical(e$comp, "betad"), flat)
  expect_length(sg, 1L)
  expect_identical(sg[[1L]]$dist$kind, "flat")
  expect_false(is.null(sg[[1L]]$innov))
})

test_that("brms's (flat) on delta is flat on delta, by default too", {
  # frmtmb holds log(delta) for an ordered family; brms declares
  # real<lower=0> delta and leaves it (flat) by default
  d <- fp_data()
  fit <- frm(bf(ord ~ x), data = d,
             family = cumulative(threshold = "equidistant"),
             dry_run = "objective")
  for (pr in list(NULL, "flat")) {
    es <- fp_rp(fit, pr)$ri$entries
    de <- Filter(function(e) identical(e$dist$kind, "flat") &&
                   identical(e$scale, "sd"), es)
    expect_length(de, 1L)
    v <- c(-2, 0.3, 1.7)
    res <- vapply(v, function(ld) {
      fp_nlp(fit, de, list(tau_raw = c(0, ld))) + ld
    }, 0)
    expect_lt(max(abs(res)), 1e-12 * max(abs(v)))
  }
  tab <- default_prior(bf(ord ~ x), data = d,
                       family = cumulative(threshold = "equidistant"),
                       route = "sample")
  expect_identical(tab$prior[tab$class == "delta"], "(flat)")
  expect_identical(tab$lb[tab$class == "delta"], 0)
})

test_that("a bounds-only specification is brms's (flat) on the slot", {
  d <- fp_data()
  fit <- frm(y ~ x + (1 | g), data = d, dry_run = "objective")
  ri <- fp_rp(fit, set_prior("", class = "sd", lb = 0.1))$ri
  sd_e <- Filter(function(e) identical(e$comp, "theta"), ri$entries)
  expect_length(sd_e, 1L)
  expect_identical(sd_e[[1L]]$dist$kind, "flat")
  expect_identical(sd_e[[1L]]$scale, "sd")
  expect_equal(unname(ri$lower[["theta_1"]]), log(0.1))
})

test_that("a flat run non-centers, its sd tail being integrable", {
  # flat on the sd is sd in theta's coordinates, which vanishes at
  # theta = -Inf; flat on theta did not, which kept such a run centered
  d <- fp_data()
  fit <- frm(y ~ x + (1 | g), data = d, dry_run = "objective")
  rp <- fp_rp(fit, "flat")
  plan <- frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, rp$ri$entries)
  expect_identical(plan$idx, 1L)
})

test_that("class theta on a random-effect parameter is refused here too", {
  # through 0.76.0 one theta_4 prior retired the default lkj(1) on all
  # three correlations of the block
  d <- fp_data()
  fit <- frm(y ~ x + (1 + x + z | g), data = d, dry_run = "objective")
  expect_error(fp_rp(fit, set_prior("normal(0, 1)", class = "theta",
                                    coef = "theta_4")),
               "has a class of its own", class = "frmtmb_error")
})

test_that("the sampling table's ar row passes back as the same density", {
  d <- fp_data()
  fo <- y ~ x + ar(t, g, cov = TRUE)
  tab <- default_prior(fo, data = d, route = "sample")
  row <- tab[tab$class == "ar", ]
  expect_identical(row$prior, "(flat)")
  expect_identical(c(row$lb, row$ub), c(-1, 1))
  fit <- frm(fo, data = d)
  r0 <- fp_rp(fit, NULL)$ri
  r1 <- fp_rp(fit, tab)$ri
  o0 <- frmtmb.sample:::prior_augmented_obj(fit, r0$entries)
  o1 <- frmtmb.sample:::prior_augmented_obj(fit, r1$entries)
  set.seed(3)
  for (j in 1:25) {
    p <- o0$par + stats::rnorm(length(o0$par), 0, 0.5)
    f0 <- o0$fn(p)
    expect_lt(abs(o1$fn(p) - f0), 1e-12 * abs(f0))
    g0 <- o0$gr(p)
    expect_lt(max(abs(o1$gr(p) - g0)), 1e-12 * max(abs(g0)))
  }
})

# the round trip of a flat run, through the objective at many points,
# never through draws: tmbstan holds the latest fit's objective in a
# process-wide pointer
fp_same_objective <- function(fit, e0, e1, n = 40) {
  o0 <- frmtmb.sample:::prior_augmented_obj(fit, e0)
  o1 <- frmtmb.sample:::prior_augmented_obj(fit, e1)
  set.seed(4)
  worst <- c(fn = 0, gr = 0)
  for (j in seq_len(n)) {
    p <- o0$par + stats::rnorm(length(o0$par), 0, 0.5)
    f0 <- o0$fn(p)
    g0 <- o0$gr(p)
    worst["fn"] <- max(worst["fn"], abs(o1$fn(p) - f0) / abs(f0))
    worst["gr"] <- max(worst["gr"],
                       max(abs(o1$gr(p) - g0)) / max(abs(g0)))
  }
  worst
}

test_that("a flat run's prior_summary() passes back as the same density", {
  # the intercept slots are in the summary too, as brms's prior_summary()
  # of a model fitted with "" lists them, so the defaults do not return
  d <- fp_data()
  fit <- allow_warnings(frm(y ~ x + (1 + x + z | g), data = d),
                        "Standard errors are not available")
  rp <- fp_rp(fit, "flat")
  pe <- as.data.frame(rp$effective)
  expect_setequal(pe$class, c("Intercept", "sigma", "sd", "cor"))
  back <- fp_rp(fit, rp$effective)
  w <- fp_same_objective(fit, rp$ri$entries, back$ri$entries)
  expect_lt(w[["fn"]], 1e-12)
  expect_lt(w[["gr"]], 1e-12)
  # and so does its table, read back with as.brmsprior()
  tb <- fp_rp(fit, as.brmsprior(pe))
  w <- fp_same_objective(fit, rp$ri$entries, tb$ri$entries)
  expect_lt(w[["fn"]], 1e-12)
  expect_lt(w[["gr"]], 1e-12)
})

test_that("\"\" on one slot is brms's (flat) there, an LKJ's too", {
  d <- fp_data()
  fit <- frm(y ~ x + (1 + x + z | g), data = d, dry_run = "objective")
  ri <- fp_rp(fit, set_prior("", class = "cor"))$ri
  cor_e <- Filter(function(e) identical(e$dist$kind, "lkj"), ri$entries)
  expect_length(cor_e, 1L)
  expect_true(isTRUE(cor_e[[1L]]$dist$flat))
  # the same entry prior = "flat" puts there
  fl <- Filter(function(e) identical(e$dist$kind, "lkj"),
               fp_rp(fit, "flat")$ri$entries)
  expect_identical(cor_e, fl)
  # the other defaults stay
  sd_e <- Filter(function(e) identical(e$comp, "theta") &&
                   identical(e$scale, "sd"), ri$entries)
  expect_true(all(vapply(sd_e, function(e) e$dist$kind, "") == "t"))
  # and "" on the intercept replaces its student_t by nothing
  ii <- fp_rp(fit, set_prior("", class = "Intercept"))$ri$entries
  bi <- Filter(function(e) identical(e$comp, "beta"), ii)
  expect_identical(vapply(bi, function(e) e$dist$kind, ""), "flat")
})

test_that("frm()'s prior table is refused, however it was edited", {
  # default_prior()'s default route is frm()'s, whose (flat) rows are
  # frm()'s defaults; read as brms's (flat) they replaced brms's
  # defaults on the Intercept, sigma, sd and cor in silence. brms keeps
  # its defaults on the rows of its own table nobody edited
  d <- fp_data()
  f <- y ~ x + (1 + x | g)
  fit <- frm(f, data = d, dry_run = "objective")
  p <- default_prior(f, data = d)
  p$prior[p$class == "b" & p$coef == "x"] <- "normal(0, 1)"
  msg <- "default_prior() or validate_prior() table of route = \"fit\""
  sp <- set_prior("normal(0, 1)", class = "b")
  for (pr in list(p, p[p$class != "cor", ], as.brmsprior(p),
                  sp + as.brmsprior(p), c(sp, as.brmsprior(p)),
                  validate_prior(sp, f, data = d))) {
    expect_error(fp_rp(fit, pr), msg, fixed = TRUE)
  }
  # frm_sample() itself, before anything is sampled
  expect_error(suppressMessages(frm_sample(fit, prior = p)), msg,
               fixed = TRUE)
  # the sample route's table is brms's defaults, edited in place
  ps <- default_prior(f, data = d, route = "sample")
  ps$prior[ps$class == "b" & ps$coef == "x"] <- "normal(0, 1)"
  e1 <- fp_rp(fit, ps)$ri$entries
  e0 <- fp_rp(fit, set_prior("normal(0, 1)", class = "b",
                             coef = "x"))$ri$entries
  key <- function(es) sort(vapply(es, function(e) {
    paste(e$comp, paste(e$idx, collapse = ","), e$dist$kind)
  }, ""))
  expect_identical(key(e1), key(e0))
  # a table with no route (built by hand) is read as written: its
  # (flat) rows are brms's (flat)
  hand <- as.brmsprior(as.data.frame(p))
  expect_null(attr(hand, "route"))
  eh <- fp_rp(fit, hand)$ri$entries
  expect_true("theta 1 flat" %in% key(eh))
})

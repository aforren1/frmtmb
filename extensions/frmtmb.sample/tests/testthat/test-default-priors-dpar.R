# frm_sample()'s defaults on every distributional parameter, against
# brms 2.23.0's default_prior() run on the same data in this test. The
# full table, family by family, is dev/sampleprior-table.R and the
# findings' section 1 (dev/sampleprior-findings.md). Through
# frmtmb.sample 0.19.0 every dpar but sigma was flat (shape, phi, nu,
# zi, hu, ...), as were the ordinal thresholds, and `sigma ~ 1` lost
# sigma's default.

# The rows brms gives a density, each compared with the row frmtmb's
# sampling route lists for the same slot. A group or coefficient row
# brms leaves empty under a class default is left out, as are the
# per-coefficient sd rows. brms's symbol min_Y is the smallest response,
# and a density is compared by its numbers, so normal(5.0, 0.8) is
# normal(5, 0.8).
dpd_mismatch <- function(form, data, fam_b, fam_f, rename = NULL,
                         miny = NA_real_) {
  rows <- function(tab) {
    tab <- as.data.frame(tab)
    pr <- ifelse(is.na(tab$prior) | !nzchar(tab$prior), "(flat)",
                 tab$prior)
    stats::setNames(pr, paste(tab$class, tab$coef, tab$group, tab$resp,
                              tab$dpar, sep = "|"))
  }
  fb <- lapply(form, stats::as.formula)
  b <- rows(suppressWarnings(suppressMessages(brms::default_prior(
    do.call(brms::bf, fb), data = data, family = fam_b))))
  ff <- if (is.null(rename)) form else {
    vapply(form, function(s) sub(rename[1], rename[2], s), "",
           USE.NAMES = FALSE)
  }
  f <- rows(suppressWarnings(suppressMessages(default_prior(
    do.call(bf, lapply(ff, stats::as.formula)), data = data,
    family = fam_f, route = "sample"))))
  if (!is.null(rename)) {
    names(f) <- vapply(strsplit(names(f), "|", fixed = TRUE), function(k) {
      k <- c(k, rep("", 5L - length(k)))
      k[k == rename[2]] <- rename[1]
      paste(k, collapse = "|")
    }, "")
  }
  num <- function(s) {
    s <- gsub("min_Y", sprintf("%.17g", miny), s, fixed = TRUE)
    m <- regmatches(s, regexec("^([a-z_]+)[(](.*)[)]$", s))[[1]]
    if (length(m) != 3L) return(list(s))
    list(m[2], as.numeric(strsplit(m[3], ",")[[1]]))
  }
  same <- function(x, y) {
    a <- num(x)
    b <- num(y)
    length(a) == length(b) && identical(a[[1]], b[[1]]) &&
      (length(a) == 1L || isTRUE(all.equal(a[[2]], b[[2]],
                                           tolerance = 1e-14)))
  }
  k <- names(b)[b != "(flat)"]
  k <- k[!grepl("^sd[|][^|]+[|]", k)]
  bad <- k[vapply(k, function(z) {
    is.na(f[z]) || !same(b[[z]], f[[z]])
  }, NA)]
  if (!length(bad)) return(character(0))
  paste0(bad, ": brms ", b[bad], ", frmtmb ", f[bad])
}

dpd_data <- function() {
  set.seed(20261014)
  n <- 60
  x <- stats::rnorm(n)
  d <- data.frame(x = x, g = factor(rep(1:6, length.out = n)))
  mu <- stats::plogis(0.2 + 0.5 * x)
  list(
    cnt = transform(d, y = stats::rpois(n, exp(0.7 + 0.3 * x))),
    zcnt = transform(d, y = ifelse(stats::runif(n) < 0.2, 0,
                                   stats::rpois(n, exp(0.7 + 0.3 * x)))),
    pos = transform(d, y = stats::rgamma(n, 4, 4 / exp(0.5 + 0.3 * x))),
    unit = transform(d, y = stats::rbeta(n, 10 * mu, 10 * (1 - mu))),
    real = transform(d, y = 1 + x + stats::rt(n, 5)),
    ord = transform(d, y = cut(0.5 * x + stats::rlogis(n),
                               c(-Inf, -1, 0, 1, Inf), labels = FALSE)),
    circ = transform(d, y = atan2(sin(0.3 * x + stats::rnorm(n, 0, 0.5)),
                                  cos(0.3 * x + stats::rnorm(n, 0, 0.5)))),
    shift = transform(d, y = 0.2 + exp(stats::rnorm(n, 0.3 * x))))
}

test_that("every dispersion and shape dpar takes brms's default", {
  skip_unless_brms()
  D <- dpd_data()
  cases <- list(
    list("y ~ x + (1 | g)", "negbinomial", D$cnt),
    list(c("y ~ x + (1 | g)", "shape ~ x"), "negbinomial", D$cnt),
    list("y ~ x + (1 | g)", "Gamma", D$pos),
    list(c("y ~ x", "shape ~ 1"), "weibull", D$pos),
    list("y ~ x", "beta", D$unit),
    list(c("y ~ x", "phi ~ 1"), "beta", D$unit),
    list("y ~ x", "student", D$real),
    list(c("y ~ x", "nu ~ x"), "student", D$real),
    list("y ~ x", "skew_normal", D$real),
    list(c("y ~ x", "alpha ~ x"), "skew_normal", D$real),
    list("y ~ x", "exgaussian", D$real),
    list(c("y ~ x", "beta ~ x"), "exgaussian", D$real),
    list("y ~ x", "asym_laplace", D$real),
    list("y ~ x", "frechet", D$pos),
    list("y ~ x", "xbeta", D$unit),
    list(c("y ~ x", "kappa ~ x"), "xbeta", D$unit),
    list("y ~ x", "von_mises", D$circ),
    list(c("y ~ x", "kappa ~ 1"), "von_mises", D$circ),
    list("y ~ x", "zero_inflated_poisson", D$zcnt),
    list(c("y ~ x", "zi ~ x"), "zero_inflated_negbinomial", D$zcnt),
    list("y ~ x", "hurdle_gamma", transform(D$pos, y = ifelse(x > 1, 0, y))),
    list("y ~ x", "zero_one_inflated_beta",
         transform(D$unit, y = ifelse(x > 1.5, 1, ifelse(x < -1.5, 0, y)))),
    list("y ~ x", "discrete_weibull", D$cnt))
  for (cs in cases) {
    expect_identical(
      dpd_mismatch(cs[[1]], cs[[3]], brms::brmsfamily(cs[[2]]),
                   frmtmb::brmsfamily(cs[[2]])),
      character(0), label = paste(cs[[2]], paste(cs[[1]], collapse = ", ")))
  }
  # com_poisson's dispersion is frmtmb's nu and brms's shape
  for (f in list("y ~ x", c("y ~ x", "shape ~ x"))) {
    expect_identical(
      dpd_mismatch(f, D$cnt, brms::brmsfamily("com_poisson"),
                   frmtmb::brmsfamily("com_poisson"),
                   rename = c("shape", "nu")),
      character(0))
  }
})

test_that("sigma ~ 1 puts the default on sigma's intercept, as brms does", {
  skip_unless_brms()
  # through 0.19.0 this was the natural-scale prior of a sigma with no
  # formula, which set_prior() reads on a slot this model does not have
  D <- dpd_data()
  for (fam in c("gaussian", "student", "lognormal")) {
    d <- if (identical(fam, "lognormal")) D$pos else D$real
    expect_identical(
      dpd_mismatch(c("y ~ x", "sigma ~ 1"), d, brms::brmsfamily(fam),
                   frmtmb::brmsfamily(fam)), character(0), label = fam)
  }
  # and sigma's scale is the response's on sigma's own scale: the
  # identity without a formula (here 47.5 where log(y) gives 2.5), the
  # log of a wide positive response with one
  set.seed(20261018)
  w <- data.frame(x = stats::rnorm(60))
  w$y <- 50 * stats::rgamma(60, 4, 4 / exp(0.5 + 0.3 * w$x))
  expect_identical(
    dpd_mismatch("y ~ x", w, brms::brmsfamily("gaussian", link = "log"),
                 frmtmb::brmsfamily("gaussian", link = "log")),
    character(0))
  w$y <- exp(stats::rnorm(60, 0, 5))
  expect_identical(
    dpd_mismatch(c("y ~ x", "sigma ~ x"), w, brms::brmsfamily("student"),
                 frmtmb::brmsfamily("student")), character(0))
})

test_that("ordinal thresholds take brms's student_t(3, 0, 2.5)", {
  skip_unless_brms()
  D <- dpd_data()
  for (fam in c("cumulative", "sratio", "cratio", "acat")) {
    expect_identical(
      dpd_mismatch("y ~ x + (1 | g)", D$ord, brms::brmsfamily(fam),
                   frmtmb::brmsfamily(fam)), character(0), label = fam)
  }
  expect_identical(
    dpd_mismatch("y ~ x", D$ord, brms::cumulative(threshold = "equidistant"),
                 frmtmb::cumulative(threshold = "equidistant")),
    character(0))
  # the density reaches the thresholds: one entry, on the threshold block
  uf <- frm(bf(y ~ x) + cumulative(), data = D$ord, dry_run = "objective")
  ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(uf, NULL))$ri
  th <- Filter(function(e) identical(e$scale, "ordthres"), ri$entries)
  expect_length(th, 1L)
  expect_identical(th[[1]]$dist$scale, 2.5)
})

test_that("a non-decision time takes uniform(0, min_Y)", {
  skip_unless_brms()
  D <- dpd_data()
  expect_identical(
    dpd_mismatch("y ~ x", D$shift, brms::brmsfamily("shifted_lognormal"),
                 frmtmb::brmsfamily("shifted_lognormal"),
                 miny = min(D$shift$y)), character(0))
  # the density is uniform below the response's minimum and NaN past
  # it, on the tape, where a guard with a constant operand would vanish
  # (dev/rtmb-pitfalls.md, item 22)
  uf <- frm(bf(y ~ x) + shifted_lognormal(), data = D$shift,
            dry_run = "objective")
  ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(uf, NULL))$ri
  e <- Filter(function(e) identical(e$dist$kind, "uniform"), ri$entries)
  expect_length(e, 1L)
  expect_identical(e[[1]]$dist$upper, min(D$shift$y))
  nlp <- neg_log_prior_fn(e)
  pars <- uf$frame$par_template
  at <- function(ndt) {
    p <- pars
    p$betad[e[[1]]$idx] <- log(ndt)
    nlp(p)
  }
  inside <- min(D$shift$y) / 2
  # -log density of uniform(0, u) at a log-linked value: log(u) - log(ndt)
  expect_equal(at(inside), log(min(D$shift$y)) - log(inside),
               tolerance = 1e-12)
  expect_true(is.nan(at(2 * min(D$shift$y))))
  tp <- RTMB::MakeTape(function(v) {
    p <- pars
    "[<-" <- RTMB::ADoverload("[<-")
    p$betad[e[[1]]$idx] <- v
    nlp(p)
  }, 0)
  expect_true(is.nan(tp(log(2 * min(D$shift$y)))))
  expect_equal(tp(log(inside)), at(inside), tolerance = 1e-12)
})

test_that("an unstr() correlation takes lkj(1), class cortime", {
  skip_unless_brms()
  D <- dpd_data()
  d <- transform(D$real, t = rep(1:10, each = 6))
  expect_identical(
    dpd_mismatch("y ~ x + unstr(time = t, gr = g)", d,
                 brms::brmsfamily("gaussian"), frmtmb::brmsfamily("gaussian")),
    character(0))
})

test_that("a mixture's sigmas take brms's default, its weights dirichlet(1)", {
  skip_unless_brms()
  set.seed(20261019)
  d <- data.frame(x = stats::rnorm(60), y = c(stats::rnorm(30, -1),
                                              stats::rnorm(30, 2)))
  expect_identical(
    suppressMessages(dpd_mismatch(
      "y ~ x", d, brms::mixture("gaussian", "gaussian"),
      frmtmb::mixture(stats::gaussian(), stats::gaussian()))),
    "theta||||: brms dirichlet(1), frmtmb NA")
  # no class names the weights here, so the Dirichlet is a density the
  # sampler adds: lgamma(K) + sum(log w) on the log ratios it samples
  fit <- frm(bf(y ~ 1), data = d,
             family = mixture(stats::gaussian(), stats::gaussian(),
                              stats::gaussian()), dry_run = "objective")
  rp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))
  expect_true(rp$mix_theta)
  fit$prior <- rp$effective
  f <- frmtmb.sample:::mix_theta_nlp(fit)
  tm <- frmtmb.sample:::mix_theta_terms(fit)[[1]]
  pars <- fit$estimates
  pars$betad[tm$idx] <- c(0.4, -1.1)
  w <- exp(c(0.4, -1.1, 0)) / sum(exp(c(0.4, -1.1, 0)))
  expect_equal(-as.numeric(f(pars)), lgamma(3) + sum(log(w)),
               tolerance = 1e-12)
  # prior = "flat" drops it with the other defaults
  rf <- suppressWarnings(suppressMessages(
    frmtmb.sample:::sample_resolve_priors(fit, "flat")))
  expect_false(isTRUE(rf$mix_theta))
})

test_that("a natural-scale default is brms's density with the link Jacobian", {
  # the log density the tape adds for a dpar without a formula is the
  # density of the parameter times |d parameter / d coefficient|
  D <- dpd_data()
  chk <- function(form, fam, data, dp, dens) {
    uf <- frm(form, family = fam, data = data, dry_run = "objective")
    ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(uf,
                                                                 NULL))$ri
    e <- Filter(function(e) e$scale %in% c("sd", "natural"), ri$entries)
    e <- Filter(function(z) {
      identical(names(uf$frame$par_template$betad)[z$idx],
                paste0(dp, "_(Intercept)"))
    }, e)
    expect_length(e, 1L)
    lp <- Filter(function(l) identical(l$dpar, dp), uf$frame$linpreds)[[1]]
    eta <- 0.7
    p <- uf$frame$par_template
    p$betad[e[[1]]$idx] <- eta
    v <- lp$link$linkinv(eta)
    want <- dens(v) + log(abs(lp$link$mu_eta(eta)))
    expect_equal(-neg_log_prior_fn(e)(p), want, tolerance = 1e-12,
                 label = paste(fam$family, dp))
  }
  dinvg <- function(v, a, b) a * log(b) - lgamma(a) - (a + 1) * log(v) - b / v
  chk(bf(y ~ x), negbinomial(), D$cnt, "shape",
      function(v) dinvg(v, 0.4, 0.3))
  chk(bf(y ~ x), Beta(), D$unit, "phi",
      function(v) stats::dgamma(v, 0.01, 0.01, log = TRUE))
  chk(bf(y ~ x), student(), D$real, "nu",
      function(v) stats::dgamma(v, 2, 0.1, log = TRUE))
  chk(bf(y ~ x), zero_inflated_poisson(), D$zcnt, "zi",
      function(v) stats::dbeta(v, 1, 1, log = TRUE))
})

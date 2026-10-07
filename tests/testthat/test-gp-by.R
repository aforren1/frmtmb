# brms's gp(..., by = ): one Gaussian process per design column of a
# factor by-variable (per level under cmc = TRUE), a numeric by-variable
# multiplying one GP, and brms's names for every sub-GP's parameters.
# The closed-form references are the gaussian marginal likelihood, which
# integrates every latent field exactly. brms's own data and Stan program
# are compared in the gated blocks at the end.

gpby_data <- function(n = 60, seed = 5) {
  set.seed(seed)
  d <- data.frame(x = round(stats::runif(n, 0, 6), 1),
                  f = factor(rep(c("a", "b", "c"), length.out = n)),
                  w = stats::runif(n, 0.5, 2))
  d$y <- 0.5 + ifelse(d$f == "a", sin(d$x),
                      ifelse(d$f == "b", cos(d$x), 0.2 * d$x)) +
    stats::rnorm(n, 0, 0.3)
  d
}

# The marginal covariance of y under independent exact GPs, sub-GP j on
# rows `rows[[j]]` multiplied by `mult[[j]]`, with frmtmb's nugget (1e-6
# on the correlation) and one latent value per distinct position.
gpby_ref_cov <- function(x, sub, sd, ell, sigma) {
  n <- length(x)
  S <- diag(sigma^2, n)
  for (j in seq_along(sub)) {
    r <- sub[[j]]$rows
    pos <- sort(unique(x[r]))
    Z <- matrix(0, n, length(pos))
    Z[cbind(r, match(x[r], pos))] <- sub[[j]]$mult
    K <- sd[j]^2 * (exp(-outer(pos, pos, "-")^2 / (2 * ell[j]^2)) +
                      diag(1e-6, length(pos)))
    S <- S + Z %*% K %*% t(Z)
  }
  S
}

gpby_ref_nll <- function(y, mu, S) {
  r <- y - mu
  0.5 * (as.numeric(determinant(S)$modulus) + sum(r * solve(S, r)) +
           length(y) * log(2 * pi))
}

# frmtmb's estimates for the reference: the sub-GP blocks in term order
gpby_est <- function(fit) {
  bks <- Filter(function(b) b$covstruct == "gp", fit$frame$re_blocks)
  th <- fit$estimates$theta
  list(sd = vapply(bks, function(b) exp(th[b$theta_idx][1]), 0),
       ell = vapply(bks, function(b) exp(th[b$theta_idx][2]), 0),
       mu = unname(fit$estimates$beta[1]),
       sigma = unname(sigma(fit))[1])
}

gpby_sub <- function(d, by) {
  if (is.numeric(d[[by]])) {
    return(list(list(rows = seq_len(nrow(d)), mult = d[[by]])))
  }
  lapply(levels(d[[by]]), function(l) {
    list(rows = which(d[[by]] == l), mult = rep(1, sum(d[[by]] == l)))
  })
}

test_that("gp(x, by = f) is one GP per level, as the closed form says", {
  d <- gpby_data()
  fit <- frm(bf(y ~ gp(x, by = f)), data = d)
  bks <- Filter(function(b) b$covstruct == "gp", fit$frame$re_blocks)
  expect_length(bks, 3L)
  # each level's GP lives on that level's own distinct positions
  for (j in 1:3) {
    lev <- levels(d$f)[j]
    expect_identical(bks[[j]]$dim,
                     length(unique(d$x[d$f == lev])))
  }
  e <- gpby_est(fit)
  S <- gpby_ref_cov(d$x, gpby_sub(d, "f"), e$sd, e$ell, e$sigma)
  ref <- gpby_ref_nll(d$y, e$mu, S)
  # the same density at the same parameters: an identity, so the gap is
  # rounding against the size of the value
  expect_lt(abs(-as.numeric(logLik(fit)) - ref), 1e-8 * abs(ref))
  # and frmtmb is at the optimum of that density
  p0 <- c(e$mu, log(e$sd), log(e$ell), log(e$sigma)) + 0.05
  nll <- function(p) {
    gpby_ref_nll(d$y, p[1], gpby_ref_cov(d$x, gpby_sub(d, "f"),
                                         exp(p[2:4]), exp(p[5:7]),
                                         exp(p[8])))
  }
  op <- stats::optim(p0, nll, method = "BFGS",
                     control = list(reltol = 1e-12, maxit = 2000))
  expect_lt(-as.numeric(logLik(fit)) - op$value, 1e-6 * abs(op$value))
})

test_that("a numeric by multiplies one GP, and cmc = FALSE fits contrasts", {
  d <- gpby_data()
  fw <- frm(bf(y ~ gp(x, by = w)), data = d)
  expect_length(Filter(function(b) b$covstruct == "gp",
                       fw$frame$re_blocks), 1L)
  e <- gpby_est(fw)
  ref <- gpby_ref_nll(d$y, e$mu, gpby_ref_cov(d$x, gpby_sub(d, "w"),
                                              e$sd, e$ell, e$sigma))
  expect_lt(abs(-as.numeric(logLik(fw)) - ref), 1e-8 * abs(ref))
  expect_identical(variables(fw),
                   c("b_Intercept", "sdgp_gpxw", "lscale_gpxw", "sigma"))

  # cmc = FALSE: an intercept GP over every row and a contrast GP over
  # the rows of each other level, brms's model.matrix(~ 1 + byval)
  fc <- frm(bf(y ~ gp(x, by = f, cmc = FALSE)), data = d)
  e <- gpby_est(fc)
  sub <- list(list(rows = seq_len(nrow(d)), mult = rep(1, nrow(d))),
              list(rows = which(d$f == "b"), mult = rep(1, sum(d$f == "b"))),
              list(rows = which(d$f == "c"), mult = rep(1, sum(d$f == "c"))))
  ref <- gpby_ref_nll(d$y, e$mu, gpby_ref_cov(d$x, sub, e$sd, e$ell,
                                              e$sigma))
  expect_lt(abs(-as.numeric(logLik(fc)) - ref), 1e-8 * abs(ref))
  expect_true(all(c("sdgp_gpxfIntercept", "sdgp_gpxfb", "sdgp_gpxfc",
                    "lscale_gpxfIntercept") %in% variables(fc)))
})

test_that("a numeric by scales the kriging covariance entry by entry", {
  # gp_krig_cov() scales by w w' one column block at a time, where
  # `S * outer(w, w)` held two more n x n (dev/ciharden-findings.md).
  # It must be that matrix to the bit, symmetric, with extra_var on its
  # diagonal. 850 rows make two blocks of 2^19 / 850 = 616 columns, and
  # the repeated rows take the expansion over distinct positions.
  d <- gpby_data()
  fw <- frm(bf(y ~ gp(x, by = w)), data = d)
  set.seed(3)
  xs <- seq(-0.5, 6.5, length.out = 800) + 1e-3
  nd <- data.frame(x = c(xs, xs[1:50]), w = stats::runif(850, 0.5, 2))
  ed <- lp_eta_design(fw, fw$frame$linpreds[["y.mu"]], nd, FALSE, FALSE)
  sp <- Filter(function(s) !is.null(s$krig), ed$sm_parts)[[1L]]
  kg <- sp$krig
  expect_length(kg$rows, 850L)
  S <- gp_krig_cov(kg)
  k1 <- kg
  k1$w <- rep(1, length(kg$w))
  expect_identical(S, gp_krig_cov(k1) * outer(kg$w, kg$w))
  expect_identical(S, t(S))
  expect_identical(diag(S), sp$extra_var[kg$rows])
})

test_that("gr = FALSE keeps one latent per row, the same model", {
  # distinct positions, so the two parameterizations are one density
  d <- gpby_data()
  d$x <- d$x + seq_len(nrow(d)) * 1e-3
  f1 <- frm(bf(y ~ gp(x, by = f)), data = d)
  f2 <- frm(bf(y ~ gp(x, by = f, gr = FALSE)), data = d)
  l1 <- as.numeric(logLik(f1))
  expect_lt(abs(as.numeric(logLik(f2)) - l1), 1e-6 * abs(l1))
  # with ties, gr = FALSE gives the duplicates separate latent values
  d2 <- gpby_data()
  f3 <- frm(bf(y ~ gp(x, by = f, gr = FALSE)), data = d2)
  bk <- Filter(function(b) b$covstruct == "gp", f3$frame$re_blocks)
  expect_identical(vapply(bk, `[[`, 0L, "dim"),
                   as.integer(table(d2$f)))
})

test_that("variables() lists a gp() term's sdgp_ and lscale_", {
  # brms's variables(fit6) lists lscale_volume_gpAgeTrt0 (brmsfit-methods
  # test line 991); before this a gp() fit listed neither parameter
  d <- gpby_data()
  f0 <- frm(bf(y ~ gp(x)), data = d)
  expect_identical(variables(f0),
                   c("b_Intercept", "sdgp_gpx", "lscale_gpx", "sigma"))
  bk <- f0$frame$re_blocks[[1]]
  th <- f0$estimates$theta[bk$theta_idx]
  # brms's length scale is on inputs divided by their largest distance
  h <- hypothesis(f0, "lscale_gpx = 0", class = NULL)
  dmax <- max(d$x) - min(d$x)
  expect_lt(abs(h$hypothesis$Estimate * dmax / exp(th[2]) - 1), 1e-12)
})

test_that("brms's names: sdgp_ and lscale_ per by-level, in brms's order", {
  d <- gpby_data()
  fit <- frm(bf(y ~ gp(x, by = f)), data = d)
  v <- variables(fit)
  expect_identical(v, c("b_Intercept", "sdgp_gpxfa", "sdgp_gpxfb",
                        "sdgp_gpxfc", "lscale_gpxfa", "lscale_gpxfb",
                        "lscale_gpxfc", "sigma"))
  s <- summary(fit)
  expect_identical(rownames(s$gp),
                   c("sdgp(gpxfa)", "sdgp(gpxfb)", "sdgp(gpxfc)",
                     "lscale(gpxfa)", "lscale(gpxfb)", "lscale(gpxfc)"))
  # brms's lscale is on each level's inputs divided by the largest
  # distance between two of that level's positions
  e <- gpby_est(fit)
  for (j in 1:3) {
    xr <- d$x[d$f == levels(d$f)[j]]
    lsj <- e$ell[j] / (max(xr) - min(xr))
    expect_lt(abs(s$gp[3 + j, "Estimate"] / lsj - 1), 1e-12)
    expect_lt(abs(s$gp[j, "Estimate"] / e$sd[j] - 1), 1e-12)
  }
  h <- hypothesis(fit, "sdgp_gpxfa > sdgp_gpxfc", class = NULL)
  expect_lt(abs(h$hypothesis$Estimate - (e$sd[1] - e$sd[3])),
            1e-12 * e$sd[1])
  # a non-isotropic term names one length scale per covariate, the
  # levels varying fastest, as brms's as.vector(sfx2)
  d$z <- stats::runif(nrow(d), 0, 3)
  f2 <- frm(bf(y ~ gp(x, z, by = f, iso = FALSE, k = 5)), data = d)
  expect_identical(rownames(summary(f2)$gp)[4:9],
                   paste0("lscale(gpxzf", rep(c("a", "b", "c"), 2),
                          rep(c("x", "z"), each = 3), ")"))
})

test_that("prediction reads each row's own level, and refuses a new one", {
  d <- gpby_data()
  fit <- frm(bf(y ~ gp(x, by = f)), data = d)
  # in sample through newdata reproduces the fit
  p <- frm_linpred(fit, newdata = d)
  expect_lt(max(abs(p - frm_linpred(fit))),
            1e-10 * max(abs(frm_linpred(fit))))
  # an unseen position kriges each level's own field, and a row of
  # level b moves only with b's coefficients
  nd <- data.frame(x = c(2.55, 2.55, 7), f = factor(c("a", "b", "b")))
  lb <- frm_lp_basis(fit, newdata = nd)
  bks <- Filter(function(b) b$covstruct == "gp", fit$frame$re_blocks)
  jc <- frm_joint_cov(fit)
  bpos <- which(jc$names == "b")
  col_a <- lb$coef_pos %in% bpos[bks[[1]]$b_idx]
  col_b <- lb$coef_pos %in% bpos[bks[[2]]$b_idx]
  expect_true(all(lb$A[1, col_b] == 0))
  expect_true(all(lb$A[2:3, col_a] == 0))
  expect_true(any(lb$A[1, col_a] != 0))
  # by-level kriging is the closed form of that level alone
  e <- gpby_est(fit)
  xa <- sort(unique(d$x[d$f == "a"]))
  Ka <- e$sd[1]^2 * (exp(-outer(xa, xa, "-")^2 / (2 * e$ell[1]^2)) +
                       diag(1e-6, length(xa)))
  ks <- e$sd[1]^2 * exp(-(2.55 - xa)^2 / (2 * e$ell[1]^2))
  ba <- fit$estimates$b[bks[[1]]$b_idx]
  mu_a <- e$mu + sum(solve(Ka, ks) * ba)
  expect_lt(abs(lb$eta[1] - mu_a), 1e-8 * abs(mu_a))

  expect_error(frm_linpred(fit, newdata = data.frame(x = 1, f = "d")),
               "New factor levels are not allowed")
  expect_error(frm_linpred(fit, newdata = data.frame(x = 1, f = "d"),
                           allow_new_levels = TRUE),
               "New factor levels are not allowed")
  # a variable missing from newdata is looked up from the formula's
  # environment, by R's rule, so any `f` a session or a test runner
  # defines would answer; this fit's formula sees only base R, so the
  # refusal does not depend on what the global environment holds
  fo <- y ~ gp(x, by = f)
  environment(fo) <- new.env(parent = baseenv())
  fit_b <- frm(bf(fo), data = d)
  expect_error(frm_linpred(fit_b, newdata = data.frame(x = 1)),
               "by variable 'f'")
})

test_that("conditional_effects() draws a gp() term's covariates", {
  # brms's get_all_effects_type(x, "gp") reads them off the term; before
  # this a y ~ gp(x) fit had "No plottable predictors"
  d <- gpby_data()
  f0 <- frm(bf(y ~ gp(x)), data = d)
  expect_identical(names(conditional_effects(f0)), "x")
  fit <- frm(bf(y ~ gp(x, by = f)), data = d)
  ce <- conditional_effects(fit)
  expect_identical(names(ce), c("x", "f", "x:f"))
  # the interaction display is each level's own curve
  g <- ce[["x:f"]]
  nd <- data.frame(x = g$x, f = g$f)
  expect_lt(max(abs(g$estimate__ - frm_linpred(fit, newdata = nd))),
            1e-10 * max(abs(g$estimate__)))
})

test_that("gp() takes brms's arguments and refuses what it does not have", {
  d <- gpby_data()
  expect_error(frm(bf(y ~ gp(x, by = factor(f))), data = d),
               "must name one column")
  expect_error(frm(bf(y ~ gp(x, cov = "matern32")), data = d),
               "not implemented here yet")
  expect_error(frm(bf(y ~ gp(x, cov = "nope")), data = d),
               "not a valid GP covariance kernel")
  expect_error(frm(bf(y ~ gp(x, bogus = 1)), data = d), "unknown: bogus")
  # by = NA and k = NA are brms's spellings of "none"
  l0 <- as.numeric(logLik(frm(bf(y ~ gp(x)), data = d)))
  l1 <- as.numeric(logLik(frm(bf(y ~ gp(x, by = NA, k = NA)), data = d)))
  expect_identical(l0, l1)
  # scale = FALSE: brms's length scale is then in data units
  fs <- frm(bf(y ~ gp(x, scale = FALSE)), data = d)
  bk <- fs$frame$re_blocks[[1]]
  expect_identical(summary(fs)$gp["lscale(gpx)", "Estimate"],
                   exp(fs$estimates$theta[bk$theta_idx][2]))
  expect_identical(as.numeric(logLik(fs)), l0)
})

test_that("brms's data for gp(x, by = f): each level's rows, scale, basis", {
  skip_unless_brms()
  d <- gpby_data()
  # exact: each sub-GP sits on its level's distinct positions, divided
  # by that level's own largest distance
  sdat <- brms_standata(brms::bf(y ~ gp(x, by = f)), data = d,
                        family = gaussian())
  fr <- frm(bf(y ~ gp(x, by = f)), data = d, dry_run = "frame")
  expect_identical(as.integer(sdat$Kgp_1), 3L)
  gis <- fr$linpreds[["y.mu"]]$gps
  Zx <- as.matrix(fr$linpreds[["y.mu"]]$Z)
  for (j in 1:3) {
    bk <- fr$re_blocks[[gis[[j]]$block_id]]
    rows <- which(rowSums(abs(Zx[, bk$c_idx, drop = FALSE])) > 0)
    expect_identical(rows, as.integer(sdat[[paste0("Igp_1_", j)]]))
    xb <- as.numeric(sdat[[paste0("Xgp_1_", j)]])
    # 1 by construction, so the ratio to it is the rounding of a division
    expect_lt(abs(max(xb) - min(xb) - 1), 8 * .Machine$double.eps)
    expect_equal(sort(xb) * bk$gp_lscale_div, gis[[j]]$positions[, 1],
                 tolerance = 1e-12)
  }
  # Hilbert-space: the same per-level scaling, center and boundary, so
  # the basis rows are brms's, row for row
  sdh <- brms_standata(brms::bf(y ~ gp(x, by = f, k = 8)), data = d,
                       family = gaussian())
  frh <- frm(bf(y ~ gp(x, by = f, k = 8)), data = d, dry_run = "frame")
  Z <- as.matrix(frh$linpreds[["y.mu"]]$Z)
  for (j in 1:3) {
    gi <- frh$linpreds[["y.mu"]]$gps[[j]]
    bk <- frh$re_blocks[[gi$block_id]]
    ig <- as.integer(sdh[[paste0("Igp_1_", j)]])
    xb <- sdh[[paste0("Xgp_1_", j)]][sdh[[paste0("Jgp_1_", j)]], ,
                                     drop = FALSE]
    expect_equal(unname(Z[ig, bk$c_idx]), unname(xb), tolerance = 1e-12)
    expect_true(all(Z[-ig, bk$c_idx] == 0))
    expect_equal(as.vector(gi$omega),
                 as.vector(sdh[[paste0("slambda_1_", j)]]),
                 tolerance = 1e-12)
  }
})

test_that("check C: gp(by = ) Hilbert-space terms are brms's joint density", {
  skip_unless_brms_fit()
  d <- gpby_data()
  for (form in c("y ~ gp(x, by = f, k = 8)", "y ~ gp(x, by = w, k = 8)")) {
    bform <- brms::bf(stats::as.formula(form))
    fit <- frm(bf(stats::as.formula(form)), data = d)
    brms_lp_check(bform, gaussian(), d, fit, joint = TRUE)
  }
})

test_that("check C: a gp(by = ) in sigma and in a multivariate model", {
  # brms suffixes a GP's data and latent names with its dpar and its
  # response (Kgp_sigma_1, zgp_y_1_1), so these reach the translator's
  # suffixed branches that a mu-only univariate model never does.
  skip_unless_brms_fit()
  d <- gpby_data()
  set.seed(31)
  d$y2 <- 1 + cos(d$x) * (d$f == "a") + stats::rnorm(nrow(d), 0, 0.4)
  d$ysig <- 0.5 + sin(d$x) +
    stats::rnorm(nrow(d), 0, exp(-1 + 0.3 * (d$f == "b") * sin(d$x)))
  # two by-level GP sds run to 0 in sigma, so they have no standard
  # error, and frm() says so; the test is about the joint density
  lost <- "Standard errors are not available"
  fs <- allow_warnings(frm(bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)),
                           data = d), lost)
  brms_lp_check(brms::bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)),
                gaussian(), d, fs, joint = TRUE)
  fm <- frm(mvbf(bf(y ~ gp(x, by = f, k = 8)), bf(y2 ~ gp(x, k = 6))) +
              set_rescor(FALSE), data = d, family = gaussian())
  brms_lp_check(brms::mvbf(brms::bf(y ~ gp(x, by = f, k = 8)),
                           brms::bf(y2 ~ gp(x, k = 6))) +
                  brms::set_rescor(FALSE), gaussian(), d, fm, joint = TRUE)
  fb <- allow_warnings(frm(mvbf(bf(y ~ gp(x, by = f, k = 8)),
                                bf(y2 ~ gp(x, by = f, k = 6))) +
                             set_rescor(FALSE),
                           data = d, family = gaussian()), lost)
  brms_lp_check(brms::mvbf(brms::bf(y ~ gp(x, by = f, k = 8)),
                           brms::bf(y2 ~ gp(x, by = f, k = 6))) +
                  brms::set_rescor(FALSE), gaussian(), d, fb, joint = TRUE)
})

test_that("priors on sdgp and lscale land where brms puts them", {
  d <- gpby_data()
  f0 <- frm(bf(y ~ gp(x, by = f)), data = d)
  pr <- set_prior("normal(0.3, 0.1)", class = "lscale", coef = "gpxfb") +
    set_prior("exponential(2)", class = "sdgp", coef = "gpxfa")
  fp <- frm(bf(y ~ gp(x, by = f)), data = d, prior = pr)
  # the MAP objective is the likelihood less brms's log densities, on
  # brms's scales with the log-Jacobian of theta: lscale is exp(theta)
  # over the level's largest distance, sdgp is exp(theta)
  bks <- Filter(function(b) b$covstruct == "gp", fp$frame$re_blocks)
  par <- fp$obj$env$last.par.best
  th <- fp$estimates$theta
  ls_b <- exp(th[bks[[2]]$theta_idx][2]) / bks[[2]]$gp_lscale_div
  sd_a <- exp(th[bks[[1]]$theta_idx][1])
  lp <- stats::dnorm(ls_b, 0.3, 0.1, log = TRUE) + log(ls_b) +
    stats::dexp(sd_a, 2, log = TRUE) + log(sd_a)
  nll0 <- f0$obj$env$f(par)
  pen <- fp$obj$env$f(par)
  expect_lt(abs((pen - nll0) / (-lp) - 1), 1e-8)
  # a bound on the length scale is on brms's scale too
  fb <- frm(bf(y ~ gp(x, by = f)), data = d,
            prior = set_prior("", class = "lscale", coef = "gpxfc",
                              lb = 0.5))
  bc <- Filter(function(b) b$covstruct == "gp", fb$frame$re_blocks)[[3]]
  expect_gte(exp(fb$estimates$theta[bc$theta_idx][2]) / bc$gp_lscale_div,
             0.5 * (1 - 1e-8))
  # class "sd" is a group's, as in brms, and does not reach a gp()
  expect_error(frm(bf(y ~ gp(x)), data = d,
                   prior = set_prior("exponential(1)", class = "sd")),
               "class \"sdgp\"")
  expect_error(frm(bf(y ~ gp(x)), data = d,
                   prior = set_prior("exponential(1)", class = "sdgp",
                                     coef = "gpz")),
               "coef = \"gpx\"")
  # brms's own table carries over row for row
  skip_if_not_installed("brms")
  bp <- brms::default_prior(brms::bf(y ~ gp(x, by = f)), data = d)
  fq <- frm(bf(y ~ gp(x, by = f)), data = d, prior = bp)
  expect_true(is.finite(as.numeric(logLik(fq))))
})

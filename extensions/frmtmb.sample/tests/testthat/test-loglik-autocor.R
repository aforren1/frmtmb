# log_lik() of draws for a response with a residual covariance term
# (ar(cov = TRUE), ma(), arma(), cosy(), unstr(), fcor(), a registered
# term): brms 2.23.0's leave-one-out conditional density of each row
# given the other rows of its group. Three references: the conditional
# of a dense multivariate normal (or t, group by group) by the partition
# formula, at hand-set draws; core's autocor_row_loglik() at each draw;
# and brms's own log_lik() at the same parameter vectors (gated: Stan's
# Fixed_param, one chain per vector).

skip_on_cran()

lla_data <- function(seed = 4, gap = FALSE) {
  set.seed(seed)
  sizes <- c(7, 5, 6, 4, 7, 3, 6, 5)
  d <- do.call(rbind, lapply(seq_along(sizes), function(g) {
    data.frame(g = g, t = seq_len(sizes[g]))
  }))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(nrow(d)) + stats::rnorm(8)[d$g]
  d <- d[sample(nrow(d)), ]
  if (gap) d <- d[!(d$g == 2 & d$t == 3) & !(d$g == 5 & d$t == 1), ]
  rownames(d) <- NULL
  d
}

# draws at chosen parameter vectors, with no sampler behind them; the
# columns are brms's scale where brms has the parameter (b_*, sigma,
# nu) and frmtmb's internal one for thetaac_k
lla_draws <- function(fit, vals) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  k <- max(lengths(vals))
  m <- matrix(0, k, length(lab), dimnames = list(NULL, lab))
  for (nm in names(vals)) {
    stopifnot(nm %in% lab)
    m[, nm] <- vals[[nm]]
  }
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}

# frmtmb's chart for a partial autocorrelation, inverted
lla_pinv <- function(r) r / sqrt(1 - r^2)

lla_cond <- function(y, m, S, nu = NULL) {
  n <- length(y)
  vapply(seq_len(n), function(i) {
    o <- seq_len(n)[-i]
    A <- solve(S[o, o, drop = FALSE])
    r <- y[o] - m[o]
    cm <- m[i] + as.numeric(S[i, o, drop = FALSE] %*% A %*% r)
    cv <- S[i, i] - as.numeric(S[i, o, drop = FALSE] %*% A %*% S[o, i])
    if (is.null(nu)) return(stats::dnorm(y[i], cm, sqrt(cv), log = TRUE))
    df <- nu + n - 1
    sc <- sqrt((nu + as.numeric(t(r) %*% A %*% r)) / df * cv)
    stats::dt((y[i] - cm) / sc, df, log = TRUE) - log(sc)
  }, 0)
}

# the dense reference for groups g whose correlation over time levels t
# is corf(lag matrix): computed here from the parameters, not read from
# frmtmb
lla_grouped <- function(d, m, s, corf, nu = NULL) {
  out <- numeric(nrow(d))
  for (g in unique(d$g)) {
    r <- which(d$g == g)
    S <- outer(s[r], s[r]) * corf(abs(outer(d$t[r], d$t[r], "-")))
    out[r] <- lla_cond(d$y[r], m[r], S, nu)
  }
  out
}

b0 <- c(1.1, 0.9, 1.3)
b1 <- c(0.4, 0.6, 0.2)
sg <- c(0.8, 1.2, 1.0)
nus <- c(5, 12, 3.5)

test_that("log_lik() of an ar() family covariance block is the dense
           conditional of each group, at each draw", {
  d <- lla_data()
  # the autocovariance at unit innovation variance: the correlation of
  # stats::ARMAacf() times the process variance from the psi weights,
  # since sigma is the innovation sd, as in brms
  acf_of <- function(a, m) function(L) {
    v <- stats::ARMAacf(ar = a[a != 0], ma = m[m != 0], lag.max = max(L))
    psi <- c(1, stats::ARMAtoMA(ar = a[a != 0], ma = m[m != 0],
                                lag.max = 5000))
    matrix(v[L + 1] * sum(psi^2), nrow(L))
  }
  cases <- list(
    ar = list(f = y ~ x + ar(t, g, cov = TRUE), a = c(0.5, -0.3, 0.8),
              m = c(0, 0, 0)),
    ma = list(f = y ~ x + ma(t, g, cov = TRUE), a = c(0, 0, 0),
              m = c(0.4, -0.6, 0.2)),
    arma = list(f = y ~ x + arma(t, g, cov = TRUE), a = c(0.5, -0.2, 0.3),
                m = c(0.3, 0.5, -0.4)))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    for (fam in c("gaussian", "student")) {
      fit <- frm(cs$f, family = get(fam)(), data = d, dry_run = "objective")
      vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
      th <- list()
      if (nm != "ma") th$a <- lla_pinv(cs$a)
      if (nm != "ar") th$m <- lla_pinv(cs$m)
      for (k in seq_along(th)) vals[[paste0("thetaac_", k)]] <- th[[k]]
      if (fam == "student") vals$nu <- nus
      ds <- lla_draws(fit, vals)
      ll <- log_lik(ds)
      expect_equal(dim(ll), c(3L, nrow(d)))
      idx <- frmtmb.sample:::draws_par_index(ds$fit)
      for (s in 1:3) {
        ref <- lla_grouped(d, b0[s] + b1[s] * d$x, rep(sg[s], nrow(d)),
                           acf_of(cs$a[s], cs$m[s]),
                           if (fam == "student") nus[s])
        expect_lt(max(abs(ll[s, ] - ref)), 1e-10 * max(abs(ref)))
        # and core's function at the same draw is what the column holds
        sh <- frmtmb.sample:::draws_fit_at(ds, s, idx)
        expect_identical(ll[s, ], frmtmb::autocor_row_loglik(
          sh, frmtmb::eval_dpars(sh), "y"))
      }
    }
  }
})

test_that("log_lik() of cosy(), unstr() and fcor() is the dense
           conditional, with gaps and sigma ~ x", {
  d <- lla_data(gap = TRUE)
  n <- nrow(d)
  for (fam in c("gaussian", "student")) {
    nu <- if (fam == "student") nus
    fit <- frm(y ~ x + cosy(t, g), family = get(fam)(), data = d,
               dry_run = "objective")
    a <- 1 / (fit$frame$autocor$y$d - 1)
    rho <- c(0.4, -0.05, 0.7)
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg,
                 thetaac_1 = stats::qlogis((rho + a) / (1 + a)), nu = nu)
    ll <- log_lik(lla_draws(fit, vals[lengths(vals) > 0]))
    for (s in 1:3) {
      ref <- lla_grouped(d, b0[s] + b1[s] * d$x, rep(sg[s], n),
                         function(L) ifelse(L == 0, 1, rho[s]), nu[s])
      expect_lt(max(abs(ll[s, ] - ref)), 1e-10 * max(abs(ref)))
    }
  }
  # unstr(): the reference correlation is frmtmb's own matrix at the
  # draw, which the rows read by time level
  fit <- frm(y ~ x + unstr(t, g), data = d, dry_run = "objective")
  ac <- fit$frame$autocor$y
  set.seed(40)
  th <- matrix(stats::rnorm(3 * ac$npar, 0, 0.5), 3)
  vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
  for (k in seq_len(ac$npar)) vals[[paste0("thetaac_", k)]] <- th[, k]
  ds <- lla_draws(fit, vals)
  ll <- log_lik(ds)
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  for (s in 1:3) {
    R <- frmtmb::autocor_matrix(frmtmb.sample:::draws_fit_at(ds, s, idx))
    ref <- numeric(n)
    for (g in unique(d$g)) {
      r <- which(d$g == g)
      ti <- match(as.character(d$t[r]), rownames(R))
      ref[r] <- lla_cond(d$y[r], b0[s] + b1[s] * d$x[r],
                         sg[s]^2 * R[ti, ti])
    }
    expect_lt(max(abs(ll[s, ] - ref)), 1e-10 * max(abs(ref)))
  }
  # fcor() with sigma ~ x: one group of every row, D V D
  V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)
  fit <- frm(bf(y ~ x + fcor(V), sigma ~ x), data = d, data2 = list(V = V),
             dry_run = "objective")
  bs <- c(-0.2, 0.1, 0.3)
  vals <- list(b_Intercept = b0, b_x = b1, b_sigma_Intercept = bs,
               b_sigma_x = c(0.3, -0.2, 0))
  ll <- log_lik(lla_draws(fit, vals))
  for (s in 1:3) {
    sv <- exp(bs[s] + vals$b_sigma_x[s] * d$x)
    ref <- lla_cond(d$y, b0[s] + b1[s] * d$x, outer(sv, sv) * V)
    expect_lt(max(abs(ll[s, ] - ref)), 1e-10 * max(abs(ref)))
  }
})

test_that("a registered term without a jacobian slot gets the same
           columns", {
  frmtmb::frmtmb_local_registries()
  # a lag SAR over a path graph, registered through the exported seam
  # alone and with no jacobian: log_lik() inverts unwhiten() instead
  frmtmb::frmtmb_register_autocor(
    "lltest", args = "W",
    build = function(spec, rows, n_all, lookup) {
      W <- lookup(spec$W, "lltest()")[rows, rows, drop = FALSE]
      ev <- Re(eigen(W, only.values = TRUE)$values)
      list(W = W, ev = ev, lo = 1 / min(ev), hi = 1 / max(ev))
    },
    npar = 1L,
    start = function(aux) stats::qlogis(-aux$lo / (aux$hi - aux$lo)),
    natural = function(theta, aux) {
      c(lagsar = aux$lo + (aux$hi - aux$lo) / (1 + exp(-theta[1])))
    },
    bounds = function(aux) cbind(lower = aux$lo, upper = aux$hi),
    whiten = function(y, mu, sigma, theta, aux) {
      rho <- aux$lo + (aux$hi - aux$lo) / (1 + exp(-theta[1]))
      list(u = (y - rho * as.vector(aux$W %*% y) - mu) / sigma,
           logjac = sum(log1p(-rho * aux$ev)) - length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      rho <- aux$lo + (aux$hi - aux$lo) / (1 + exp(-theta[1]))
      solve(diag(nrow(aux$W)) - rho * aux$W, mu + sigma * u)
    },
    sigma_varies = FALSE)
  set.seed(6)
  n <- 20
  P <- matrix(0, n, n)
  P[cbind(1:(n - 1), 2:n)] <- 1
  P <- P + t(P)
  W <- P / rowSums(P)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- as.vector(solve(diag(n) - 0.5 * W, 1 + 0.5 * d$x + stats::rnorm(n)))
  fit <- frm(y ~ x + lltest(W), data = d, data2 = list(W = W),
             dry_run = "objective")
  aux <- fit$frame$autocor$y$aux
  rho <- c(0.4, -0.3, 0.6)
  ll <- log_lik(lla_draws(fit, list(
    b_Intercept = b0, b_x = b1, sigma = sg,
    thetaac_1 = stats::qlogis((rho - aux$lo) / (aux$hi - aux$lo)))))
  for (s in 1:3) {
    A <- diag(n) - rho[s] * W
    ref <- lla_cond(d$y, as.vector(solve(A, b0[s] + b1[s] * d$x)),
                    sg[s]^2 * solve(crossprod(A)))
    expect_lt(max(abs(ll[s, ] - ref)), 1e-10 * max(abs(ref)))
  }
})

test_that("loo(), waic() and the loo_*() functions on a covariance block
           are the loo package's on the same log_lik() matrix", {
  skip_if_not_installed("loo")
  d <- lla_data()
  n <- nrow(d)
  set.seed(5)
  k <- 40L
  base <- list(
    b_Intercept = stats::rnorm(k, 1, 0.1), b_x = stats::rnorm(k, 0.5, 0.1),
    sigma = exp(stats::rnorm(k, log(1.2), 0.05)))
  fa <- frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = "objective")
  fc <- frm(y ~ x + cosy(t, g), data = d, dry_run = "objective")
  a <- 1 / (fc$frame$autocor$y$d - 1)
  cases <- list(
    ar = lla_draws(fa, c(base, list(
      thetaac_1 = lla_pinv(stats::rnorm(k, 0.5, 0.05))))),
    cosy = lla_draws(fc, c(base, list(
      thetaac_1 = stats::qlogis((stats::rnorm(k, 0.3, 0.03) + a) /
                                  (1 + a))))))
  for (ds in cases) {
    ll <- log_lik(ds)
    # the reference: the loo package's own estimators on the same matrix,
    # with the relative efficiency of one chain, as frmtmb.sample takes it
    # for draws with no sampler behind them; loo warns that it capped an
    # effective sample size
    cap <- "The ESS has been capped"
    allow_warnings(
      reff <- loo::relative_eff(exp(ll), chain_id = rep(1L, k)), cap)
    lo <- suppressWarnings(loo(ds))
    ref <- suppressWarnings(loo::loo(ll, r_eff = reff))
    expect_equal(lo$estimates, ref$estimates,
                 tolerance = 1e3 * .Machine$double.eps)
    expect_equal(lo$pointwise[, "elpd_loo"], ref$pointwise[, "elpd_loo"],
                 tolerance = 1e3 * .Machine$double.eps)
    # a column is an observation, so loo() names no other unit
    expect_no_message(suppressWarnings(loo(ds)))
    wa <- suppressWarnings(waic(ds))
    rw <- suppressWarnings(loo::waic(ll))
    expect_equal(wa$estimates, rw$estimates,
                 tolerance = 1e3 * .Machine$double.eps)
    # loo_epred() and loo_predict() are E_loo() of the draws under the
    # same smoothed weights; posterior_predict() draws, so both sides
    # start from one random state
    ps <- suppressWarnings(loo::psis(-ll, r_eff = reff))
    le <- suppressWarnings(loo_epred(ds))
    allow_warnings(
      re <- loo::E_loo(posterior_epred(ds), ps, type = "mean")$value, cap)
    expect_equal(as.numeric(le), as.numeric(re),
                 tolerance = 1e3 * .Machine$double.eps)
    # and E_loo() is the weighted mean, written out
    w <- exp(stats::weights(ps, log = TRUE, normalize = TRUE))
    expect_lt(max(abs(as.numeric(le) - colSums(w * posterior_epred(ds)))),
              1e-10 * max(abs(re)))
    set.seed(9)
    lp <- suppressWarnings(loo_predict(ds))
    set.seed(9)
    allow_warnings(
      rp <- loo::E_loo(posterior_predict(ds), ps, type = "mean")$value, cap)
    expect_equal(as.numeric(lp), as.numeric(rp),
                 tolerance = 1e3 * .Machine$double.eps)
    expect_length(lp, n)
  }
})

## ---- brms 2.23.0's log_lik() at the same parameter vectors -----------

lla_stan_model <- local({
  mods <- new.env(parent = emptyenv())
  function(code) {
    f <- tempfile(fileext = ".stan")
    on.exit(unlink(f), add = TRUE)
    writeLines(c(code, paste0("// rstan ", utils::packageVersion("rstan"))),
               f)
    key <- unname(tools::md5sum(f))
    if (!is.null(mods[[key]])) return(mods[[key]])
    dir <- Sys.getenv("FRMTMB_STAN_CACHE", "")
    if (!nzchar(dir)) {
      dir <- testthat::test_path("..", "..", "..", "..", "dev", "stan-cache")
    }
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    path <- file.path(dir, paste0(key, ".rds"))
    mod <- if (file.exists(path)) try(readRDS(path), silent = TRUE)
    if (is.null(mod) || inherits(mod, "try-error")) {
      mod <- rstan::stan_model(model_code = code, save_dso = TRUE)
      saveRDS(mod, path)
    }
    mods[[key]] <- mod
    mod
  }
})

# a brmsfit whose k-th draw is inits[[k]]: one Fixed_param chain each
lla_brms_at <- function(bform, family, data, inits, ...) {
  prior <- brms::get_prior(bform, data = data, family = family, ...)
  prior$prior <- ""
  code <- brms::make_stancode(bform, data = data, family = family,
                              prior = prior, ...)
  sdat <- suppressMessages(brms::make_standata(bform, data = data,
                                               family = family,
                                               prior = prior, ...))
  bfit <- suppressMessages(brms::brm(bform, data = data, family = family,
                                     prior = prior, empty = TRUE, ...))
  bfit$fit <- suppressMessages(rstan::sampling(
    lla_stan_model(code), data = sdat, chains = length(inits), iter = 1,
    warmup = 0, algorithm = "Fixed_param", init = inits, refresh = 0,
    seed = 1))
  brms::log_lik(brms::rename_pars(bfit))
}

lla_expect_brms <- function(llf, llb) {
  expect_equal(dim(llf), dim(llb))
  expect_lt(max(abs(llf - llb)), 1e-10 * max(abs(llb)))
}

test_that("log_lik() is brms's under ar(), ma() and arma() covariance
           terms", {
  skip_unless_brms_fit()
  d <- lla_data()
  # brms's sigma is the innovation sd, and so is frmtmb's: the draws'
  # sigma goes to brms unchanged
  cases <- list(
    ar = list(f = "y ~ x + ar(t, g, cov = TRUE)", a = c(0.5, -0.3, 0.8),
              m = c(0, 0, 0)),
    ma = list(f = "y ~ x + ma(t, g, cov = TRUE)", a = c(0, 0, 0),
              m = c(0.4, -0.6, 0.2)),
    arma = list(f = "y ~ x + arma(t, g, cov = TRUE)", a = c(0.5, -0.2, 0.3),
                m = c(0.3, 0.5, -0.4)))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    for (fam in c("gaussian", "student")) {
      fit <- frm(stats::as.formula(cs$f), family = get(fam)(), data = d,
                 dry_run = "objective")
      vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
      th <- list()
      if (nm != "ma") th$a <- lla_pinv(cs$a)
      if (nm != "ar") th$m <- lla_pinv(cs$m)
      for (k in seq_along(th)) vals[[paste0("thetaac_", k)]] <- th[[k]]
      if (fam == "student") vals$nu <- nus
      inits <- lapply(1:3, function(s) {
        p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                  sigma = sg[s])
        if (nm != "ma") p$ar <- array(cs$a[s], 1)
        if (nm != "ar") p$ma <- array(cs$m[s], 1)
        if (fam == "student") p$nu <- nus[s]
        p
      })
      lla_expect_brms(log_lik(lla_draws(fit, vals)),
                      lla_brms_at(brms::bf(stats::as.formula(cs$f)),
                                  brms::brmsfamily(fam), d, inits))
    }
  }
  # sigma ~ x: brms's log sigma, the log innovation sd, is frmtmb's
  for (fam in c("gaussian", "student")) {
    fit <- frm(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x),
               family = get(fam)(), data = d, dry_run = "objective")
    a <- c(0.5, -0.4, 0.2)
    bs0 <- c(-0.2, 0.1, 0)
    bs1 <- c(0.3, -0.2, 0.1)
    vals <- list(b_Intercept = b0, b_x = b1, b_sigma_Intercept = bs0,
                 b_sigma_x = bs1, thetaac_1 = lla_pinv(a))
    if (fam == "student") vals$nu <- nus
    inits <- lapply(1:3, function(s) {
      p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                ar = array(a[s], 1), b_sigma = array(bs1[s], 1),
                Intercept_sigma = bs0[s] + bs1[s] * mean(d$x))
      if (fam == "student") p$nu <- nus[s]
      p
    })
    lla_expect_brms(log_lik(lla_draws(fit, vals)),
                    lla_brms_at(brms::bf(y ~ x + ar(t, g, cov = TRUE),
                                         sigma ~ x),
                                brms::brmsfamily(fam), d, inits))
  }
})

test_that("log_lik() is brms's under cosy(), unstr() and fcor()", {
  skip_unless_brms_fit()
  # cosy() and unstr() on groups with gaps: every pair shares one
  # correlation, and brms's unstr() carries the time index, so the two
  # readings of a gap agree
  d <- lla_data(gap = TRUE)
  for (fam in c("gaussian", "student")) {
    fit <- frm(y ~ x + cosy(t, g), family = get(fam)(), data = d,
               dry_run = "objective")
    a <- 1 / (fit$frame$autocor$y$d - 1)
    rho <- c(0.4, 0.1, 0.7)
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg,
                 thetaac_1 = stats::qlogis((rho + a) / (1 + a)))
    if (fam == "student") vals$nu <- nus
    inits <- lapply(1:3, function(s) {
      p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                sigma = sg[s], cosy = rho[s])
      if (fam == "student") p$nu <- nus[s]
      p
    })
    lla_expect_brms(log_lik(lla_draws(fit, vals)),
                    lla_brms_at(brms::bf(y ~ x + cosy(t, g)),
                                brms::brmsfamily(fam), d, inits))
    fit <- frm(y ~ x + unstr(t, g), family = get(fam)(), data = d,
               dry_run = "objective")
    ac <- fit$frame$autocor$y
    set.seed(40)
    th <- matrix(stats::rnorm(3 * ac$npar, 0, 0.5), 3)
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
    for (k in seq_len(ac$npar)) vals[[paste0("thetaac_", k)]] <- th[, k]
    if (fam == "student") vals$nu <- nus
    ds <- lla_draws(fit, vals)
    idx <- frmtmb.sample:::draws_par_index(ds$fit)
    inits <- lapply(1:3, function(s) {
      # brms's Lcortime is the Cholesky factor of frmtmb's matrix
      R <- frmtmb::autocor_matrix(frmtmb.sample:::draws_fit_at(ds, s, idx))
      p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                sigma = sg[s], Lcortime = t(chol(unname(R))))
      if (fam == "student") p$nu <- nus[s]
      p
    })
    lla_expect_brms(log_lik(ds),
                    lla_brms_at(brms::bf(y ~ x + unstr(t, g)),
                                brms::brmsfamily(fam), d, inits))
  }
  d <- lla_data()
  n <- nrow(d)
  V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)
  for (fam in c("gaussian", "student")) {
    fit <- frm(y ~ x + fcor(V), family = get(fam)(), data = d,
               data2 = list(V = V), dry_run = "objective")
    vals <- list(b_Intercept = b0, b_x = b1, sigma = sg)
    if (fam == "student") vals$nu <- nus
    inits <- lapply(1:3, function(s) {
      p <- list(b = array(b1[s], 1), Intercept = b0[s] + b1[s] * mean(d$x),
                sigma = sg[s])
      if (fam == "student") p$nu <- nus[s]
      p
    })
    lla_expect_brms(log_lik(lla_draws(fit, vals)),
                    lla_brms_at(brms::bf(y ~ x + fcor(V)),
                                brms::brmsfamily(fam), d, inits,
                                data2 = list(V = V)))
  }
})

test_that("log_lik() says once per call that a term takes the dense
           route", {
  frmtmb::frmtmb_local_registries()
  # above 1000 rows a term with no jacobian slot warns from core at
  # every draw; log_lik() passes the first warning on and muffles the
  # rest
  frmtmb::frmtmb_register_autocor(
    "lldense", args = "W",
    build = function(spec, rows, n_all, lookup) {
      list(W = lookup(spec$W, "lldense()")[rows, rows, drop = FALSE])
    },
    npar = 1L, start = function(aux) 0,
    natural = function(theta, aux) c(rho = tanh(theta[1])),
    whiten = function(y, mu, sigma, theta, aux) {
      A <- diag(length(y)) - tanh(theta[1]) * aux$W
      list(u = as.vector(A %*% (y - mu)) / sigma,
           logjac = as.numeric(determinant(A)$modulus) -
             length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      mu + sigma * solve(diag(nrow(aux$W)) - tanh(theta[1]) * aux$W, u)
    },
    sigma_varies = FALSE)
  set.seed(3)
  n <- 1001
  P <- matrix(0, n, n)
  P[cbind(1:(n - 1), 2:n)] <- 1
  P <- P + t(P)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- d$x + stats::rnorm(n)
  fit <- frm(y ~ x + lldense(W), data = d,
             data2 = list(W = P / rowSums(P)), dry_run = "objective")
  ds <- lla_draws(fit, list(b_Intercept = c(0, 0.1, -0.1), b_x = 1,
                            sigma = 1, thetaac_1 = c(0.2, 0.3, 0.4)))
  seen <- 0L
  ll <- withCallingHandlers(log_lik(ds),
                            frmtmb_autocor_dense_fallback = function(w) {
                              seen <<- seen + 1L
                              invokeRestart("muffleWarning")
                            })
  expect_identical(seen, 1L)
  expect_equal(dim(ll), c(3L, n))
  expect_true(all(is.finite(ll)))
})

## ---- the 0.76.0 consolidation: the round's three changes together -----
# mo()'s coding from the minimum, the residual-correlation priors (cosy
# sampled on brms's [0, 1] through one logistic, sigma's prior read at
# the innovation sd under a first-order cov = TRUE ar()) and log_lik()
# of a residual covariance term were each tested without the others.
# These blocks sample real chains and read their draws: the mo() mean is
# written from brms's definition, D * sum(simo[1:(k - min)]), and the
# reference is the dense conditional of each group at each draw, from
# the draws' own columns (dev/rel076-overlap.R is the probe).

test_that("cosy() and mo() together: log_lik() of the draws is the dense
           conditional, with brms's mo() coding and cosy on [0, 1]", {
  skip_sampler()
  set.seed(20261027)
  d <- expand.grid(t = 1:4, g = 1:30)
  d$x <- stats::rnorm(nrow(d))
  d$k <- sample(1:4, nrow(d), replace = TRUE)
  d$y <- 1 + 0.5 * d$x + c(0, 0.6, 0.9, 1.5)[d$k] +
    stats::rnorm(30, 0, 0.4)[d$g] + stats::rnorm(nrow(d), 0, 0.9)
  ds <- allow_warnings(suppressMessages(frm_sample(
    y ~ x + mo(k) + cosy(t, g), data = d, chains = 1, iter = 400,
    seed = 3, refresh = 0)), c("Effective Samples Size", "R-hat"))
  m <- as.matrix(ds$draws)
  sx <- grep("^simo_mok1\\[", colnames(m))
  # values 1 to 4: three weights, as brms codes them
  expect_length(sx, 3L)
  expect_true(all(m[, "cosy"] >= 0 & m[, "cosy"] <= 1))
  ll <- log_lik(ds)
  expect_equal(dim(ll), c(nrow(m), nrow(d)))
  for (s in seq(1, nrow(m), by = 20)) {
    cs <- c(0, cumsum(m[s, sx]))
    mu <- m[s, "b_Intercept"] + m[s, "b_x"] * d$x +
      m[s, "bsp_mok"] * length(sx) * cs[d$k - min(d$k) + 1]
    rho <- m[s, "cosy"]
    ref <- lla_grouped(d, mu, rep(m[s, "sigma"], nrow(d)),
                       function(L) ifelse(L == 0, 1, rho))
    expect_lt(max(abs(ll[s, ] - ref)), 1e-10 * max(abs(ref)))
  }
  skip_if_not_installed("loo")
  lo <- suppressWarnings(loo(ds))
  expect_true(is.finite(lo$estimates["elpd_loo", "Estimate"]))
  expect_length(lo$pointwise[, "elpd_loo"], nrow(d))
})

test_that("ar(cov = TRUE) under the default priors: log_lik() of the
           draws reads sigma as the innovation sd", {
  skip_sampler()
  # the draws' sigma column is the innovation sd, as brms's is, and the
  # conditional is of the covariance sigma^2 ar^|lag| / (1 - ar^2).
  # frmtmb 0.76.0's sigma was the marginal sd, with covariance
  # sigma^2 ar^|lag|; this block pinned that reading
  set.seed(20261027)
  d <- expand.grid(t = 1:6, g = 1:8)
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + unlist(lapply(1:8, function(g) {
    as.numeric(stats::arima.sim(list(ar = 0.9), 6, sd = 0.45))
  }))
  ds <- allow_warnings(suppressMessages(frm_sample(
    y ~ x + ar(t, g, cov = TRUE), data = d, chains = 1, iter = 400,
    seed = 3, refresh = 0)), c("Effective Samples Size", "R-hat"))
  m <- as.matrix(ds$draws)
  ll <- log_lik(ds)
  expect_equal(dim(ll), c(nrow(m), nrow(d)))
  for (s in seq(1, nrow(m), by = 20)) {
    a1 <- m[s, "ar[1]"]
    ref <- lla_grouped(d, m[s, "b_Intercept"] + m[s, "b_x"] * d$x,
                       rep(m[s, "sigma"], nrow(d)),
                       function(L) a1^L / (1 - a1^2))
    expect_lt(max(abs(ll[s, ] - ref)), 1e-10 * max(abs(ref)))
  }
  skip_if_not_installed("loo")
  lo <- suppressWarnings(loo(ds))
  expect_true(is.finite(lo$estimates["elpd_loo", "Estimate"]))
})

test_that("conditional_effects(method = 'posterior_predict') on draws draws
           at brms's sigma", {
  skip_sampler()
  # brms drops autocorrelation terms in conditional_effects()
  # (incl_autocor = FALSE), so its predictive band is iid draws at its
  # sigma, the innovation sd; frmtmb.sample's band reads the draws'
  # sigma the same way. Through frmtmb 0.76.0 the draws' sigma was the
  # marginal sd, 1 / sqrt(1 - 0.8^2) = 1.67 times the innovation sd of
  # these data (0.5)
  set.seed(20261028)
  ng <- 20
  k <- 8
  d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + unlist(lapply(seq_len(ng), function(i) {
    as.numeric(stats::arima.sim(list(ar = 0.8), k, sd = 0.5))
  }))
  ds <- allow_warnings(suppressMessages(frm_sample(
    y ~ x + ar(t, g, cov = TRUE), data = d, chains = 1, iter = 600,
    seed = 3, refresh = 0)), c("Effective Samples Size", "R-hat"))
  ce <- conditional_effects(ds, effects = "x", method = "posterior_predict",
                            seed = 1)[[1L]]
  half <- mean((ce$upper__ - ce$lower__) / (2 * stats::qnorm(0.975)))
  # the innovation sd within a quarter of the factor's distance from 1
  expect_lt(abs(half / 0.5 - 1), 0.25 * (1 / sqrt(1 - 0.8^2) - 1))
})

test_that("log_lik(newdata =) names only an alternative that answers", {
  # the 0.76.0 release review, finding 5: the refusal named
  # posterior_predict(newdata =), which a correlated residual refuses
  d <- lla_data()
  fit <- frm(y ~ x + ar(t, g, cov = TRUE), data = d, dry_run = "objective")
  ds <- lla_draws(fit, list(b_Intercept = b0, b_x = b1, sigma = sg,
                            thetaac_1 = lla_pinv(c(0.5, -0.3, 0.8))))
  e <- tryCatch(log_lik(ds, newdata = d), error = identity)
  expect_s3_class(e, "frmtmb_error")
  expect_match(conditionMessage(e), "does not take newdata", fixed = TRUE)
  expect_match(conditionMessage(e), "Use posterior_epred(newdata =) for",
               fixed = TRUE)
  expect_match(conditionMessage(e), "is refused for this model too",
               fixed = TRUE)
  # the alternative it rules out is refused indeed
  pe <- tryCatch(posterior_predict(ds, newdata = d), error = identity)
  expect_s3_class(pe, "error")
  expect_match(conditionMessage(pe), "its draws are structured",
               fixed = TRUE)
  # the guard absent: without the residual term both are named
  f2 <- frm(y ~ x, data = d, dry_run = "objective")
  ds2 <- lla_draws(f2, list(b_Intercept = b0, b_x = b1, sigma = sg))
  e2 <- tryCatch(log_lik(ds2, newdata = d), error = identity)
  expect_s3_class(e2, "frmtmb_error")
  expect_match(conditionMessage(e2),
               "Use posterior_epred(newdata =) or posterior_predict(",
               fixed = TRUE)
})

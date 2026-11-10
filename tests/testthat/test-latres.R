# brms's latent residuals: ar(), ma(), arma(), cosy() and unstr() on a family
# with no residual y - mu. test-latres-brms.R holds the identities with brms's
# Stan program; this file holds what needs no brms: the law of the latents, an
# exact marginal likelihood, the refusals, and the readers.

lr_data <- function(seed = 167, G = 12L, K = 5L, sd = 1) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  e <- as.vector(replicate(G, stats::arima.sim(list(ar = 0.6), K, sd = sd)))
  eta <- 0.6 + 0.4 * d$x + e
  d$cnt <- stats::rpois(n, exp(eta))
  d$bin <- stats::rbinom(n, 1, stats::plogis(eta - 0.6))
  d$pos <- exp(eta + stats::rnorm(n, 0, 0.4))
  d$gam <- stats::rgamma(n, 4, 4 / exp(eta))
  d
}

lr_terms <- c(ar1 = "ar(t, g, cov = TRUE)", ma1 = "ma(t, g, cov = TRUE)",
              arma1 = "arma(t, g, cov = TRUE)", condar = "ar(t, g, p = 2)",
              cosy = "cosy(t, g)", unstr = "unstr(t, g)")

lr_fit <- function(resp, term, family, data, ...) {
  allow_warnings(
    frm(bf(stats::as.formula(paste(resp, "~ x +", term))), family = family,
        data = data, ...),
    c("Optimizer did not report convergence",
      "Standard errors are not available",
      "Large maximum absolute gradient",
      "non-finite gradient",
      "The Laplace approximation of",
      "The dispersion `"))
}

# A form frm() refuses and frm_sample() reads: the objective and the methods
# hold it, so their tests build it past the refusal.
lr_fit_unread <- function(...) {
  local_mocked_bindings(latres_refuse_laplace = function(spec, frame) NULL,
                        latres_refuse_margin = function(fit) NULL,
                        .package = "frmtmb")
  lr_fit(...)
}

# C at unit scale, from the definitions of brms's fun_cholesky_cor_*.stan,
# restated here so the chain below is checked against the model and not
# against itself
lr_target_cov <- function(kind, k, ar = NA, ma = NA, rho = NA) {
  lag <- abs(outer(seq_len(k), seq_len(k), "-"))
  switch(kind,
    ar1 = ar^lag / (1 - ar^2),
    ma1 = (1 + ma^2) * (lag == 0) + ma * (lag == 1),
    arma1 = {
      g1 <- (1 + ar * ma) * (ar + ma)
      m <- g1 * ar^pmax(lag - 1, 0)
      diag(m) <- 1 + ma^2 + 2 * ar * ma
      m / (1 - ar^2)
    },
    cosy = rho + (1 - rho) * (lag == 0))
}

# the block's own map and density as numbers: Cov(err) = A Q^-1 A', at a theta
# the test chooses, where neither factor is near singular
lr_block_cov <- function(fit, th) {
  bk <- frame_latres(fit$frame)[[1L]]
  pj <- latres_jacobian(bk, th)
  A <- Matrix::sparseMatrix(i = pj$i, j = pj$j, x = pj$x,
                            dims = c(bk$aux$n, bk$aux$m))
  # the precision read back from the density, by its value at unit vectors:
  # a gaussian density has no other way to answer
  m <- bk$aux$m
  f0 <- latres_lpdf(numeric(m), th, bk)
  Q <- matrix(0, m, m)
  for (i in seq_len(m)) {
    ei <- numeric(m)
    ei[i] <- 1
    Q[i, i] <- -2 * (latres_lpdf(ei, th, bk) - f0)
  }
  for (i in seq_len(m)) {
    for (j in seq_len(i - 1L)) {
      ei <- numeric(m)
      ei[c(i, j)] <- 1
      Q[i, j] <- Q[j, i] <-
        -(latres_lpdf(ei, th, bk) - f0) - (Q[i, i] + Q[j, j]) / 2
    }
  }
  as.matrix(A %*% solve(Q) %*% Matrix::t(A))
}

test_that("each chain has brms's covariance, sderr squared times C", {
  d <- lr_data(G = 4L, K = 5L)
  # ragged groups, so a leading block of C is what a short group reads
  d <- d[-c(2L, 8L, 9L), ]
  for (kind in c("ar1", "ma1", "arma1", "cosy")) {
    fit <- lr_fit("cnt", lr_terms[[kind]], poisson(), d)
    bk <- frame_latres(fit$frame)[[1L]]
    th <- c(log(0.7), c(0.6, -0.3)[seq_len(length(bk$theta_idx) - 1L)])
    cf <- latres_coefs(th[-1L], bk$aux)
    V <- lr_block_cov(fit, th)
    for (g in levels(d$g)) {
      rows <- which(d$g == g)
      tgt <- exp(2 * th[1L]) * lr_target_cov(
        kind, length(rows), ar = cf$ar %||% NA, ma = cf$ma %||% NA,
        rho = cf$rho %||% NA)
      expect_lt(max(abs(V[rows, rows] - tgt)), 1e-10 * max(abs(tgt)))
    }
    # and groups do not covary
    other <- outer(d$g, d$g, "!=")
    expect_lt(max(abs(V[other])), 1e-12 * max(abs(V)))
  }
})

test_that("unstr() reads the time levels of each row", {
  d <- lr_data(G = 6L, K = 4L)
  d <- d[-c(2L, 7L, 8L), ]
  fit <- lr_fit("cnt", lr_terms[["unstr"]], poisson(), d)
  bk <- frame_latres(fit$frame)[[1L]]
  th <- c(log(0.7), 0.1 * seq_len(length(bk$theta_idx) - 1L) - 0.3)
  C <- us_chol_cor(th[-1L], 4L)
  V <- lr_block_cov(fit, th)
  rows <- which(d$g == "1")
  tgt <- exp(2 * th[1L]) * C[d$t[rows], d$t[rows]]
  expect_lt(max(abs(V[rows, rows] - tgt)), 1e-10 * max(abs(tgt)))
})

test_that("the cov = FALSE chain is brms's moving average of iid latents", {
  d <- lr_data(G = 3L, K = 6L)
  fit <- lr_fit("cnt", "ar(t, g, p = 2)", poisson(), d)
  bk <- frame_latres(fit$frame)[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  ar <- latres_coefs(th[-1L], bk$aux)$ar
  expect_true(all(abs(ar) < 1))
  set.seed(1)
  z <- stats::rnorm(bk$aux$m)
  err <- latres_expand(z, bk, th)
  # brms's model block, written out row by row
  s <- exp(th[1L])
  by_hand <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    by_hand[i] <- s * z[i]
    if (d$t[i] > 1L) by_hand[i] <- by_hand[i] + s * ar[1L] * z[i - 1L]
    if (d$t[i] > 2L) by_hand[i] <- by_hand[i] + s * ar[2L] * z[i - 2L]
  }
  expect_equal(err, by_hand, tolerance = 1e-14)
})

test_that("on a gaussian inner problem the fit is the exact marginal", {
  # lognormal: log y | err is normal, so the Laplace value IS the marginal
  # likelihood of log y ~ N(mu, sigma^2 I + sderr^2 C), which a dense density
  # gives with no latent at all. At a chosen point of the outer parameters and
  # not at the fit's, which on these few rows runs off toward a bound where
  # 1 - ar^2 in the reference has no digits left
  d <- lr_data(G = 6L, K = 4L)
  for (kind in c("ar1", "ma1", "arma1", "cosy")) {
    fit <- lr_fit("pos", lr_terms[[kind]], lognormal(), d)
    bk <- frame_latres(fit$frame)[[1L]]
    par <- fit$opt$par
    th <- c(log(0.7), c(0.6, -0.3)[seq_len(length(bk$theta_idx) - 1L)])
    par[names(par) == "theta"] <- th
    par[names(par) == "betad"] <- log(0.5)
    tape <- -fit$obj$fn(par)
    cf <- latres_coefs(th[-1L], bk$aux)
    s2 <- exp(2 * th[1L])
    sigma <- 0.5
    mu <- as.numeric(cbind(1, d$x) %*% par[names(par) == "beta"])
    ll <- 0
    C <- s2 * lr_target_cov(kind, 4L, ar = cf$ar %||% NA, ma = cf$ma %||% NA,
                            rho = cf$rho %||% NA) + diag(sigma^2, 4L)
    R <- chol(C)
    for (g in levels(d$g)) {
      r <- which(d$g == g)
      res <- backsolve(R, log(d$pos[r]) - mu[r], transpose = TRUE)
      ll <- ll - sum(log(diag(R))) - 2 * log(2 * pi) - 0.5 * sum(res^2) -
        sum(log(d$pos[r]))
    }
    expect_lt(abs(tape - ll), 1e-8 * abs(ll))
  }
})

test_that("the readers agree with the tape: rows plus latents is the joint", {
  d <- lr_data(G = 12L, K = 5L, sd = 0.5)
  cases <- list(c("cnt", "ar(t, g, cov = TRUE)", "poisson"),
                c("cnt", "ar(t, g)", "poisson"),
                c("gam", "ar(t, g)", "Gamma_log"),
                c("gam", "ar(t, g, p = 2)", "Gamma_log"),
                c("bin", "ar(t, g)", "probit"),
                c("cnt", "cosy(t, g)", "poisson"))
  for (cs in cases) {
    fam <- switch(cs[[3]], poisson = poisson(),
                  Gamma_log = Gamma(link = "log"),
                  probit = bernoulli(link = "probit"))
    fit <- lr_fit_unread(cs[[1]], cs[[2]], fam, d)
    bk <- frame_latres(fit$frame)[[1L]]
    th <- fit$estimates$theta[bk$theta_idx]
    ep <- arma_cond_dpars(fit, eval_dpars(fit))
    y <- fit$frame$y[[1L]]
    rows <- row_lpdf(fit$spec$responses[[1L]]$family, y, y, ep[[1L]],
                     fit$frame$aterm_values[[1L]], fit_extras(fit))
    joint <- sum(rows) + latres_lpdf(fit$estimates$b[bk$b_idx], th, bk)
    env <- fit$obj$env
    tape <- -env$f(env$last.par.best)
    expect_lt(abs(joint - tape), 1e-8 * max(1, abs(tape)))
  }
})

test_that("the latents leave the inner Hessian sparse", {
  # a count, not a clock: stored entries per row of the upper triangle, which
  # a dense group of k rows would make k(k + 1) / 2 per k rows
  d <- lr_data(G = 10L, K = 20L)
  n <- nrow(d)
  nnz <- function(term) {
    fit <- lr_fit("cnt", term, poisson(), d)
    env <- fit$obj$env
    H <- env$spHess(env$last.par.best, random = TRUE)
    length(H@x)
  }
  dense <- 10L * 20L * 21L / 2L
  # a chain's tridiagonal is at most two entries a row in the upper triangle,
  # which is itself at most twice the stored entries of a lower-only store
  for (term in lr_terms[c("ar1", "ma1", "arma1", "cosy")]) {
    expect_lte(nnz(term), 3L * (n + 10L))
    expect_lt(nnz(term), dense / 3L)
  }
  expect_lte(nnz("ar(t, g, p = 2)"), 5L * n)
})

test_that("brms's refusals are frm()'s, in brms's words", {
  d <- lr_data(G = 6L, K = 4L)
  fe <- function(term, why, fam = poisson()) {
    expect_error(frm(bf(stats::as.formula(paste("cnt ~ x +", term))),
                     family = fam, data = d), why)
  }
  ma_words <- "Please set cov = TRUE when modeling MA structures for this"
  fe("ma(t, g)", ma_words)
  fe("arma(t, g)", ma_words)
  fe("ma(t, g, q = 2)", ma_words)
  ord_words <- paste0("Covariance formulation of ARMA structures is only ",
                      "possible for effects of maximal order one")
  fe("ar(t, g, p = 2, cov = TRUE)", ord_words)
  fe("ma(t, g, q = 2, cov = TRUE)", ord_words)
  fe("arma(t, g, p = 2, q = 1, cov = TRUE)", ord_words)
  V <- diag(nrow(d))
  expect_error(frm(bf(cnt ~ x + fcor(V)), family = poisson(), data = d,
                   data2 = list(V = V)),
               "FCOR terms are not implemented for this family")
  expect_error(
    frm(bf(cnt ~ a * x, a ~ 1, nl = TRUE, autocor = ~ ar(t, g, cov = TRUE)),
        family = poisson(), data = d),
    "non-linear|nonlinear")
  d$cat <- factor(c("a", "b", "c")[(d$cnt %% 3L) + 1L])
  expect_error(frm(bf(cat ~ x + ar(t, g, cov = TRUE)), family = categorical(),
                   data = d), "can only be specified on 'mu'")
  expect_error(frm(bf(cnt ~ x + ar(t, g, cov = TRUE)),
                   family = mixture(poisson(), poisson()), data = d),
               "can only be specified on 'mu'")
})

test_that("the guard's absent case: gaussian and student keep their forms", {
  d <- lr_data(G = 6L, K = 4L)
  for (fam in list(gaussian(), student())) {
    fit <- suppressWarnings(frm(bf(pos ~ x + ar(t, g, cov = TRUE)),
                                family = fam, data = d))
    expect_length(frame_latres(fit$frame), 0L)
    expect_false(is.null(fit$frame$autocor[[1L]]))
    expect_null(fit$frame$autocor[[1L]][["latent"]])
    expect_true("sigma" %in% names(fit$frame$linpreds[[2L]][["dpar"]]) ||
                  length(fit$estimates$thetaac) == 1L)
  }
  # a higher order stays frmtmb's documented divergence on these two
  fit <- suppressWarnings(frm(bf(pos ~ x + ar(t, g, p = 2, cov = TRUE)),
                              family = gaussian(), data = d))
  expect_length(fit$estimates$thetaac, 2L)
  # and a latent term no longer needs the time levels to be dense
  fit <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  expect_length(fit$estimates$thetaac, 0L)
  expect_length(frame_latres(fit$frame), 1L)
})

test_that("fitted() carries the latent residual, under re_formula = NA too", {
  d <- lr_data()
  fit <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  bk <- frame_latres(fit$frame)[[1L]]
  err <- as.numeric(coef_b(fit)[bk$c_idx])
  X <- cbind(1, d$x)
  eta <- as.numeric(X %*% fit$estimates$beta) + err
  expect_equal(as.numeric(fitted(fit)[, 1L]), exp(eta), tolerance = 1e-12)
  expect_equal(as.numeric(fitted(fit, re_formula = NA)[, 1L]), exp(eta),
               tolerance = 1e-12)
  expect_equal(as.numeric(frm_linpred(fit)), eta, tolerance = 1e-12)
  expect_equal(as.numeric(residuals(fit, type = "response")[, 1L]),
               d$cnt - exp(eta), tolerance = 1e-12)
  # a standard error reads the latent's own uncertainty, through the chain
  se <- frm_linpred(fit, se.fit = TRUE)$se.fit
  expect_true(all(is.finite(se) & se > 0))
  # it is an autocorrelation term and not a group-level one
  expect_length(ranef(fit), 0L)
  expect_null(ngrps(fit))
})

test_that("the post-link AR part of a cov = FALSE term reaches the mean", {
  d <- lr_data()
  fit <- lr_fit_unread("gam", "ar(t, g)", Gamma(link = "log"), d)
  bk <- frame_latres(fit$frame)[[1L]]
  expect_true(isTRUE(bk$aux$on_response))
  th <- fit$estimates$theta[bk$theta_idx]
  b <- fit$estimates$b[bk$b_idx]
  ar <- latres_coefs(th[-1L], bk$aux)$ar
  s <- exp(th[1L])
  lag <- c(NA, b[-length(b)])
  lag[d$t == 1L] <- 0
  mu <- exp(as.numeric(cbind(1, d$x) %*% fit$estimates$beta) + s * b) +
    s * ar * lag
  expect_equal(as.numeric(fitted(fit)[, 1L]), mu, tolerance = 1e-10)
  expect_equal(as.numeric(frm_linpred(fit, type = "response")), mu,
               tolerance = 1e-10)
  # a family whose density reads the link scale keeps it there
  fit2 <- lr_fit("cnt", "ar(t, g)", poisson(), d)
  expect_false(isTRUE(frame_latres(fit2$frame)[[1L]]$aux$on_response))
})

test_that("new rows draw latents under cov = TRUE, and cov = FALSE refuses", {
  d <- lr_data()
  fit <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  nd <- d[d$g == "1", ]
  nd$g <- factor("new")
  eta <- as.numeric(cbind(1, nd$x) %*% fit$estimates$beta)
  bk <- frame_latres(fit$frame)[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  ar <- latres_coefs(th[-1L], bk$aux)$ar
  v <- exp(2 * th[1L]) / (1 - ar^2)
  # the expected value on new rows is over the latent, the predictor is at 0
  expect_equal(as.numeric(fitted(fit, newdata = nd)[, 1L]), exp(eta + v / 2),
               tolerance = 1e-10)
  expect_equal(as.numeric(frm_linpred(fit, newdata = nd, type = "response")),
               exp(eta), tolerance = 1e-12)
  # a predictive draw reads it too, so its log mean carries sderr^2 C / 2
  set.seed(3)
  sims <- simulate(fit, nsim = 4000L, newdata = nd)
  m <- rowMeans(as.matrix(sims))
  se <- sqrt(apply(as.matrix(sims), 1L, stats::var) / 4000)
  expect_lt(max(abs(m - exp(eta + v / 2)) / se), 4.5)
  expect_gt(max(exp(eta + v / 2) / exp(eta)), 1.0001)
  set.seed(3)
  pp <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000L,
                propagate_error = FALSE)
  se2 <- sqrt(apply(pp, 2L, stats::var) / 4000)
  expect_lt(max(abs(colMeans(pp) - exp(eta + v / 2)) / se2), 4.5)
  # brms's words
  fitf <- lr_fit("cnt", "ar(t, g)", poisson(), d)
  w <- "Cannot predict new latent residuals when using cov = FALSE"
  expect_error(fitted(fitf, newdata = nd), w)
  expect_error(predict(fitf, newdata = nd), w)
  expect_error(simulate(fitf, newdata = nd), w)
  expect_error(frm_linpred(fitf, newdata = nd), w)
  # a grid drops the term, as brms's incl_autocor = FALSE does
  ce <- conditional_effects(fitf, "x")
  expect_true(is.data.frame(ce[[1L]]))
})

test_that("autocor_matrix() gives the latent correlation", {
  d <- lr_data()
  fit <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  bk <- frame_latres(fit$frame)[[1L]]
  ar <- latres_coefs(fit$estimates$theta[bk$theta_idx][-1L], bk$aux)$ar
  R <- autocor_matrix(fit)
  expect_equal(unname(R[1L, 3L]), ar^2, tolerance = 1e-10)
  expect_equal(unname(diag(R)), rep(1, 5L), tolerance = 1e-10)
  fc <- lr_fit("cnt", "cosy(t, g)", poisson(), d)
  bc <- frame_latres(fc$frame)[[1L]]
  rho <- latres_coefs(fc$estimates$theta[bc$theta_idx][-1L], bc$aux)$rho
  expect_equal(unname(autocor_matrix(fc)[2L, 5L]), rho, tolerance = 1e-12)
  # a moving average of independent latents, by lag
  fu <- lr_fit("cnt", "ar(t, g)", poisson(), d)
  bu <- frame_latres(fu$frame)[[1L]]
  a1 <- latres_coefs(fu$estimates$theta[bu$theta_idx][-1L], bu$aux)$ar
  Ru <- autocor_matrix(fu)
  expect_equal(unname(Ru[1L, 2L]), a1 / (1 + a1^2), tolerance = 1e-10)
  expect_equal(unname(Ru[1L, 3L]), 0, tolerance = 1e-12)
})

test_that("variables, the summary and the intervals carry brms's names", {
  d <- lr_data()
  fit <- lr_fit("cnt", "arma(t, g, cov = TRUE)", poisson(), d)
  # the row the equivalence survey's pair p167 read refused through 0.78.0
  expect_identical(frm_compat("poisson", "arma()")$status, "conditional")
  expect_true(all(c("ar[1]", "ma[1]", "sderr") %in% variables(fit)))
  s <- summary(fit)
  expect_setequal(rownames(s$cor_pars), c("ar[1]", "ma[1]", "sderr"))
  ci <- confint_varcorr(fit)
  expect_true(all(c("sderr", "ar[1]", "ma[1]") %in% ci$term))
  expect_equal(ci$estimate[ci$term == "sderr"],
               exp(fit$estimates$theta[1L]), tolerance = 1e-10)
})

test_that("the prior classes are brms's, on their natural scales", {
  d <- lr_data()
  p <- default_prior(bf(cnt ~ x + ar(t, g, cov = TRUE)), family = poisson(),
                     data = d)
  expect_true(all(c("ar", "sderr") %in% p$class))
  expect_equal(p$lb[p$class == "sderr"], 0)
  expect_false("sd" %in% p$class)
  pc <- default_prior(bf(cnt ~ x + cosy(t, g)), family = poisson(), data = d)
  expect_equal(pc$lb[pc$class == "cosy"], 0)
  expect_equal(pc$ub[pc$class == "cosy"], 1)
  # a prior on sderr is that density at the natural parameter, carried onto
  # its log with the log map's Jacobian, as class sd is
  f0 <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  f1 <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d,
               prior = set_prior("normal(0, 1)", class = "sderr"))
  par <- f0$opt$par
  s <- exp(par[names(par) == "theta"][1L])
  gap <- f1$obj$fn(par) - f0$obj$fn(par)
  want <- -(stats::dnorm(s, 0, 1, log = TRUE) + log(s))
  expect_lt(abs(gap - want), 1e-8 * max(1, abs(want)))
  f2 <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d,
               prior = set_prior("normal(0, 0.1)", class = "ar"))
  expect_lt(abs(f2$opt$objective - f0$opt$objective), abs(f0$opt$objective))
  # a class that names another thing is refused, not applied to this block
  expect_error(lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d,
                      prior = set_prior("normal(0, 1)", class = "sd")),
               "No random-effect|no random-effect|standard deviations")
  expect_error(frm(bf(pos ~ x + ar(t, g, cov = TRUE)), family = gaussian(),
                   data = d,
                   prior = set_prior("normal(0, 1)", class = "sderr")),
               "No latent residual")
  # one parameter takes one prior
  expect_error(lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d,
                      prior = set_prior("normal(0, 1)", class = "theta",
                                        coef = "theta_2")),
               "class = \"ar\"")
})

test_that("bounds on a cov = FALSE coefficient keep it in brms's box", {
  d <- lr_data()
  fit <- lr_fit("cnt", "ar(t, g, p = 2)", poisson(), d,
                prior = set_prior("", class = "ar", lb = 0, ub = 0.5))
  bk <- frame_latres(fit$frame)[[1L]]
  ar <- latres_coefs(fit$estimates$theta[bk$theta_idx][-1L], bk$aux)$ar
  expect_true(all(ar >= 0 & ar <= 0.5))
})

test_that("a multivariate model takes one term per response", {
  d <- lr_data()
  d$cnt2 <- rev(d$cnt)
  fit <- allow_warnings(
    frm(bf(cnt ~ x + ar(t, g, cov = TRUE)) + bf(cnt2 ~ x + cosy(t, g)),
        family = poisson(), data = d),
    c("The Laplace approximation",
      "Optimizer did not report convergence",
      "Standard errors are not available",
      "non-finite gradient"))
  expect_named(frame_latres(fit$frame), c("cnt", "cnt2"), ignore.order = TRUE)
  expect_true(all(c("ar_cnt[1]", "sderr_cnt", "cosy_cnt2", "sderr_cnt2") %in%
                    variables(fit)))
  ep <- fitted(fit)
  expect_equal(dim(ep)[1L], nrow(d))
  p <- default_prior(bf(cnt ~ x + ar(t, g, cov = TRUE)) +
                       bf(cnt2 ~ x + cosy(t, g)), family = poisson(), data = d)
  expect_true(all(c("cnt", "cnt2") %in% p$resp[p$class == "sderr"]))
})

test_that("the structure's name is no formula spelling", {
  d <- lr_data(G = 6L, K = 4L)
  expect_error(frm(bf(cnt ~ x + latres(1 | g)), family = poisson(), data = d))
})

test_that("the rows may come in any order", {
  d <- lr_data()
  fit <- lr_fit("cnt", "arma(t, g, cov = TRUE)", poisson(), d)
  set.seed(9)
  shuf <- d[sample(nrow(d)), ]
  fit2 <- lr_fit("cnt", "arma(t, g, cov = TRUE)", poisson(), shuf)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(fit2))),
            1e-6 * abs(as.numeric(logLik(fit))))
  # each fitted row is its own row's
  o <- order(as.integer(shuf$g), shuf$t)
  expect_lt(max(abs(fitted(fit)[, 1L] - fitted(fit2)[o, 1L]) /
                  fitted(fit)[, 1L]), 1e-3)
})

test_that("the delta method reads the chain's own uncertainty", {
  d <- lr_data()
  for (term in c("ar(t, g, cov = TRUE)", "ma(t, g, cov = TRUE)",
                 "cosy(t, g)", "ar(t, g)")) {
    fit <- lr_fit("cnt", term, poisson(), d)
    lp <- fit$frame$linpreds[[1L]]
    ana <- frm_linpred(fit, se.fit = TRUE)$se.fit
    # central differences of the same predictor over every coefficient
    f <- function(x) frm_linpred(x)
    fd <- autocor_cond_fd_se(fit, f, TRUE)
    ok <- is.finite(ana) & is.finite(fd) & fd > 0
    expect_gt(sum(ok), nrow(d) / 2)
    expect_lt(max(abs(ana[ok] / fd[ok] - 1)), 0.02)
  }
})

test_that("what a latent residual cannot be combined with says so", {
  d <- lr_data()
  d$w2 <- 2
  # one latent per row, under weights above one: no maximum, as (1 | obs)
  expect_error(frm(bf(gam | weights(w2) ~ x + ar(t, g, cov = TRUE)),
                   family = Gamma(link = "log"), data = d), "weights")
  expect_error(lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d,
                      quadrature = TRUE), "quadrature")
  expect_error(lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d,
                      importance = 100), "cannot be combined with a latent")
  # and a term the response has two of is one term too many
  expect_error(frm(bf(cnt ~ x + ar(t, g, cov = TRUE) + cosy(t, g)),
                   family = poisson(), data = d), "residual correlation terms")
})

test_that("a prior on cortime and cosy reach their parameters", {
  d <- lr_data(G = 12L, K = 3L)
  f0 <- lr_fit("cnt", "unstr(t, g)", poisson(), d)
  f1 <- lr_fit("cnt", "unstr(t, g)", poisson(), d,
               prior = set_prior("lkj(4)", class = "cortime"))
  gap <- f1$obj$fn(f0$opt$par) - f0$obj$fn(f0$opt$par)
  expect_true(is.finite(gap) && abs(gap) > 0)
  g0 <- lr_fit("cnt", "cosy(t, g)", poisson(), d)
  g1 <- lr_fit("cnt", "cosy(t, g)", poisson(), d,
               prior = set_prior("beta(2, 2)", class = "cosy"))
  par <- g0$opt$par
  th <- par[names(par) == "theta"][2L]
  rho <- stats::plogis(th)
  # the density at the natural parameter, with the logistic map's Jacobian
  want <- -(stats::dbeta(rho, 2, 2, log = TRUE) + log(rho) + log(1 - rho))
  gapc <- g1$obj$fn(par) - g0$obj$fn(par)
  expect_lt(abs(gapc - want), 1e-8 * max(1, abs(want)))
})

test_that("a mean a latent pushes off the family's support has no density", {
  # brms's Stan program rejects it; on the tape a zero count folds log(mu) away,
  # so without the guard a row with y = 0 was a finite density at a negative
  # rate, and the likelihood had no maximum (rtmb-pitfalls item 22)
  d <- lr_data(sd = 0.5)
  fit <- allow_warnings(
    lr_fit_unread("cnt", "ar(t, g, cov = TRUE)", poisson(link = "identity"), d),
    "NaNs produced")
  p <- fit$obj$env$last.par.best
  bi <- fit$obj$env$random
  zero <- which(d$cnt == 0L)
  skip_if(!length(zero))
  bk <- frame_latres(fit$frame)[[1L]]
  p[bi[bk$b_idx[zero[1L]]]] <- -1e3
  expect_true(is.nan(fit$obj$env$f(p)))
})

test_that("ar(cov = TRUE) is the same model as an ar1() random effect", {
  # two routes through frmtmb, one per row and one over the time factor, whose
  # charts differ and whose marginal likelihood cannot
  d <- lr_data(sd = 0.5)
  d$ft <- factor(d$t)
  a <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  b <- allow_warnings(frm(bf(cnt ~ x + ar1(ft + 0 | g)), family = poisson(),
                          data = d),
                      c("Optimizer did not report convergence",
                        "Standard errors are not available",
                        "Large maximum absolute gradient"))
  la <- as.numeric(logLik(a))
  expect_lt(abs(la - as.numeric(logLik(b))), 1e-6 * abs(la))
})

lr_series <- function(seed, G, K, b0, s = 1, ar = 0.6, tr = 1L, y = "pois") {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  d$x <- stats::rnorm(nrow(d))
  e <- as.vector(replicate(G, stats::arima.sim(list(ar = ar), K, sd = s)))
  eta <- b0 + 0.5 * d$x + e
  d$tr <- tr
  d$y <- switch(y, pois = stats::rpois(nrow(d), exp(eta)),
                bin = stats::rbinom(nrow(d), tr, stats::plogis(eta)))
  d
}

lr_laplace_warned <- function(expr) {
  said <- character()
  withCallingHandlers(
    expr,
    frmtmb_latres_laplace = function(w) {
      said <<- c(said, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
  said
}

lr_warned_on <- function(formula, family, data, ...) {
  length(lr_laplace_warned(frm(bf(formula), family = family, data = data, ...)))
}

test_that("frm() warns where a row tells its latent too little", {
  f <- y ~ x + ar(t, g, cov = TRUE)
  ft <- y | trials(tr) ~ x + ar(t, g, cov = TRUE)
  bern <- lr_series(3L, 20L, 10L, 0.2, y = "bin")
  said <- lr_laplace_warned(frm(bf(f), family = bernoulli(), data = bern))
  expect_length(said, 1L)
  expect_match(said, "frm_sample()", fixed = TRUE)
  expect_match(said, "ar(t, g, cov = TRUE)", fixed = TRUE)
  # above the row limit: three trials has a row's share of a poisson of mean 1
  # that reads true, so the condition reads the trials
  b3 <- lr_series(3L, 15L, 10L, 0.2, 0.8, tr = 3L, y = "bin")
  said <- withr::with_options(list(frmtmb.latres_se_rows = 10L),
    lr_laplace_warned(frm(bf(ft), family = binomial(), data = b3)))
  expect_length(said, 1L)
  expect_match(said, "3 trials", fixed = TRUE)
  # a low-count poisson overstates sderr, and is the same condition
  lo <- lr_series(3L, 20L, 15L, -2.5, 0.8)
  expect_equal(lr_warned_on(f, poisson(), lo), 1L)
  # the fine designs stay quiet: a poisson of mean 1, 10 trials
  p1 <- lr_series(3L, 10L, 15L, -0.2, 0.5)
  expect_equal(lr_warned_on(f, poisson(), p1), 0L)
  b10 <- lr_series(3L, 15L, 10L, 0.2, 0.8, tr = 10L, y = "bin")
  expect_equal(lr_warned_on(ft, binomial(), b10), 0L)
  # the row's information is the second derivative the family's density has:
  # exp(eta) on poisson, at the fit's own mean
  fit <- lr_fit("y", "ar(t, g, cov = TRUE)", poisson(), p1)
  info <- latres_row_info(fit, frame_latres(fit$frame)[[1L]])
  expect_equal(info, as.numeric(eval_dpars(fit)[[1L]][["mu"]]),
               tolerance = 1e-6)
})

test_that("the warning does not read units, and silences by its class", {
  f <- y ~ x + ar(t, g, cov = TRUE)
  set.seed(3L)
  d <- expand.grid(t = 1:8, g = factor(1:25))
  d <- d[order(d$g, d$t), ]
  d$x <- stats::rnorm(200L)
  e <- as.vector(replicate(25, stats::arima.sim(list(ar = 0.6), 8, sd = 0.5)))
  # Laplace is exact here, whatever the rows tell the latent
  d$y <- stats::rlnorm(200L, e, 2)
  fl <- suppressMessages(suppressWarnings(
    frm(bf(f), family = lognormal(), data = d)))
  expect_equal(lr_warned_on(f, lognormal(), d), 0L)
  # the stand-aside by construction, above the row limit: a thin row with no
  # third derivative reads no reason, and the same row with one does, so no
  # other condition answers for it
  withr::local_options(frmtmb.latres_se_rows = 10L)
  n <- nrow(d)
  bk <- frame_latres(fl$frame)[[1L]]
  local({
    local_mocked_bindings(latres_row_derivs = function(fit, bk) {
      list(d2 = rep(0.01, n), d3 = numeric(n), var = 1)
    }, .package = "frmtmb")
    expect_null(latres_laplace_reason(fl, bk))
  })
  local({
    local_mocked_bindings(latres_row_derivs = function(fit, bk) {
      list(d2 = rep(0.01, n), d3 = rep(0.01, n), var = 1)
    }, .package = "frmtmb")
    expect_match(latres_laplace_reason(fl, bk), "of its latent", fixed = TRUE)
  })
  # one model in two units reads one way
  z <- 0.5 * stats::rnorm(200L)
  d$y <- e + z
  w1 <- lr_warned_on(f, skew_normal(), d)
  d$y <- 10 * (e + z)
  expect_equal(lr_warned_on(f, skew_normal(), d), w1)
  # it is about the estimates, so the SE check's setting does not reach it
  bern <- lr_series(3L, 20L, 10L, 0.2, y = "bin")
  expect_equal(lr_warned_on(f, bernoulli(), bern,
                            control = frmtmb_control(check_se = "ignore")), 1L)
  seen <- character()
  withCallingHandlers(
    suppressWarnings(frm(bf(f), family = bernoulli(), data = bern),
                     classes = "frmtmb_latres_laplace"),
    frmtmb_latres_laplace = function(w) seen <<- c(seen, conditionMessage(w)),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
  expect_length(seen, 0L)
})

test_that("a dispersion that trades off against the latent warns", {
  f <- y ~ x + ar(t, g, cov = TRUE)
  ft <- y | trials(tr) ~ x + ar(t, g, cov = TRUE)
  dat <- function(seed, K = 8L) {
    set.seed(seed)
    d <- expand.grid(t = seq_len(K), g = factor(1:25))
    d <- d[order(d$g, d$t), ]
    n <- nrow(d)
    d$x <- stats::rnorm(n)
    d$e <- as.vector(replicate(25, stats::arima.sim(list(ar = 0.6), K,
                                                    sd = 0.5)))
    d$m <- stats::plogis(0.4 * d$x + d$e)
    d$tr <- 10L
    d
  }
  d <- dat(3001L)
  d$y <- pmin(pmax(stats::rbeta(nrow(d), d$m * 4, (1 - d$m) * 4), 1e-6),
              1 - 1e-6)
  expect_equal(lr_warned_on(f, Beta(), d), 1L)
  d <- dat(2003L, 6L)
  p <- stats::rbeta(nrow(d), d$m * 10, (1 - d$m) * 10)
  d$y <- stats::rbinom(nrow(d), 10L, p)
  expect_equal(lr_warned_on(ft, beta_binomial(), d), 1L)
  d <- dat(2003L, 6L)
  v <- stats::rbeta(nrow(d), d$m * 10, (1 - d$m) * 10)
  u <- stats::runif(nrow(d))
  d$y <- ifelse(u < 0.1, stats::rbinom(nrow(d), 1L, 0.5), v)
  expect_equal(lr_warned_on(f, zero_one_inflated_beta(), d), 1L)
  d <- dat(3002L)
  d$y <- stats::rnbinom(nrow(d), size = 2, mu = exp(1.6 + 0.4 * d$x + d$e))
  expect_equal(withr::with_options(list(frmtmb.latres_se_rows = 10L),
                                  lr_warned_on(f, negbinomial(), d)), 1L)
  # hurdle and zero-inflated rows have no dispersion to trade
  d <- dat(3001L)
  d$y <- stats::rpois(nrow(d), exp(1 + 0.4 * d$x + d$e)) *
    stats::rbinom(nrow(d), 1L, 0.7)
  expect_equal(lr_warned_on(f, zero_inflated_poisson(), d), 0L)
})

lr_review_ctrl <- function(seed) {
  set.seed(seed)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  e <- as.vector(replicate(10, stats::arima.sim(list(ar = 0.4), 6, sd = 0.25)))
  d$y <- stats::rpois(n, exp(0.5 + 0.3 * d$x + e))
  d
}

test_that("predict() holds a latent block's parameters at their estimates", {
  # log(sderr) near a boundary has a standard error of 10, and a draw of it
  # multiplied the latents by exp(30) before the inverse link
  d <- lr_review_ctrl(7L)
  fit <- allow_boundary(lr_fit("y", "ar(t, g, cov = TRUE)", poisson(), d))
  set.seed(7L)
  p <- predict(fit, ndraws = 500L)[, 1L]
  expect_lt(max(p) / max(d$y), 3)
  # the absent guard, on the same seed: the draws leave the double range
  local({
    local_mocked_bindings(latres_held_pos = function(fit, map) integer(0),
                          .package = "frmtmb")
    set.seed(7L)
    wild <- suppressWarnings(predict(fit, summary = FALSE, ndraws = 500L))
    expect_true(!all(is.finite(wild)) || max(wild) / max(d$y) > 1e3)
  })
})

test_that("a fit at a limit of its latent residual says so", {
  set.seed(1L)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  d$x <- stats::rnorm(60L)
  d$y <- stats::rpois(60L, exp(0.5 + 0.3 * d$x +
                                 rep(stats::rnorm(10L, 0, 0.5), each = 6L)))
  said <- character()
  hold <- function(m) {
    said <<- c(said, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
  fit <- withCallingHandlers(lr_fit("y", "ar(t, g, cov = TRUE)", poisson(), d),
                             frmtmb_boundary_fit = hold)
  expect_length(said, 1L)
  expect_match(said, "ar[1] = 1", fixed = TRUE)
  expect_match(said, "cosy() or (1 | group)", fixed = TRUE)
  expect_no_match(said, "theta_", fixed = TRUE)
  cp <- summary(fit)$cor_pars
  expect_true(all(is.na(cp[c("ar[1]", "sderr"), "Est.Error"])))
  expect_true(all(is.na(cp[c("ar[1]", "sderr"), "l-95% CI"])))
  sg <- suppressWarnings(diagnose(fit, quiet = TRUE))$singular
  expect_setequal(sg$term, c("ar[1]", "sderr"))
  # the absent case: a fit inside the range keeps its numbers
  d2 <- lr_data(sd = 0.5)
  f2 <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d2)
  a2 <- frame_latres(f2$frame)[[1L]]$aux
  expect_false(any(latres_at_limit(c(log(0.5), 0.5), a2)$at))
  # each limit by kind, on the chart each kind reads, with no number to say the
  # end is determined
  expect_identical(latres_at_limit(c(log(0.5), 1e6), a2)$at, c(TRUE, TRUE))
  expect_identical(latres_at_limit(c(log(1e-6), 0.5), a2)$at, c(TRUE, FALSE))
  expect_identical(latres_at_limit(c(log(1e-6), 0.5), a2)$na, c(TRUE, TRUE))
  cz <- frame_latres(lr_fit("cnt", "cosy(t, g)", poisson(), d2)$frame)[[1L]]$aux
  expect_identical(latres_at_limit(c(log(0.5), -50), cz)$at, c(FALSE, TRUE))
  # the fit's numbers decide: a small SE of log sderr at ar 0.9992 is no limit,
  # and a small sderr against a small response is none either
  expect_false(any(latres_at_limit(c(log(0.05), 25), a2, c(0.07, 10))$at))
  expect_false(any(latres_at_limit(c(log(5e-5), 0.5), a2, c(0.07, 0.2),
                                   scale = 1e-4)$at))
  # a theta name a user cannot read is brms's name in the message
  expect_identical(latres_theta_names(frame_latres(f2$frame)[[1L]]),
                   c("sderr", "ar[1]"))
  tag <- " (in ar(t, g, cov = TRUE))"
  expect_identical(latres_par_text(f2, "theta_1, theta_2, theta_12"),
                   paste0("sderr", tag, ", ar[1]", tag, ", theta_12"))
  expect_identical(flat_par_display(f2, "theta_1"), paste0("sderr", tag))
  expect_true(all(is.finite(summary(f2)$cor_pars[, "Est.Error"])))
})

lr_said <- function(expr) {
  said <- character()
  val <- withCallingHandlers(
    expr,
    message = function(m) {
      said <<- c(said, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      said <<- c(said, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(val = val, said = paste(said, collapse = " | "))
}

test_that("the limit's message is true of the fit it names", {
  # ar at -1: an alternating group effect is no cosy(), which it fits worse by
  # up to 32 nats (seed 35)
  set.seed(35L)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  d$x <- stats::rnorm(60L)
  re <- as.vector(replicate(10L, stats::rnorm(1L, 0, 0.6) *
                              rep(c(1, -1), length.out = 6L)))
  d$y <- stats::rpois(60L, exp(0.5 + 0.3 * d$x + re))
  r <- lr_said(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = poisson(),
                   data = d))
  expect_no_match(r$said, "cosy()", fixed = TRUE)
  expect_match(r$said, "alternates", fixed = TRUE)
  # a strong AR on a long series is well determined
  set.seed(2L)
  n <- 2000L
  e <- numeric(n)
  e[1L] <- stats::rnorm(1L, 0, 0.05 / sqrt(1 - 0.9995^2))
  for (i in 2:n) e[i] <- 0.9995 * e[i - 1L] + stats::rnorm(1L, 0, 0.05)
  ds <- data.frame(t = seq_len(n), g = factor(1), x = stats::rnorm(n))
  ds$y <- stats::rpois(n, exp(1 + 0.3 * ds$x + e))
  r <- lr_said(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = poisson(),
                   data = ds))
  expect_no_match(r$said, "Boundary", fixed = TRUE)
  expect_true(all(is.finite(summary(r$val)$cor_pars[, "Est.Error"])))
  # one fit in small units is at no limit
  set.seed(7L)
  dd <- expand.grid(t = 1:8, g = factor(1:25))
  dd <- dd[order(dd$g, dd$t), ]
  dd$x <- stats::rnorm(nrow(dd))
  e <- as.vector(replicate(25, stats::arima.sim(list(ar = 0.6), 8, sd = 0.5)))
  dd$y <- 1e-4 * (0.3 * dd$x + e + stats::rnorm(nrow(dd), 0, 0.3))
  r <- lr_said(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = skew_normal(),
                   data = dd))
  expect_no_match(r$said, "no latent residual", fixed = TRUE)
  expect_true(is.finite(summary(r$val)$cor_pars["sderr", "Est.Error"]))
  # cancelling roots: white noise, no shared latent, and both rows named
  set.seed(2L)
  da <- expand.grid(t = 1:6, g = factor(1:10))
  da$x <- stats::rnorm(60L)
  da$y <- stats::rpois(60L, exp(0.5 + 0.3 * da$x +
                                  rep(stats::rnorm(10L, 0, 0.5), each = 6L)))
  r <- lr_said(frm(bf(y ~ x + arma(t, g, cov = TRUE)), family = poisson(),
                   data = da))
  expect_match(r$said, "ma[1]", fixed = TRUE)
  expect_match(r$said, "cancel", fixed = TRUE)
  expect_no_match(r$said, "cosy()", fixed = TRUE)
  # sderr at 0 names sderr alone; ar is not at a limit (seed 5)
  set.seed(5L)
  d5 <- expand.grid(t = 1:6, g = factor(1:10))
  d5$x <- stats::rnorm(60L)
  d5$y <- stats::rpois(60L, exp(0.5 + 0.3 * d5$x +
                                  rep(stats::rnorm(10L, 0, 0.5), each = 6L)))
  r <- lr_said(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = poisson(),
                   data = d5))
  expect_match(r$said, "no latent residual", fixed = TRUE)
  expect_no_match(r$said, "ar[1] =", fixed = TRUE)
  expect_identical(diagnose(r$val, quiet = TRUE)$singular$term, "sderr")
  # the generic check stands aside for the rows said, not the block's others
  fx <- r$val
  bk <- frame_latres(fx$frame)[[1L]]
  local_mocked_bindings(
    latres_block_limit = function(fit, bk, V = NULL) {
      list(at = c(FALSE, TRUE), na = c(FALSE, TRUE), se = c(1, NaN))
    }, .package = "frmtmb")
  expect_identical(latres_limit_par(fx, paste0("theta_", bk$theta_idx)),
                   c(FALSE, TRUE))
})

test_that("an unbounded objective is never reported beside a latent", {
  sim <- function(seed) {
    set.seed(seed)
    d <- expand.grid(t = 1:6, g = factor(1:25))
    d <- d[order(d$g, d$t), ]
    n <- nrow(d)
    d$x <- stats::rnorm(n)
    d$e <- as.vector(replicate(25, stats::arima.sim(list(ar = 0.5), 6,
                                                    sd = 0.5)))
    d$tr <- 10L
    d
  }
  # +6.5e17 on seed 3, and the binomial is the limit it reads
  d <- sim(3L)
  d$y <- stats::rbinom(nrow(d), 10L, stats::plogis(0.3 * d$x + d$e))
  r <- lr_said(frm(bf(y | trials(tr) ~ x + ar(t, g, cov = TRUE)),
                   family = beta_binomial(), data = d))
  lim <- lr_fit("y | trials(tr)", "ar(t, g, cov = TRUE)", binomial(), d)
  expect_lt(abs(as.numeric(logLik(r$val)) - as.numeric(logLik(lim))), 5)
  # +1.9e30 on seed 10, where the dispersion stops at the box and says so
  d <- sim(10L)
  m <- stats::plogis(0.3 * d$x + d$e)
  d$y <- pmin(pmax(stats::rbeta(nrow(d), 1000 * m, 1000 * (1 - m)), 1e-4),
              1 - 1e-4)
  r <- lr_said(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = Beta(), data = d))
  expect_lt(as.numeric(logLik(r$val)), 1e3)
  expect_match(r$said, "is a limit and not an estimate", fixed = TRUE)
  expect_match(r$said, "ar(t, g, cov = TRUE)", fixed = TRUE)
  # a dispersion with a predictor has no box: the guard reads NaN past the cap
  expect_true(is.nan(latres_disp_guard(exp(latres_disp_cap + 1))))
  expect_identical(latres_disp_guard(5), 5)
  # the box reaches the constant dispersion and leaves the rest of the vector
  fr <- r$val
  om <- outer_par_map(fr)
  up <- latres_disp_upper(fr$frame, fr$spec, rep(Inf, length(om$names)), om)
  expect_identical(up[om$comp == "betad"], latres_disp_box)
  expect_true(all(is.infinite(up[om$comp != "betad"])))
})

test_that("frm() refuses a post-link form before the fit", {
  d <- lr_data(sd = 0.5)
  w <- "no number from the fit says"
  refused <- list(
    c("gam", "ar(t, g)", "Gamma(link = 'log')"),
    c("gam", "ar(t, g)", "Gamma(link = 'inverse')"),
    c("bin", "ar(t, g)", "bernoulli(link = 'probit')"),
    c("cnt", "ar(t, g)", "poisson(link = 'sqrt')"))
  for (cs in refused) {
    f <- stats::as.formula(paste(cs[[1L]], "~ x +", cs[[2L]]))
    fam <- eval(str2lang(cs[[3L]]))
    expect_error(frm(bf(f), family = fam, data = d), w, fixed = TRUE,
                 info = cs[[3L]])
    expect_error(frm(bf(f), family = fam, data = d), "frm_sample()",
                 fixed = TRUE, info = cs[[3L]])
    # frm_sample() builds the same objective and is not refused
    expect_true(is.list(frm(bf(f), family = fam, data = d,
                            dry_run = "objective")), info = cs[[3L]])
  }
  # the same terms on the link scale are not refused before a fit
  fit <- lr_fit("gam", "ar(t, g, cov = TRUE)", Gamma(link = "log"), d)
  expect_length(frame_latres(fit$frame), 1L)
  fit <- lr_fit("cnt", "ar(t, g)", poisson(), d)
  expect_length(frame_latres(fit$frame), 1L)
  expect_null(latres_refuse_laplace(fit$spec, fit$frame))
})

test_that("a link that does not hold the mean is read by its margin", {
  sim <- function(seed, eta, y, K = 6L) {
    set.seed(seed)
    d <- expand.grid(t = seq_len(K), g = factor(1:25))
    d <- d[order(d$g, d$t), ]
    d$x <- stats::rnorm(nrow(d))
    d$y <- y(eta(d))
    d
  }
  ar1 <- function(s) {
    as.vector(replicate(25, stats::arima.sim(list(ar = 0.5), 6, sd = s)))
  }
  # Gamma() at shape 50, inverse link: far from the edge, and it fits
  d <- sim(2001L, function(d) 2 + 0.1 * d$x + ar1(0.3),
           function(e) stats::rgamma(length(e), 50, 50 * e))
  fit <- lr_fit("y", "ar(t, g, cov = TRUE)", Gamma(), d)
  m <- latres_margin(fit, frame_latres(fit$frame)[[1L]])
  expect_gte(m, latres_margin_floor)
  # poisson("identity") at a mean of 50 fits
  d <- sim(2001L, function(d) 50 + 3 * d$x + ar1(4),
           function(e) stats::rpois(length(e), pmax(e, 0.1)))
  expect_length(frame_latres(lr_fit("y", "ar(t, g, cov = TRUE)",
                                    poisson("identity"), d)$frame), 1L)
  # at a mean of 2 it refuses after the fit, with the margin in its sentence
  d <- sim(2001L, function(d) 2 + 0.5 * d$x + ar1(1),
           function(e) stats::rpois(length(e), pmax(e, 0.1)))
  e <- tryCatch(lr_fit("y", "ar(t, g, cov = TRUE)", poisson("identity"), d),
                error = conditionMessage)
  expect_match(e, "marginal sds of the latent from the edge", fixed = TRUE)
  expect_match(e, "frm_sample()", fixed = TRUE)
  # the frame of that fit builds, so frm_sample() reaches the model
  expect_true(is.list(frm(bf(y ~ x + ar(t, g, cov = TRUE)),
                          family = poisson("identity"), data = d,
                          dry_run = "objective")))
  # a unit range has two edges, and a held link has none
  a <- frame_latres(lr_fit("y", "ar(t, g, cov = TRUE)", poisson(),
                           lr_series(3L, 10L, 15L, 0, 0.5))$frame)[[1L]]
  expect_null(latres_margin(
    lr_fit("y", "ar(t, g, cov = TRUE)", poisson(),
           lr_series(3L, 10L, 15L, 0, 0.5)), a))
})

test_that("mi() beside a latent residual: refused only where a row is read", {
  # a missing value no other density reads leaves the objective
  # (mi_ml_frame()), and its row's latent then has no observation: the model
  # of a weight of 0 on that row. 0.79.1 refused it with a reason that is
  # false there; a read value still shares the inner problem with the
  # latents, which nobody has measured (dev/fitedge-findings.md, item 4)
  d <- lr_data(sd = 0.5)
  fit <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  expect_null(latres_refuse_laplace(fit$spec, fit$frame))
  set.seed(12301)
  dg <- expand.grid(t = 1:8, g = factor(1:10))
  dg <- dg[order(dg$g, dg$t), ]
  dg$x <- stats::rnorm(nrow(dg))
  dg$y <- stats::rgamma(nrow(dg), 5, 5 / exp(0.3 * dg$x))
  dg$y[c(3L, 20L, 41L)] <- NA
  fm <- lr_fit("y | mi()", "ar(t, g, cov = TRUE)", Gamma(link = "log"), dg)
  dw <- dg
  dw$w <- ifelse(is.na(dg$y), 0, 1)
  dw$y[is.na(dw$y)] <- 1
  fw <- lr_fit("y | weights(w)", "ar(t, g, cov = TRUE)", Gamma(link = "log"),
               dw)
  l0 <- as.numeric(logLik(fw))
  expect_lt(abs(as.numeric(logLik(fm)) - l0), 1e-8 * abs(l0))
  e0 <- unlist(fw$estimates[c("beta", "betad", "theta")])
  e1 <- unlist(fm$estimates[c("beta", "betad", "theta")])
  expect_lt(max(abs(e1 - e0)), 1e-6 * max(abs(e0)))
  # lognormal: the value is the exact marginal of log y over the observed
  # rows, sigma^2 I + sderr^2 C restricted to them, at a point of the test's
  # choosing
  dl <- lr_data(G = 6L, K = 4L)
  miss <- c(2L, 7L, 16L)
  dl$pos[miss] <- NA
  fl <- lr_fit("pos | mi()", lr_terms[["ar1"]], lognormal(), dl)
  par <- fl$opt$par
  th <- c(log(0.7), 0.6)
  par[names(par) == "theta"] <- th
  par[names(par) == "betad"] <- log(0.5)
  bk <- frame_latres(fl$frame)[[1L]]
  ar <- latres_coefs(th[-1L], bk$aux)$ar
  mu <- as.numeric(cbind(1, dl$x) %*% par[names(par) == "beta"])
  ll <- 0
  for (g in levels(dl$g)) {
    r <- which(dl$g == g & !is.na(dl$pos))
    k <- which(!is.na(dl$pos[dl$g == g]))
    C <- exp(2 * th[1L]) * lr_target_cov("ar1", 4L, ar = ar)[k, k] +
      diag(0.25, length(k))
    R <- chol(C)
    res <- backsolve(R, log(dl$pos[r]) - mu[r], transpose = TRUE)
    ll <- ll - sum(log(diag(R))) - 0.5 * length(k) * log(2 * pi) -
      0.5 * sum(res^2) - sum(log(dl$pos[r]))
  }
  expect_lt(abs(-fl$obj$fn(par) - ll), 1e-8 * abs(ll))
  # a value another response reads is still refused
  dg$w2 <- stats::rnorm(nrow(dg))
  expect_error(frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE),
                      family = Gamma(link = "log")) +
                     bf(w2 ~ mi(y), family = gaussian()) + set_rescor(FALSE),
                   data = dg),
               "another density reads its missing values", fixed = TRUE)
})

test_that("a read mi() predictor beside a latent residual is exact", {
  # lognormal reading a gaussian predictor's missing values: log y and xm are
  # jointly gaussian, so the missing xm integrate out in closed form, and the
  # Laplace value of the shared inner problem is that marginal, here at a
  # point of the test's choosing (dev/fitedge-findings.md, fix round, m3)
  set.seed(73001)
  G <- 6L
  K <- 4L
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$xm <- stats::rnorm(n, 0.5, 1)
  e <- as.vector(replicate(G, stats::arima.sim(list(ar = 0.6), K, sd = 0.5)))
  d$ln <- exp(0.6 + 0.4 * d$xm + e + stats::rnorm(n, 0, 0.3))
  d$xm[c(2L, 7L, 16L)] <- NA
  fit <- allow_warnings(
    frm(bf(ln ~ mi(xm) + ar(t, g, cov = TRUE), family = lognormal()) +
          bf(xm | mi() ~ 1, family = gaussian()) + set_rescor(FALSE),
        data = d),
    c("Optimizer did not report", "Standard errors are not available",
      "Large maximum absolute gradient", "The Laplace approximation of"))
  par <- fit$opt$par
  b <- c(0.6, 0.4, 0.5)
  sg <- c(0.3, 0.9)
  th <- c(log(0.5), 0.6)
  expect_identical(names(fit$estimates$beta),
                   c("ln_(Intercept)", "ln_mixm", "xm_(Intercept)"))
  par[names(par) == "beta"] <- b
  par[names(par) == "betad"] <- log(sg)
  par[names(par) == "theta"] <- th
  bk <- frame_latres(fit$frame)[[1L]]
  ar <- latres_coefs(th[-1L], bk$aux)$ar
  C <- exp(2 * th[1L]) * lr_target_cov("ar1", K, ar = ar)
  o <- !is.na(d$xm)
  m <- ifelse(o, d$xm, b[3L])
  ll <- sum(stats::dnorm(d$xm[o], b[3L], sg[2L], log = TRUE))
  for (g in levels(d$g)) {
    i <- which(d$g == g)
    V <- C + diag(sg[1L]^2 + b[2L]^2 * sg[2L]^2 * (!o[i]), K)
    R <- chol(V)
    res <- backsolve(R, log(d$ln[i]) - b[1L] - b[2L] * m[i], transpose = TRUE)
    ll <- ll - sum(log(diag(R))) - 0.5 * K * log(2 * pi) - 0.5 * sum(res^2) -
      sum(log(d$ln[i]))
  }
  expect_lt(abs(-fit$obj$fn(par) - ll), 1e-8 * abs(ll))
})

test_that("a nonlinear dpar beside a linear mu keeps its latent residual", {
  d <- lr_data(sd = 0.5)
  fit <- allow_boundary(allow_warnings(
    frm(bf(cnt ~ x + ar(t, g, cov = TRUE)) + nlf(shape ~ exp(a)) + lf(a ~ 1),
        family = negbinomial(), data = d),
    c("Optimizer did not report convergence", "The dispersion `",
      "Standard errors are not available", "The Laplace approximation of")))
  expect_length(frame_latres(fit$frame), 1L)
  fit <- allow_boundary(allow_warnings(
    frm(bf(cnt ~ x + ar(t, g, cov = TRUE)) + nlf(zi ~ exp(a) / (1 + exp(a))) +
          lf(a ~ 1),
        family = zero_inflated_poisson(), data = d),
    c("Optimizer did not report convergence",
      "Standard errors are not available", "The Laplace approximation of")))
  expect_length(frame_latres(fit$frame), 1L)
  # mu's own body is what has no linear predictor to add to, wherever the
  # parser's way to it goes
  ac <- list(fn = "ar", struct = "ar", cov = TRUE, p = 1L, q = 0L)
  expect_error(check_autocor_latent(
    list(family = poisson(), nl_dpars = "mu"), ac),
    "Latent residuals are not implemented for non-linear models", fixed = TRUE)
  expect_true(check_autocor_latent(
    list(family = negbinomial(), nl_dpars = "shape"), ac)[["latent"]])
})

test_that("fitted() on new rows is the mean over the latent", {
  d <- lr_data(sd = 0.5)
  nd <- data.frame(t = 1:5, g = factor(rep("new", 5L)), x = 0)
  fp <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  for (s in c(1, 2)) {
    fp$estimates$theta[frame_latres(fp$frame)[[1L]]$theta_idx] <-
      c(log(s), 0.5 / sqrt(0.75))
    ms <- s / sqrt(0.75)
    b0 <- unname(fp$estimates$beta[1L])
    v <- as.numeric(fitted(fp, newdata = nd)[, 1L])
    expect_equal(v, rep(exp(b0 + ms^2 / 2), 5L), tolerance = 1e-8)
    expect_equal(as.numeric(frm_linpred(fp, newdata = nd, type = "response")),
                 rep(exp(b0), 5L), tolerance = 1e-12)
  }
  # logit, against a quadrature that shares nothing with the method
  fb <- lr_fit("bin", "ar(t, g, cov = TRUE)", bernoulli(), d)
  ib <- frame_latres(fb$frame)[[1L]]$theta_idx
  fb$estimates$theta[ib] <- c(log(2), 0.5 / sqrt(0.75))
  ms <- 2 / sqrt(0.75)
  b0 <- unname(fb$estimates$beta[1L])
  ex <- stats::integrate(function(u) {
    stats::plogis(b0 + ms * u) * stats::dnorm(u)
  }, -Inf, Inf, rel.tol = 1e-10)$value
  # a logistic has poles near the real line, so this is the rule's harder case
  expect_equal(as.numeric(fitted(fb, newdata = nd)[, 1L]), rep(ex, 5L),
               tolerance = 1e-6)
  # marginal sd by kind: Var(err) of the kind's own covariance
  a <- frame_latres(fp$frame)[[1L]]$aux
  for (kd in c("ma", "arma", "cosy", "unstr")) {
    tm <- switch(kd, ma = "ma(t, g, cov = TRUE)",
                 arma = "arma(t, g, cov = TRUE)", cosy = "cosy(t, g)",
                 unstr = "unstr(t, g)")
    fk <- allow_boundary(lr_fit("cnt", tm, poisson(), d))
    bk <- frame_latres(fk$frame)[[1L]]
    th <- fk$estimates$theta[bk$theta_idx]
    set.seed(1L)
    draws <- replicate(
      4000L, as.numeric(latres_expand(latres_draw(th, bk), bk, th)))
    sdn <- stats::sd(draws)
    expect_lt(abs(latres_marginal_sd(th, bk$aux) / sdn - 1), 0.05, label = kd)
  }
  # the standard error is differenced where no design row exists
  se <- fitted(lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d),
               newdata = nd)[, "Est.Error"]
  expect_true(all(is.finite(se) & se > 0))
  # the absent case, a mean that is mu on an identity link
  ds <- d
  ds$yr <- 1 + 0.4 * ds$x + stats::rnorm(nrow(ds))
  fi <- allow_boundary(lr_fit("yr", "ar(t, g, cov = TRUE)", skew_normal(), ds))
  expect_equal(as.numeric(fitted(fi, newdata = nd)[, 1L]),
               rep(unname(fi$estimates$beta[1L]), 5L), tolerance = 1e-10)
  # a family whose fitted() is a matrix of probabilities averages cell by cell;
  # simulate() on new rows draws the same latent, so it is the check
  cats <- findInterval(0.5 * d$cnt + stats::rlogis(nrow(d)), c(1, 2, 3))
  d$ord <- factor(cats + 1L, ordered = TRUE)
  fo <- allow_boundary(lr_fit("ord", "ar(t, g, cov = TRUE)", cumulative(), d))
  bo <- frame_latres(fo$frame)[[1L]]
  fo$estimates$theta[bo$theta_idx] <- c(log(1.5), 0.5 / sqrt(0.75))
  po <- suppressWarnings(fitted(fo, newdata = nd))[, "Estimate", ]
  expect_equal(unname(rowSums(po)), rep(1, 5L), tolerance = 1e-10)
  set.seed(1L)
  sims <- as.matrix(simulate(fo, nsim = 4000L, newdata = nd))
  emp <- vapply(seq_len(ncol(po)), function(k) mean(sims[1L, ] == k), 0)
  se <- sqrt(po[1L, ] * (1 - po[1L, ]) / 4000)
  expect_lt(max(abs(emp - po[1L, ]) / se), 4.5)
  at0 <- frm_linpred(fo, newdata = nd, type = "response")
  expect_gt(max(abs(po[1L, ] - at0[1L, ]) / se), 4.5)
  # Gauss-Hermite's own accuracy, on the moment it is held to
  gh <- gh_std_normal(41L)
  for (s in c(1, 2.31, 4)) {
    expect_lt(abs(sum(gh$w * exp(s * gh$x)) / exp(s^2 / 2) - 1), 1e-8)
  }
})

test_that("response_mean_gh() is one rule for a row's total latent variance", {
  # by row: a total of 0 is the plain mean, and the variances add, so two
  # sources read as one call. The family objects are the fit's own.
  d <- lr_data(sd = 0.5)
  fp <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  fam <- fp$spec$responses[[1L]]$family
  eta <- c(0.5, 0.5, -1)
  v <- c(0, 1.3, 0.2 + 0.5)
  got <- response_mean_gh(eta, v, fam, list(mu = exp(eta)))
  expect_equal(got, exp(eta + v / 2), tolerance = 1e-10)
  # a mean that reads another dpar integrates whole: zi is not touched
  fz <- lr_fit("cnt", "ar(t, g, cov = TRUE)", zero_inflated_poisson(), d)
  famz <- fz$spec$responses[[1L]]$family
  zi <- c(0.1, 0.3, 0.5)
  got <- response_mean_gh(eta, v, famz, list(mu = exp(eta), zi = zi))
  expect_equal(got, (1 - zi) * exp(eta + v / 2), tolerance = 1e-10)
  # a caller after mu alone reads it through mean_fn
  got <- response_mean_gh(eta, v, famz, list(mu = exp(eta), zi = zi),
                          mean_fn = function(fam, dpars, aterms) dpars[["mu"]])
  expect_equal(got, exp(eta + v / 2), tolerance = 1e-10)
  # a scalar variance is every row's
  expect_equal(response_mean_gh(eta, 0.5, fam, list(mu = exp(eta))),
               exp(eta + 0.25), tolerance = 1e-10)
  # the standard error of that mean moves with sderr and ar as well as beta: by
  # the closed form exp(b0 + v / 2), with v = sderr^2 (1 + t^2)
  nd <- data.frame(t = 1:2, g = factor(rep("new", 2L)), x = 0)
  bk <- frame_latres(fp$frame)[[1L]]
  th <- fp$estimates$theta[bk$theta_idx]
  vv <- exp(2 * th[1L]) * (1 + th[2L]^2)
  m <- exp(unname(fp$estimates$beta[1L]) + vv / 2)
  V <- as.matrix(vcov(fp, full = TRUE))
  g <- stats::setNames(numeric(nrow(V)), rownames(V))
  g[rownames(V) == "theta_1"] <- m * vv
  g[rownames(V) == "theta_2"] <- m * exp(2 * th[1L]) * th[2L]
  g[1L] <- m
  v2 <- fitted(fp, newdata = nd)
  draw <- m * sqrt(exp(vv) - 1)
  got <- unname(v2[1L, "Est.Error"])
  expect_gte(got, 0.999 * draw)
  expect_lt(got, 1.05 * sqrt(draw^2 + drop(crossprod(g, V %*% g))))
})

test_that("the SEs of many fitted rows are read from the sparse factor", {
  # above 200 latents, so the fitted rows take the selected inverse and not
  # a dense block of the inner covariance (jcs_block_max)
  d <- lr_data(sd = 0.5, G = 30L, K = 8L)
  fit <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  # the reference is the joint covariance formed densely from the joint
  # precision (helper-latres-ref.R), not the sparse route
  ref <- latres_ref_se(fit)
  # above the row limit no dense block is formed, so nothing is skipped
  withr::local_options(frmtmb.latres_se_rows = 10L)
  said <- 0L
  count <- function(m) {
    said <<- said + 1L
    invokeRestart("muffleMessage")
  }
  v <- withCallingHandlers(fitted(fit), frmtmb_latres_se_skipped = count)
  eps <- sqrt(.Machine$double.eps)
  # poisson's log link: d mu / d eta is mu
  expect_lt(latres_ref_ratio(v[, "Est.Error"], v[, "Estimate"] * ref), eps)
  expect_true(all(is.finite(v[, "Estimate"])))
  l <- withCallingHandlers(frm_linpred(fit, se.fit = TRUE),
                           frmtmb_latres_se_skipped = count)
  expect_lt(latres_ref_ratio(l$se.fit, ref), eps)
  expect_equal(l$fit, frm_linpred(fit))
  expect_identical(said, 0L)
  # rows of a grid read beta and theta alone, above the limit too
  nd <- data.frame(t = 1:5, g = factor(rep("new", 5L)), x = 0)
  v <- withCallingHandlers(fitted(fit, newdata = nd),
                           frmtmb_latres_se_skipped = count)
  expect_true(all(is.finite(v[, "Est.Error"])))
  l <- withCallingHandlers(frm_linpred(fit, newdata = nd, se.fit = TRUE),
                           frmtmb_latres_se_skipped = count)
  expect_true(all(is.finite(l$se.fit)))
  cu <- withCallingHandlers(frm_curve(fit, newdata = nd),
                            frmtmb_latres_se_skipped = count)
  expect_true(is.data.frame(cu) || is.list(cu))
  ce <- withCallingHandlers(conditional_effects(fit, "x"),
                            frmtmb_latres_se_skipped = count)
  expect_true(all(is.finite(ce[[1L]]$lower__)))
  expect_identical(said, 0L)
})

test_that("a far node of negligible weight says nothing about the mean", {
  # an ordinary fit under trunc(ub = 100): a node's mu of 5000 reads NaN there
  d <- lr_data(seed = 167, sd = 0.5)
  fit <- lr_fit("cnt | trunc(ub = 100)", "ar(t, g, cov = TRUE)", poisson(), d)
  nd <- data.frame(t = 1:3, g = factor(rep("new", 3L)), x = c(-1, 0, 1))
  v <- fitted(fit, newdata = nd)
  expect_true(all(is.finite(v[, "Estimate"])))
  # hurdle_poisson() at a marginal sd of 4 reads what integrate() reads
  f2 <- allow_boundary(lr_fit("cnt", "ar(t, g, cov = TRUE)", hurdle_poisson(),
                              d))
  fam <- f2$spec$responses[[1L]]$family
  exact <- stats::integrate(function(u) {
    m <- exp(1 + 4 * u)
    0.7 * m / -expm1(-m) * stats::dnorm(u)
  }, -10, 10, rel.tol = 1e-8)$value
  got <- response_mean_gh(1, 16, fam, list(mu = exp(1), hu = 0.3))
  expect_true(is.finite(got))
  expect_lt(abs(got / exact - 1), 1e-3)
  # the mean itself, by its limit at a small mu
  expect_equal(response_mean(fam, list(mu = 1e-20, hu = 0.3), list()), 0.7,
               tolerance = 1e-12)
  # by construction: a node in the tail is dropped, a node in the middle is not
  pf <- fp <- list(links = list(mu = list(linkinv = identity)))
  far <- function(fam, dpars, aterms) ifelse(abs(dpars$mu) > 7, NaN, dpars$mu^2)
  expect_equal(response_mean_gh(0, 1, pf, list(), mean_fn = far), 1,
               tolerance = 1e-6)
  near <- function(fam, dpars, aterms) {
    ifelse(abs(dpars$mu) < 1, NaN, dpars$mu^2)
  }
  expect_true(is.nan(response_mean_gh(0, 1, fp, list(), mean_fn = near)))
})

test_that("two latent blocks read as two sets of states", {
  d <- lr_data(sd = 0.5)
  fit <- allow_boundary(suppressWarnings(
    frm(bf(cnt ~ x + ar(t, g, cov = TRUE)) + bf(gam ~ x + cosy(t, g)) +
          set_rescor(FALSE),
        family = list(poisson(), Gamma(link = "log")), data = d)))
  lab <- brms_r_labels(fit)
  expect_true(all(c("latent_state_cnt[1]", "latent_state_gam[1]") %in% lab))
  expect_false(anyDuplicated(lab) > 0L)
})

test_that("coef() lists no group for the latents", {
  d <- lr_data(sd = 0.5)
  fit <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  # as a fit with no group-level effect at all
  plain <- allow_boundary(frm(bf(cnt ~ x), family = poisson(), data = d))
  expect_identical(class(coef(fit)), class(coef(plain)))
  expect_identical(names(coef(fit)), names(coef(plain)))
  d$h <- factor(rep(1:6, length.out = nrow(d)))
  both <- allow_boundary(lr_fit("cnt", "ar(t, g, cov = TRUE) + (1 | h)",
                                poisson(), d))
  expect_identical(names(coef(both)), "h")
})

test_that("new rows name the column they need", {
  d <- lr_data(sd = 0.5)
  fit <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  expect_error(simulate(fit, newdata = d[, c("x", "g")]), "newdata needs `t`",
               fixed = TRUE)
  expect_error(predict(fit, newdata = d[, c("x", "t")]), "newdata needs `g`",
               fixed = TRUE)
})

test_that("asym_laplace() is refused in the term's own words", {
  d <- lr_data(sd = 0.5)
  d$yr <- stats::rnorm(nrow(d))
  e <- tryCatch(frm(bf(yr ~ x + ar(t, g, cov = TRUE)), family = asym_laplace(),
                    data = d), error = conditionMessage)
  expect_match(e, "ar(t, g, cov = TRUE)", fixed = TRUE)
  expect_no_match(e, "latres(", fixed = TRUE)
})

test_that("every group of one row is refused under cov = TRUE", {
  d <- data.frame(t = 1:40, g = factor(1:40), x = stats::rnorm(40))
  d$y <- stats::rpois(40L, 2)
  for (tm in c("ar(t, g, cov = TRUE)", "cosy(t, g)", "ma(t, g, cov = TRUE)")) {
    expect_error(
      frm(bf(stats::as.formula(paste("y ~ x +", tm))), family = poisson(),
          data = d, dry_run = "frame"),
      "every group has one row", fixed = TRUE, info = tm)
  }
  # one group of two lets the correlation show, and the check stands aside
  d2 <- d
  d2$g[2L] <- d2$g[1L]
  expect_true(is.list(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = poisson(),
                          data = d2, dry_run = "frame")))
})

test_that("a multi-category or simplex family reads no latent block", {
  d <- lr_data(sd = 0.5)
  d$un <- stats::rbeta(nrow(d), 2, 2)
  expect_error(frm(bf(cbind(un, 1 - un) ~ x + ar(t, g, cov = TRUE)),
                   family = dirichlet(), data = d, dry_run = "frame"),
               "Explicit covariance terms", fixed = TRUE)
  # the check by hand, on the others the parser reaches first today
  ac <- list(fn = "ar", struct = "ar", cov = TRUE, p = 1L, q = 0L)
  for (fam in list(categorical(), multinomial(), dirichlet_multinomial())) {
    expect_error(check_autocor_latent(list(family = fam), ac),
                 "Explicit covariance terms", fixed = TRUE)
  }
})

test_that("a chain's generations are its positions, in linear time", {
  # the rows of each generation by hand, unbalanced groups in any order
  gidx <- c(1L, 2L, 1L, 3L, 2L, 1L, 1L, 2L)
  tidx <- c(1L, 1L, 2L, 1L, 2L, 3L, 4L, 3L)
  ac <- list(struct = "ar", cov = TRUE, p = 1L, q = 0L, fn = "ar",
             label = "ar(t, g)")
  a <- latres_aux(ac, gidx, tidx, 8L, FALSE)
  st <- latres_steps(a)
  got <- lapply(st, function(r) sort(a$group[r]))
  expect_identical(got, list(c(1L, 2L), c(1L, 2L), 1L))
  # no generation of the chain reads a row its predecessor has not filled
  done <- logical(a$m)
  done[a$start] <- TRUE
  for (r in st) {
    expect_true(all(done[a$prev[r]]))
    done[r] <- TRUE
  }
  expect_true(all(done))
  # one series is one generation per row and the call does not scan them
  n <- 1e5L
  one <- latres_aux(ac, rep(1L, n), seq_len(n), n, FALSE)
  expect_length(latres_steps(one), n - 1L)
})

test_that("the refusal on new rows surfaces through marginaleffects", {
  skip_if_not_installed("marginaleffects")
  d <- lr_data(sd = 0.5)
  fu <- lr_fit("cnt", "ar(t, g)", poisson(), d)
  e <- tryCatch(marginaleffects::predictions(fu, newdata = d[1:3, ]),
                error = conditionMessage)
  expect_match(paste(e, collapse = " "),
               "Cannot predict new latent residuals", fixed = TRUE)
  # on cov = TRUE, the mean fitted() reads
  ft <- lr_fit("cnt", "ar(t, g, cov = TRUE)", poisson(), d)
  nd <- data.frame(t = 1:3, g = factor(rep("new", 3L)), x = 0)
  pr <- suppressWarnings(marginaleffects::predictions(
    ft, newdata = nd, allow_new_levels = TRUE, vcov = FALSE))
  expect_equal(as.numeric(pr$estimate),
               as.numeric(fitted(ft, newdata = nd)[, 1L]), tolerance = 1e-10)
})

test_that("a cs() term above the row limit gives a matrix se.fit", {
  d <- lr_series(3L, 15L, 10L, 0, 0.8)
  d$o <- ordered(cut(d$x + stats::rlogis(nrow(d)), c(-Inf, -1, 0, 1, Inf),
                     labels = 1:4))
  fit <- suppressWarnings(frm(bf(o ~ cs(x) + ar(t, g, cov = TRUE)),
                              family = acat(), data = d))
  v <- frm_linpred(fit)
  expect_length(dim(v), 2L)
  withr::local_options(frmtmb.latres_se_rows = 10L)
  l <- suppressMessages(frm_linpred(fit, se.fit = TRUE))
  expect_identical(dim(l$se.fit), dim(v))
  expect_true(all(is.na(l$se.fit)))
})

test_that("the correction reads the factor of its own point", {
  # Matrix keeps a factor in the slot of the template TMB fills, so a Cholesky()
  # of it could hand back an earlier point's: C read 0.97 for 1.10
  d <- lr_series(3L, 15L, 10L, 0.2, 0.8)
  fit <- lr_fit("y", "ar(t, g, cov = TRUE)", poisson(), d)
  bk <- frame_latres(fit$frame)[[1L]]
  par <- fit$opt$par
  map <- outer_par_map(fit)
  # another point first, so that the template's slot holds its factor
  got <- latres_tk(fit, par + 0.2, bk, map)
  got <- latres_tk(fit, par, bk, map)
  env <- fit$obj$env
  full <- env$last.par
  H <- as.matrix(env$spHess(full, random = TRUE))
  pj <- latres_jacobian(bk, fit$estimates$theta[bk$theta_idx])
  A <- matrix(0, bk$aux$n, length(env$random))
  A[cbind(pj$i, pj$j)] <- pj$x
  v <- A %*% solve(H, t(A))
  dv <- diag(v)
  y <- d$y
  eta <- as.numeric(log(eval_dpars(fit)[[1L]]$mu))
  # a poisson's -log density has h3 = h4 = exp(eta) in eta
  m <- exp(eta)
  ref <- -sum(m * dv^2) / 8 + sum(outer(m * dv, m * dv) * v) / 8 +
    sum(outer(m, m) * v^3) / 12
  # a step of 0.02 sd in h3 and h4 reads to about 0.2 percent; the stale factor,
  # 12 percent off
  expect_lt(abs(got / ref - 1), 0.01)
})

test_that("the Laplace warning goes by the correction's shift", {
  f <- y ~ x + ar(t, g, cov = TRUE)
  # binary rows: the correction moves log(sderr) by far more than 0.15
  set.seed(1L)
  d <- expand.grid(t = 1:10, g = factor(1:20))
  d <- d[order(d$g, d$t), ]
  d$x <- stats::rnorm(200L)
  e <- as.vector(replicate(20, stats::arima.sim(list(ar = 0.6), 10, sd = 1)))
  d$y <- stats::rbinom(200L, 1L, stats::plogis(0.2 + 0.5 * d$x + e))
  expect_equal(lr_warned_on(f, bernoulli(), d), 1L)
  # Beta at phi 4, a row with little to say
  m <- stats::plogis(0.5 * e)
  d$y <- pmin(pmax(stats::rbeta(200L, 4 * m, 4 * (1 - m)), 1e-5), 1 - 1e-5)
  expect_equal(lr_warned_on(f, Beta(), d), 1L)
  # a skew_normal() has C of 0, at both units
  z <- 0.5 * stats::rnorm(200L)
  d$y <- e + z
  expect_equal(lr_warned_on(f, skew_normal(), d), 0L)
  d$y <- 10 * (e + z)
  expect_equal(lr_warned_on(f, skew_normal(), d), 0L)
})

test_that("a joint covariance past the limit is refused, by name", {
  d <- lr_series(3L, 10L, 15L, -0.2, 0.5)
  fit <- lr_fit("y", "ar(t, g, cov = TRUE)", poisson(), d)
  withr::local_options(frmtmb.latres_dense_rows = 10L)
  nd <- d[1:5, ]
  # what returns a dense block over every latent
  expect_error(frm_lp_basis(fit), "frmtmb.latres_dense_rows", fixed = TRUE)
  expect_error(frm_joint_cov(fit), "frmtmb.latres_dense_rows", fixed = TRUE)
  # new rows load no fitted latent, and read the estimates alone
  expect_true(all(is.finite(frm_curve(fit, newdata = nd)$.se)))
  # what names no joint covariance still runs
  expect_true(all(is.finite(frm_linpred(fit, newdata = nd))))
  expect_true(all(is.finite(fitted(fit, newdata = nd)[, "Estimate"])))
})

test_that("a coefficient that lost its value beside a chain is named", {
  # ma[1] reads NA with a finite and large standard error (190 to 730), and no
  # sentence named it (dev/review-latres-m3d.R, seed 1)
  set.seed(1L)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rpois(nrow(d),
                      exp(0.5 + 0.3 * d$x + rep(stats::rnorm(10, 0, 0.5), each = 6)))
  ms <- character()
  withCallingHandlers(
    frm(y ~ x + arma(t, g, cov = TRUE), family = poisson(), data = d),
    message = function(m) {
      ms <<- c(ms, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) invokeRestart("muffleWarning"))
  expect_true(any(grepl("ma[1]", ms, fixed = TRUE)))
})

test_that("Est.Error on a new row carries the latent's draw", {
  # M3's rule on a latent residual, as on a me() row: 0.133 against the draws'
  # sd of 1.005 on poisson before it
  d <- lr_series(3L, 25L, 12L, 0.5, 0.7)
  fit <- lr_fit("y", "ar(t, g, cov = TRUE)", poisson(), d)
  nd <- data.frame(t = 1:3, g = factor(rep("new", 3L)), x = c(-0.5, 0, 0.5))
  bk <- frame_latres(fit$frame)[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  s <- latres_marginal_sd(th, bk$aux)
  f <- fitted(fit, newdata = nd)
  eta <- as.numeric(frm_linpred(fit, newdata = nd, type = "link"))
  draws <- exp(outer(s * stats::qnorm(stats::ppoints(4000L)), eta, `+`))
  ref <- apply(draws, 2L, stats::sd)
  expect_true(all(f[, "Est.Error"] >= 0.95 * as.numeric(ref)))
  expect_lt(max(f[, "Est.Error"] / ref), 1.3)
  # the interval holds the draw's quantiles, and stays in the support
  expect_true(all(f[, "Q2.5"] > 0))
  # the mean is the mean of that draw
  expect_lt(max(abs(f[, "Estimate"] / colMeans(draws) - 1)), 1e-3)
})

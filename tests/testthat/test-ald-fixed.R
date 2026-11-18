# asym_laplace() fits without random effects (backlog items 4982 and
# 5490): the Hessian standard error is refused by name, the fit stops at
# the exact quantile regression solution and says nothing about its
# convergence there, a stop the proof does not cover still warns, and a
# response whose every row is censored on one side has no estimate.

# every warning and message of `expr`, by class
ald_conds <- function(expr) {
  w <- list()
  m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w[[length(w) + 1L]] <<- x
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x
    invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m,
       wclass = vapply(w, function(x) class(x)[1L], ""))
}

# the asymmetric Laplace log-likelihood at the fit's sigma and quantile,
# with the location coefficients `b`
ald_ll_at <- function(fit, b, X, y) {
  dp <- frmtmb:::eval_dpars(fit)[[1L]]
  r <- y - drop(X %*% b)
  p <- dp$quantile
  sum(log(p) + log1p(-p) - log(dp$sigma) - r * (p - (r < 0)) / dp$sigma)
}

ald_data <- function(seed, n = 100) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n)
  d
}

test_that("the Hessian standard error is refused by name", {
  d <- ald_data(5)
  r <- ald_conds(frm(y ~ x, family = asym_laplace(), data = d))
  f <- r$value
  # one warning, the refusal; no convergence warning at the exact
  # solution
  expect_identical(r$wclass, "frmtmb_ald_no_se")
  expect_match(conditionMessage(r$warnings[[1L]]),
               "linear in mu on each side of its mode", fixed = TRUE)
  expect_match(conditionMessage(r$warnings[[1L]]), "frm_bootstrap()",
               fixed = TRUE)
  lost <- sdr_of(f)$se_lost
  expect_identical(sort(names(lost)), sort(frmtmb:::outer_par_names(f)))
  expect_true(all(lost == "kink"))
  # every reader agrees: NaN, and the reason in words
  s <- summary(f)
  expect_true(all(is.nan(s$fixed[, "Est.Error"])))
  expect_match(paste(s$se_lost, collapse = " "),
               "no Hessian standard error, by design", fixed = TRUE)
  expect_true(all(is.nan(vcov(f))))
  ci <- confint(f)
  expect_true(all(is.nan(ci[, "lwr"])))
  ft <- allow_warnings(stats::fitted(f), "no Hessian standard error by design",
                       require = "no Hessian standard error by design")
  expect_true(all(is.nan(ft[, "Est.Error"])))
  out <- gsub("\\s+", " ", paste(utils::capture.output(diagnose(f)),
                                 collapse = " "))
  expect_match(out, "exact quantile regression solution", fixed = TRUE)
  expect_match(out, "no Hessian standard error, by design", fixed = TRUE)
  # under check_se the refusal is silent or an error, by its class
  r2 <- ald_conds(frm(y ~ x, family = asym_laplace(), data = d,
                      control = frmtmb_control(check_se = "ignore")))
  expect_length(r2$warnings, 0L)
  expect_error(frm(y ~ x, family = asym_laplace(), data = d,
                   control = frmtmb_control(check_se = "stop")),
               class = "frmtmb_ald_no_se")
})

test_that("the fit is the exact quantile regression solution", {
  skip_if_not_installed("quantreg")
  # the 4982 cases: gaussian data, a fixed quantile, asymmetric Laplace
  # data; quantreg::rq() at the fit's own quantile
  d <- ald_data(5)
  set.seed(43)
  da <- data.frame(x = stats::rnorm(160))
  da$y <- brms::rasym_laplace(160, 1 + 0.5 * da$x, 1, 0.3)
  cases <- list(list(fo = y ~ x, d = d),
                list(fo = bf(y ~ x, quantile = 0.5), d = d),
                list(fo = y ~ x, d = da))
  for (k in seq_along(cases)) {
    cs <- cases[[k]]
    r <- ald_conds(frm(cs$fo, family = asym_laplace(), data = cs$d))
    f <- r$value
    expect_identical(r$wclass, "frmtmb_ald_no_se", info = paste(k))
    expect_identical(f$opt$convergence, 0L, info = paste(k))
    expect_true(isTRUE(f$opt$ald_cert$ok), info = paste(k))
    p <- frmtmb:::eval_dpars(f)[[1L]]$quantile[1L]
    rq <- quantreg::rq(y ~ x, tau = p, data = cs$d)
    X <- stats::model.matrix(~x, cs$d)
    b <- frmtmb:::fixef_by_dpar(f)$mu
    l_fit <- ald_ll_at(f, b, X, cs$d$y)
    l_rq <- ald_ll_at(f, stats::coef(rq), X, cs$d$y)
    # rq's objective is not above the fit's, to rounding of the sum
    expect_lte((l_rq - l_fit) / abs(l_fit), 1e-13, label = paste(k))
    expect_lte(max(abs(b - stats::coef(rq))) / max(abs(stats::coef(rq))),
               1e-8, label = paste(k))
  }
})

test_that("the proof holds on every rounding of the data", {
  skip_on_cran()
  d0 <- ald_data(5, 60)
  set.seed(109)
  dt <- data.frame(x = sample(0:2, 60, TRUE))
  dt$y0 <- stats::rpois(60, 3 + dt$x)
  for (k in 0:20) {
    d <- d0
    d$y <- d0$y * (1 + 2e-15 * k)
    r <- ald_conds(frm(y ~ x, family = asym_laplace(), data = d))
    expect_identical(r$wclass, "frmtmb_ald_no_se", info = paste(k))
    expect_true(isTRUE(r$value$opt$ald_cert$global), info = paste(k))
    # tied data, where the vertex is degenerate
    dt$y <- dt$y0 * (1 + 2e-15 * k)
    r <- ald_conds(frm(bf(y ~ x, quantile = 0.5), family = asym_laplace(),
                       data = dt))
    expect_identical(r$wclass, "frmtmb_ald_no_se", info = paste("ties", k))
    expect_true(isTRUE(r$value$opt$ald_cert$global), info = paste("ties", k))
  }
})

test_that("an estimated quantile is searched for the highest maximum", {
  skip_if_not_installed("quantreg")
  # the review's design B, seed 109 (dev/reviews/2026-11-17-alapse.md,
  # finding 1): the first proved maximum is at quantile 0.4197, logLik
  # -137.159; the highest, from the quantile process, is at 0.7247,
  # -134.506, and the fit used to stop at the first with code 0
  set.seed(2109)
  d <- data.frame(x = stats::rnorm(100))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(100)
  r <- ald_conds(frm(y ~ x, family = asym_laplace(), data = d))
  f <- r$value
  expect_identical(r$wclass, "frmtmb_ald_no_se")
  expect_identical(f$opt$convergence, 0L)
  expect_true(isTRUE(f$opt$ald_cert$global))
  # the highest maximum over every solution of the quantile process
  pr <- suppressWarnings(quantreg::rq(y ~ x, tau = -1, data = d))
  taus <- pr$sol[1L, ]
  B <- pr$sol[-(1:3), , drop = FALSE]
  X <- stats::model.matrix(~x, d)
  top <- -Inf
  for (j in seq_len(ncol(B))) {
    lo <- max(taus[j], 1e-9)
    hi <- min(if (j < ncol(B)) taus[j + 1L] else 1, 1 - 1e-9)
    if (hi <= lo) next
    rr <- d$y - drop(X %*% B[, j])
    prof <- function(p) {
      100 * (log(p) + log1p(-p) - log(sum(rr * (p - (rr < 0))) / 100) - 1)
    }
    top <- max(top, stats::optimize(prof, c(lo, hi), maximum = TRUE,
                                    tol = 1e-12)$objective)
  }
  ll <- as.numeric(logLik(f))
  expect_lte((top - ll) / abs(ll), 1e-10)
  expect_gt(frmtmb:::eval_dpars(f)[[1L]]$quantile[1L], 0.7)
})

test_that("a quantile or scale model is a local maximum, and says so", {
  # the review's designs D (sigma ~ z, seed 7) and E (quantile ~ z, seed
  # 18): a multi-start maximum was 0.313 and 1.885 higher than the
  # first proved one, which then had code 0 and no warning
  set.seed(4007)
  d <- data.frame(x = stats::rnorm(200), z = stats::rnorm(200))
  d$y <- 1 + 0.5 * d$x + exp(0.4 * d$z) * stats::rnorm(200)
  r <- ald_conds(frm(bf(y ~ x, sigma ~ z), family = asym_laplace(),
                     data = d))
  expect_setequal(r$wclass, c("frmtmb_ald_local", "frmtmb_ald_no_se"))
  expect_false(identical(r$value$opt$convergence, 0L))
  expect_true(isTRUE(r$value$opt$ald_cert$ok))
  expect_false(isTRUE(r$value$opt$ald_cert$global))
  w <- r$warnings[r$wclass == "frmtmb_ald_local"][[1L]]
  expect_match(conditionMessage(w), "Fix the quantile", fixed = TRUE)
  set.seed(5018)
  d <- data.frame(x = stats::rnorm(200), z = stats::rnorm(200))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(200)
  r <- ald_conds(frm(bf(y ~ x, quantile ~ z), family = asym_laplace(),
                     data = d))
  expect_setequal(r$wclass, c("frmtmb_ald_local", "frmtmb_ald_no_se"))
  expect_false(identical(r$value$opt$convergence, 0L))
  expect_false(isTRUE(r$value$opt$ald_cert$global))
})

test_that("a profile interval is refused for a refused parameter", {
  d <- ald_data(5)
  f <- allow_warnings(frm(y ~ x, family = asym_laplace(), data = d),
                      "no Hessian standard error")
  expect_error(confint(f, parm = "x", method = "profile"),
               class = "frmtmb_ald_no_se")
  expect_error(confint(f, parm = "x", method = "uniroot"),
               class = "frmtmb_ald_no_se")
})

test_that("distributional and zero-inflated fits are covered", {
  d <- ald_data(7, 150)
  d$y <- d$y * exp(0.3 * d$z)
  r <- ald_conds(frm(bf(y ~ x, sigma ~ z, quantile ~ z),
                     family = asym_laplace(), data = d))
  expect_setequal(r$wclass, c("frmtmb_ald_local", "frmtmb_ald_no_se"))
  expect_true(isTRUE(r$value$opt$ald_cert$ok))
  expect_true(all(sdr_of(r$value)$se_lost == "kink"))
  expect_length(sdr_of(r$value)$se_lost, 6L)
  # zi reads only the zero rows: its coefficients keep the standard
  # errors of their own block, a logistic regression on the zeros
  set.seed(9)
  dz <- data.frame(x = stats::rnorm(200))
  dz$y <- ifelse(stats::runif(200) < stats::plogis(-1.5 + dz$x), 0,
                 brms::rasym_laplace(200, 1 + 0.5 * dz$x, 1, 0.5))
  r <- ald_conds(frm(bf(y ~ x, zi ~ x), family = zero_inflated_asym_laplace(),
                     data = dz))
  f <- r$value
  expect_identical(r$wclass, "frmtmb_ald_no_se")
  lost <- sdr_of(f)$se_lost
  expect_setequal(names(lost), c("(Intercept)", "x", "sigma_(Intercept)",
                                 "quantile_(Intercept)"))
  se <- sqrt(diag(sdr_of(f)$cov.fixed))
  zi <- grepl("^zi_", frmtmb:::outer_par_names(f))
  expect_true(all(is.finite(se[zi]) & se[zi] > 0))
  g <- stats::glm(I(y == 0) ~ x, family = stats::binomial(), data = dz)
  ratio <- se[zi] / sqrt(diag(stats::vcov(g)))
  expect_lt(max(abs(ratio - 1)), 1e-4)
})

test_that("a stop the proof does not cover still warns", {
  d <- ald_data(5)
  # a stop held short (three iterations), from which the exact solution
  # is higher
  ctl <- frmtmb_control(optCtrl = list(iter.max = 3), restarts = 0)
  ok <- ald_conds(frm(y ~ x, family = asym_laplace(), data = d,
                      control = ctl))
  expect_identical(ok$wclass, "frmtmb_ald_no_se")
  # the vertex search taken away: the optimizer's own stop keeps its
  # convergence warnings beside the refusal
  testthat::local_mocked_bindings(ald_vertex = function(...) NULL)
  r <- ald_conds(frm(y ~ x, family = asym_laplace(), data = d,
                     control = ctl))
  f <- r$value
  expect_false(isTRUE(f$opt$ald_cert$ok))
  expect_gt((as.numeric(logLik(ok$value)) - as.numeric(logLik(f))) /
              abs(as.numeric(logLik(f))), 1e-6)
  expect_true(any(r$wclass %in% c("frmtmb_conv_code",
                                  "frmtmb_conv_gradient")))
  expect_true("frmtmb_ald_no_se" %in% r$wclass)
  code <- r$warnings[r$wclass == "frmtmb_conv_code"]
  for (w in code) {
    expect_match(conditionMessage(w), "quantile regression's optimum failed",
                 fixed = TRUE)
  }
})

test_that("a censored row, another link or a random effect is not proved", {
  d <- ald_data(5)
  d$cy <- ifelse(d$y > 2, "right", "none")
  d$yc <- pmin(d$y, 2)
  r <- ald_conds(frm(yc | cens(cy) ~ x, family = asym_laplace(), data = d))
  # censored rows read curved tail probabilities: no proof, the
  # refusal still applies (the exact rows' kinks)
  expect_false(isTRUE(r$value$opt$ald_cert$ok))
  expect_match(r$value$opt$ald_cert$why, "censored", fixed = TRUE)
  expect_true("frmtmb_ald_no_se" %in% r$wclass)
  # a random effect outside mu: the outer Hessian is the Laplace
  # marginal's, which the refusal by name does not cover
  d$g <- factor(rep(1:10, 10))
  f <- allow_warnings(frm(bf(y ~ x, sigma ~ (1 | g)),
                          family = asym_laplace(), data = d),
                      c("did not report convergence", "Large maximum",
                        "Standard errors are not", "Boundary"))
  expect_null(f$opt$ald_cert)
  expect_false(any(sdr_of(f)$se_lost == "kink"))
  # the guard's absent case: a gaussian fit is untouched
  r <- ald_conds(frm(y ~ x, family = gaussian(), data = d))
  expect_length(r$warnings, 0L)
  expect_null(r$value$opt$ald_cert)
})

test_that("every row censored on one side has no estimate", {
  # backlog item 5490 (dev/wcheck2-allcens.R): 80 rows right censored at
  # their minimum, with and without a group term
  set.seed(20261021)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), g4 = factor(rep(1:20, each = 4)))
  d$ya <- 1 + 0.5 * d$x + stats::rnorm(n)
  d$ca <- "right"
  d$yac <- rep(min(d$ya), n)
  for (fo in list(yac | cens(ca) ~ x, yac | cens(ca) ~ x + (1 | g4))) {
    r <- ald_conds(frm(fo, family = asym_laplace(), data = d))
    sep <- r$warnings[r$wclass == "frmtmb_separation"]
    expect_length(sep, 1L)
    expect_match(conditionMessage(sep[[1L]]),
                 "maximum likelihood estimates do not exist", fixed = TRUE)
    expect_false(any(r$wclass %in% c("frmtmb_conv_code",
                                     "frmtmb_ald_no_se")))
    expect_true(all(sdr_of(r$value)$se_lost == "separation"))
  }
  # the guard's absent case: 72 of 80 censored, a finite estimate, and
  # the exact rows' refusal by name
  d$ca[1:8] <- "none"
  d$yac[1:8] <- d$ya[1:8] + abs(d$ya[1:8] - min(d$ya)) + 0.1
  r <- ald_conds(frm(yac | cens(ca) ~ x, family = asym_laplace(), data = d))
  expect_false("frmtmb_separation" %in% r$wclass)
  expect_true("frmtmb_ald_no_se" %in% r$wclass)
  expect_null(frmtmb:::ald_cens_sep(r$value))
})

test_that("a bound keeps the fit inside it, and the proof reads it", {
  # the short review of the fix round (finding 1): on design B seed 109,
  # whose highest maximum is at quantile 0.7247 (logit 0.968, log sigma
  # -1.267), a bound on the quantile's or sigma's intercept was crossed,
  # and the fit called "the highest maximum" outside it
  set.seed(2109)
  d <- data.frame(x = stats::rnorm(100))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(100)
  inside_box <- function(f) {
    bx <- f$cache$bounds
    np <- length(f$opt$par)
    all(f$opt$par >= rep_len(bx$lower %||% -Inf, np) &
          f$opt$par <= rep_len(bx$upper %||% Inf, np))
  }
  for (pr in list(set_prior("", class = "Intercept", dpar = "quantile",
                            lb = -2, ub = 0.3),
                  set_prior("", class = "Intercept", dpar = "sigma",
                            lb = -1.1, ub = 2))) {
    fo <- if (identical(pr$dpar, "quantile")) {
      bf(y ~ x, quantile ~ 1)
    } else bf(y ~ x, sigma ~ 1)
    r <- ald_conds(frm(fo, family = asym_laplace(), data = d, prior = pr))
    f <- r$value
    expect_true(inside_box(f), info = pr$dpar)
    expect_false(isTRUE(f$opt$ald_cert$global), info = pr$dpar)
    expect_false(identical(f$opt$convergence, 0L), info = pr$dpar)
    expect_true("frmtmb_ald_local" %in% r$wclass, info = pr$dpar)
  }
  # a box on sigma with a fixed quantile keeps the concave proof
  r <- ald_conds(frm(bf(y ~ x, sigma ~ 1, quantile = 0.5),
                     family = asym_laplace(), data = d,
                     prior = set_prior("", class = "Intercept",
                                       dpar = "sigma", lb = -1.1, ub = 2)))
  expect_true(inside_box(r$value))
  expect_true(isTRUE(r$value$opt$ald_cert$global))
  expect_identical(r$wclass, "frmtmb_ald_no_se")
})

test_that("a constant sigma on another link keeps the proof", {
  d <- ald_data(5)
  fl <- allow_warnings(frm(bf(y ~ x, quantile = 0.5),
                           family = asym_laplace(), data = d),
                       "no Hessian standard error")
  r <- ald_conds(frm(bf(y ~ x, quantile = 0.5),
                     family = asym_laplace(link_sigma = "softplus"),
                     data = d))
  expect_identical(r$wclass, "frmtmb_ald_no_se")
  expect_identical(r$value$opt$convergence, 0L)
  expect_true(isTRUE(r$value$opt$ald_cert$global))
  l1 <- as.numeric(logLik(r$value))
  l0 <- as.numeric(logLik(fl))
  expect_lte(abs(l1 - l0) / abs(l0), 1e-10)
  # and an estimated quantile on the identity link is searched
  r <- ald_conds(frm(y ~ x, family = asym_laplace(link_quantile = "identity"),
                     data = d))
  expect_identical(r$wclass, "frmtmb_ald_no_se")
  expect_true(isTRUE(r$value$opt$ald_cert$global))
})

test_that("the local warning's advice fits the model", {
  d <- ald_data(5)
  d$y2 <- d$y + stats::rnorm(nrow(d))
  r <- ald_conds(frm(bf(y ~ x) + bf(y2 ~ x) + set_rescor(FALSE),
                     family = asym_laplace(), data = d))
  w <- r$warnings[r$wclass == "frmtmb_ald_local"]
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]), "Each response fitted alone",
               fixed = TRUE)
  expect_false(grepl("drop the scale model", conditionMessage(w[[1L]]),
                     fixed = TRUE))
})

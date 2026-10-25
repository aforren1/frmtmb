# sar(), brms's spatial simultaneous autoregressive term, registered
# through frmtmb_register_autocor(). The references are dense
# multivariate densities written here; the brms rows at the end compare
# with brms 2.23.0's own Stan program.

sar_lattice <- function(k, queen = FALSE) {
  g <- expand.grid(r = seq_len(k), c = seq_len(k))
  n <- nrow(g)
  outer(seq_len(n), seq_len(n), function(i, j) {
    dr <- abs(g$r[i] - g$r[j])
    dc <- abs(g$c[i] - g$c[j])
    as.numeric(if (queen) pmax(dr, dc) == 1 else dr + dc == 1)
  })
}

sar_data <- function(k = 7, rho = 0.5, seed = 2, queen = FALSE) {
  set.seed(seed)
  A <- sar_lattice(k, queen)
  W <- A / rowSums(A)
  n <- nrow(W)
  d <- data.frame(x = rnorm(n))
  IA <- diag(n) - rho * W
  d$ylag <- as.vector(solve(IA, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
  d$yerr <- 1 + 0.5 * d$x + as.vector(solve(IA, rnorm(n, 0, 0.7)))
  list(d = d, W = W, A = A, n = n)
}

dense_mvn <- function(y, m, S) {
  R <- chol(S)
  z <- backsolve(R, y - m, transpose = TRUE)
  -0.5 * length(y) * log(2 * pi) - sum(log(diag(R))) - 0.5 * sum(z^2)
}

sar_rho_of <- function(fit) {
  ac <- fit$frame[["autocor"]][[1L]]
  unname(frmtmb:::autocor_natural(
    fit$estimates[["thetaac"]][ac[["theta_idx"]]], ac))
}

test_that("the lag SAR is the dense gaussian of A y - mu, with mean
           A^-1 mu", {
  s <- sar_data()
  d <- s$d
  W <- s$W
  fit <- frm(ylag ~ x + sar(W), data = d, data2 = list(W = W))
  rho <- sar_rho_of(fit)
  A <- diag(s$n) - rho * W
  eta <- frm_linpred(fit)
  m <- as.vector(solve(A, eta))
  S <- sigma(fit)^2 * solve(crossprod(A))
  ll <- as.numeric(logLik(fit))
  # an identity at the estimates, so rounding only
  expect_equal(ll, dense_mvn(d$ylag, m, S), tolerance = 1e-10)
  # and a maximum: the same density in its change-of-variables form,
  # log|A| from a dense determinant, optimized on its own
  X <- cbind(1, d$x)
  ev <- eigen(W, only.values = TRUE)$values
  lo <- 1 / min(Re(ev))
  hi <- 1 / max(Re(ev))
  nll <- function(p) {
    r <- lo + (hi - lo) / (1 + exp(-p[4]))
    Ap <- diag(s$n) - r * W
    e <- as.vector(Ap %*% d$ylag - X %*% p[1:2])
    -(sum(stats::dnorm(e, 0, exp(p[3]), log = TRUE)) +
        as.numeric(determinant(Ap)$modulus))
  }
  op <- stats::optim(c(0, 0, 0, 0), nll, method = "BFGS",
                     control = list(reltol = 1e-12, maxit = 1000))
  expect_lte(-op$value, ll + 1e-7 * abs(ll))
  # brms's posterior_epred(): the mean of the response
  expect_equal(unname(fitted(fit)[, "Estimate"]), m, tolerance = 1e-10)
  r <- residuals(fit)[, "Estimate"]
  expect_equal(unname(r), d$ylag - m, tolerance = 1e-10)
  rp <- residuals(fit, type = "pearson")[, "Estimate"]
  expect_equal(unname(rp), (d$ylag - m) / sqrt(diag(S)), tolerance = 1e-9)
  expect_equal(unname(autocor_matrix(fit)), solve(crossprod(A)),
               tolerance = 1e-9)
  # brms's names
  expect_true("lagsar" %in% variables(fit))
  cv <- confint_varcorr(fit)
  expect_identical(cv$term, "lagsar")
  expect_equal(cv$estimate, rho)
  expect_true(cv$lwr < rho && rho < cv$upr)
})

test_that("the error SAR is the dense gaussian of A (y - mu), with mean
           mu", {
  s <- sar_data()
  d <- s$d
  W <- s$W
  fit <- frm(yerr ~ x + sar(W, type = "error"), data = d,
             data2 = list(W = W))
  rho <- sar_rho_of(fit)
  A <- diag(s$n) - rho * W
  eta <- frm_linpred(fit)
  expect_equal(as.numeric(logLik(fit)),
               dense_mvn(d$yerr, eta, sigma(fit)^2 * solve(crossprod(A))),
               tolerance = 1e-10)
  expect_equal(unname(fitted(fit)[, "Estimate"]), unname(eta),
               tolerance = 1e-12)
  expect_true("errorsar" %in% variables(fit))
  # brms reads a unique abbreviation of the type
  f2 <- frm(yerr ~ x + sar(W, type = "e"), data = d, data2 = list(W = W))
  expect_identical(as.numeric(logLik(f2)), as.numeric(logLik(fit)))
})

test_that("sar() under student() is the multivariate t of the
           innovations", {
  s <- sar_data(seed = 3)
  d <- s$d
  W <- s$W
  fit <- frm(ylag ~ x + sar(W), family = student(), data = d,
             data2 = list(W = W))
  # one response is one multivariate-t draw, so nu is checked at a
  # point: nu = 5, rho = 0.3 (nu's link is log(nu - 1))
  p <- fit$obj$env$last.par.best
  bd <- which(names(p) == "betad")
  p[bd[2L]] <- log(5 - 1)
  aux <- fit$frame[["autocor"]][[1L]][["aux"]]
  p[names(p) == "thetaac"] <- stats::qlogis((0.3 - aux$lo) /
                                              (aux$hi - aux$lo))
  b <- p[names(p) == "beta"]
  sg <- exp(p[bd[1L]])
  A <- diag(s$n) - 0.3 * W
  r <- as.vector(A %*% d$ylag - cbind(1, d$x) %*% b)
  n <- s$n
  nu <- 5
  ref <- lgamma((nu + n) / 2) - lgamma(nu / 2) - n / 2 * log(nu * pi) -
    n * log(sg) + as.numeric(determinant(A)$modulus) -
    (nu + n) / 2 * log1p(sum(r^2) / (sg^2 * nu))
  expect_equal(-fit$obj$fn(p), unname(ref), tolerance = 1e-10)
})

test_that("rho lives on (1 / min eigenvalue, 1 / max eigenvalue)", {
  # a queen lattice, whose row-standardized W has its smallest
  # eigenvalue well above -1: brms's interval (min, max eigenvalue)
  # cuts off rho below it, and the data are drawn there
  s <- sar_data(k = 8, rho = -0.8, seed = 7, queen = TRUE)
  W <- s$W
  ev <- range(Re(eigen(W, only.values = TRUE)$values))
  expect_gt(ev[1L], -0.8)
  fit <- frm(ylag ~ x + sar(W), data = s$d, data2 = list(W = W))
  rho <- sar_rho_of(fit)
  expect_lt(rho, ev[1L])
  expect_gt(rho, 1 / ev[1L])
  aux <- fit$frame[["autocor"]][[1L]][["aux"]]
  expect_equal(c(aux$lo, aux$hi), 1 / ev, tolerance = 1e-10)
  # an unstandardized W: the interval is (1 / min, 1 / max) of ITS
  # eigenvalues, so rho * max eigenvalue stays below 1
  A <- s$A
  fa <- frm(ylag ~ x + sar(A), data = s$d, data2 = list(A = A))
  ra <- sar_rho_of(fa)
  expect_lt(ra * max(eigen(A, only.values = TRUE)$values), 1)
})

test_that("sar() takes a sparse Matrix, an nb and a listw object, and
           drops the rows the model drops", {
  s <- sar_data()
  d <- s$d
  W <- s$W
  f0 <- frm(ylag ~ x + sar(W), data = d, data2 = list(W = W))
  fs <- frm(ylag ~ x + sar(Ws), data = d,
            data2 = list(Ws = Matrix::Matrix(W, sparse = TRUE)))
  expect_equal(as.numeric(logLik(fs)), as.numeric(logLik(f0)),
               tolerance = 1e-12)
  # spdep's classes, built by hand: an nb object is a list of neighbor
  # indices, and nb2mat() row-standardizes it, which brms relies on
  nb <- structure(lapply(seq_len(s$n), function(i) which(s$A[i, ] > 0)),
                  class = "nb")
  fn <- frm(ylag ~ x + sar(nb), data = d, data2 = list(nb = nb))
  expect_equal(as.numeric(logLik(fn)), as.numeric(logLik(f0)),
               tolerance = 1e-12)
  lw <- structure(list(neighbours = nb,
                       weights = lapply(nb, function(v) {
                         rep(1 / length(v), length(v))
                       })), class = c("listw", "nb"))
  fl <- frm(ylag ~ x + sar(lw), data = d, data2 = list(lw = lw))
  expect_equal(as.numeric(logLik(fl)), as.numeric(logLik(f0)),
               tolerance = 1e-12)
  nb0 <- nb
  nb0[[3]] <- 0L
  expect_error(frm(ylag ~ x + sar(nb0), data = d, data2 = list(nb0 = nb0)),
               "without neighbors", class = "frmtmb_error")
  # brms drops the rows of W the model drops, and does not
  # re-standardize
  dn <- d
  dn$x[c(4, 30)] <- NA
  fna <- suppressMessages(frm(ylag ~ x + sar(W), data = dn,
                              data2 = list(W = W)))
  keep <- setdiff(seq_len(s$n), c(4, 30))
  fk <- frm(ylag ~ x + sar(Wk), data = d[keep, ],
            data2 = list(Wk = W[keep, keep]))
  expect_equal(as.numeric(logLik(fna)), as.numeric(logLik(fk)),
               tolerance = 1e-12)
})

test_that("sar() draws carry the fitted covariance, and newdata on the
           fitted rows is a counterfactual", {
  s <- sar_data(k = 4)
  d <- s$d
  W <- s$W
  fit <- frm(ylag ~ x + sar(W), data = d, data2 = list(W = W))
  A <- diag(s$n) - sar_rho_of(fit) * W
  S <- sigma(fit)^2 * solve(crossprod(A))
  m <- as.vector(solve(A, frm_linpred(fit)))
  sims <- as.matrix(simulate(fit, nsim = 4000, seed = 9))
  # moment checks against the sampling error at 4000 draws
  expect_lt(max(abs(rowMeans(sims) - m)), 6 * sqrt(max(diag(S)) / 4000))
  expect_lt(max(abs(stats::cov(t(sims)) - S)) / max(diag(S)),
            6 * sqrt(2 / 4000))
  # newdata on the fitted rows, with x moved: the mean is
  # A^-1 X_new b, as brms's posterior_epred() gives it, and the draws
  # carry the same covariance
  nd <- d
  nd$x <- nd$x + 1
  eta <- frm_linpred(fit, newdata = nd)
  mn <- as.vector(solve(A, eta))
  expect_equal(unname(fitted(fit, newdata = nd)[, "Estimate"]), mn,
               tolerance = 1e-10)
  sn <- as.matrix(simulate(fit, newdata = nd, nsim = 4000, seed = 9))
  expect_lt(max(abs(rowMeans(sn) - mn)), 6 * sqrt(max(diag(S)) / 4000))
  expect_identical(dim(predict(fit, newdata = nd)), c(s$n, 4L))
  # the error SAR moves no row's mean
  fe <- frm(yerr ~ x + sar(W, type = "error"), data = d,
            data2 = list(W = W))
  expect_equal(unname(fitted(fe, newdata = nd)[, "Estimate"]),
               unname(frm_linpred(fe, newdata = nd)), tolerance = 1e-12)
  # rows that differ from the fitted ones are refused, as brms refuses
  # them without newdata2
  n3 <- d[1:3, ]
  expect_error(fitted(fit, newdata = n3), "newdata2",
               class = "frmtmb_error")
  expect_error(simulate(fit, newdata = n3), "newdata2",
               class = "frmtmb_error")
  expect_error(predict(fit, newdata = n3), "one row per row of the term",
               class = "frmtmb_error")
  expect_equal(unname(frm_linpred(fit, newdata = n3)),
               unname(frm_linpred(fit)[1:3]), tolerance = 1e-12)
})

test_that("after rows are dropped, newdata covers every row of W, as in
           brms", {
  # brms reads newdata against W from data2, so it takes nrow(W) rows,
  # the two the fit dropped among them, and refuses the kept rows alone
  # (dev/sarfcor-fr2-rows.R)
  s <- sar_data(k = 5)
  d <- s$d
  W <- s$W
  d$x[c(4, 17)] <- NA
  fit <- suppressMessages(frm(ylag ~ x + sar(W), data = d,
                              data2 = list(W = W)))
  keep <- setdiff(seq_len(s$n), c(4, 17))
  nd <- d
  nd$x <- ifelse(is.na(nd$x), 0, nd$x) + 1
  rho <- sar_rho_of(fit)
  b <- unname(coef(fit))
  # the mean over all 25 rows of W, at the fitted rho
  want <- as.vector(solve(diag(s$n) - rho * W, b[1L] + b[2L] * nd$x))
  expect_equal(unname(fitted(fit, newdata = nd)[, "Estimate"]), want,
               tolerance = 1e-10)
  expect_identical(dim(predict(fit, newdata = nd)), c(as.integer(s$n), 4L))
  expect_identical(dim(as.matrix(simulate(fit, newdata = nd, nsim = 2,
                                          seed = 1))),
                   c(as.integer(s$n), 2L))
  # a missing predictor gives NA, which the lag spreads to every row,
  # as in brms
  expect_true(all(is.na(fitted(fit, newdata = transform(nd, x = d$x + 1))[
    , "Estimate"])))
  # the kept rows alone are refused, as brms refuses them
  expect_error(fitted(fit, newdata = nd[keep, ]),
               "one row per row of the term's matrix: 25 rows",
               fixed = TRUE, class = "frmtmb_error")
  expect_error(simulate(fit, newdata = nd[keep, ]), "dropped 2 of those",
               class = "frmtmb_error")
  # the absent case: without newdata, the fitted rows, as in brms
  expect_length(fitted(fit)[, "Estimate"], length(keep))
})

test_that("sar() is refused where brms refuses it, and where W gives no
           log-determinant", {
  s <- sar_data()
  d <- s$d
  W <- s$W
  d$cnt <- stats::rpois(s$n, 3)
  d$w <- 1.5
  fe <- function(...) expect_error(..., class = "frmtmb_error")
  fe(frm(cnt ~ x + sar(W), family = poisson(), data = d,
         data2 = list(W = W)), "SAR terms are not implemented")
  fe(frm(bf(ylag ~ x + sar(W), sigma ~ x), data = d, data2 = list(W = W)),
     "SAR models are not implemented when predicting")
  fe(frm(bf(ylag ~ x + sar(W), nu ~ x), family = student(), data = d,
         data2 = list(W = W)), "constant nu")
  fe(frm(ylag | weights(w) ~ x + sar(W), data = d, data2 = list(W = W)),
     "weights[(][)] cannot be combined")
  fe(frm(ylag ~ x + sar(W) + fcor(V), data = d,
         data2 = list(W = W, V = diag(s$n))), "residual correlation terms")
  fe(frm(ylag ~ x + sar(W, type = "mixed"), data = d, data2 = list(W = W)),
     "type must be")
  fe(frm(ylag ~ x + sar(Wb), data = d, data2 = list(Wb = diag(3))),
     "must be equal to the number of observations")
  fe(frm(ylag ~ x + sar(Wd), data = d,
         data2 = list(Wd = as.data.frame(W))), "must be of class")
  # a directed cycle: complex eigenvalues, so log|I - rho W| has no real
  # form over them
  Cc <- matrix(0, s$n, s$n)
  Cc[cbind(seq_len(s$n), c(2:s$n, 1))] <- 1
  fe(frm(ylag ~ x + sar(Cc), data = d, data2 = list(Cc = Cc)),
     "eigenvalues of 'M' are complex")
  fe(frm(ylag ~ x + sar(I), data = d, data2 = list(I = diag(s$n))),
     "needs a negative and a positive eigenvalue")
  # the absent case: a row-standardized W, whose eigen() values carry a
  # rounding imaginary part, is not refused, and sar()'s map passes
  # core's check of whiten() against unwhiten() at the frame, dense or
  # sparse, lag or error
  expect_true(is.complex(eigen(W, only.values = TRUE)$values))
  expect_s3_class(frm(ylag ~ x + sar(W), data = d, data2 = list(W = W)),
                  "frmtmb_fit")
  Wm <- Matrix::Matrix(W, sparse = TRUE)
  expect_s3_class(frm(yerr ~ x + sar(Wm, type = "error"), data = d,
                      data2 = list(Wm = Wm)), "frmtmb_fit")
})

test_that("a W whose extreme eigenvalue repeats warns that rho runs to
           the singular end", {
  # the complete graph: eigenvalues 1 and -1 / (n - 1), n - 1 times, so
  # I - rho W is singular at rho = 1 / lambda_min and the intercept and
  # rho are not separately identified there
  set.seed(4)
  K <- (matrix(1, 9, 9) - diag(9)) / 8
  dk <- data.frame(x = stats::rnorm(9), y = stats::rnorm(9))
  fk <- allow_warnings(
    frm(y ~ x + sar(K), data = dk, data2 = list(K = K)),
    c("is repeated", "Standard errors are not available",
      "did not report convergence", "Large maximum absolute gradient"),
    require = "is repeated")
  expect_s3_class(fk, "frmtmb_fit")
  # the absent case: a lattice W does not warn
  s <- sar_data(k = 4)
  expect_no_warning(frm(ylag ~ x + sar(W), data = s$d,
                        data2 = list(W = s$W)))
})

test_that("sar() works with a random effect beside it", {
  s <- sar_data(seed = 8)
  d <- s$d
  W <- s$W
  d$g <- factor(rep(seq_len(7), length.out = s$n))
  d$ylag <- d$ylag + rnorm(7, 0, 1.5)[d$g]
  fit <- frm(ylag ~ x + (1 | g) + sar(W), data = d, data2 = list(W = W))
  A <- diag(s$n) - sar_rho_of(fit) * W
  tau <- VarCorr(fit)$g$sd[1, "Estimate"]
  Z <- stats::model.matrix(~ 0 + g, d)
  # y = A^-1 (X b + Z u + e): the marginal of the lag model
  Ai <- solve(A)
  S <- Ai %*% (sigma(fit)^2 * diag(s$n) + tau^2 * tcrossprod(Z)) %*% t(Ai)
  m <- as.vector(Ai %*% frm_linpred(fit, re_formula = NA))
  expect_equal(as.numeric(logLik(fit)), dense_mvn(d$ylag, m, S),
               tolerance = 1e-8)
})

test_that("frm_compat() describes sar()", {
  expect_identical(frm_compat("sar()", "gaussian")$status, "works")
  expect_identical(frm_compat("sar()", "poisson")$status, "refused")
  expect_identical(frm_compat("sar()", "rescor")$status, "refused")
  expect_identical(frm_compat("sar()", "fitted")$status, "conditional")
})

# brms 2.23.0's Stan program at frmtmb's estimates (helper-brms.R,
# brms_lp_check()). brms's normal_lagsar_lpdf and normal_errorsar_lpdf
# leave out -N/2 log(2 pi), and student_t_lagsar_lpdf leaves out
# -N/2 log(pi), so those are the constants. brms passes eigen()'s values
# to Stan as they come, and a row-standardized W gives them a rounding
# imaginary part that rstan does not take, so the rows use the
# symmetric normalization D^-1/2 A D^-1/2, whose eigenvalues are real.
test_that("brms row: sar() log density equals brms's Stan program", {
  skip_unless_brms_fit()
  s <- sar_data()
  dg <- 1 / sqrt(rowSums(s$A))
  Ws <- s$A * outer(dg, dg)
  d <- s$d
  n <- s$n
  fl <- frm(ylag ~ x + sar(Ws), data = d, data2 = list(Ws = Ws))
  brms_lp_check(brms::bf(ylag ~ x + sar(Ws)), gaussian(), d, fl,
                const = 0.5 * n * log(2 * pi), data2 = list(Ws = Ws))
  fe <- frm(yerr ~ x + sar(Ws, type = "error"), data = d,
            data2 = list(Ws = Ws))
  brms_lp_check(brms::bf(yerr ~ x + sar(Ws, type = "error")), gaussian(),
                d, fe, const = 0.5 * n * log(2 * pi),
                data2 = list(Ws = Ws))
})

# brms's Stan program under student() at hand-set parameters, nu = 5:
# one response is one multivariate-t draw, so nu runs off at the
# estimates, and the density there is close to the gaussian one
sar_brms_lp_at <- function(bform, fit, data, data2, pn, nu) {
  fam <- brms::student()
  prior <- brms_flat_prior(bform, data = data, family = fam, data2 = data2)
  code <- brms::make_stancode(bform, data = data, family = fam,
                              prior = prior, data2 = data2)
  sdat <- brms_standata(bform, data = data, family = fam, prior = prior,
                        data2 = data2)
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code), data = sdat,
                                         chains = 0))
  b0 <- 0.8
  b1 <- 0.4
  sg <- 0.9
  rho <- 0.35
  pars <- list(b = array(b1, 1), Intercept = b0 + b1 * mean(data$x),
               sigma = sg, nu = nu)
  pars[[pn]] <- rho
  lp <- rstan::log_prob(sf, rstan::unconstrain_pars(sf, pars),
                        adjust_transform = FALSE, gradient = FALSE)
  p <- fit$obj$env$last.par.best
  p[names(p) == "beta"] <- c(b0, b1)
  bd <- which(names(p) == "betad")
  p[bd[1L]] <- log(sg)
  p[bd[2L]] <- log(nu - 1)
  aux <- fit$frame[["autocor"]][[1L]][["aux"]]
  p[names(p) == "thetaac"] <- stats::qlogis((rho - aux$lo) /
                                              (aux$hi - aux$lo))
  c(brms = lp, frmtmb = -fit$obj$fn(p))
}

test_that("brms row: sar() under student() equals brms's Stan program at
           nu = 5", {
  skip_unless_brms_fit()
  s <- sar_data()
  dg <- 1 / sqrt(rowSums(s$A))
  Ws <- s$A * outer(dg, dg)
  d <- s$d
  # student_t_lagsar_lpdf and student_t_errorsar_lpdf leave out
  # -N/2 log(pi)
  cst <- 0.5 * s$n * log(pi)
  for (type in c("lag", "error")) {
    y <- if (type == "lag") "ylag" else "yerr"
    fo <- paste0(y, " ~ x + sar(Ws, type = '", type, "')")
    ft <- frm(stats::as.formula(fo), family = student(), data = d,
              data2 = list(Ws = Ws))
    v <- sar_brms_lp_at(brms::bf(stats::as.formula(fo)), ft, d,
                        list(Ws = Ws), paste0(type, "sar"), 5)
    expect_lt(abs(v[["brms"]] - v[["frmtmb"]] - cst) / abs(v[["frmtmb"]]),
              1e-10)
  }
})

test_that("brms row: fitted() on the fitted rows with new x is brms's
           posterior_epred()", {
  skip_unless_brms_fit()
  # brms's draws held at frmtmb's estimates (Stan's Fixed_param), so a
  # difference is a difference in the method. The lag SAR only: brms
  # 2.23.0's posterior_epred() fails on every error SAR fit
  # (dev/upstream-bugs.md)
  s <- sar_data()
  dg <- 1 / sqrt(rowSums(s$A))
  Ws <- s$A * outer(dg, dg)
  d <- s$d
  fit <- frm(ylag ~ x + sar(Ws), data = d, data2 = list(Ws = Ws))
  b <- unname(coef(fit))
  pars <- list(b = array(b[2L], 1), Intercept = b[1L] + b[2L] * mean(d$x),
               sigma = sigma(fit), lagsar = sar_rho_of(fit))
  bform <- brms::bf(ylag ~ x + sar(Ws))
  fam <- stats::gaussian()
  prior <- brms_flat_prior(bform, data = d, family = fam,
                           data2 = list(Ws = Ws))
  code <- brms::make_stancode(bform, data = d, family = fam, prior = prior,
                              data2 = list(Ws = Ws))
  sdat <- brms_standata(bform, data = d, family = fam, prior = prior,
                        data2 = list(Ws = Ws))
  bfit <- suppressMessages(brms::brm(bform, data = d, family = fam,
                                     prior = prior, data2 = list(Ws = Ws),
                                     empty = TRUE))
  bfit$fit <- suppressMessages(rstan::sampling(
    brms_stan_model(code), data = sdat, chains = 1, iter = 2, warmup = 0,
    algorithm = "Fixed_param", init = list(pars), refresh = 0, seed = 1))
  bfit <- brms::rename_pars(bfit)
  nd <- d
  nd$x <- nd$x + 1
  be <- brms::posterior_epred(bfit, newdata = nd)[1L, ]
  fe <- unname(fitted(fit, newdata = nd)[, "Estimate"])
  expect_lt(max(abs(be - fe)) / max(abs(fe)), 1e-10)
  # and brms refuses rows that differ, as frmtmb does
  expect_error(brms::posterior_epred(bfit, newdata = nd[1:5, ]))
  # after the fit drops rows: brms takes newdata over every row of W and
  # refuses the kept rows alone, and so does frmtmb
  dd <- d
  dd$x[c(5, 20)] <- NA
  keep <- setdiff(seq_len(s$n), c(5, 20))
  fd <- suppressMessages(frm(ylag ~ x + sar(Ws), data = dd,
                             data2 = list(Ws = Ws)))
  b <- unname(coef(fd))
  pd <- list(b = array(b[2L], 1),
             Intercept = b[1L] + b[2L] * mean(dd$x[keep]),
             sigma = sigma(fd), lagsar = sar_rho_of(fd))
  # brms says so each time it drops the two rows
  na_note <- "Rows containing NAs were excluded from the model"
  prior <- allow_warnings(brms_flat_prior(bform, data = dd, family = fam,
                                          data2 = list(Ws = Ws)), na_note)
  code <- allow_warnings(brms::make_stancode(bform, data = dd,
                                             family = fam, prior = prior,
                                             data2 = list(Ws = Ws)),
                         na_note)
  sdat <- allow_warnings(brms_standata(bform, data = dd, family = fam,
                                       prior = prior,
                                       data2 = list(Ws = Ws)), na_note)
  bd <- allow_warnings(suppressMessages(
    brms::brm(bform, data = dd, family = fam, prior = prior,
              data2 = list(Ws = Ws), empty = TRUE)), na_note)
  bd$fit <- suppressMessages(rstan::sampling(
    brms_stan_model(code), data = sdat, chains = 1, iter = 2, warmup = 0,
    algorithm = "Fixed_param", init = list(pd), refresh = 0, seed = 1))
  bd <- brms::rename_pars(bd)
  nf <- dd
  nf$x <- ifelse(is.na(nf$x), 0, nf$x) + 1
  be <- brms::posterior_epred(bd, newdata = nf)[1L, ]
  fe <- unname(fitted(fd, newdata = nf)[, "Estimate"])
  expect_length(be, s$n)
  expect_lt(max(abs(be - fe)) / max(abs(fe)), 1e-10)
  expect_error(brms::posterior_epred(bd, newdata = nf[keep, ]),
               "must be equal to the number of observations")
  expect_error(fitted(fd, newdata = nf[keep, ]), "one row per row",
               class = "frmtmb_error")
})

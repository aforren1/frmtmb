# brms's fcor(M): a fixed residual correlation (or covariance) matrix
# over all the observations, core's first term of the residual
# correlation registry (R/autocor-seam.R). The references are dense
# multivariate densities written here; the brms rows at the end compare
# with brms 2.23.0's own Stan program.

fcor_data <- function(n = 30, seed = 1) {
  set.seed(seed)
  V <- 0.6^abs(outer(seq_len(n), seq_len(n), "-"))
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / 5), each = 5)))
  d$y <- 1 + 0.5 * d$x + as.vector(t(chol(V)) %*% rnorm(n)) * 0.8
  list(d = d, V = V, n = n)
}

dense_mvn <- function(y, m, S) {
  R <- chol(S)
  z <- backsolve(R, y - m, transpose = TRUE)
  -0.5 * length(y) * log(2 * pi) - sum(log(diag(R))) - 0.5 * sum(z^2)
}

test_that("fcor() is the dense gaussian with covariance sigma^2 M", {
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  fit <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  mu <- frm_linpred(fit)
  ll <- as.numeric(logLik(fit))
  # an identity at the estimates, so rounding only
  expect_equal(ll, dense_mvn(d$y, mu, sigma(fit)^2 * V), tolerance = 1e-10)
  # and a maximum: the dense likelihood optimized on its own is not
  # higher
  X <- cbind(1, d$x)
  nll <- function(p) -dense_mvn(d$y, X %*% p[1:2], exp(2 * p[3]) * V)
  op <- stats::optim(c(0, 0, 0), nll, method = "BFGS",
                     control = list(reltol = 1e-12))
  expect_lte(-op$value, ll + 1e-7 * abs(ll))
  # the mean is mu: fcor() moves no row's mean
  expect_equal(unname(fitted(fit)[, "Estimate"]), unname(mu),
               tolerance = 1e-12)
  expect_equal(unname(autocor_matrix(fit)), V, tolerance = 1e-12)
  # no parameter, so no row in the reports of the residual terms
  expect_false(any(grepl("fcor", variables(fit))))
})

test_that("an fcor() covariance is scaled by sigma, and a vector is its
           diagonal, as in brms", {
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  f1 <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  f4 <- frm(y ~ x + fcor(V4), data = d, data2 = list(V4 = 4 * V))
  expect_equal(sigma(f4), sigma(f1) / 2, tolerance = 1e-6)
  expect_equal(as.numeric(logLik(f4)), as.numeric(logLik(f1)),
               tolerance = 1e-9)
  v <- seq(0.5, 2, length.out = fd$n)
  fv <- frm(y ~ x + fcor(v), data = d, data2 = list(v = v))
  fm <- frm(y ~ x + fcor(Dv), data = d, data2 = list(Dv = diag(v)))
  expect_identical(as.numeric(logLik(fv)), as.numeric(logLik(fm)))
})

test_that("a predicted sigma gives fcor() the covariance D M D", {
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  fit <- frm(bf(y ~ x + fcor(V), sigma ~ x), data = d, data2 = list(V = V))
  s <- frm_linpred(fit, dpar = "sigma", type = "response")
  expect_equal(as.numeric(logLik(fit)),
               dense_mvn(d$y, frm_linpred(fit), diag(s) %*% V %*% diag(s)),
               tolerance = 1e-10)
  # pearson residuals divide by the marginal sd, sigma_i sqrt(M_ii)
  r <- residuals(fit, type = "pearson")[, "Estimate"]
  expect_equal(unname(r), unname((d$y - frm_linpred(fit)) /
                                   (s * sqrt(diag(V)))), tolerance = 1e-10)
})

test_that("fcor() under student() is the multivariate t", {
  skip_if_not_installed("mvtnorm")
  fd <- fcor_data(seed = 4)
  d <- fd$d
  V <- fd$V
  d$y <- 1 + 0.5 * d$x +
    as.vector(t(chol(V)) %*% stats::rt(fd$n, 4)) * 0.8
  fit <- frm(y ~ x + fcor(V), family = student(), data = d,
             data2 = list(V = V))
  # one response is ONE multivariate-t draw, so its nu is not
  # identified and runs off; the density is checked at a point instead,
  # nu = 5 (nu's link is log(nu - 1))
  p <- fit$obj$env$last.par.best
  bd <- which(names(p) == "betad")
  p[bd[2L]] <- log(5 - 1)
  b <- p[names(p) == "beta"]
  ref <- mvtnorm::dmvt(d$y, delta = as.vector(cbind(1, d$x) %*% b),
                       sigma = exp(2 * p[bd[1L]]) * V, df = 5, log = TRUE)
  expect_equal(-fit$obj$fn(p), ref, tolerance = 1e-10)
})

test_that("fcor() drops the rows the model drops from M, as brms does", {
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  dn <- d
  dn$x[c(3, 17)] <- NA
  fn <- suppressMessages(frm(y ~ x + fcor(V), data = dn,
                             data2 = list(V = V)))
  keep <- setdiff(seq_len(fd$n), c(3, 17))
  fk <- frm(y ~ x + fcor(Vk), data = d[keep, ],
            data2 = list(Vk = V[keep, keep]))
  expect_identical(as.numeric(logLik(fn)), as.numeric(logLik(fk)))
  expect_equal(unname(autocor_matrix(fn)), V[keep, keep])
  # M over the kept rows only is refused: brms reads M over the data
  expect_error(suppressMessages(frm(y ~ x + fcor(Vk), data = dn,
                                    data2 = list(Vk = V[keep, keep]))),
               "must be equal to the number of observations",
               class = "frmtmb_error")
})

test_that("fcor() with random effects is the dense marginal", {
  fd <- fcor_data(n = 40, seed = 6)
  d <- fd$d
  V <- fd$V
  d$y <- d$y + rnorm(nlevels(d$g))[d$g]
  fit <- frm(y ~ x + (1 | g) + fcor(V), data = d, data2 = list(V = V))
  tau <- VarCorr(fit)$g$sd[1, "Estimate"]
  Z <- stats::model.matrix(~ 0 + g, d)
  S <- sigma(fit)^2 * V + tau^2 * tcrossprod(Z)
  mu <- frm_linpred(fit, re_formula = NA)
  # the Laplace approximation is exact for a gaussian response
  expect_equal(as.numeric(logLik(fit)), dense_mvn(d$y, mu, S),
               tolerance = 1e-8)
})

test_that("fcor() under REML is generalized least squares", {
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  ml <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  re <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V), REML = TRUE)
  # with M known, sigma is the only covariance parameter: the GLS
  # coefficients are the same, and sigma^2 moves by n / (n - p)
  expect_equal(unname(coef(re)), unname(coef(ml)), tolerance = 1e-6)
  expect_equal(sigma(re)^2 / sigma(ml)^2, fd$n / (fd$n - 2),
               tolerance = 1e-6)
  pr <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V),
            control = frmtmb_control(profile = TRUE))
  expect_equal(as.numeric(logLik(pr)), as.numeric(logLik(ml)),
               tolerance = 1e-8)
})

test_that("fcor() draws carry sigma^2 M, and take newdata on the
           fitted rows", {
  fd <- fcor_data(n = 10, seed = 2)
  d <- fd$d
  V <- fd$V
  fit <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  sims <- as.matrix(simulate(fit, nsim = 4000, seed = 9))
  C <- stats::cov(t(sims))
  S <- sigma(fit)^2 * V
  # a moment check, against the sampling error of a covariance at 4000
  # draws: about sqrt(2 / 4000) of the diagonal
  expect_lt(max(abs(C - S)) / max(diag(S)), 6 * sqrt(2 / 4000))
  expect_lt(max(abs(rowMeans(sims) - frm_linpred(fit))),
            6 * sqrt(max(diag(S)) / 4000))
  # newdata on the fitted rows, with x moved: a counterfactual whose
  # mean is the new linear predictor, as in brms
  nf <- d
  nf$x <- nf$x + 1
  expect_equal(unname(fitted(fit, newdata = nf)[, "Estimate"]),
               unname(frm_linpred(fit, newdata = nf)), tolerance = 1e-12)
  sn <- as.matrix(simulate(fit, newdata = nf, nsim = 4000, seed = 9))
  expect_lt(max(abs(rowMeans(sn) - frm_linpred(fit, newdata = nf))),
            6 * sqrt(max(diag(S)) / 4000))
  expect_lt(max(abs(stats::cov(t(sn)) - S)) / max(diag(S)),
            6 * sqrt(2 / 4000))
  expect_identical(dim(predict(fit, newdata = nf)), c(as.integer(fd$n), 4L))
  # rows that differ from the fitted ones are refused, as brms refuses
  # them without newdata2
  nd <- d[1:4, ]
  expect_error(fitted(fit, newdata = nd), "newdata2",
               class = "frmtmb_error")
  expect_error(simulate(fit, newdata = nd), "newdata2",
               class = "frmtmb_error")
  expect_error(predict(fit, newdata = nd), "one row per row of the term",
               class = "frmtmb_error")
  # after rows are dropped, newdata covers every row of M, as brms reads
  # it, and the kept rows alone are refused, as brms refuses them
  dn <- d
  dn$x[c(3, 7)] <- NA
  fn <- suppressMessages(frm(y ~ x + fcor(V), data = dn,
                             data2 = list(V = V)))
  nf <- dn
  nf$x <- ifelse(is.na(nf$x), 0, nf$x) + 1
  expect_equal(unname(fitted(fn, newdata = nf)[, "Estimate"]),
               unname(frm_linpred(fn, newdata = nf)), tolerance = 1e-12)
  expect_length(fitted(fn, newdata = nf)[, "Estimate"], fd$n)
  sn2 <- as.matrix(simulate(fn, newdata = nf, nsim = 2, seed = 1))
  expect_identical(dim(sn2), c(as.integer(fd$n), 2L))
  expect_error(fitted(fn, newdata = nf[-c(3, 7), ]),
               "one row per row of the term's matrix: 10 rows",
               fixed = TRUE, class = "frmtmb_error")
  # the linear predictor does not use the term
  expect_equal(unname(frm_linpred(fit, newdata = nd)),
               unname(frm_linpred(fit)[1:4]), tolerance = 1e-12)
})

test_that("fcor() is refused where brms refuses it, and where M is not
           a covariance", {
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  d$cnt <- stats::rpois(fd$n, 3)
  d$w <- 1.5
  fe <- function(...) expect_error(..., class = "frmtmb_error")
  fe(frm(cnt ~ x + fcor(V), family = poisson(), data = d,
         data2 = list(V = V)),
     "FCOR terms are not implemented for this family")
  fe(frm(y | weights(w) ~ x + fcor(V), data = d, data2 = list(V = V)),
     "weights[(][)] cannot be combined")
  fe(frm(bf(y ~ x + fcor(V), nu ~ x), family = student(), data = d,
         data2 = list(V = V)), "constant nu")
  fe(frm(bf(y ~ x, sigma ~ fcor(V)), data = d, data2 = list(V = V)),
     "can only be written on 'mu'")
  fe(frm(y ~ x + fcor(V) + cosy(gr = g), data = d, data2 = list(V = V)),
     "residual correlation terms")
  fe(frm(y ~ x * fcor(V), data = d, data2 = list(V = V)),
     "cannot be crossed|whole term")
  Va <- V
  Va[1, 2] <- 0.9
  fe(frm(y ~ x + fcor(Va), data = d, data2 = list(Va = Va)),
     "must be symmetric")
  # one negative eigenvalue: brms 2.23.0's check passes it
  Vi <- V
  Vi[1, 2] <- Vi[2, 1] <- 1.2
  expect_lt(min(eigen(Vi, only.values = TRUE)$values), 0)
  fe(frm(y ~ x + fcor(Vi), data = d, data2 = list(Vi = Vi)),
     "must be positive definite")
  fe(frm(y ~ x + fcor(Vb), data = d, data2 = list(Vb = diag(3))),
     "must be equal to the number of observations")
  fe(frm(y ~ x + fcor(), data = d), "Argument 'M' is missing")
  # the absent case: a symmetric positive definite M of the right size
  expect_s3_class(frm(y ~ x + fcor(V), data = d, data2 = list(V = V)),
                  "frmtmb_fit")
})

test_that("each response of a multivariate model takes its own fcor()", {
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  d$y2 <- 2 - 0.3 * d$x + stats::rnorm(fd$n)
  f <- frm(bf(y ~ x + fcor(V)) + bf(y2 ~ x) + set_rescor(FALSE),
           data = d, data2 = list(V = V))
  f1 <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  f2 <- frm(y2 ~ x, data = d)
  expect_equal(as.numeric(logLik(f)),
               as.numeric(logLik(f1)) + as.numeric(logLik(f2)),
               tolerance = 1e-8)
  expect_error(frm(bf(y ~ x + fcor(V)) + bf(y2 ~ x) + set_rescor(TRUE),
                   data = d, data2 = list(V = V)),
               "rescor", class = "frmtmb_error")
})

test_that("frm_compat() describes fcor()", {
  expect_identical(frm_compat("fcor()", "gaussian")$status, "works")
  expect_identical(frm_compat("fcor()", "poisson")$status, "refused")
  expect_identical(frm_compat("fcor()", "weights()")$status, "refused")
  expect_identical(frm_compat("fcor()", "fitted")$status, "conditional")
})

# brms 2.23.0's Stan program for the same model, at frmtmb's estimates
# (helper-brms.R, brms_lp_check()). normal_fcor_hom_lpdf and
# student_t_fcor_hom_lpdf are multi_normal_cholesky and multi_student_t,
# full densities, so the constant is 0.
test_that("brms row: fcor() log density equals brms's Stan program", {
  skip_unless_brms_fit()
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  fit <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  brms_lp_check(brms::bf(y ~ x + fcor(V)), gaussian(), d, fit,
                data2 = list(V = V))
  fh <- frm(bf(y ~ x + fcor(V), sigma ~ x), data = d, data2 = list(V = V))
  brms_lp_check(brms::bf(y ~ x + fcor(V), sigma ~ x), gaussian(), d, fh,
                data2 = list(V = V))
})

test_that("brms row: fcor() under student() equals brms's Stan program at
           nu = 5", {
  skip_unless_brms_fit()
  # one response is one multivariate-t draw, so nu runs off at the
  # estimates; the density is compared at hand-set parameters instead
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  bform <- brms::bf(y ~ x + fcor(V))
  fam <- brms::student()
  prior <- brms_flat_prior(bform, data = d, family = fam,
                           data2 = list(V = V))
  code <- brms::make_stancode(bform, data = d, family = fam, prior = prior,
                              data2 = list(V = V))
  sdat <- brms_standata(bform, data = d, family = fam, prior = prior,
                        data2 = list(V = V))
  sf <- suppressMessages(rstan::sampling(brms_stan_model(code), data = sdat,
                                         chains = 0))
  ft <- frm(y ~ x + fcor(V), family = student(), data = d,
            data2 = list(V = V))
  b0 <- 0.8
  b1 <- 0.4
  sg <- 0.9
  nu <- 5
  lp <- rstan::log_prob(
    sf, rstan::unconstrain_pars(sf, list(
      b = array(b1, 1), Intercept = b0 + b1 * mean(d$x), sigma = sg,
      nu = nu)), adjust_transform = FALSE, gradient = FALSE)
  p <- ft$obj$env$last.par.best
  p[names(p) == "beta"] <- c(b0, b1)
  bd <- which(names(p) == "betad")
  p[bd[1L]] <- log(sg)
  p[bd[2L]] <- log(nu - 1)
  # multi_student_t is the full density, so the constant is 0
  ours <- -ft$obj$fn(p)
  expect_lt(abs(lp - ours) / abs(ours), 1e-10)
})

test_that("brms row: fitted() on the fitted rows with new x is brms's
           posterior_epred()", {
  skip_unless_brms_fit()
  # brms's draws held at frmtmb's estimates (helper-brms-methods.R)
  fd <- fcor_data()
  d <- fd$d
  V <- fd$V
  fit <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  bfit <- brms_fixed_fit(brms::bf(y ~ x + fcor(V)), stats::gaussian(), d,
                         fit, ndraws = 2, data2 = list(V = V))
  nd <- d
  nd$x <- nd$x + 1
  be <- brms::posterior_epred(bfit, newdata = nd)[1L, ]
  fe <- unname(fitted(fit, newdata = nd)[, "Estimate"])
  expect_lt(max(abs(be - fe)) / max(abs(fe)), 1e-10)
  expect_error(brms::posterior_epred(bfit, newdata = nd[1:5, ]),
               "must be equal to the number of observations")
})

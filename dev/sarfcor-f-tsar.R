# Fix round 1: test-sar.R (F1 absent case, F3 same-row newdata and the
# brms posterior_epred() pin, F6 student at nu = 5, F8 complete graph).
# Rscript <this> <root>
root <- commandArgs(TRUE)[1]
p <- file.path(root, "extensions/frmtmb.spatial/tests/testthat/test-sar.R")
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1('test_that("sar() draws carry the fitted covariance, in sample only", {',
     'test_that("sar() draws carry the fitted covariance, and newdata on the
           fitted rows is a counterfactual", {')
rep1('  nd <- d[1:3, ]
  expect_error(fitted(fit, newdata = nd), "newdata2",
               class = "frmtmb_error")
  expect_error(simulate(fit, newdata = nd), "newdata2",
               class = "frmtmb_error")
  expect_equal(unname(frm_linpred(fit, newdata = nd)),
               unname(frm_linpred(fit)[1:3]), tolerance = 1e-12)
})',
'  # newdata on the fitted rows, with x moved: the mean is
  # A^-1 X_new b, as brms\'s posterior_epred() gives it, and the draws
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
  # the error SAR moves no row\'s mean
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
  expect_error(predict(fit, newdata = n3), "with the fitted rows",
               class = "frmtmb_error")
  expect_equal(unname(frm_linpred(fit, newdata = n3)),
               unname(frm_linpred(fit)[1:3]), tolerance = 1e-12)
})')
rep1('  fe(frm(ylag ~ x + sar(I), data = d, data2 = list(I = diag(s$n))),
     "needs a negative and a positive eigenvalue")
  # the absent case: a row-standardized W, whose eigen() values carry a
  # rounding imaginary part, is not refused
  expect_true(is.complex(eigen(W, only.values = TRUE)$values))
  expect_s3_class(frm(ylag ~ x + sar(W), data = d, data2 = list(W = W)),
                  "frmtmb_fit")
})',
'  fe(frm(ylag ~ x + sar(I), data = d, data2 = list(I = diag(s$n))),
     "needs a negative and a positive eigenvalue")
  # the absent case: a row-standardized W, whose eigen() values carry a
  # rounding imaginary part, is not refused, and sar()\'s map passes
  # core\'s check of whiten() against unwhiten() at the frame, dense or
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
})')
rep1('  ft <- frm(ylag ~ x + sar(Ws), family = student(), data = d,
            data2 = list(Ws = Ws))
  brms_lp_check(brms::bf(ylag ~ x + sar(Ws)), brms::student(), d, ft,
                const = 0.5 * n * log(pi), data2 = list(Ws = Ws))
})',
'})

# brms\'s Stan program under student() at hand-set parameters, nu = 5:
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

test_that("brms row: sar() under student() equals brms\'s Stan program at
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
    fo <- paste0(y, " ~ x + sar(Ws, type = \'", type, "\')")
    ft <- frm(stats::as.formula(fo), family = student(), data = d,
              data2 = list(Ws = Ws))
    v <- sar_brms_lp_at(brms::bf(stats::as.formula(fo)), ft, d,
                        list(Ws = Ws), paste0(type, "sar"), 5)
    expect_lt(abs(v[["brms"]] - v[["frmtmb"]] - cst) / abs(v[["frmtmb"]]),
              1e-10)
  }
})

test_that("brms row: fitted() on the fitted rows with new x is brms\'s
           posterior_epred()", {
  skip_unless_brms_fit()
  # brms\'s draws held at frmtmb\'s estimates (Stan\'s Fixed_param), so a
  # difference is a difference in the method. The lag SAR only: brms
  # 2.23.0\'s posterior_epred() fails on every error SAR fit
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
})')
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("edited", p, "\n")

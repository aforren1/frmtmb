# Fix round 1: test-fcor.R (F3 same-row newdata and the brms
# posterior_epred() pin, F6 student at nu = 5). Rscript <this> <root>
root <- commandArgs(TRUE)[1]
p <- file.path(root, "tests/testthat/test-fcor.R")
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1('test_that("fcor() draws carry sigma^2 M, in sample only", {',
     'test_that("fcor() draws carry sigma^2 M, and take newdata on the
           fitted rows", {')
rep1('  nd <- d[1:4, ]
  expect_error(fitted(fit, newdata = nd), "newdata2",
               class = "frmtmb_error")
  expect_error(simulate(fit, newdata = nd), "newdata2",
               class = "frmtmb_error")
  expect_error(predict(fit, newdata = nd), "newdata",
               class = "frmtmb_error")',
'  # newdata on the fitted rows, with x moved: a counterfactual whose
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
  expect_error(predict(fit, newdata = nd), "with the fitted rows",
               class = "frmtmb_error")')
rep1('  ft <- frm(y ~ x + fcor(V), family = student(), data = d,
            data2 = list(V = V))
  brms_lp_check(brms::bf(y ~ x + fcor(V)), brms::student(), d, ft,
                data2 = list(V = V))
})',
'})

test_that("brms row: fcor() under student() equals brms\'s Stan program at
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

test_that("brms row: fitted() on the fitted rows with new x is brms\'s
           posterior_epred()", {
  skip_unless_brms_fit()
  # brms\'s draws held at frmtmb\'s estimates (helper-brms-methods.R)
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
})')
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
cat("edited", p, "\n")

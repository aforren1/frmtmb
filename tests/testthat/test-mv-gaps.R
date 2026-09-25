# Three multivariate gaps closed together (lane mv, 2026-09-25): any
# number of bf() formulas summed with `+`, student() with rescor = TRUE,
# and the ordinal families inside a multivariate model. The references
# and their numbers are in dev/mv-findings.md.

mv_gap_data <- function(seed = 5, n = 150) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  data.frame(x = x, z = z,
             y1 = 1 + 0.5 * x + rnorm(n), y2 = -1 + rnorm(n),
             y3 = 0.3 * x + rnorm(n),
             o = cut(x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE),
             o2 = cut(-x + rlogis(n), c(-Inf, 0, 1, Inf), labels = FALSE))
}

# ---------------------------------------------------------------- `+`

test_that("three or more bf() formulas add into one multivariate formula", {
  # R's Ops dispatch fell back to the internal `+` when the two operands
  # had different methods, so the third bf() stopped the sum
  f <- bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + gaussian()
  expect_s3_class(f, "frmtmb_mvformula")
  expect_identical(vapply(f$forms, function(g) deparse1(g$formula), ""),
                   c("y1 ~ x", "y2 ~ x", "y3 ~ x"))
  f4 <- bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + bf(y4 ~ x)
  expect_length(f4$forms, 4L)
  # a parenthesized sum flattens the same way
  expect_length((bf(y1 ~ x) + (bf(y2 ~ x) + bf(y3 ~ x)))$forms, 3L)
  # set_rescor() anywhere in the chain
  f <- bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE) + bf(y3 ~ x)
  expect_true(f$rescor)
  expect_length(f$forms, 3L)
  # a family after a bf() is that response's; a family added to the sum
  # fills only the responses without one
  f <- bf(o ~ x) + cumulative() + bf(y1 ~ x) + gaussian() + bf(y2 ~ x) +
    poisson()
  expect_identical(vapply(f$forms, function(g) g$family[["family"]], ""),
                   c("cumulative", "gaussian", "poisson"))
  # lf() / nlf() name their response with resp =
  f <- bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + lf(sigma ~ z, resp = "y3") +
    nlf(sigma ~ a * z, a ~ 1, resp = "y1")
  expect_named(f$forms[[3]]$pforms, "sigma")
  expect_named(f$forms[[1]]$nlforms, "sigma")
  expect_length(f$forms[[2]]$pforms, 0L)
})

test_that("the `+` refusals name what to write instead", {
  expect_error(bf(y1 ~ x) + bf(y2 ~ x) + lf(sigma ~ z),
               "resp = \"y2\"")
  expect_error(bf(y1 ~ x) + bf(y2 ~ x) + nlf(sigma ~ a * z, a ~ 1),
               "does not say which response")
  expect_error(bf(y1 ~ x) + bf(y2 ~ x) + lf(sigma ~ z, resp = "y9"),
               "not one of the responses: y1, y2")
  expect_error(bf(y1 ~ x) + lf(sigma ~ z, resp = "y9"),
               "models 'y1'")
  expect_error(gaussian() + bf(y1 ~ x), "Start the sum with the bf")
  expect_error(lf(sigma ~ z, resp = c("a", "b")), "single response name")
})

test_that("a three-response gaussian rescor fit prints and names brms's way", {
  dd <- mv_gap_data()
  fit <- frm(bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + set_rescor(TRUE),
             data = dd)
  expect_identical(dim(rescor_matrix(fit)), c(3L, 3L))
  expect_true(all(c("rescor__y1__y2", "rescor__y1__y3", "rescor__y2__y3") %in%
                    variables(fit)))
  out <- capture.output(print(fit))
  expect_true(any(grepl("Family: MV(gaussian, gaussian, gaussian)", out,
                        fixed = TRUE)))
  expect_true(any(grepl("^ +y3 ~ x", out)))
})

# ------------------------------------------------------- student rescor

test_that("student() rescor is the multivariate t, with one shared nu", {
  testthat::skip_if_not_installed("mvtnorm")
  set.seed(9)
  n <- 120
  x <- rnorm(n)
  z <- rnorm(n)
  E <- (matrix(rnorm(2 * n), n) %*% chol(matrix(c(1, 0.5, 0.5, 1), 2))) /
    sqrt(rchisq(n, 4) / 4)
  dd <- data.frame(x = x, z = z, y1 = 1 + x + E[, 1],
                   y2 = exp(0.3 * z) * E[, 2])
  fit <- frm(bf(y1 ~ x) + bf(y2 ~ x, sigma ~ z) + set_rescor(TRUE) +
               student(), data = dd)
  expect_true("nu" %in% variables(fit))
  expect_false(any(grepl("^nu_", variables(fit))))

  # the taped objective at a point off the optimum, against dmvt
  set.seed(1)
  p <- fit$opt$par + rnorm(length(fit$opt$par), 0, 0.05)
  pl <- fit$obj$env$parList(p)
  bn <- names(fit$frame$par_template$betad)
  bd <- function(nm) pl$betad[match(nm, bn)]
  mu1 <- pl$beta[1] + pl$beta[2] * x
  mu2 <- pl$beta[3] + pl$beta[4] * x
  s1 <- exp(bd("y1_sigma_(Intercept)"))
  s2 <- exp(bd("y2_sigma_(Intercept)") + bd("y2_sigma_z") * z)
  nu <- 1 + exp(bd("nu_(Intercept)"))
  rho <- us_chol_cor(pl$thetar, 2L)[1, 2]
  ref <- sum(vapply(seq_len(n), function(i) {
    S <- matrix(c(s1^2, rho * s1 * s2[i], rho * s1 * s2[i], s2[i]^2), 2)
    mvtnorm::dmvt(c(dd$y1[i], dd$y2[i]), delta = c(mu1[i], mu2[i]),
                  sigma = S, df = nu, log = TRUE)
  }, 0))
  expect_lt(abs(-fit$obj$fn(p) - ref), 1e3 * .Machine$double.eps * abs(ref))

  # the pointwise joint density a log_lik() reads sums to logLik()
  rl <- rescor_row_loglik(fit, eval_dpars(fit))
  expect_length(rl, n)
  ll <- as.numeric(logLik(fit))
  expect_lt(abs(sum(rl) - ll), 1e3 * .Machine$double.eps * abs(ll))

  # the post-fit surface is the gaussian rescor model's
  expect_identical(dim(fitted(fit)), c(as.integer(n), 4L, 2L))
  expect_identical(dim(predict(fit, ndraws = 50)), c(as.integer(n), 4L, 2L))
  expect_error(simulate(fit), "multivariate")
  expect_error(residuals(fit), "multivariate")
  gp <- get_prior(bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE) + student(),
                  data = dd)
  expect_identical(gp$resp[gp$class == "nu"], "")
})

test_that("student() rescor refuses what brms refuses", {
  dd <- mv_gap_data()
  expect_error(frm(bf(y1 ~ x, nu ~ z) + bf(y2 ~ x) + set_rescor(TRUE) +
                     student(), data = dd, dry_run = "spec"),
               "Cannot predict or fix 'nu'")
  expect_error(frm(bf(y1 ~ x, nu = 5) + bf(y2 ~ x) + set_rescor(TRUE) +
                     student(), data = dd, dry_run = "spec"),
               "Cannot predict or fix 'nu'")
  expect_error(frm(bf(y1 ~ x) + student() + bf(y2 ~ x) + gaussian() +
                     set_rescor(TRUE), data = dd, dry_run = "spec"),
               "all to be student")
  expect_error(frm(bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE) + student(),
                   data = dd,
                   prior = set_prior("gamma(2, 0.1)", class = "nu",
                                     resp = "y1")),
               "takes no resp")
})

# ------------------------------------------------ ordinal in a mv model

test_that("an ordinal response in a multivariate model is its own factor", {
  dd <- mv_gap_data()
  mv <- frm(bf(o ~ x) + cumulative() + bf(o2 ~ x) + sratio() +
              bf(y1 ~ x) + gaussian(), data = dd)
  u <- list(o = frm(bf(o ~ x) + cumulative(), data = dd),
            o2 = frm(bf(o2 ~ x) + sratio(), data = dd),
            y1 = frm(bf(y1 ~ x) + gaussian(), data = dd))
  # IDENTITY: no parameter is shared, so the joint objective at any
  # point is the sum of the univariate ones at that point's pieces
  set.seed(2)
  pm <- mv$opt$par + rnorm(length(mv$opt$par), 0, 0.1)
  lm_ <- mv$obj$env$parList(pm)
  tot <- 0
  for (r in names(u)) {
    f <- u[[r]]
    pl <- f$obj$env$parList(f$opt$par)
    for (nm in names(pl)) {
      pl[[nm]] <- if (nm %in% c("beta", "betad")) {
        src <- names(mv$frame$par_template[[nm]])
        lm_[[nm]][match(paste0(r, "_", names(f$frame$par_template[[nm]])),
                        src)]
      } else {
        lm_[[paste0(r, "_", nm)]]
      }
    }
    tot <- tot + f$obj$fn(unlist(pl))
  }
  expect_lt(abs(mv$obj$fn(pm) - tot),
            64 * .Machine$double.eps * abs(tot))

  # brms's names, one threshold block per ordinal response
  v <- variables(mv)
  expect_true(all(c("b_o_Intercept[1]", "b_o_Intercept[3]",
                    "b_o2_Intercept[2]") %in% v))
  expect_true(all(c("o_tau_raw_1", "o2_tau_raw_2") %in%
                    rownames(vcov(mv, full = TRUE))))
  fx <- fixef(mv)
  expect_true(all(c("o_Intercept[1]", "o2_Intercept[1]", "o_x") %in%
                    rownames(fx)))
  hy <- hypothesis(mv, "o_Intercept[2] > o_Intercept[1]")
  expect_equal(nrow(hy$hypothesis), 1L)

  # fitted() reads the response's own thresholds: the cumulative logit
  # probabilities written out from fixef()
  th <- fx[c("o_Intercept[1]", "o_Intercept[2]", "o_Intercept[3]"), 1]
  eta <- fx["o_x", 1] * dd$x
  cdf <- cbind(0, plogis(outer(-eta, th, "+")), 1)
  P <- cdf[, -1] - cdf[, -5]
  fo <- fitted(mv, resp = "o")
  expect_lt(max(abs(fo[, "Estimate", ] - P)), 64 * .Machine$double.eps)
  expect_identical(dim(fitted(mv)), c(nrow(dd), 4L, 4L + 3L + 1L))
  expect_identical(dim(predict(mv, resp = "o", ndraws = 20)),
                   c(nrow(dd), 4L))

  # a threshold prior names its response, as in brms
  gp <- get_prior(bf(o ~ x) + cumulative() + bf(y1 ~ x) + gaussian(),
                  data = dd)
  expect_true(any(gp$class == "Intercept" & gp$resp == "o"))
  mp <- frm(bf(o ~ x) + cumulative() + bf(y1 ~ x) + gaussian(), data = dd,
            prior = set_prior("normal(0, 0.1)", class = "Intercept",
                              resp = "o"))
  th_mp <- fixef(mp)[c("o_Intercept[1]", "o_Intercept[3]"), 1]
  expect_lt(diff(th_mp), diff(th[c(1, 3)]))
})

test_that("a shared |ID| effect across ordinal and gaussian is Laplace-exact", {
  skip_on_cran()
  set.seed(31)
  n_g <- 30
  m <- 8
  g <- factor(rep(seq_len(n_g), each = m))
  x <- rnorm(n_g * m)
  U <- matrix(rnorm(n_g * 2), n_g) %*%
    chol(matrix(c(0.64, 0.3, 0.3, 0.49), 2))
  o <- cut(0.8 * x + U[g, 1] + rlogis(n_g * m), c(-Inf, -1, 0.5, Inf),
           labels = FALSE)
  y1 <- 0.5 * x + U[g, 2] + rnorm(n_g * m, 0, 0.8)
  dd <- data.frame(x, g, o, y1)
  fit <- frm(bf(o ~ x + (1 | p | g)) + cumulative() +
               bf(y1 ~ x + (1 | p | g)) + gaussian(), data = dd)
  gi <- as.integer(g)
  nll_ref <- function(q) {
    "c" <- RTMB::ADoverload("c")
    S <- frmtmb:::us_sigma(q$theta, 2L)
    Um <- RTMB::matrix(q$u, n_g, 2)
    tau2 <- q$tr[1] + exp(q$tr[2])
    eta <- q$bo * x + Um[gi, 1]
    pr <- (o == 1) * RTMB::plogis(q$tr[1] - eta) +
      (o == 2) * (RTMB::plogis(tau2 - eta) - RTMB::plogis(q$tr[1] - eta)) +
      (o == 3) * (1 - RTMB::plogis(tau2 - eta))
    mu <- q$b1[1] + q$b1[2] * x + Um[gi, 2]
    -sum(RTMB::dmvnorm(Um, 0, S, log = TRUE)) - sum(log(pr)) -
      sum(RTMB::dnorm(y1, mu, exp(q$ls), log = TRUE))
  }
  ob <- RTMB::MakeADFun(nll_ref, list(bo = 0, tr = c(0, 0), b1 = c(0, 0),
                                      ls = 0, theta = numeric(3),
                                      u = numeric(2 * n_g)),
                        random = "u", silent = TRUE)
  pf <- fit$obj$env$parList(fit$opt$par)
  bn <- names(fit$frame$par_template$beta)
  qref <- c(pf$beta[match("o_x", bn)], pf[["o_tau_raw"]],
            pf$beta[match(c("y1_(Intercept)", "y1_x"), bn)], pf$betad,
            pf$theta)
  ref <- ob$fn(qref)
  expect_lt(abs(fit$opt$objective - ref), 1e4 * .Machine$double.eps * abs(ref))
})

test_that("a family whose extras are not per-response blocks stays refused", {
  set.seed(4)
  n <- 40
  dd <- data.frame(t = rexp(n), cc = rbinom(n, 1, 0.3), x = rnorm(n),
                   y = rnorm(n))
  expect_error(frm(bf(t | cens(cc) ~ x) + cox() + bf(y ~ x) + gaussian(),
                   data = dd, dry_run = "frame"),
               "not supported in multivariate fits yet. The ordinal")
})

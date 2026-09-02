# mixture_mvn(): multivariate gaussian mixture components
# (mclust-style model-based clustering with linear-predictor means).

sim_mvnmix_data <- function(seed = 42, n = 400, p1 = 0.4,
                            m1 = c(0, 0), m2 = c(3, 4)) {
  set.seed(seed)
  cl <- rbinom(n, 1, p1)
  L1 <- t(chol(matrix(c(1, 0.5, 0.5, 1), 2)))
  L2 <- t(chol(matrix(c(0.5, -0.2, -0.2, 0.8), 2)))
  E <- matrix(rnorm(2 * n), 2)
  Y <- t(ifelse(matrix(cl == 1, 2, n, byrow = TRUE),
                m1 + L1 %*% E, m2 + L2 %*% E))
  dd <- data.frame(row = seq_len(n))
  dd$Y <- Y
  list(dd = dd, cl = cl)
}

test_that("mixture_mvn matches direct ML", {
  sim <- sim_mvnmix_data()
  fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = sim$dd)

  Y <- sim$dd$Y
  # hand-coded bivariate normal over Cholesky-parameterized covariances:
  # p = (m1, m2, log diag / off chol 1, same for 2, logit pi1)
  ldmvn2 <- function(Y, m, S) {
    det_ <- S[1, 1] * S[2, 2] - S[1, 2]^2
    d1 <- Y[, 1] - m[1]
    d2 <- Y[, 2] - m[2]
    q <- (S[2, 2] * d1^2 - 2 * S[1, 2] * d1 * d2 + S[1, 1] * d2^2) / det_
    -log(2 * pi) - 0.5 * log(det_) - 0.5 * q
  }
  nll <- function(p) {
    L1 <- matrix(c(exp(p[5]), p[7], 0, exp(p[6])), 2)
    L2 <- matrix(c(exp(p[8]), p[10], 0, exp(p[9])), 2)
    pi1 <- stats::plogis(p[11])
    -sum(log(pi1 * exp(ldmvn2(Y, p[1:2], L1 %*% t(L1))) +
               (1 - pi1) * exp(ldmvn2(Y, p[3:4], L2 %*% t(L2)))))
  }
  # BFGS probes regions where the hand-rolled density overflows to NaN
  # before recovering; only the converged value matters here
  op <- suppressWarnings(
    stats::optim(c(0, 0, 3, 4, 0, 0, 0.5, log(0.7), log(0.8), -0.3,
                   stats::qlogis(0.4)),
                 nll, method = "BFGS",
                 control = list(reltol = 1e-13, maxit = 5000))
  )
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-5)

  # posterior class probabilities recover the simulated classes (up to
  # label swap)
  P <- mixture_probs(fit)
  expect_equal(dim(P), c(nrow(Y), 2L))
  expect_equal(unname(rowSums(P)), rep(1, nrow(Y)), tolerance = 1e-10)
  acc <- mean((P[, 1] > 0.5) == (sim$cl == 1))
  expect_gt(max(acc, 1 - acc), 0.95)

  # fitted() is the n x D mixture-mean matrix; its column means equal
  # the response column means at the gaussian-mixture ML optimum
  fv <- fitted(fit)
  expect_equal(dim(fv), dim(Y))
  expect_vector_equal(colMeans(fv), colMeans(Y), tol = 1e-4)
  expect_equal(dim(residuals(fit)), dim(Y))
})

test_that("mixture_mvn recovers the faithful clusters", {
  Y <- data.matrix(datasets::faithful)   # eruptions, waiting
  dd <- data.frame(row = seq_len(nrow(Y)))
  dd$Y <- Y
  fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)

  fx <- fixef(fit)
  m <- rbind(c(fx$mu1d1, fx$mu1d2), c(fx$mu2d1, fx$mu2d2))
  m <- m[order(m[, 1]), ]   # label-swap invariant: order by eruptions
  # reference means from the standard 2-cluster EM solution
  expect_lt(abs(m[1, 1] - 2.04), 0.3)
  expect_lt(abs(m[1, 2] - 54.5), 2)
  expect_lt(abs(m[2, 1] - 4.29), 0.3)
  expect_lt(abs(m[2, 2] - 80.0), 2)

  # class assignments agree with the eruptions < 3 threshold split
  P <- mixture_probs(fit)
  short <- Y[, 1] < 3
  agree <- mean((P[, 1] > 0.5) == short)
  expect_gt(max(agree, 1 - agree), 0.95)
})

test_that("mixture_mvn class means can depend on covariates", {
  set.seed(43)
  n <- 300
  x <- rnorm(n)
  cl <- rbinom(n, 1, 0.5)
  L <- t(chol(matrix(c(0.49, 0.15, 0.15, 0.49), 2)))
  E <- L %*% matrix(rnorm(2 * n), 2)
  Y <- cbind(ifelse(cl == 1, 0, 4) + 1.2 * x + E[1, ],
             ifelse(cl == 1, 1, 5) - 0.8 * x + E[2, ])
  dd <- data.frame(x = x)
  dd$Y <- Y

  f_x <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), data = dd)
  f_0 <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
  expect_gt(as.numeric(logLik(f_x)), as.numeric(logLik(f_0)) + 20)
  # both classes share the true slopes (loose recovery)
  fx <- fixef(f_x)
  for (nm in c("mu1d1", "mu2d1")) {
    expect_lt(abs(fx[[nm]][["x"]] - 1.2), 0.3)
  }
  for (nm in c("mu1d2", "mu2d2")) {
    expect_lt(abs(fx[[nm]][["x"]] + 0.8), 0.3)
  }
})

test_that("mixture_mvn dpar overrides and covariate gating work", {
  sim <- sim_mvnmix_data(seed = 44, n = 150)
  dd <- sim$dd
  dd$x <- rnorm(nrow(dd))

  # per-class per-dimension override
  f_ov <- frm(bf(Y ~ x, mu2d1 ~ 1) + mixture_mvn(K = 2, D = 2),
              data = dd)
  fx <- fixef(f_ov)
  expect_named(fx$mu1d1, c("(Intercept)", "x"))
  expect_named(fx$mu2d1, "(Intercept)")

  # gating: mixing weights take a full linear predictor
  f_gate <- frm(bf(Y ~ 1, theta1 ~ x) + mixture_mvn(K = 2, D = 2),
                data = dd)
  f_flat <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
  expect_named(fixef(f_gate)$theta1, c("(Intercept)", "x"))
  expect_gte(as.numeric(logLik(f_gate)),
             as.numeric(logLik(f_flat)) - 1e-6)
})

test_that("mixture_mvn validation and guards", {
  expect_error(mixture_mvn(2), "K >= 2")
  expect_error(mixture_mvn(1, 2), "K >= 2")
  expect_error(mixture_mvn(2, 1), "mixture\\(gaussian")

  sim <- sim_mvnmix_data(seed = 45, n = 80)
  dd <- sim$dd
  expect_error(frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 3), data = dd),
               "n x 3")
  expect_error(frm(bf(row ~ 1) + mixture_mvn(K = 2, D = 2), data = dd),
               "numeric matrix")
  # extra-parameter families refuse multivariate specs
  dd$z <- rnorm(nrow(dd))
  expect_error(
    frm(mvbf(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2),
             bf(z ~ 1) + gaussian()), data = dd),
    "multivariate"
  )
  # cens()/trunc() need a CDF, which a mixture density does not carry
  dd$cc <- 0
  expect_error(frm(bf(Y | cens(cc) ~ 1) + mixture_mvn(K = 2, D = 2),
                   data = dd),
               "CDF")

  fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
  expect_error(residuals(fit, type = "pearson"), "variance function")
  expect_error(simulate(fit), "no simulator")

  expect_error(mixture_mvn(2, 2, model = "VEV"),
               "unknown covariance model 'VEV'")
  expect_error(mixture_mvn(2, 2, model = "vvv"), "Supported models")
  expect_error(mixture_mvn(2, 2, model = c("EII", "VVV")),
               "unknown covariance model")
})

## mclust's volume-shape-orientation taxonomy -----------------------------

mvn_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VVV")

test_that("mixture_mvn covariance models size the template correctly", {
  sim <- sim_mvnmix_data(seed = 46, n = 120)
  # covariance parameters per model at K = 2, D = 2, plus the 4 class
  # means and the one gating intercept
  n_cov <- c(EII = 1, VII = 2, EEI = 2, VEI = 3, EVI = 3, VVI = 4,
             EEE = 3, VVV = 6)
  want <- list(
    EII = c(sigmaraw = 1), VII = c(sigmaraw1 = 1, sigmaraw2 = 1),
    EEI = c(sigmaraw = 2), VEI = c(sigmavol1 = 1, sigmavol2 = 1,
                                   sigmashape = 1),
    EVI = c(sigmavol = 1, sigmashape1 = 1, sigmashape2 = 1),
    VVI = c(sigmaraw1 = 2, sigmaraw2 = 2), EEE = c(sigmaraw = 3),
    VVV = c(sigmaraw1 = 3, sigmaraw2 = 3)
  )
  for (m in mvn_models) {
    fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = m),
               data = sim$dd)
    tpl <- fit$frame$par_template
    got <- vapply(fit$frame$extra_names, function(nm) length(tpl[[nm]]),
                  numeric(1))
    expect_equal(got, want[[m]], info = m)
    expect_equal(attr(logLik(fit), "df"), 5 + n_cov[[m]], info = m)
    # confint() labels follow the shrunken template
    cn <- rownames(confint(fit))
    expect_true(all(paste0(rep(names(want[[m]]), want[[m]]), "_",
                           unlist(lapply(want[[m]], seq_len))) %in% cn),
                info = m)
    # every class covariance is symmetric positive definite
    ex <- fit$estimates[fit$frame$extra_names]
    sg <- fit$spec$responses[[1]]$family$mix$sigma
    for (k in 1:2) {
      S <- sg(ex, k)
      expect_equal(S, t(S), info = m)
      expect_gt(min(eigen(S, only.values = TRUE)$values), 0)
      # the diagonal models carry no off-diagonal covariance
      if (m %in% c("EII", "VII", "EEI", "VEI", "EVI", "VVI")) {
        expect_equal(S[1, 2], 0, info = m)
      }
      # the equal-volume/shape models share one covariance across classes
      if (m %in% c("EII", "EEI", "EEE")) {
        expect_equal(S, sg(ex, 1), info = m)
      }
    }
  }
})

test_that("mixture_mvn covariance models agree at a shared parameter point", {
  # Every model contains the spherical-equal point Sigma_k = s^2 I, so
  # all eight likelihoods must be bit-identical there. This is the
  # nesting identity; the fitted maxima only have to be ORDERED by it.
  set.seed(3)
  n <- 200
  g <- rbinom(n, 1, 0.5)
  Y <- cbind(ifelse(g == 1, 0, 3), ifelse(g == 1, 0, 3)) +
    matrix(rnorm(2 * n, sd = 0.9), n, 2)
  dd <- data.frame(row = seq_len(n))
  dd$Y <- Y
  ls <- log(0.9)
  point <- list(
    EII = list(sigmaraw = ls),
    VII = list(sigmaraw1 = ls, sigmaraw2 = ls),
    EEI = list(sigmaraw = c(ls, ls)),
    VEI = list(sigmavol1 = ls, sigmavol2 = ls, sigmashape = 0),
    EVI = list(sigmavol = ls, sigmashape1 = 0, sigmashape2 = 0),
    VVI = list(sigmaraw1 = c(ls, ls), sigmaraw2 = c(ls, ls)),
    EEE = list(sigmaraw = c(ls, ls, 0)),
    VVV = list(sigmaraw1 = c(ls, ls, 0), sigmaraw2 = c(ls, ls, 0))
  )
  ll_point <- numeric(0)
  ll_fit <- numeric(0)
  for (m in mvn_models) {
    st <- c(list(beta = c(0, 0, 3, 3), betad = 0.2), point[[m]])
    fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = m),
               data = dd, start = st)
    v <- unlist(make_start(fit$frame, st), use.names = FALSE)
    ll_point[m] <- -fit$obj$fn(v)
    ll_fit[m] <- as.numeric(logLik(fit))
  }
  expect_vector_equal(ll_point, rep(ll_point[["EII"]], length(mvn_models)),
                      tol = 1e-10)

  # nesting chains: a submodel can never beat the model containing it
  chains <- list(c("EII", "VII", "VEI", "VVI", "VVV"),
                 c("EII", "EEI", "VEI", "VVI", "VVV"),
                 c("EII", "EEI", "EVI", "VVI", "VVV"),
                 c("EII", "EEI", "EEE", "VVV"))
  for (ch in chains) {
    for (i in seq_len(length(ch) - 1L)) {
      expect_gte(ll_fit[[ch[i + 1L]]], ll_fit[[ch[i]]] - 1e-6)
    }
  }
})

test_that("mixture_mvn class means take covariates under every model", {
  # mclust cannot do this: the class means are full linear predictors,
  # so the covariance taxonomy must not restrict them.
  set.seed(47)
  n <- 300
  x <- rnorm(n)
  cl <- rbinom(n, 1, 0.5)
  Y <- cbind(ifelse(cl == 1, 0, 4) + 1.2 * x + rnorm(n, 0, 0.7),
             ifelse(cl == 1, 1, 5) - 0.8 * x + rnorm(n, 0, 0.7))
  dd <- data.frame(x = x)
  dd$Y <- Y
  # mixture surfaces are flat near the optimum, so a couple of extra
  # optimizer restarts keep the gradient check quiet
  ctl <- frmtmb_control(restarts = 3)
  for (m in mvn_models) {
    f_x <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2, model = m),
               data = dd, control = ctl)
    f_0 <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = m),
               data = dd, control = ctl)
    expect_gt(as.numeric(logLik(f_x)), as.numeric(logLik(f_0)) + 20)
    fx <- fixef(f_x)
    for (nm in c("mu1d1", "mu2d1")) {
      expect_lt(abs(fx[[nm]][["x"]] - 1.2), 0.3)
    }
    for (nm in c("mu1d2", "mu2d2")) {
      expect_lt(abs(fx[[nm]][["x"]] + 0.8), 0.3)
    }
  }
})

# mclust's fitted covariances, encoded in our extra parameters. Used to
# start our ML at mclust's EM solution so the two optimizers cannot be
# confused for a likelihood disagreement.
mclust_extras <- function(model, Sig, K, D) {
  lsd <- vapply(Sig, function(S) 0.5 * log(diag(S)), numeric(D))
  us <- function(S) {
    c(log(sqrt(diag(S))), us_theta_cor(stats::cov2cor(S)))
  }
  per <- function(pre, v) stats::setNames(v, paste0(pre, seq_len(K)))
  free <- function(k) (lsd[, k] - mean(lsd[, k]))[seq_len(D - 1L)]
  switch(model,
    EII = list(sigmaraw = mean(lsd)),
    VII = per("sigmaraw", as.list(colMeans(lsd))),
    EEI = list(sigmaraw = rowMeans(lsd)),
    VEI = c(per("sigmavol", as.list(colMeans(lsd))),
            list(sigmashape = free(1))),
    EVI = c(list(sigmavol = mean(lsd)),
            per("sigmashape", lapply(seq_len(K), free))),
    VVI = per("sigmaraw", lapply(seq_len(K), function(k) lsd[, k])),
    EEE = list(sigmaraw = us(Sig[[1]])),
    VVV = per("sigmaraw", lapply(Sig, us))
  )
}

expect_mclust_agreement <- function(Y, K, D) {
  dd <- data.frame(row = seq_len(nrow(Y)))
  dd$Y <- Y
  # Mclust() evaluates a constructed mclustBIC() call in ITS CALLER's
  # frame, so the name has to be visible here even though we qualify
  # every call with mclust::
  mclustBIC <- mclust::mclustBIC
  for (m in mvn_models) {
    # a tight EM tolerance: mclust's default stops on a loglik change
    # of 1e-5, which is looser than the agreement we are checking
    mc <- mclust::Mclust(Y, G = K, modelNames = m, verbose = FALSE,
                         control = mclust::emControl(tol = 1e-13,
                                                     itmax = 1e6))
    Sig <- lapply(seq_len(K),
                  function(k) mc$parameters$variance$sigma[, , k])
    st <- c(list(beta = as.vector(mc$parameters$mean),
                 betad = log(mc$parameters$pro[-K] / mc$parameters$pro[K])),
            mclust_extras(m, Sig, K, D))
    fit <- frm(bf(Y ~ 1) + mixture_mvn(K = K, D = D, model = m),
               data = dd, start = st)
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(mc$loglik)), 1e-4)

    # label order is arbitrary in both fits; sort by the first mean
    fx <- fixef(fit)
    ourM <- t(vapply(seq_len(K), function(k) {
      vapply(seq_len(D), function(j) {
        unname(fx[[paste0("mu", k, "d", j)]])
      }, numeric(1))
    }, numeric(D)))
    sg <- fit$spec$responses[[1]]$family$mix$sigma
    ex <- fit$estimates[fit$frame$extra_names]
    ourS <- lapply(seq_len(K), function(k) sg(ex, k))
    mcM <- t(mc$parameters$mean)
    ord <- order(ourM[, 1])
    mo <- order(mcM[, 1])
    expect_vector_equal(ourM[ord, ], mcM[mo, ], tol = 1e-4)
    expect_vector_equal(unlist(ourS[ord]), unlist(Sig[mo]), tol = 1e-4)
    expect_vector_equal(mixture_probs(fit)[, ord], mc$z[, mo], tol = 1e-4)
  }
}

test_that("mixture_mvn matches mclust for every covariance model", {
  skip_if_not_installed("mclust")
  # intercept-only means: our model and mclust's are then the same
  # likelihood, so the fits must agree parameter for parameter
  Y <- unname(data.matrix(datasets::faithful))
  Y[, 2] <- Y[, 2] / 10       # put the columns on a common scale
  expect_mclust_agreement(Y, 2L, 2L)

  set.seed(11)
  n <- 600
  g <- sample(1:3, n, TRUE)
  mus <- rbind(c(0, 0, 0), c(3, 1, -2), c(-2, 3, 2))
  Y3 <- mus[g, ] + matrix(rnorm(3 * n, sd = 0.8), n, 3) %*%
    chol(matrix(c(1, .3, .1, .3, 1, -.2, .1, -.2, 1), 3))
  expect_mclust_agreement(Y3, 3L, 3L)
})

# The pointwise log-density of a response with a residual covariance
# term, autocor_row_loglik() (R/autocor-loglik.R): brms's leave-one-out
# conditional density of each row given the other rows of its group.
# The reference is the conditional of a dense multivariate normal (or t,
# group by group) by the partition formula, built from the fitted
# model's own covariance, with no precision matrix. Agreement with brms
# 2.23.0's log_lik() at the same parameter vectors is in frmtmb.sample's
# test-loglik-autocor.R (the ar() family and fcor()) and frmtmb.spatial's
# test-sar.R (sar()).

# the conditional of each row of y ~ N(m, S), or the multivariate t with
# scale S and nu, given every other row
ll_dense_cond <- function(y, m, S, nu = NULL) {
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

# the joint density of the same groups
ll_dense_joint <- function(y, m, S, grp, nu = NULL) {
  tot <- 0
  for (g in unique(grp)) {
    r <- which(grp == g)
    e <- y[r] - m[r]
    k <- length(r)
    Sg <- S[r, r, drop = FALSE]
    q <- sum(e * solve(Sg, e))
    ld <- as.numeric(determinant(Sg)$modulus)
    tot <- tot + if (is.null(nu)) {
      -0.5 * k * log(2 * pi) - 0.5 * ld - 0.5 * q
    } else {
      lgamma((nu + k) / 2) - lgamma(nu / 2) - 0.5 * k * log(nu * pi) -
        0.5 * ld - 0.5 * (nu + k) * log1p(q / nu)
    }
  }
  tot
}

# a model assembled for sampling, read at parameters set by hand, as
# frmtmb.sample reads one draw (its draws_base_fit())
ll_unfit <- function(f, ...) {
  fit <- frm(f, ..., dry_run = "objective")
  class(fit) <- setdiff(class(fit), "frmtmb_unfitted")
  fit
}

# ragged groups, rows shuffled, and one group missing an interior time
# point: the conditional is then the model's own, lags from the levels
ll_ac_data <- function(seed = 3) {
  set.seed(seed)
  sizes <- c(7, 5, 6, 4, 7, 3, 6, 5)
  d <- do.call(rbind, lapply(seq_along(sizes), function(g) {
    data.frame(g = factor(g), t = seq_len(sizes[g]))
  }))
  d <- d[!(d$g == 2 & d$t == 3), ]
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(nrow(d)) + stats::rnorm(8)[d$g]
  d <- d[sample(nrow(d)), ]
  rownames(d) <- NULL
  d
}

# the fit's covariance over its rows, group by group, from its own
# correlation matrix over the time levels and its sigma
ll_ar_model <- function(fit, d) {
  dp <- eval_dpars(fit)[[1]]
  n <- nrow(d)
  s <- rep_len(as.numeric(dp$sigma), n)
  R <- autocor_matrix(fit)
  ti <- match(as.character(d$t), rownames(R))
  S <- matrix(0, n, n)
  for (g in unique(d$g)) {
    r <- which(d$g == g)
    S[r, r] <- outer(s[r], s[r]) * R[ti[r], ti[r]]
  }
  list(m = rep_len(as.numeric(dp$mu), n), S = S,
       nu = if (!is.null(dp$nu)) as.numeric(dp$nu)[1])
}

ll_grouped_ref <- function(fit, d) {
  md <- ll_ar_model(fit, d)
  out <- numeric(nrow(d))
  for (g in unique(d$g)) {
    r <- which(d$g == g)
    out[r] <- ll_dense_cond(d$y[r], md$m[r], md$S[r, r, drop = FALSE],
                            md$nu)
  }
  out
}

test_that("the ar() family's rows are the dense conditionals of their
           groups", {
  d <- ll_ac_data()
  forms <- list(y ~ x + ar(t, g, cov = TRUE),
                y ~ x + ar(t, g, p = 2, cov = TRUE),
                y ~ x + ma(t, g, cov = TRUE),
                y ~ x + arma(t, g, cov = TRUE),
                y ~ x + cosy(t, g),
                y ~ x + ar(t, g, cov = TRUE) + (1 | g),
                bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x))
  for (f in forms) {
    for (fam in list(gaussian(), student())) {
      if (inherits(f, "brmsformula") && fam$family == "student") next
      fit <- frm(f, family = fam, data = d)
      ll <- autocor_row_loglik(fit, eval_dpars(fit), "y")
      ref <- ll_grouped_ref(fit, d)
      expect_length(ll, nrow(d))
      expect_lt(max(abs(ll - ref)), 1e-10 * max(abs(ref)))
    }
  }
})

test_that("unstr() rows are the dense conditionals at chosen
           correlations", {
  # at the estimates seven levels over eight groups are near singular,
  # so the correlations are set by hand
  d <- ll_ac_data()
  for (fam in list(gaussian(), student())) {
    fit <- ll_unfit(y ~ x + unstr(t, g), family = fam, data = d)
    ac <- fit$frame$autocor$y
    set.seed(11)
    fit$estimates$beta[] <- c(1, 0.5)
    fit$estimates$betad[] <- c(log(0.9), if (fam$family == "student")
      log(5 - 1))
    fit$estimates$thetaac[] <- stats::rnorm(ac$npar, 0, 0.5)
    ll <- autocor_row_loglik(fit, eval_dpars(fit), "y")
    ref <- ll_grouped_ref(fit, d)
    expect_lt(max(abs(ll - ref)), 1e-10 * max(abs(ref)))
    if (fam$family == "student") {
      expect_equal(as.numeric(eval_dpars(fit)$y$nu)[1], 5)
    }
  }
})

test_that("fcor() rows are the dense conditionals, with sigma ~ x", {
  d <- ll_ac_data()
  n <- nrow(d)
  V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)
  for (f in list(y ~ x + fcor(V), bf(y ~ x + fcor(V), sigma ~ x))) {
    fit <- frm(f, data = d, data2 = list(V = V))
    dp <- eval_dpars(fit)$y
    s <- rep_len(as.numeric(dp$sigma), n)
    ref <- ll_dense_cond(d$y, as.numeric(dp$mu), outer(s, s) * V)
    ll <- autocor_row_loglik(fit, eval_dpars(fit), "y")
    expect_lt(max(abs(ll - ref)), 1e-10 * max(abs(ref)))
  }
  fs <- ll_unfit(y ~ x + fcor(V), family = student(), data = d,
            data2 = list(V = V))
  fs$estimates$beta[] <- c(1, 0.5)
  fs$estimates$betad[] <- c(log(1.1), log(4 - 1))
  ll <- autocor_row_loglik(fs, eval_dpars(fs), "y")
  ref <- ll_dense_cond(d$y, 1 + 0.5 * d$x, 1.1^2 * V, nu = 4)
  expect_lt(max(abs(ll - ref)), 1e-10 * max(abs(ref)))
})

test_that("the row sum is not the joint density, except without
           correlation", {
  d <- ll_ac_data()
  d <- d[order(d$g, d$t), ]
  fit <- ll_unfit(y ~ x + ar(t, g, cov = TRUE), data = d)
  fit$estimates$beta[] <- c(1, 0.5)
  fit$estimates$betad[] <- log(1.2)
  for (rho in c(0.8, 0)) {
    fit$estimates$thetaac[] <- rho / sqrt(1 - rho^2)
    ll <- autocor_row_loglik(fit, eval_dpars(fit), "y")
    md <- ll_ar_model(fit, d)
    joint <- ll_dense_joint(d$y, md$m, md$S, d$g)
    # the objective is the joint density, as the reference says
    expect_lt(abs(-as.numeric(build_objective(fit$frame)(fit$estimates)) -
                    joint), 1e-10 * abs(joint))
    if (rho > 0) {
      # each row conditions on its group's other rows, so the columns
      # are not factors of the likelihood
      expect_gt(abs(sum(ll) - joint), 0.01 * abs(joint))
    } else {
      # no correlation: each conditional is the row's own density
      expect_lt(abs(sum(ll) - joint), 1e-10 * abs(joint))
      expect_lt(max(abs(ll - stats::dnorm(d$y, 1 + 0.5 * d$x, 1.2,
                                          log = TRUE))),
                1e-10 * abs(joint))
    }
  }
})

seam_fix_term <- function() {
  list(
    build = function(spec, rows, n_all, lookup) {
      V <- lookup(spec$V, "llfix()")[rows, rows, drop = FALSE]
      L <- t(chol(V))
      list(L = L, Linv = forwardsolve(L, diag(nrow(L))),
           logdet = sum(log(diag(L))))
    },
    whiten = function(y, mu, sigma, theta, aux) {
      list(u = as.vector(aux$Linv %*% ((y - mu) / sigma)),
           logjac = -aux$logdet - length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      mu + sigma * (aux$L %*% u)
    })
}

test_that("a registered term's jacobian slot and the fallback agree with
           fcor()", {
  local_frmtmb_registries()
  tt <- seam_fix_term()
  d <- ll_ac_data()
  n <- nrow(d)
  V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)
  frmtmb_register_autocor("llfix", args = "V", build = tt$build,
                          whiten = tt$whiten, unwhiten = tt$unwhiten,
                          sigma_varies = FALSE)
  frmtmb_register_autocor("llfixj", args = "V", build = tt$build,
                          whiten = tt$whiten, unwhiten = tt$unwhiten,
                          jacobian = function(sigma, theta, aux) {
                            aux$Linv / sigma
                          },
                          sigma_varies = FALSE)
  f0 <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
  # each fit at fcor()'s estimates, so the three read one parameter
  # vector
  f1 <- frm(y ~ x + llfix(V), data = d, data2 = list(V = V))
  f2 <- frm(y ~ x + llfixj(V), data = d, data2 = list(V = V))
  f1$estimates <- f0$estimates
  f2$estimates <- f0$estimates
  l0 <- autocor_row_loglik(f0, eval_dpars(f0), "y")
  l1 <- autocor_row_loglik(f1, eval_dpars(f1), "y")
  l2 <- autocor_row_loglik(f2, eval_dpars(f2), "y")
  # the slot is fcor()'s matrix, up to dividing by sigma or multiplying
  # by its reciprocal
  expect_lt(max(abs(l2 - l0)), 1e-12 * max(abs(l0)))
  # the fallback inverts unwhiten() of the identity: the same values
  # through another route
  expect_lt(max(abs(l1 - l0)), 1e-10 * max(abs(l0)))
  # the absent case: a fit at other parameters gives other values
  l3 <- autocor_row_loglik(f0, list(y = list(mu = eval_dpars(f0)$y$mu,
                                             sigma = 2)), "y")
  expect_gt(max(abs(l3 - l0)), 0.01 * max(abs(l0)))
})

test_that("a jacobian that is not the map's is refused by name", {
  local_frmtmb_registries()
  tt <- seam_fix_term()
  d <- ll_ac_data()
  n <- nrow(d)
  V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)
  frmtmb_register_autocor("llbad", args = "V", build = tt$build,
                          whiten = tt$whiten, unwhiten = tt$unwhiten,
                          # the transpose: right shape, wrong matrix
                          jacobian = function(sigma, theta, aux) {
                            t(aux$Linv) / sigma
                          },
                          sigma_varies = FALSE)
  expect_error(frm(y ~ x + llbad(V), data = d, data2 = list(V = V)),
               "jacobian() must give the n x n matrix K", fixed = TRUE,
               class = "frmtmb_error")
  frmtmb_register_autocor("llshape", args = "V", build = tt$build,
                          whiten = tt$whiten, unwhiten = tt$unwhiten,
                          jacobian = function(sigma, theta, aux) 1 / sigma,
                          sigma_varies = FALSE)
  expect_error(frm(y ~ x + llshape(V), data = d, data2 = list(V = V)),
               "jacobian() must give the n x n matrix K", fixed = TRUE,
               class = "frmtmb_error")
  expect_error(frmtmb_register_autocor("llarg", args = "V",
                                       build = tt$build,
                                       whiten = tt$whiten,
                                       unwhiten = tt$unwhiten,
                                       jacobian = function(sigma) sigma),
               "jacobian =) must be a function of 3 arguments",
               fixed = TRUE)
  # the guard absent: fcor()'s own slot passes the same check
  expect_s3_class(frm(y ~ x + fcor(V), data = d, data2 = list(V = V)),
                  "frmtmb_fit")
})

test_that("autocor_row_loglik() refuses a response with no covariance
           term", {
  d <- ll_ac_data()
  fit <- frm(y ~ x + ar(t, g), data = d)
  expect_error(autocor_row_loglik(fit, eval_dpars(fit), "y"),
               "no residual covariance term", class = "frmtmb_error")
})

# A jacobian slot that is right at the start (sigma = 1, the start
# parameter) and wrong elsewhere: the four slips of the fix round's
# review (dev/reviews/2026-10-27-loglik.md, item 1), each accepted by
# the frame check of the first version and wrong in log_lik() by 9 to 44
# percent. The check now also runs at a generic point.
slip_fix <- function() {
  list(
    build = function(spec, rows, n_all, lookup) {
      V <- lookup(spec$V, "slip()")[rows, rows, drop = FALSE]
      L <- t(chol(V))
      list(L = L, Linv = forwardsolve(L, diag(nrow(L))),
           logdet = sum(log(diag(L))))
    },
    # sigma_varies = TRUE: one sigma per row
    whiten = function(y, mu, sigma, theta, aux) {
      list(u = as.vector(aux$Linv %*% ((y - mu) / sigma)),
           logjac = -aux$logdet - sum(log(rep_len(sigma, length(y)))))
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      mu + sigma * (aux$L %*% u)
    })
}
# a map with a parameter, rho = tanh(theta), started at rho = 0
slip_rho <- function() {
  list(
    build = function(spec, rows, n_all, lookup) {
      list(W = lookup(spec$W, "slip()")[rows, rows, drop = FALSE])
    },
    whiten = function(y, mu, sigma, theta, aux) {
      A <- diag(length(y)) - tanh(theta[1]) * aux$W
      list(u = as.vector(A %*% (y - mu)) / sigma,
           logjac = as.numeric(determinant(A)$modulus) -
             length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      A <- diag(nrow(aux$W)) - tanh(theta[1]) * aux$W
      mu + sigma * solve(A, u)
    })
}

test_that("a jacobian slot wrong away from the start point is refused", {
  local_frmtmb_registries()
  set.seed(81)
  n <- 40
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n)
  V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)
  P <- matrix(0, n, n)
  P[cbind(1:(n - 1), 2:n)] <- 1
  P <- P + t(P)
  # row-standardized: asymmetric, so K and its transpose differ
  W <- P / rowSums(P)
  tf <- slip_fix()
  tr <- slip_rho()
  reg_f <- function(name, jac) {
    frmtmb_register_autocor(name, args = "V", build = tf$build,
                            whiten = tf$whiten, unwhiten = tf$unwhiten,
                            sigma_varies = TRUE, jacobian = jac)
  }
  reg_r <- function(name, jac) {
    frmtmb_register_autocor(name, args = "W", build = tr$build,
                            whiten = tr$whiten, unwhiten = tr$unwhiten,
                            npar = 1L, start = function(aux) 0,
                            natural = function(theta, aux) {
                              c(rho = tanh(theta[1]))
                            },
                            sigma_varies = FALSE, jacobian = jac)
  }
  # forgets to divide by sigma
  reg_f("slipnosig", function(sigma, theta, aux) aux$Linv)
  # divides the rows by sigma, where a sigma per row divides the columns
  reg_f("sliprows", function(sigma, theta, aux) aux$Linv / sigma)
  # ignores rho
  reg_r("slipnorho", function(sigma, theta, aux) diag(nrow(aux$W)) / sigma)
  # the transpose of the right matrix
  reg_r("sliptrans", function(sigma, theta, aux) {
    t(diag(nrow(aux$W)) - tanh(theta[1]) * aux$W) / sigma
  })
  # the right ones: the guard's absent cases
  reg_f("slipokf", function(sigma, theta, aux) {
    aux$Linv * rep(1 / rep_len(sigma, n), each = n)
  })
  reg_r("slipokr", function(sigma, theta, aux) {
    (diag(nrow(aux$W)) - tanh(theta[1]) * aux$W) / sigma
  })
  d2 <- list(V = V, W = W)
  for (nm in c("slipnosig(V)", "sliprows(V)", "slipnorho(W)",
               "sliptrans(W)")) {
    expect_error(frm(stats::as.formula(paste("y ~ x +", nm)), data = d,
                     data2 = d2),
                 "jacobian() must give the n x n matrix K", fixed = TRUE,
                 class = "frmtmb_error")
  }
  ff <- frm(bf(y ~ x + slipokf(V), sigma ~ x), data = d, data2 = d2)
  fr <- frm(y ~ x + slipokr(W), data = d, data2 = d2)
  expect_s3_class(ff, "frmtmb_fit")
  expect_s3_class(fr, "frmtmb_fit")
  # and the accepted slots give the dense conditionals
  dp <- eval_dpars(ff)$y
  s <- as.numeric(dp$sigma)
  ref <- ll_dense_cond(d$y, as.numeric(dp$mu), outer(s, s) * V)
  ll <- autocor_row_loglik(ff, eval_dpars(ff), "y")
  expect_lt(max(abs(ll - ref)), 1e-10 * max(abs(ref)))
  rho <- tanh(fr$estimates$thetaac)
  A <- diag(n) - rho * W
  dp <- eval_dpars(fr)$y
  ref <- ll_dense_cond(d$y, as.numeric(dp$mu),
                       as.numeric(dp$sigma)[1]^2 * solve(crossprod(A)))
  ll <- autocor_row_loglik(fr, eval_dpars(fr), "y")
  expect_lt(max(abs(ll - ref)), 1e-10 * max(abs(ref)))
})

test_that("a term without a jacobian slot warns above 1000 rows", {
  local_frmtmb_registries()
  tr <- slip_rho()
  frmtmb_register_autocor("nojac", args = "W", build = tr$build,
                          whiten = tr$whiten, unwhiten = tr$unwhiten,
                          npar = 1L, start = function(aux) 0,
                          natural = function(theta, aux) {
                            c(rho = tanh(theta[1]))
                          },
                          sigma_varies = FALSE)
  mk <- function(n) {
    set.seed(2)
    P <- matrix(0, n, n)
    P[cbind(1:(n - 1), 2:n)] <- 1
    P <- P + t(P)
    d <- data.frame(x = stats::rnorm(n))
    d$y <- d$x + stats::rnorm(n)
    fit <- ll_unfit(y ~ x + nojac(W), data = d,
                    data2 = list(W = P / rowSums(P)))
    fit$estimates$beta[] <- c(0, 1)
    fit$estimates$betad[] <- 0
    fit$estimates$thetaac[] <- 0.3
    fit
  }
  big <- mk(1001)
  allow_warnings(
    v <- autocor_row_loglik(big, eval_dpars(big), "y"),
    "has no jacobian slot", require = "has no jacobian slot")
  expect_true(all(is.finite(v)))
  w <- tryCatch(autocor_row_loglik(big, eval_dpars(big), "y"),
                warning = identity)
  expect_s3_class(w, "frmtmb_autocor_dense_fallback")
  expect_match(conditionMessage(w), "jacobian =", fixed = TRUE)
  expect_match(conditionMessage(w), "nojac(W)", fixed = TRUE)
  # the absent case: 1000 rows take the same route without a word
  small <- mk(1000)
  expect_no_warning(autocor_row_loglik(small, eval_dpars(small), "y"))
})

# a parameter that can only move down from its start: rho = theta on
# (-1, 0], started at 0, a map that stops for theta > 0 on numbers. The
# fix round 2 review (dev/reviews/2026-10-27-loglik.md, "Fix round 1
# check") found that the generic point stepped upward only, so a slot
# that ignores rho passed and log_lik() was 35 percent off.
slip_oneside <- function(jac) {
  rho1 <- function(theta) {
    if (!inherits(theta, "advector") && theta[1] > 0) {
      stop("rho must be <= 0")
    }
    theta[1]
  }
  frmtmb_register_autocor(
    "oneside", args = "W",
    build = function(spec, rows, n_all, lookup) {
      W <- lookup(spec$W, "oneside()")[rows, rows, drop = FALSE]
      list(W = W, ev = Re(eigen(W, only.values = TRUE)$values))
    },
    whiten = function(y, mu, sigma, theta, aux) {
      A <- diag(length(y)) - rho1(theta) * aux$W
      list(u = as.vector(A %*% (y - mu)) / sigma,
           logjac = sum(log1p(-rho1(theta) * aux$ev)) -
             length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      mu + sigma * solve(diag(nrow(aux$W)) - rho1(theta) * aux$W, u)
    },
    npar = 1L, start = function(aux) 0,
    natural = function(theta, aux) c(rho = theta[1]),
    sigma_varies = FALSE, jacobian = jac, replace = TRUE)
}

test_that("the generic point also steps a parameter down from its start", {
  local_frmtmb_registries()
  set.seed(101)
  n <- 40
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n)
  P <- matrix(0, n, n)
  P[cbind(1:(n - 1), 2:n)] <- 1
  P <- P + t(P)
  W <- P / rowSums(P)
  # the slot ignores rho: right at the start only
  slip_oneside(function(sigma, theta, aux) diag(nrow(aux$W)) / sigma)
  expect_error(ll_unfit(y ~ x + oneside(W), data = d, data2 = list(W = W)),
               "jacobian() must give the n x n matrix K", fixed = TRUE,
               class = "frmtmb_error")
  # the absent case: the same map with the right slot is accepted, and
  # gives the dense conditional at rho = -0.5
  slip_oneside(function(sigma, theta, aux) {
    (diag(nrow(aux$W)) - theta[1] * aux$W) / sigma
  })
  fit <- ll_unfit(y ~ x + oneside(W), data = d, data2 = list(W = W))
  fit$estimates$beta[] <- c(1, 0.5)
  fit$estimates$betad[] <- log(0.9)
  fit$estimates$thetaac[] <- -0.5
  A <- diag(n) + 0.5 * W
  ref <- ll_dense_cond(d$y, 1 + 0.5 * d$x, 0.81 * solve(crossprod(A)))
  ll <- autocor_row_loglik(fit, eval_dpars(fit), "y")
  expect_lt(max(abs(ll - ref)), 1e-10 * max(abs(ref)))
})

test_that("a slot whose map cannot move off its start is checked there
           only, and says so", {
  local_frmtmb_registries()
  set.seed(102)
  n <- 30
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n)
  P <- matrix(0, n, n)
  P[cbind(1:(n - 1), 2:n)] <- 1
  P <- P + t(P)
  W <- P / rowSums(P)
  # a map defined at its start alone: every moved point stops
  rho0 <- function(theta) {
    if (!inherits(theta, "advector") && theta[1] != 0) stop("rho is 0")
    theta[1]
  }
  frmtmb_register_autocor(
    "pinned", args = "W",
    build = function(spec, rows, n_all, lookup) {
      list(W = lookup(spec$W, "pinned()")[rows, rows, drop = FALSE])
    },
    whiten = function(y, mu, sigma, theta, aux) {
      A <- diag(length(y)) - rho0(theta) * aux$W
      list(u = as.vector(A %*% (y - mu)) / sigma,
           logjac = as.numeric(determinant(A)$modulus) -
             length(y) * log(sigma))
    },
    unwhiten = function(u, mu, sigma, theta, aux) {
      mu + sigma * solve(diag(nrow(aux$W)) - rho0(theta) * aux$W, u)
    },
    npar = 1L, start = function(aux) 0,
    natural = function(theta, aux) c(rho = theta[1]),
    sigma_varies = FALSE,
    jacobian = function(sigma, theta, aux) {
      (diag(nrow(aux$W)) - theta[1] * aux$W) / sigma
    })
  w <- tryCatch(ll_unfit(y ~ x + pinned(W), data = d, data2 = list(W = W)),
                warning = identity)
  expect_s3_class(w, "frmtmb_autocor_jacobian_unchecked")
  expect_match(conditionMessage(w), "only at the start", fixed = TRUE)
  expect_match(conditionMessage(w), "pinned(W)", fixed = TRUE)
})

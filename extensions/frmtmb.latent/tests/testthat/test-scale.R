## Phase 0 of dev/extension-gaps-plan.md: the hmm and lca rows.
##
## hmm: 50 sequences x 500 steps, K = 3 gaussian, `tr12 ~ (1 | id)`.
## 25,000 rows. Tape-build growth past T = 20,000 was measured in the
## probe, so what this row adds is the POST-FIT passes, hmm_probs() and
## hmm_viterbi(), which are R loops over sequences and are timed
## separately here.
##
## lca: n = 2000, 10 binary items, K = 4, two covariates on membership.
## The plan calls this row cheap and includes it for completeness.
##
## Neither row asserts recovery of a state-indexed parameter, because a
## finite mixture is only identified up to a relabeling and Phase 0 is a
## cost measurement. The estimates against the truth are RECORDED, and
## Phase 2's items 2.3 and 2.4 are where they are validated against
## depmixS4, hmmTMB and poLCA.
##
## See dev/scale-findings.md for the numbers this produced.

latent_truth <- list(
  mu = c(-2, 0, 3), sigma = 0.7,
  # a multinomial logit per row of the transition matrix, with state 1
  # as the reference destination
  eta = rbind(c(-1.5, -2.5), c(2.0, -1.0), c(-1.0, 2.0)),
  sd_tr12 = 0.6,
  # lca: 10 binary items, 4 classes, two covariates on membership
  lca_n = 2000L, lca_J = 10L, lca_K = 4L)

latent_ci <- function(ci, key) {
  j <- grep(key, rownames(ci), fixed = TRUE)
  if (!length(j)) return(c(NA_real_, NA_real_))
  as.numeric(ci[j[1L], 1:2])
}

# 50 sequences of 500 steps. Sequence `i` has its own first-row
# transition, which is the random effect the row exists to price.
latent_hmm_data <- function(seed = 20260908L,
                            ns = if (scale_small()) 5L else 50L,
                            tl = if (scale_small()) 40L else 500L) {
  tr <- latent_truth
  set.seed(seed)
  u <- stats::rnorm(ns, 0, tr$sd_tr12)
  out <- vector("list", ns)
  for (i in seq_len(ns)) {
    e <- tr$eta
    e[1L, 1L] <- e[1L, 1L] + u[i]
    G <- t(apply(cbind(0, e), 1L, function(z) exp(z) / sum(exp(z))))
    s <- integer(tl)
    s[1L] <- sample.int(3L, 1L)
    for (t in seq_len(tl)[-1L]) {
      s[t] <- sample.int(3L, 1L, prob = G[s[t - 1L], ])
    }
    out[[i]] <- data.frame(id = i, t = seq_len(tl), state = s,
                           y = stats::rnorm(tl, tr$mu[s], tr$sigma))
  }
  d <- do.call(rbind, out)
  d$id <- factor(d$id)
  d
}

# n = 2000, 10 binary items, 4 classes, two covariates on membership.
latent_lca_data <- function(seed = 20260908L) {
  tr <- latent_truth
  n <- if (scale_small()) 400L else tr$lca_n
  J <- tr$lca_J
  K <- tr$lca_K
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1L, 0.5)
  # a multinomial logit on membership with class 1 as the reference
  gam <- rbind(c(0, 0, 0), c(-0.4, 0.8, -0.5), c(0.2, -0.6, 0.9),
               c(-0.1, 0.3, 0.4))
  eta <- cbind(1, x1, x2) %*% t(gam)
  pr <- exp(eta) / rowSums(exp(eta))
  cl <- apply(pr, 1L, function(p) sample.int(K, 1L, prob = p))
  # each class endorses a different block of the ten items
  base <- matrix(0.2, K, J)
  for (k in seq_len(K)) {
    base[k, ((k - 1L) * 2L + 1L):((k - 1L) * 2L + 3L)] <- 0.85
  }
  Y <- matrix(0L, n, J)
  for (j in seq_len(J)) Y[, j] <- 1L + stats::rbinom(n, 1L, base[cl, j])
  dd <- data.frame(x1 = x1, x2 = factor(x2))
  dd$Y <- Y
  list(dd = dd, cl = cl, base = base)
}

# The shared measurement: tape build, one gradient at the start and one
# at the optimum, the whole frm() call including sdreport(), the peak R
# heap, and whatever post-fit pass the row is really about.
latent_scale_run <- function(row, form, fam, d, post = NULL) {
  scale_mem_reset()
  bd <- scale_build(form, family = fam, data = d)
  g0 <- scale_grad(bd$dry$obj, bd$dry$obj$par)
  ctl <- scale_control(bd$dry$obj, bd$dry$obj$par, g0$calls)
  bd$dry <- NULL

  fit <- NULL
  t_fit <- scale_elapsed(fit <- frm(form, family = fam, data = d,
                                    se = TRUE))
  g1 <- scale_grad(fit$obj, fit$opt$par)
  mem <- scale_mem_peak_mb()
  extra <- if (is.null(post)) list() else post(fit)

  args <- c(list(row, rows = nrow(d), n_par = length(fit$opt$par),
                 build_s = bd$build_s, frame_s = bd$frame_s,
                 grad_start_s = g0$seconds, grad_start_calls = g0$calls,
                 grad_opt_s = g1$seconds, control_ratio = ctl,
                 fit_s = t_fit, mem_mb = mem,
                 logLik = as.numeric(stats::logLik(fit))),
            extra, list(diag = scale_diag(fit)))
  do.call(scale_record, args)
  fit
}

test_that("the hmm scale row fits and reports its cost", {
  skip_unless_scale()
  scale_row_on("hmm")
  d <- latent_hmm_data()
  form <- bf(y ~ 1, tr12 ~ 1 + (1 | id))
  # init = "uniform": the default stationary initial distribution is
  # refused when a transition carries a predictor, and the simulator
  # above draws the first state uniformly, so uniform is also the
  # correct model rather than only the allowed one.
  fam <- hmm(K = 3, gaussian(), time = t, group = id, init = "uniform")
  post <- function(fit) {
    pf <- scale_interleave(list(probs = function() hmm_probs(fit),
                                viterbi = function() hmm_viterbi(fit)))
    b <- unlist(fixef(fit))
    ci <- suppressWarnings(stats::confint(fit))
    i12 <- latent_ci(ci, "tr12_(Intercept)")
    mu <- sort(c(unname(b["mu1.(Intercept)"]), unname(b["mu2.(Intercept)"]),
                 unname(b["mu3.(Intercept)"])))
    list(probs_first_s = pf$first[["probs"]],
         probs_s = pf$seconds[["probs"]],
         viterbi_first_s = pf$first[["viterbi"]],
         viterbi_s = pf$seconds[["viterbi"]],
         probs_spread = pf$spread[["probs"]], rounds = pf$rounds,
         mu_sorted = paste(formatC(mu, digits = 4, format = "g"),
                           collapse = ";"),
         mu_true = paste(latent_truth$mu, collapse = ";"),
         tr12 = unname(b["tr12.(Intercept)"]),
         tr12_true = latent_truth$eta[1L, 1L],
         tr12_lo = i12[1L], tr12_hi = i12[2L],
         sd_tr12 = sqrt(VarCorr(fit)[[1L]][1L, 1L]),
         sd_tr12_true = latent_truth$sd_tr12)
  }
  fit <- latent_scale_run("hmm", form, fam, d, post)
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  # NOT expect_identical(convergence, 0L): nlminb returns 1 on
  # perfectly good fits, as this tier's own ode row does, and the code
  # can differ with a different BLAS. The pair that means "converged"
  # is a positive definite Hessian and a gradient small relative to
  # something the run itself measures, so the gradient is scaled by the
  # log-likelihood the fit reached.
  d1 <- frmtmb::diagnose(fit, quiet = TRUE)
  expect_true(isTRUE(d1$pdHess))
  expect_lt(d1$max_grad / abs(as.numeric(stats::logLik(fit))), 1e-3)
})

test_that("the lca scale row fits and reports its cost", {
  skip_unless_scale()
  scale_row_on("lca")
  s <- latent_lca_data()
  form <- bf(Y ~ x1 + x2)
  fam <- lca(K = latent_truth$lca_K)
  post <- function(fit) {
    tm <- scale_interleave(list(probs = function() lca_probs(fit),
                                profiles = function() lca_profiles(fit)))
    pr <- lca_probs(fit)
    pf <- lca_profiles(fit)
    # the classes come back in whatever order the optimizer found them,
    # so the recorded number is the best matching of fitted profiles to
    # true ones, as a maximum absolute error over the 4 x 10 table
    # lca_profiles() returns one K x ncat table PER ITEM, so the
    # endorsement probabilities stack as columns into a K x J matrix
    P <- vapply(pf, function(m) m[, 2L], numeric(nrow(pf[[1L]])))
    perms <- as.matrix(expand.grid(1:4, 1:4, 1:4, 1:4))
    perms <- perms[apply(perms, 1L, function(p) length(unique(p)) == 4L), ]
    err <- min(apply(perms, 1L, function(p) {
      max(abs(P[p, , drop = FALSE] - s$base))
    }))
    list(probs_first_s = tm$first[["probs"]],
         probs_s = tm$seconds[["probs"]],
         profiles_s = tm$seconds[["profiles"]], rounds = tm$rounds,
         profile_max_err = err,
         class_share = paste(formatC(sort(colMeans(pr)), digits = 3,
                                     format = "g"), collapse = ";"))
  }
  fit <- latent_scale_run("lca", form, fam, s$dd, post)
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  # see the hmm row above for why this is not a status-code assertion
  d1 <- frmtmb::diagnose(fit, quiet = TRUE)
  expect_true(isTRUE(d1$pdHess))
  expect_lt(d1$max_grad / abs(as.numeric(stats::logLik(fit))), 1e-3)
})

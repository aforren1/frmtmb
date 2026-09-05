# Pins vignette("reinforcement-learning"). The family under test is the
# one the vignette shows, read from inst/rl/rw-delta.R by helper-rl.R,
# so a change to the vignette's code has to keep these passing.

# One dataset, built once: the fits below are the slow part and the
# design is the same for all of them.
rl_fixture <- local({
  cache <- new.env(parent = emptyenv())
  function(n_subj = 20L, n_trial = 60L, seed = 42L, drop = FALSE) {
    key <- paste(n_subj, n_trial, seed, drop, sep = "-")
    if (!is.null(cache[[key]])) return(cache[[key]])
    set.seed(seed)
    dd <- rl_simulate(rl_bandit_design(n_subj, n_trial), nsim = 1,
                      seed = seed)[[1L]]
    if (drop) {
      # unequal trial counts, so the block pads and the mask is exercised
      keep <- unlist(lapply(split(seq_len(nrow(dd)), dd$id), function(r) {
        r[seq_len(sample(seq(n_trial %/% 2L, n_trial), 1L))]
      }))
      dd <- dd[sort(keep), ]
    }
    fit <- suppressWarnings(
      frm(rl_bform, family = rl_family(), data = dd))
    cache[[key]] <- list(data = dd, fit = fit)
    cache[[key]]
  }
})

# The recursion written the way the model is stated, one subject and one
# trial at a time, over quantities rebuilt from the public accessors.
# Nothing here shares a line with the taped version.
rl_reference <- function(fit, data) {
  u <- ranef(fit)[[1L]]
  sid <- match(as.character(data$id), rownames(u))
  xa <- unname(stats::model.matrix(~ condition, data))
  a <- stats::plogis(as.vector(xa %*% unname(fixef(fit)$alpha)) + u[sid, 1L])
  b <- exp(unname(fixef(fit)$beta)[[1L]] + u[sid, 2L])
  ll <- 0
  for (s in unique(sid)) {
    r <- which(sid == s)
    r <- r[order(data$trial[r])]
    q <- c(0, 0)
    for (i in r) {
      p <- stats::plogis(b[i] * (q[1L] - q[2L]))
      ll <- ll + stats::dbinom(data$choice[i], 1L, p, log = TRUE)
      k <- if (data$choice[i] == 1) 1L else 2L
      pay <- if (k == 1L) data$pay1[i] else data$pay2[i]
      q[k] <- q[k] + a[i] * (pay - q[k])
    }
  }
  # the subject effects' own density, so the total is the JOINT one the
  # objective holds at the modes
  sig <- unname(VarCorr(fit)[[1L]])
  si <- solve(sig)
  ldet <- as.numeric(determinant(sig, logarithm = TRUE)$modulus)
  re <- sum(vapply(seq_len(nrow(u)), function(k) {
    v <- as.numeric(u[k, ])
    -0.5 * (2 * log(2 * pi) + ldet + sum(v * (si %*% v)))
  }, 0))
  list(data = ll, joint = ll + re)
}

test_that("the block turns a long data frame into subject by trial", {
  skip_unless_rl()
  fx <- rl_fixture(drop = TRUE)
  blk <- frame_block_of(fx$fit$frame, "choice")
  expect_identical(dim(blk[["idx"]]), dim(blk[["mask"]]))
  expect_identical(blk[["n_subj"]], nlevels(fx$data$id))
  expect_identical(blk[["n_trial"]], max(blk[["len"]]))
  # one cell per real trial, and every real cell points at its own row
  expect_equal(sum(blk[["mask"]]), nrow(fx$data))
  expect_true(any(blk[["mask"]] == 0))
  real <- blk[["idx"]][blk[["mask"]] == 1]
  expect_setequal(real, seq_len(nrow(fx$data)))
  # each row of idx is that subject's rows in trial order
  for (s in seq_len(blk[["n_subj"]])) {
    r <- blk[["idx"]][s, seq_len(blk[["len"]][s])]
    expect_false(is.unsorted(fx$data$trial[r]))
    expect_equal(unique(as.character(fx$data$id[r])), blk[["levels"]][s])
  }
})

test_that("one trial per subject keeps the subjects apart", {
  skip_unless_rl()
  # The regression this pins: `vapply()` drops to a plain vector when
  # the longest subject has one trial, and `t()` of a vector is
  # 1-by-n_subj, so the block came back transposed and the recursion
  # read six subjects as one learner with six trials, carrying Q across
  # subject boundaries. It fitted without complaint and gave a wrong
  # likelihood.
  set.seed(3)
  dd <- rl_simulate(rl_bandit_design(6, 1), nsim = 1, seed = 3)[[1L]]
  fr <- frm(rl_bform, family = rl_family(), data = dd, dry_run = "frame")
  blk <- frame_block_of(fr, "choice")
  expect_identical(dim(blk[["idx"]]), c(6L, 1L))
  expect_identical(dim(blk[["mask"]]), c(6L, 1L))
  # Every subject's FIRST trial starts from Q1 = Q2 = 0, so its choice
  # probability is exactly one half whatever alpha and beta are. That
  # makes the log-likelihood 6 * log(0.5) for ANY parameter value, and
  # a transposed block cannot produce it.
  ll <- rw_loglik(fr[["y"]][["choice"]], list(alpha = 0.3, beta = 2),
                  fr[["aterm_values"]][["choice"]], 1, blk, NULL)
  expect_equal(as.numeric(ll), 6 * log(0.5), tolerance = 1e-12)
})

test_that("the taped recursion equals an independent scalar reference", {
  skip_unless_rl()
  for (drop in c(FALSE, TRUE)) {
    fx <- rl_fixture(drop = drop)
    ref <- rl_reference(fx$fit, fx$data)
    joint <- -fx$fit$obj$env$f(fx$fit$obj$env$last.par.best)
    expect_equal(joint, ref$joint, tolerance = 1e-8)
  }
})

test_that("fitted() is the per-trial choice probability", {
  skip_unless_rl()
  fx <- rl_fixture(drop = TRUE)
  p <- fitted(fx$fit)
  expect_length(p, nrow(fx$data))
  expect_true(all(p > 0 & p < 1))
  # the data log-likelihood the reference computes is the one these
  # probabilities imply, so fitted() and the tape agree row by row
  ref <- rl_reference(fx$fit, fx$data)
  expect_equal(sum(stats::dbinom(fx$data$choice, 1L, p, log = TRUE)),
               ref$data, tolerance = 1e-8)
  expect_equal(unname(predict(fx$fit, type = "response")), unname(p),
               tolerance = 1e-10)
  expect_equal(unname(residuals(fx$fit, type = "response")),
               fx$data$choice - unname(p), tolerance = 1e-10)
  # pearson divides by the binomial variance of that probability
  expect_equal(unname(residuals(fx$fit, type = "pearson")),
               (fx$data$choice - unname(p)) / sqrt(p * (1 - p)),
               tolerance = 1e-8)
})

test_that("the structured simulator draws whole sequences", {
  skip_unless_rl()
  fx <- rl_fixture()
  sm <- simulate(fx$fit, nsim = 20, seed = 7)
  expect_identical(dim(sm), c(nrow(fx$data), 20L))
  expect_true(all(as.matrix(sm) %in% c(0, 1)))
  # the fit reproduces its own choice rate: 20 draws of 1200 trials put
  # the observed rate well inside the simulated spread
  rate <- colMeans(as.matrix(sm))
  expect_lt(abs(mean(rate) - mean(fx$data$choice)),
            4 * stats::sd(rate) + 0.02)
  # frm_simulate() reaches the same slot with no fit in hand
  sim2 <- rl_simulate(rl_bandit_design(6, 30), nsim = 2, seed = 3)
  expect_length(sim2, 2L)
  expect_true(all(sim2[[1L]]$choice %in% c(0, 1)))
  expect_false(identical(sim2[[1L]]$choice, sim2[[2L]]$choice))
})

test_that("a rw_delta() fit refuses what its structure declares", {
  skip_unless_rl()
  fx <- rl_fixture()
  expect_error(conditional_effects(fx$fit), "whole trial history")
  expect_error(residuals(fx$fit, type = "osa"), "registered observation")
  expect_error(residuals(fx$fit, type = "deviance"), "returns one total")
  expect_error(predict(fx$fit, newdata = fx$data, type = "response"),
               "carries no block to replay")
  expect_error(
    frm(rl_bform, family = rl_family(), data = fx$data, REML = TRUE),
    "REML")
  # the addition terms the recursion cannot honor, refused in the
  # family's own words rather than through a missing CDF
  wform <- bf(choice | reward(pay1, pay2) + weights(w) ~ condition,
              beta ~ 1)
  dw <- fx$data
  dw$w <- 1
  expect_error(frm(wform, family = rl_family(), data = dw),
               "conditional on every earlier trial")
  # loo() refuses one step earlier than the structure would: an elpd is
  # a posterior quantity and this is a maximum-likelihood fit. The
  # structure's `unit` string is read by frmtmb.sample's loo(), which is
  # where a sampled rw_delta() fit would meet the pointwise gap
  # dev/rl-findings.md records.
  expect_error(loo(fx$fit), "maximum-likelihood fit")
})

test_that("the response and its data are checked before the tape", {
  skip_unless_rl()
  fx <- rl_fixture()
  d2 <- fx$data
  d2$choice[1L] <- 2
  expect_error(frm(rl_bform, family = rl_family(), data = d2),
               "the arm chosen on each trial")
  d3 <- fx$data
  d3$trial[d3$id == levels(d3$id)[1L]][2L] <- 1
  expect_error(frm(rl_bform, family = rl_family(), data = d3),
               "unique within a subject")
  # reward() is not optional: the recursion has nothing to learn from
  expect_error(frm(bf(choice ~ condition, beta ~ 1),
                   family = rl_family(), data = fx$data),
               "reward")
})

test_that("the fixed effects recover their simulated values", {
  skip_unless_rl()
  skip_on_cran()
  set.seed(19)
  dd <- rl_simulate(rl_bandit_design(60, 150), nsim = 1, seed = 19)[[1L]]
  fit <- suppressWarnings(frm(rl_bform, family = rl_family(), data = dd))
  ci <- confint(fit)
  truth <- c(`alpha_(Intercept)` = rl_truth$alpha_Intercept,
             alpha_conditiontrt = rl_truth$alpha_conditiontrt,
             `beta_(Intercept)` = rl_truth$beta_Intercept)
  for (nm in names(truth)) {
    expect_gte(truth[[nm]], ci[nm, "lwr"])
    expect_lte(truth[[nm]], ci[nm, "upr"])
  }
})

# --- the Stan tier ----------------------------------------------------

test_that("the taped joint density equals a Stan program's log_prob", {
  skip_unless_rl_stan()
  fx <- rl_fixture(n_subj = 30L, n_trial = 100L, seed = 5L)
  # the reshape the map relies on: `b` is level-major, so a subject's
  # two coefficients are contiguous
  pars <- rl_stan_pars(fx$fit$obj$env$last.par.best)
  # ranef() carries a term attribute since 0.52.0, so compare the values
  expect_equal(as.vector(pars$u), as.vector(ranef(fx$fit)[[1L]]),
               tolerance = 1e-10)
  out <- rl_lp_check(fx$fit, fx$data)
  expect_lt(abs(out$const), 1e-6)
})

test_that("the Stan identity holds away from the optimum too", {
  skip_unless_rl_stan()
  fx <- rl_fixture(n_subj = 30L, n_trial = 100L, seed = 5L)
  par <- fx$fit$obj$env$last.par.best
  # a point that is nobody's optimum, so the agreement cannot be an
  # accident of both sides being stationary. The covariance parameters
  # are left alone: they ride into Stan as data.
  move <- names(par) %in% c("beta", "betad", "b")
  par[move] <- par[move] + 0.25 * cos(seq_len(sum(move)))
  rl_lp_check(fx$fit, fx$data, par = par, check_grad = FALSE)
})

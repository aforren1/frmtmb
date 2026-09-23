# Lane wt-reunc: coverage of predict() and fitted() at a grouping level
# the fit SAW.
#
# THE CLAIM. predict()'s 95 percent interval at a known level is a
# prediction interval for a NEW observation at that level, and its
# coverage is AVERAGED over data sets and over the group effects: in
# every replicate the truth draws fresh group effects u ~ N(0, tau^2),
# the fit is made, and new points land at levels chosen at random. It
# does NOT claim coverage conditional on one group's realized effect.
# The best linear unbiased predictor shrinks toward zero, so for a
# group whose effect is far from zero the interval sits too close to
# the population and under-covers, and for a group near zero it
# over-covers; only the average over groups is nominal. That is the
# frequentist reading of what brms's posterior interval at a known
# level covers, and it is the one Henderson's prediction error variance
# delivers exactly when the variance parameters are known. The side
# table "by |u| tercile" measures the conditional behavior so the two
# claims are not confused.
#
# fitted()'s interval claims the same average coverage for the
# CONDITIONAL MEAN at that level, x'beta + u_g, not for an observation.
#
# THE ARMS.
#   pred_ml, pred_reml, pred_prof   predict() under ML, REML = TRUE and
#                                   control(profile = TRUE)
#   pred_plug                       predict(propagate_error = FALSE), ML
#   oracle_pred                     POSITIVE CONTROL. Henderson's mixed
#                                   model equations at the TRUE sigma and
#                                   tau: yhat +- z sqrt(sigma^2 + a'C^-1 a).
#                                   The prediction error is exactly normal
#                                   with that variance, marginally over u,
#                                   so this covers 0.95 by construction.
#   fit_ml, fit_reml                fitted() covering x'beta + u_g
#   oracle_mean                     POSITIVE CONTROL for the mean:
#                                   yhat +- z sqrt(a'C^-1 a)
#   wald_obs                        NEGATIVE CONTROL: fitted()'s interval
#                                   asked to cover the new OBSERVATION. It
#                                   carries no observation noise, so it
#                                   must under-cover badly.
#
# Monte Carlo error: the m points of one replicate share a fit, so the
# standard error is the replicate-level one, sd(per-replicate coverage)
# / sqrt(R), and z is against 0.95 with it.
#
#   Rscript dev/reunc-coverage.R <lib> <design> <nrep> <outdir>
# designs: mixed12 (the dev/shapes-coverage.R mixed design: n = 60,
# 12 groups, sigma = 1, tau = 0.7) and few4 (4 groups of 10).
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]
design <- a[2]
nrep <- as.integer(a[3])
outdir <- a[4]
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages(library(frmtmb))
lane <- grepl("reunc-lib", find.package("frmtmb"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ng <- switch(design, mixed12 = 12L, few4 = 4L)
n <- switch(design, mixed12 = 60L, few4 = 40L)
m <- 10L
ndraws <- 1000L
seed0 <- switch(design, mixed12 = 20260922L, few4 = 20260923L)
TAU <- 0.7
SIGMA <- 1
BETA <- c(1, 0.5)
z <- stats::qnorm(0.975)

henderson <- function(d, nd) {
  X <- cbind(1, d$x)
  Z <- stats::model.matrix(~ 0 + g, d)
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + diag(SIGMA^2 / TAU^2,
                                                       ncol(Z)))) / SIGMA^2
  rhs <- c(crossprod(X, d$y), crossprod(Z, d$y)) / SIGMA^2
  Ci <- solve(C)
  est <- drop(Ci %*% rhs)
  An <- cbind(1, nd$x, stats::model.matrix(~ 0 + g, nd))
  list(fit = drop(An %*% est), pev = rowSums((An %*% Ci) * An))
}

arms <- c("pred_ml", "pred_reml", "pred_prof", "pred_plug", "oracle_pred",
          "fit_ml", "fit_reml", "oracle_mean", "wald_obs")
rows <- list()
t_start <- proc.time()[["elapsed"]]
for (r in seq_len(nrep)) {
  seed <- seed0 + 1000L * r
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng),
                                               length.out = n)))
  u <- rnorm(ng, 0, TAU)
  d$y <- rnorm(n, BETA[1] + BETA[2] * d$x + u[d$g], SIGMA)
  gsel <- sample(seq_len(ng), m, TRUE)
  nd <- data.frame(x = rnorm(m), g = factor(gsel, levels = levels(d$g)))
  mu_new <- BETA[1] + BETA[2] * nd$x + u[gsel]
  y_new <- rnorm(m, mu_new, SIGMA)
  fits <- list(
    ml = tryCatch(suppressWarnings(frm(bf(y ~ x + (1 | g)), data = d)),
                  error = function(e) NULL),
    reml = tryCatch(suppressWarnings(frm(bf(y ~ x + (1 | g)), data = d,
                                         REML = TRUE)),
                    error = function(e) NULL),
    prof = tryCatch(suppressWarnings(
      frm(bf(y ~ x + (1 | g)), data = d,
          control = frmtmb_control(profile = TRUE))),
      error = function(e) NULL))
  if (any(vapply(fits, is.null, NA))) {
    rows[[length(rows) + 1L]] <- data.frame(rep = r, seed = seed,
                                            failed = TRUE)
    next
  }
  pr <- function(fit, s, ...) {
    set.seed(s)
    tryCatch(suppressWarnings(predict(fit, newdata = nd, ndraws = ndraws,
                                      ...)),
             error = function(e) NULL)
  }
  ivs <- list(
    pred_ml = pr(fits$ml, seed + 1L),
    pred_reml = pr(fits$reml, seed + 2L),
    pred_prof = pr(fits$prof, seed + 3L),
    # WARNING: this arm no longer reproduces the pred_plug rows in
    # dev/reunc-findings.md section 3. Those were measured when the
    # argument was `param_uncertainty`, which held the PARAMETERS and
    # still drew the group effects. Since the 2026-09-23 rename
    # `propagate_error = FALSE` holds the group effects too, so this
    # arm is now an interval carrying the observation noise alone and
    # under-covers by construction. The recorded rows are labelled OLD
    # flag and were deliberately not rerun; re-running this script
    # produces a DIFFERENT quantity, not a check on them.
    pred_plug = pr(fits$ml, seed + 4L, propagate_error = FALSE),
    fit_ml = suppressWarnings(fitted(fits$ml, newdata = nd)),
    fit_reml = suppressWarnings(fitted(fits$reml, newdata = nd)))
  hd <- henderson(d, nd)
  ivs$oracle_pred <- cbind(hd$fit, NA, hd$fit - z * sqrt(SIGMA^2 + hd$pev),
                           hd$fit + z * sqrt(SIGMA^2 + hd$pev))
  ivs$oracle_mean <- cbind(hd$fit, NA, hd$fit - z * sqrt(hd$pev),
                           hd$fit + z * sqrt(hd$pev))
  ivs$wald_obs <- ivs$fit_ml
  target <- list(pred_ml = y_new, pred_reml = y_new, pred_prof = y_new,
                 pred_plug = y_new, oracle_pred = y_new, fit_ml = mu_new,
                 fit_reml = mu_new, oracle_mean = mu_new, wald_obs = y_new)
  if (any(vapply(ivs, is.null, NA))) {
    rows[[length(rows) + 1L]] <- data.frame(rep = r, seed = seed,
                                            failed = TRUE)
    next
  }
  for (arm in arms) {
    iv <- ivs[[arm]]
    rows[[length(rows) + 1L]] <- data.frame(
      rep = r, seed = seed, failed = FALSE, arm = arm, point = seq_len(m),
      absu = abs(u[gsel]), n_g = as.integer(table(d$g)[gsel]),
      hit = target[[arm]] >= iv[, 3L] & target[[arm]] <= iv[, 4L],
      width = iv[, 4L] - iv[, 3L])
  }
  if (r %% 20L == 0L) {
    cat(sprintf("rep %d of %d, %.0f s\n", r, nrep,
                proc.time()[["elapsed"]] - t_start))
  }
}
res <- do.call(rbind, lapply(rows, function(x) {
  for (v in c("arm", "point", "absu", "n_g", "hit", "width")) {
    if (is.null(x[[v]])) x[[v]] <- NA
  }
  x
}))
tag <- if (lane) "lane" else "base"
saveRDS(list(design = design, lib = lib, lane = lane, seed0 = seed0,
             nrep = nrep, m = m, ndraws = ndraws, res = res,
             version = format(packageVersion("frmtmb"))),
        file.path(outdir, sprintf("cov-%s-%s.rds", design, tag)))
cat("done", design, tag, "in", proc.time()[["elapsed"]] - t_start, "s\n")

# Reviewer, lane splinecurve: what the three curve functions return at an
# UNSEEN level of a BAR term, a route allow_new_levels = TRUE opens and
# the lane guards only for a contrast.
#   Rscript splinecurve-rev-02-extravar.R > splinecurve-rev-log/02-extravar.txt
#
# frm_lp_basis() returns a new level's variance as a per-ROW vector,
# extra_var, with no covariance between rows. One unseen subject is ONE
# draw u (or (u0, u1) with a slope) shared by every row of the grid, so
# the grid's covariance is C V C' + Z S Z', where Z is the new level's
# design rows and S the block covariance. Each function is compared with
# that.
# Data: the curve-inference vignette's simulation, set.seed(4).
LIB <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb.spline", format(packageVersion("frmtmb.spline")), "from",
    dirname(find.package("frmtmb.spline")), "\n\n")

set.seed(4)
n_sub <- 20; n_rep <- 12; n_t <- 30
peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
sub <- rep(seq_len(n_sub), each = n_rep * n_t)
d <- data.frame(subject = factor(sub),
                t = rep(seq(0, 1, length.out = n_t), times = n_sub * n_rep))
h_sub <- rnorm(n_sub, 1, 0.12)
s_sub <- rnorm(n_sub, 0.5, 0.04)
d$v <- peak(d$t, h_sub[sub], s_sub[sub]) + rnorm(nrow(d), 0, 0.06)
lev_u <- c(levels(d$subject), "new")
grid <- data.frame(t = seq(0.05, 0.95, length.out = 60),
                   subject = factor("new", levels = lev_u))

# simultaneous coverage of a band est +/- crit * se when the deviation
# process is N(0, S): the fraction of draws whose max |dev / se| <= crit
cover <- function(S, se, crit, nsim = 20000L, seed = 99L) {
  set.seed(seed)
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow(S))
  z <- matrix(rnorm(nrow(S) * nsim), nrow(S), nsim)
  mx <- apply(abs((L %*% z) / se), 2L, max)
  c(coverage = mean(mx <= crit),
    crit_true = unname(quantile(mx, 0.95, type = 8)))
}

run <- function(label, form, Zfun) {
  cat("====", label, "====\n")
  fit <- frm(bf(form), family = gaussian(), data = d)
  vc <- VarCorr(fit)$subject
  sds <- vc$sd[, 1]
  # VarCorr() is brms-shaped: cov[row, "Estimate", col]
  S <- if (!is.null(vc$cov)) {
    as.matrix(vc$cov[, "Estimate", ])
  } else {
    matrix(sds^2, 1, 1)
  }
  cat("block sd:", format(sds, digits = 6), "\n")
  lb <- frm_lp_basis(fit, newdata = grid, re_formula = NULL,
                     allow_new_levels = TRUE)
  C <- as.matrix(lb$A)
  Z <- Zfun(grid$t)
  zsz <- Z %*% S %*% t(Z)
  cat("extra_var vs diag(Z S Z'): max rel diff",
      format(max(abs(lb$extra_var / diag(zsz) - 1)), digits = 3), "\n")
  cv <- frm_curve(fit, newdata = grid, re_formula = NULL,
                  allow_new_levels = TRUE, nsim = 20000, seed = 1)
  full <- C %*% lb$V %*% t(C) + zsz
  cat("frm_curve pointwise .se vs sqrt(diag(full)): max rel diff",
      format(max(abs(cv$.se / sqrt(diag(full)) - 1)), digits = 3), "\n")
  Sig <- attr(cv, "Sigma")
  cat("attr Sigma: max |diag(Sigma) / .se^2 - 1|",
      format(max(abs(diag(Sig) / cv$.se^2 - 1)), digits = 4),
      " (share of .se^2 that is extra_var: median",
      format(median(lb$extra_var / cv$.se^2), digits = 4), ")\n")
  cr <- cover(full, cv$.se, cv$.crit_sim[1])
  cat("frm_curve simultaneous: crit reported", format(cv$.crit_sim[1],
                                                      digits = 5),
      " crit from the full covariance", format(cr[["crit_true"]], digits = 5),
      " qnorm(0.975) 1.96\n",
      "  coverage of the reported band under the full covariance",
      format(cr[["coverage"]], digits = 4), "(nominal 0.95, 20000 draws)\n")

  d1 <- frm_curve_deriv(fit, var = "t", newdata = grid, re_formula = NULL,
                        allow_new_levels = TRUE, simultaneous = FALSE)
  e <- attr(d1, "eps")
  Zd <- (Zfun(grid$t + e) - Zfun(grid$t - e)) / (2 * e)
  lbp <- frm_lp_basis(fit, newdata = transform(grid, t = t + e),
                      re_formula = NULL, allow_new_levels = TRUE)
  lbm <- frm_lp_basis(fit, newdata = transform(grid, t = t - e),
                      re_formula = NULL, allow_new_levels = TRUE)
  D <- (as.matrix(lbp$A) - as.matrix(lbm$A)) / (2 * e)
  se_true <- sqrt(diag(D %*% lbp$V %*% t(D)) + diag(Zd %*% S %*% t(Zd)))
  cat("frm_curve_deriv: .se / se with the new level's slope variance:",
      "range", format(range(d1$.se / se_true), digits = 4), "\n")

  pk <- frm_curve_feature(fit, var = "t", type = "crossing", at = 0.5,
                          newdata = grid, re_formula = NULL,
                          allow_new_levels = TRUE)
  at_root <- frm_curve(fit, newdata = transform(grid[seq_len(nrow(pk)), ],
                                                t = pk$.estimate),
                       re_formula = NULL, allow_new_levels = TRUE,
                       simultaneous = FALSE)
  zr <- Zfun(pk$.estimate)
  cat("frm_curve_feature(crossing at 0.5): roots",
      format(pk$.estimate, digits = 5), "\n",
      "  .value_se", format(pk$.value_se, digits = 5),
      " frm_curve .se at the same t", format(at_root$.se, digits = 5), "\n",
      "  crossing .se", format(pk$.se, digits = 5),
      " with the new level's variance",
      format(pk$.se * at_root$.se / pk$.value_se, digits = 5), "\n\n")
  invisible(NULL)
}

run("(1 | subject)", v ~ s(t, k = 12) + (1 | subject),
    function(t) matrix(1, length(t), 1))
run("(1 + t | subject)", v ~ s(t, k = 12) + (1 + t | subject),
    function(t) cbind(1, t))

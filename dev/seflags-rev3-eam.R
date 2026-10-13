# Reviewer of lane seflags, final check of m5's consequence:
# frmtmb.eam's test-defects.R:92 (wiener + lognormal mixture, weight at
# 6e-11). Which parameters are lost and why; for residuals(type =
# "response"), the gradient G of every row along each lost direction
# (natural units, so G . n is the change in the row's fitted mean per
# unit step along n), the change in the fitted mean when the lost
# parameters move by +-1 and +-5 on their own scale (the mean refit by
# eval of the mixture, not linearized), against the rows' kept SE.
#   Rscript dev/seflags-rev3-eam.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
set.seed(3)
dat <- frmtmb.eam:::ddm_simulate(300, mu = 0.9, bs = 1.4, ndt = 0.30)
k <- sample(300, 18)
dat$rt[k] <- stats::runif(18, 0.12, 2.5)
w <- character()
fit <- withCallingHandlers(frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
  family = mixture(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                   lognormal()), data = dat),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
cat("fit warnings:\n"); for (x in w) cat(" *", substr(x, 1, 400), "\n")
sdr <- suppressWarnings(ns$sdr_of(fit))
on <- ns$outer_par_names(fit)
cat("\nouter parameters and lost reasons:\n")
print(data.frame(par = on, est = signif(fit$opt$par, 4),
                 lost = ifelse(on %in% names(sdr$se_lost),
                               sdr$se_lost[on], "")), row.names = FALSE)
cat("se_lost:\n"); print(sdr$se_lost)
cap_g <- NULL
trace("jc_nonest", where = ns, print = FALSE, tracer = quote({
  cap_g <<- list(G = as.matrix(G), pos = pos, N = jc$null, flat = jc$flat,
                 units = jc$units)
}))
rr <- suppressWarnings(stats::residuals(fit, type = "response"))
untrace("jc_nonest", where = ns)
cat("\nresidual Est.Error: NaN", sum(is.nan(rr[, "Est.Error"])), "of",
    nrow(rr), "\n")
if (!is.null(cap_g)) {
  G <- cap_g$G; N <- cap_g$N; pos <- cap_g$pos
  cat("null columns", NCOL(N), "| G", paste(dim(G), collapse = " x "),
      "| pos", length(pos), "\n")
  gn <- abs(G %*% N[pos, , drop = FALSE])
  cat("max |G . n| per null column (natural units, mean change per unit step):",
      signif(apply(gn, 2, max), 3), "\n")
  cat("which outer parameters load on each null column (|n| > 0.05):\n")
  for (j in seq_len(NCOL(N))) {
    l <- which(abs(N[, j]) > 0.05)
    cat("  col", j, ":", paste(sprintf("%d(%.2f)", l, N[l, j]), collapse = " "), "\n")
  }
}
# the rows' SEs with the prediction test off: what the kept directions give
op <- ns$se_pred_tol
utils::assignInNamespace("jc_nonest", function(jc, G, pos) rep(FALSE, NROW(G)),
                         ns = "frmtmb")
r0 <- suppressWarnings(stats::residuals(fit, type = "response"))
cat("kept-directions SE of the residuals: range", signif(range(r0[, "Est.Error"]), 3), "\n")
# move the lost parameters on their own scale and recompute the mean
lost <- which(on %in% names(sdr$se_lost))
if (length(lost)) {
  mean_at <- function(par) {
    f2 <- fit
    f2$opt$par <- par
    f2$estimates <- fit$obj$env$parList(par)
    stats::fitted(f2, summary = TRUE)[, "Estimate"]
  }
  ok <- TRUE
  cat("estimates_at available:", ok, "\n")
  if (ok) {
    m0 <- mean_at(fit$opt$par)
    for (st in c(-5, -1, 1, 5)) for (j in lost) {
      p <- fit$opt$par; p[j] <- p[j] + st
      m1 <- tryCatch(mean_at(p), error = function(e) rep(NA, length(m0)))
      cat(sprintf("  %-25s %+d: max |change in fitted mean| %.3g\n", on[j], st,
                  max(abs(m1 - m0))))
    }
  }
}
# the mixing weight: a slice of the objective along theta1 (others at
# the estimates) and the change in the fitted mean there
j <- which(grepl("^theta1", on))
if (length(j) == 1L) {
  f0 <- fit$obj$fn(fit$opt$par)
  m0 <- mean_at(fit$opt$par)
  for (t in c(-23.44, -15, -10, -8, -7, -6, -5, -4.5, -4, -3.5, -3)) {
    p <- fit$opt$par; p[j] <- t
    cat(sprintf("  theta1 %6.2f (weight %.3g): objective rise %.4g | max |change in mean| %.3g | kept SE min %.3g\n",
                t, stats::plogis(t), fit$obj$fn(p) - f0,
                max(abs(mean_at(p) - m0)), min(r0[, "Est.Error"])))
  }
}

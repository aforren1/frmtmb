# rel084-kt copy of dev/sparsese-fix1.R: rellib-r22 removed from .libPaths() (lib = argument 1, then user lib). No other change.
# Lane sparsese, fix round 1 (review 2026-11-17, F1 and F2). Two designs
# of dev/reviews/sparsese-review-designs.R, each fitted twice in one
# process: the dense route (limit 1e6) and the sparse one (limit 20), the
# option in force at fit and read time.
#  F1: frm_curve() on the gaussian "base" design, and frm_lp_basis()'s
#      V_root recipe against frm_linpred(se.fit = TRUE).
#  F2: fitted() Est.Error of a cumulative fit (fit_fd_se()), new rows and
#      in-sample.
#   Rscript sparsese-fix1.R <lib: lane|<path>> [seed] [group multiplier]
a <- commandArgs(trailingOnly = TRUE)
lib <- if (a[1] == "lane") "C:/Users/adf44/source/r/wt-sparsese-lib" else a[1]
seed <- if (length(a) >= 2L) as.integer(a[2]) else 1L
# a multiplier of the group effect; 0 runs the group sd to its boundary
gm <- if (length(a) >= 3L) as.numeric(a[3]) else 1
.libPaths(c(lib,
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "seed", seed, "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
lims <- c(dense = 1e6, sparse = 20)
on <- function(arm, expr) {
  op <- options(frmtmb.joint_root_max = lims[[arm]])
  on.exit(options(op))
  expr
}
mk <- function(seed) {
  set.seed(seed)
  ng <- 30L
  n <- 300L
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(seq_len(ng), length.out = n)),
                  f = factor(rep(c("a", "b"), length.out = n)))
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  d$u <- gm * stats::rnorm(ng, 0, 0.5)[d$g]
  d$e <- stats::rnorm(n, 0, 0.3)
  d$z <- stats::rlogis(n)
  d
}
d <- mk(seed)
d$y <- sin(3 * d$x) + d$u + 0.3 * (d$f == "b") + d$e
eta <- sin(3 * d$x) + d$u + 0.3 * (d$f == "b")
d$yc <- factor(cut(0.8 * eta + d$z, c(-Inf, 0, 0.8, Inf),
                   labels = c("lo", "mid", "hi")), ordered = TRUE)
grid <- data.frame(x = seq(1, 2, by = 0.1), f = factor("a", c("a", "b")))
for (k in 1:7) grid[[paste0("x", k)]] <- grid$x^k
cat("== F1, gaussian REML, y ~ f + x1..x7 + (1 | g)\n")
res <- list()
for (arm in names(lims)) on(arm, {
  fit <- q(frm(y ~ f + x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g), data = d,
               REML = TRUE))
  invisible(q(frmtmb:::sdr_of(fit)))
  cv <- tryCatch(q(frm_curve(fit, newdata = grid, re_formula = NA,
                             simultaneous = FALSE)),
                 error = function(e) e)
  se <- q(frm_linpred(fit, newdata = grid, re_formula = NA,
                      se.fit = TRUE))$se.fit
  lb <- q(frm_lp_basis(fit, newdata = grid, re_formula = NA))
  rec <- if (is.null(lb$V_root)) NA else {
    sqrt(rowSums((lb$A %*% lb$V_root)^2) + rowSums((lb$A %*% lb$V_rest) *
                                                     lb$A))
  }
  cat(arm, ": joint", nrow(frmtmb:::sdr_of(fit)$jointPrecision),
      " route", class(frmtmb:::get_joint_cov(fit)$V)[1],
      " V_root", !is.null(lb$V_root), "\n")
  if (inherits(cv, "error")) {
    cat("  frm_curve ERROR:", substr(conditionMessage(cv), 1, 200), "\n")
  } else {
    cat("  frm_curve cov_rel_error",
        format(attr(cv, "check")$cov_rel_error, digits = 3),
        "; max |.se / se.fit - 1|", format(max(abs(cv$.se / se - 1)),
                                           digits = 3), "\n")
  }
  cat("  V_root recipe: max |recipe / se.fit - 1|",
      format(max(abs(rec / se - 1)), digits = 3), "\n")
  res[[arm]] <- list(se = se, cv = if (!inherits(cv, "error")) cv$.se)
})
cat("se.fit sparse / dense - 1:",
    format(max(abs(res$sparse$se / res$dense$se - 1)), digits = 3), "\n")
if (!is.null(res$sparse$cv)) {
  cat("frm_curve .se sparse / dense - 1:",
      format(max(abs(res$sparse$cv / res$dense$cv - 1)), digits = 3), "\n")
}
cat("== F2, cumulative REML, yc ~ x1..x5 + (1 | g)\n")
nd <- d[c(1:6, 151:156), ]
out <- list()
for (arm in names(lims)) on(arm, {
  fit <- q(frm(yc ~ x1 + x2 + x3 + x4 + x5 + (1 | g), data = d,
               family = cumulative(), REML = TRUE))
  invisible(q(frmtmb:::sdr_of(fit)))
  cat(arm, ": joint", nrow(frmtmb:::sdr_of(fit)$jointPrecision),
      " whitened", !is.null(fit$beta_whiten), " lost [",
      paste(frmtmb:::sdr_of(fit)$se_lost), "]\n")
  t0 <- proc.time()[["elapsed"]]
  fn <- q(fitted(fit, newdata = nd))
  fi <- q(fitted(fit))
  cat(sprintf("  fitted() new + in-sample: %.2f s\n",
              proc.time()[["elapsed"]] - t0))
  out[[arm]] <- list(n = fn[, "Est.Error", ], i = fi[, "Est.Error", ],
                     e = fi[, "Estimate", ])
})
cat("Estimate identical:", identical(out$dense$e, out$sparse$e), "\n")
for (w in c("n", "i")) {
  x <- out$dense[[w]]
  y <- out$sparse[[w]]
  cat(sprintf("%s Est.Error: max |sparse / dense - 1| %.3e; max sparse %.4f",
              if (w == "n") "new rows" else "in-sample",
              max(abs(y / x - 1)), max(y)),
      sprintf("; max dense %.4f\n", max(x)))
}

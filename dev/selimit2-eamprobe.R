# Lane selimit2: why frmtmb.eam's test-defects.R "max_ndt above min(rt)"
# mixture gets NaN residual SEs from jc_limit_extrap() (item 5433's
# non-local move test). Prints the limit info, the level ranges
# se_limit_ranges() finds, and per direction the largest move of the
# mixture mean over the grid against the row SEs.
#   Rscript dev/selimit2-eamprobe.R
.libPaths(c("C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
ns <- asNamespace("frmtmb")
set.seed(3)
dat <- ddm_simulate(300, mu = 0.9, bs = 1.4, ndt = 0.30)
k <- sample(300, 18)
dat$rt[k] <- stats::runif(18, 0.12, 2.5)
fit <- suppressWarnings(frm(bf(rt | dec(upper) ~ 1, bias1 = 0.5),
  family = mixture(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                   lognormal()), data = dat))
li <- ns$se_limit_info(fit)
print(li[c("kind", "est", "dir", "gain", "joint", "sep")])
sdr <- ns$sdr_of(fit)
print(sdr$se_lost)
print(fit$opt$par)
u <- ns$se_limit_units(fit, ns$outer_par_map(fit))
print(u)
lr <- ns$se_limit_ranges_of(fit, sdr)
for (e in lr) {
  cat("pos", e$pos, "w", signif(e$w, 4), "t range", signif(range(e$t), 4),
      "\n")
}
r <- suppressWarnings(stats::residuals(fit, type = "response"))
cat("NaN SE rows:", sum(is.nan(r[, "Est.Error"])), "of", nrow(r), "\n")

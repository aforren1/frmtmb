# Reviewer of lane surface, claim 1: pp_mixture() on a fit against
# brms 2.23.0's pp_mixture() on the same data (dev/surface-rev-brms.R:
# 4 chains x 1500 post-warmup draws, seed 1, brms's default priors).
#   Rscript dev/surface-rev-ppbrms.R > dev/surface-rev-out/ppbrms.txt
b <- readRDS("dev/surface-rev-out/brms-ppmix.rds")
f <- readRDS("dev/surface-rev-out/ppmix.rds")
pb <- b$pm; pf <- f$pf
cat("brms dim", dim(pb), "| frmtmb dim", dim(pf), "\n")
cat("dimnames identical:", identical(dimnames(pb), dimnames(pf)), "\n")
str(dimnames(pb)); str(dimnames(pf))
r <- cor(pb[, "Estimate", 1], pf[, "Estimate", 1])
cat("component correlation (1 vs 1):", r, "\n")
k <- if (r < 0) c(2, 1) else c(1, 2)
d <- (pf[, "Estimate", 1] - pb[, "Estimate", k[1]]) / pb[, "Est.Error", k[1]]
cat("|frmtmb Estimate - brms mean| / brms Est.Error: median",
    signif(median(abs(d)), 3), "| 95th pct", signif(quantile(abs(d), 0.95), 3),
    "| max", signif(max(abs(d)), 3), "\n")
inside <- pf[, "Estimate", 1] >= pb[, "Q2.5", k[1]] &
  pf[, "Estimate", 1] <= pb[, "Q97.5", k[1]]
cat("frmtmb Estimate inside brms 95% interval:", sum(inside), "of",
    length(inside), "\n")
p <- pf[, "Estimate", 1]
bins <- cut(pmin(p, 1 - p), c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.5),
            right = FALSE)
rat <- pf[, "Est.Error", 1] / pb[, "Est.Error", k[1]]
ln <- f$ln_sd / pb[, "Est.Error", k[1]]
cat("\nEst.Error / brms posterior SD by min(p, 1 - p):\n")
print(data.frame(n = tapply(rat, bins, length),
                 delta_median = signif(tapply(rat, bins, median), 3),
                 logitnormal_median = signif(tapply(ln, bins, median), 3)))
cat("all rows: delta median", signif(median(rat), 3),
    "| logit-normal median", signif(median(ln), 3), "\n")
wf <- pf[, "Q97.5", 1] - pf[, "Q2.5", 1]
wb <- pb[, "Q97.5", k[1]] - pb[, "Q2.5", k[1]]
cat("interval width frmtmb / brms: median", signif(median(wf / wb), 3),
    "| rows with brms width > 0.02:",
    signif(median((wf / wb)[wb > 0.02]), 3), "\n")
cat("brms posterior summary:\n"); print(round(b$summ[1:5, ], 3))

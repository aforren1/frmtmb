d <- read.csv("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/calib.csv")
d$gap <- d$ll_sn - d$ll
d$gain <- d$ll_refit - d$ll
d$best <- pmax(d$ll_sn, d$ll_refit)
d$short <- d$best - d$ll
cat("rows", nrow(d), "\n\n")
cat("arm            n  short>1e-6  short>0.01  short>1  maxshort  ",
    "refit>sn\n", sep = "")
for (a in unique(d$arm)) {
  s <- d[d$arm == a, ]
  cat(sprintf("%-13s %2d %10d %11d %8d %9.4f %8d\n", a, nrow(s),
              sum(s$short > 1e-6), sum(s$short > 0.01),
              sum(s$short > 1), max(s$short),
              sum(s$ll_refit > s$ll_sn + 1e-6)))
}
cat("\n-- |alpha| against shortfall, all 240 rows\n")
br <- c(-1, 1e-4, 1e-3, 1e-2, 0.05, 0.2, 1, Inf)
g <- cut(abs(d$alpha), br)
print(table(bin = g, short_gt_0.01 = d$short > 0.01))
cat("\n-- shortfall by |alpha| bin\n")
print(tapply(d$short, g, function(v) c(n = length(v), max = max(v))))
cat("\n-- rows with short > 0.01 and |alpha| > 0.05\n")
sel <- d$short > 0.01 & abs(d$alpha) > 0.05
print(d[sel, c("arm", "seed", "alpha", "alpha_sn", "ll", "ll_sn",
               "ll_refit", "short", "se_alpha")])
cat("\n-- rows with |alpha| < 1e-2 : how often is a refit worth it\n")
z <- d[abs(d$alpha) < 1e-2, ]
cat("n =", nrow(z), " short>1e-6:", sum(z$short > 1e-6),
    " short>0.01:", sum(z$short > 0.01), " short>1:", sum(z$short > 1), "\n")
print(table(arm = z$arm, worth = z$short > 0.01))
cat("\n-- se_alpha on the |alpha| < 1e-2 rows\n")
print(summary(z$se_alpha))
cat("NA se:", sum(is.na(z$se_alpha)), "\n")
cat("\n-- probe on the |alpha| < 1e-2 rows, split by worth\n")
for (p in c("probe_0.25", "probe_0.5", "probe_1", "probe_2")) {
  cat(sprintf("%-11s worth: %s\n", p,
              paste(sprintf("%.4f", quantile(z[[p]][z$short > 0.01],
                                             c(0, .05, .5, 1))),
                    collapse = " ")))
  cat(sprintf("%-11s not  : %s\n", "",
              paste(sprintf("%.4f", quantile(z[[p]][z$short <= 0.01],
                                             c(0, .5, .95, 1))),
                    collapse = " ")))
}

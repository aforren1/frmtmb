# M3: the convergence-code and condition class rate over both 40-seed
# streams, base against the lane. logLik alone was not enough: a fit
# that reaches the optimum and warns is still a defect.
d <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/"
a <- read.csv(paste0(d, "accept-base.csv"))
b <- read.csv(paste0(d, "accept-fixed.csv"))
stopifnot(nrow(a) == nrow(b), identical(a$seed, b$seed),
          identical(a$stream, b$stream))
cat("fits:", nrow(a), "\n\n")
cat(sprintf("%-10s %-6s %8s %8s %8s %8s\n", "stream", "build",
            "conv!=0", "warns>0", "msgs>0", "worst|g|"))
for (st in unique(a$stream)) {
  for (nm in c("base", "lane")) {
    s <- if (nm == "base") a[a$stream == st, ] else b[b$stream == st, ]
    cat(sprintf("%-10s %-6s %8d %8d %8d %8.2g\n", st, nm,
                sum(s$conv != 0), sum(s$warns > 0), sum(s$msgs > 0),
                max(s$grad)))
  }
}
cat("\n-- per-fit transitions, all 80\n")
tb <- table(base = ifelse(a$conv == 0, "conv0", "conv!=0"),
            lane = ifelse(b$conv == 0, "conv0", "conv!=0"))
print(tb)
new_bad <- which(a$conv == 0 & b$conv != 0)
fixed_bad <- which(a$conv != 0 & b$conv == 0)
cat("\nlane introduces a nonzero code on:", length(new_bad),
    if (length(new_bad)) paste0(a$stream[new_bad], " seed ", a$seed[new_bad],
                                collapse = "; "), "\n")
cat("lane removes a nonzero code on:", length(fixed_bad),
    if (length(fixed_bad)) paste0(a$stream[fixed_bad], " seed ",
                                  a$seed[fixed_bad], collapse = "; "), "\n")
cat("\nwarnings: base", sum(a$warns), " lane", sum(b$warns), "\n")
cat("messages: base", sum(a$msgs), " lane", sum(b$msgs), "\n")
cat("\n-- on the fits where BOTH converge cleanly, did logLik move?\n")
k <- a$conv == 0 & b$conv == 0
cat("n =", sum(k), " max |d logLik| =",
    format(max(abs(b$ll[k] - a$ll[k])), digits = 6), "\n")
cat("\n-- shortfall against sn by build\n")
for (nm in c("base", "lane")) {
  s <- if (nm == "base") a else b
  cat(sprintf("  %-5s fits short of sn by > 1e-6: %d of %d; worst %.6f\n",
              nm, sum(s$gap < -1e-6), nrow(s), -min(s$gap)))
}

# The convergence class rate over the 200-fit spread, base vs lane.
d <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/"
a <- read.csv(paste0(d, "falsealarm-base.csv"))
b <- read.csv(paste0(d, "falsealarm-fixed.csv"))
stopifnot(nrow(a) == nrow(b), identical(a$arm, b$arm),
          identical(a$seed, b$seed))
cat("fits:", nrow(a), "\n")
cat("conv != 0   base", sum(a$conv != 0), "  lane", sum(b$conv != 0), "\n")
cat("warnings>0  base", sum(a$warns > 0), "  lane", sum(b$warns > 0), "\n")
print(table(base = ifelse(a$conv == 0, "conv0", "conv!=0"),
            lane = ifelse(b$conv == 0, "conv0", "conv!=0")))
i <- which(a$conv == 0 & b$conv != 0)
cat("lane-only nonzero codes:", length(i),
    if (length(i)) paste(a$arm[i], a$seed[i], collapse = "; "), "\n")
j <- which(a$conv != 0 & b$conv == 0)
cat("base-only nonzero codes:", length(j),
    if (length(j)) paste(a$arm[j], a$seed[j], collapse = "; "), "\n")
cat("\nby arm, conv != 0 (base -> lane):\n")
for (nm in unique(a$arm)) {
  cat(sprintf("  %-9s %d -> %d\n", nm, sum(a$conv[a$arm == nm] != 0),
              sum(b$conv[b$arm == nm] != 0)))
}

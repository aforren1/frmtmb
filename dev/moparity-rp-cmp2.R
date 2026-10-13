# Compares dev/moparity-rp-cmp.R's base and lane runs.
a <- readRDS("dev/moparity-log/rp-cmp-base.rds")
b <- readRDS("dev/moparity-log/rp-cmp-lane.rds")
for (cp in names(a$est)) {
  cat("est", cp, "identical:", identical(a$est[[cp]], b$est[[cp]]),
      "max diff", max(abs(a$est[[cp]] - b$est[[cp]])), "\n")
}
ia <- unlist(a$init); ib <- unlist(b$init)
cat("inits identical:", identical(ia, ib), "max diff", max(abs(ia - ib)), "\n")
print(a$step[, c("stepsize__", "treedepth__", "n_leapfrog__", "energy__")])
print(b$step[, c("stepsize__", "treedepth__", "n_leapfrog__", "energy__")])
for (k in c("cen", "ncp")) {
  cat(k, "draws identical:", identical(a[[k]], b[[k]]), "; max |diff|",
      max(abs(a[[k]] - b[[k]])), "\n")
}

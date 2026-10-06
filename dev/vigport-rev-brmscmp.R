# Reviewer: the reviewer's brms rerun of brms_distreg (same runner,
# same budget, same port_seed) against the lane's brms records.
#
#   Rscript dev/vigport-rev-brmscmp.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
a <- readRDS(file.path(root, "vigport-port-out/r5/brms/brms_distreg.rds"))
b <- readRDS(file.path(root, "vigport-rev-out/brmsrerun/brms/brms_distreg.rds"))
for (id in names(b)) {
  if (is.null(b[[id]]$fit)) next
  fa <- a[[id]]$fit$fixef; fb <- b[[id]]$fit$fixef
  k <- intersect(names(fa), names(fb))
  cat(sprintf("%-20s identical fixef: %s; max |diff| %.3g; rhat %.3f / %.3f\n",
              id, identical(fa, fb), max(abs(fa[k] - fb[k])),
              a[[id]]$fit$rhat_max, b[[id]]$fit$rhat_max))
}

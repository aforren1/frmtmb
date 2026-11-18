# Lane squant, fix round 1: every frm_compat() row of quantile_elf next
# to asym_laplace's, with every extension loaded; prints the rows whose
# status differs. Run: Rscript dev/squant-compat2.R > dev/squant-out/compat2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  for (p in c("frmtmb.sample", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
              "frmtmb.ode", "frmtmb.survival", "frmtmb.spatial",
              "frmtmb.coupling")) library(p, character.only = TRUE)
})
other <- function(x, fam) ifelse(x$feature_a == fam, x$feature_b, x$feature_a)
a <- frm_compat("asym_laplace")
q <- frm_compat("quantile_elf")
A <- stats::setNames(a$status, other(a, "asym_laplace"))
Q <- stats::setNames(q$status, other(q, "quantile_elf"))
cat("rows: asym_laplace", length(A), "quantile_elf", length(Q),
    "same features", setequal(names(A), names(Q)), "\n")
k <- names(Q)[Q != A[names(Q)]]
for (f in k) cat(sprintf("%-14s asym_laplace %-12s quantile_elf %-12s | %s\n",
                         f, A[[f]], Q[[f]],
                         substr(q$note[other(q, "quantile_elf") == f], 1, 110)))
cat("quantile_elf status counts:\n")
print(table(Q))

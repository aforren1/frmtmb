# Compare two bitwise-battery recordings with identical(), element by
# element, and describe every fit that moved.
#   Rscript dev/predfix-bitwise-cmp.R base1 base2   (the control)
#   Rscript dev/predfix-bitwise-cmp.R base1 lane
a <- commandArgs(trailingOnly = TRUE)
dir <- "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/"
A <- readRDS(paste0(dir, "bitwise-", a[1], ".rds"))
B <- readRDS(paste0(dir, "bitwise-", a[2], ".rds"))
stopifnot(identical(names(A), names(B)))
els <- c("error", "ll", "par", "conv", "se_full", "se_fix", "sd_re")
same <- vapply(names(A), function(k) {
  identical(A[[k]][els], B[[k]][els])
}, NA)
tgt <- grepl("^TARGET", names(A))
err <- vapply(A, function(r) !is.null(r$error), NA)
fam <- sub("[|].*", "", names(A)[!tgt])
cat(a[1], "against", a[2], "\n")
cat("battery fits:", sum(!tgt), "in", length(unique(fam)), "families,",
    sum(err[!tgt]), "of them an identical error on both arms;",
    "bitwise identical:", sum(same[!tgt]), "of", sum(!tgt), "\n")
for (m in c("ML", "REML", "profile")) {
  k <- !tgt & grepl(paste0("[|]", m, "$"), names(A))
  cat(sprintf("  %-8s %3d of %3d identical (%d answered, %d errored)\n", m,
              sum(same[k]), sum(k), sum(k & !err), sum(k & err)))
}
nd <- grepl("[|]noint[|]|[|]noint_f[|]", names(A)) & !tgt
cat("  no-intercept designs:", sum(same[nd]), "of", sum(nd), "identical\n")
if (any(!same[!tgt])) {
  cat("battery fits that MOVED:\n")
  print(names(A)[!tgt][!same[!tgt]])
}
cat("\ntargeted fits (a column spread below 1e-3):\n")
cat(sprintf("%-46s %9s %18s %18s %12s %5s %5s\n", "fit", "identical",
            "logLik a", "logLik b", "b - a", "unitA", "unitB"))
for (k in names(A)[tgt]) {
  la <- A[[k]]$ll
  lb <- B[[k]]$ll
  num <- is.numeric(la) && is.numeric(lb)
  cat(sprintf("%-46s %9s %18.9f %18.9f %12.3e %5s %5s\n",
              sub("^TARGET[|]", "", k), same[[k]],
              if (num) la else NA, if (num) lb else NA,
              if (num) lb - la else NA,
              !is.null(A[[k]]$units), !is.null(B[[k]]$units)))
}

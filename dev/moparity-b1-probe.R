# Punch round 1, B1: the prior spellings update(prior = ) must merge.
# Usage: Rscript dev/moparity-b1-probe.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else
  "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(11)
dd <- data.frame(x = rnorm(60), z = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
tbl <- default_prior(y ~ x + z, data = dd, family = gaussian())
print(class(tbl))
print(as.data.frame(tbl))
sub <- tbl[tbl$class == "sigma", ]
print(class(sub))
sub$prior <- "exponential(2)"
print(sub)
print(unclass(frmtmb:::as_priorlist(sub)))
flat <- tbl
flat$prior <- ""
print(frmtmb:::as_priorlist(flat))
if (requireNamespace("brms", quietly = TRUE)) {
  cat("brms", format(packageVersion("brms")), "\n")
  print(brms::empty_prior())
  print(frmtmb:::as_priorlist(brms::empty_prior()))
}
print(unclass(empty_prior()))

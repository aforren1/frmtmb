# Reviewer, lane dblbar: re_formula with a || term on a fit whose blocks
# are written out, (1 | g) + (0 + x | g). Data seed 20261013.
# Usage: Rscript ... base|lane. Output: stdout
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
dd <- dblbar_data()
tc <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", substr(conditionMessage(e), 1, 140))
})
W <- frm(bf(y ~ x + (1 | g) + (0 + x | g)), family = gaussian(), data = dd)
nd <- dd[1:6, ]
full <- fitted(W, newdata = nd)
for (rf in c("~ (1 + x || g)", "~ (1 + x | g)", "~ (x || g)",
             "~ (1 | g)", "~ (0 + x | g)")) {
  r <- tc(fitted(W, newdata = nd, re_formula = as.formula(rf)))
  cat(sprintf("%-18s %s\n", rf, if (is.character(r)) r else
    sprintf("max |diff from full| %.3e", max(abs(r[, 1] - full[, 1])))))
}

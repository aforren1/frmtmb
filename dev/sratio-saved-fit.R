# Lane sratio: a fit saved by a build that held sratio()'s thresholds
# ordered, read back by a build that holds them unordered. The fit
# carries its family's closures (the log-density, the simulator and the
# threshold map), so the old fit should keep reading its own
# (first threshold, log increments) storage. This checks that nothing
# after the fit re-derives the storage from the family's name.
#
#   Rscript dev/sratio-saved-fit.R save  base
#   Rscript dev/sratio-saved-fit.R check /opt/rlib/lane-sratio
#     > dev/sratio-saved-fit-log.txt 2>&1   (both runs appended)
#
# The save run uses the round's base build (sratio ordered).
args <- commandArgs(trailingOnly = TRUE)
lib <- args[2]
.libPaths(c(if (lib != "base") lib, "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
cat("arm", args[1], "frmtmb from", find.package("frmtmb"), "\n")
path <- "/tmp/lanes/sratio/sratio-saved-fit.rds"
data("inhaler", package = "brms")
report <- function(fit) {
  fe <- fixef(fit)
  cat("fixef:\n")
  print(signif(fe[, 1:2], 10))
  P <- fitted(fit)[1:3, "Estimate", ]
  cat("fitted, rows 1 to 3:\n")
  print(signif(P, 10))
  cat("variables:", paste(head(variables(fit), 3), collapse = " "), "\n")
  h <- hypothesis(fit, "Intercept[3] - Intercept[2] = 0")$hypothesis
  cat("hypothesis Intercept[3] - Intercept[2]:", signif(h$Estimate, 10),
      "\n")
  set.seed(1)
  cat("simulate, first 10:", as.integer(simulate(fit)[1:10, 1]), "\n")
}
if (identical(args[1], "save")) {
  fit <- frm(bf(rating ~ period + carry + cs(treat)) + sratio(),
             data = inhaler)
  saveRDS(fit, path)
  report(fit)
} else {
  fit <- readRDS(path)
  report(fit)
}

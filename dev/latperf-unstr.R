# Lane latperf, fix round (review m6): the fit whose Schur complement base
# solve() refuses (unstr(t, g) with (1 | id), dev/review-latperf-why.R) leaves
# the sparse route for the dense one. The dense route must read what the base
# build reads there. Usage: Rscript dev/latperf-unstr.R lane|base
# Data: dev/review-latperf-algebra.R's mk(), seed 20261106.
arm <- commandArgs(TRUE)[1]
base <- c("C:/Users/adf44/source/r/rellib-r18",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lane") c("C:/Users/adf44/source/r/wt-latperf-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat(sprintf("lib: %s frmtmb %s\n", dirname(find.package("frmtmb")),
            format(packageVersion("frmtmb"))))
src <- readLines(
  "C:/Users/adf44/source/r/frmtmb-wt-latperf/dev/review-latperf-algebra.R")
eval(parse(text = src[grep("^mk <- function", src):
                        (grep("^designs <- list", src) - 1L)]))
d <- mk()
q <- function(e) suppressWarnings(suppressMessages(e))
fit <- q(frm(cnt0 ~ x + (1 | id) + unstr(t, g), family = poisson(), data = d))
if (arm == "lane") {
  cat("sparse route:", !is.null(frmtmb:::jcs_joint_cov(fit)), "\n")
}
se <- q(frm_linpred(fit, se.fit = TRUE))$se.fit
cat(sprintf("RESULT %s finite %d of %d, first %.10g, median %.10g, sum %.10g\n",
            arm, sum(is.finite(se)), length(se), se[1],
            stats::median(se, na.rm = TRUE), sum(se[is.finite(se)])))

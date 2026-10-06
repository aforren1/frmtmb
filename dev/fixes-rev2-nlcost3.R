# Reviewer of lane fixes, re-check: frm() time at n = 32000 for the
# model of dev/fixes-rev2-nlcost2.R, with the guard replaced by a no-op
# (the control) where the build has one.
#   Rscript dev/fixes-rev2-nlcost3.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
has <- exists("check_nl_identified", ns)
if (has) assignInNamespace("check_nl_identified",
                           function(spec, frame, prior) invisible(NULL),
                           ns = "frmtmb")
fo <- bf(y ~ a * exp(b * x) + c0, a ~ 1 + z, b ~ 1 + z, c0 ~ 1, nl = TRUE)
n <- 32000
set.seed(31)
d <- data.frame(x = runif(n), z = rnorm(n))
d$y <- 2 * exp((0.5 + 0.2 * d$z) * d$x) + 1 + rnorm(n, 0, 0.3)
t <- sapply(1:2, function(i) system.time(suppressMessages(suppressWarnings(
  frm(fo, data = d, start = list(beta = c(2, 0, 0.5, 0, 1))))))[["elapsed"]])
cat(find.package("frmtmb"), "guard replaced by no-op:", has, " fit s:", t, "\n")

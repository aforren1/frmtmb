# Review of lane parity2, item 3919: smooth terms the frame now parses
# with str2lang(): a non-syntactic column, an expression with a by
# variable, t2() of an expression. Arm base|lane.
#   Rscript dev/reviews/parity2-review-sexp.R base|lane
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-parity2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(4)
d <- data.frame(x = runif(80, 0.1, 2), z = rnorm(80),
                f = factor(rep(c("a", "b"), 40)))
d$y <- sin(2 * d$x) + 0.3 * d$z + rnorm(80, sd = 0.3)
d[["my x"]] <- d$x
tr <- function(label, e) {
  r <- tryCatch(suppressMessages(suppressWarnings(e)),
                error = function(err) err)
  if (inherits(r, "error")) {
    cat(label, ": ERROR ", conditionMessage(r), "\n", sep = "")
  } else {
    cat(label, ": logLik ", format(as.numeric(logLik(r)), digits = 8), "  ",
        paste(grep("^(bs|sds)_", variables(r), value = TRUE),
              collapse = " "), "\n", sep = "")
  }
  invisible(r)
}
tr("s(`my x`)", frm(y ~ s(`my x`), data = d))
tr("s(x)", frm(y ~ s(x), data = d))
tr("s(log(x), by = f)", frm(y ~ s(log(x), by = f) + f, data = d))
tr("t2(exp(x), z)", frm(y ~ t2(exp(x), z), data = d))
tr("s(x, z)", frm(y ~ s(x, z), data = d))
tr("sigma ~ s(exp(x))", frm(bf(y ~ z, sigma ~ s(exp(x))), data = d))

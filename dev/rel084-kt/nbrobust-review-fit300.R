# rel084-kt copy of dev/reviews/nbrobust-review-fit300.R: arm lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Review of lane nbrobust, check 3b: item 5224's data (seed 4, n = 300),
# negbinomial() and brmsfamily("negbinomial2"), at a few roundings of x,
# check_se "warning".
#   Rscript dev/reviews/nbrobust-review-fit300.R <base|lib>
args <- commandArgs(TRUE)
arm <- args[1]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/rellib-r23", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), " CORETYPE:",
    Sys.getenv("OPENBLAS_CORETYPE"), " R:", R.home(), "\n")
catch <- function(expr) {
  w <- character(0)
  mc <- character(0)
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    mc <<- c(mc, class(x)[1])
    invokeRestart("muffleMessage")
  })
  list(v = v, w = w, mc = mc)
}
set.seed(4)
d0 <- data.frame(x = stats::rnorm(300))
d0$y <- stats::rpois(300, exp(1 + 0.3 * d0$x))
for (k in c(0, 5, 9, 13, 17)) {
  d <- d0
  d$x <- d$x * (1 + 2e-15 * k)
  for (fn in c("negbinomial", "negbinomial2")) {
    fam <- if (fn == "negbinomial") negbinomial() else
      brmsfamily("negbinomial2")
    r <- catch(frm(y ~ x, family = fam, data = d))
    f <- r$v
    cat(sprintf("k %2d %-12s code %d betad %8.4f logLik %.12f nwarn %d msg %s\n",
                k, fn, f$opt$convergence, f$opt$par[["betad"]],
                as.numeric(logLik(f)), length(r$w),
                paste(unique(r$mc), collapse = ",")))
  }
}

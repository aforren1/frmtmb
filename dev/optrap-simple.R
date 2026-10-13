# Lane optrap: test-bcm-simple.R's SIMPLE_1 fit, its optimizer stages
# and its gradient at the optimum, on a library (run under each BLAS).
#   Rscript dev/optrap-simple.R <lib or "">
a <- commandArgs(TRUE)
.libPaths(c(if (length(a) && nzchar(a[1])) a[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib", dirname(find.package("frmtmb")), "| R", R.home(), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
for (h in Sys.glob("tests/testthat/helper-bcm*.R")) sys.source(h, envir = env)
src <- readLines("tests/testthat/test-bcm-simple.R")
i <- grep("^bcm_simple_sets <- function", src)
j <- grep("^test_that", src)[1] - 1L
eval(parse(text = src[i:j]), envir = env)
d <- env$bcm_simple_data()
f <- suppressWarnings(withCallingHandlers(
  frm(bf(k | trials(n) ~ 0 + set, s ~ 0 + set, t ~ 0 + set),
      family = env$bcm_simple_family(), data = d,
      start = env$bcm_simple_start(d), control = env$bcm_simple_control(),
      verbose = TRUE),
  message = function(m) {
    if (grepl("optimize|restart|search|escape|done", conditionMessage(m)))
      cat("  ", conditionMessage(m))
    invokeRestart("muffleMessage")
  }))
cat(sprintf("code %d obj %.10f max|grad| %.3g\n", f$opt$convergence,
            f$opt$objective, max(abs(f$obj$gr(f$opt$par)))))

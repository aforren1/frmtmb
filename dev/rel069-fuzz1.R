# The one fuzz spec that failed test-fuzz.R under OpenBLAS 0.3.32 (BLAS
# and LAPACK, 4 threads) on the 0.69.0 release tree: Gamma, REML, ar1,
# mo_int, dpar_x, seed 20379118. Fits that spec alone, in one process,
# and prints the outcome and the optimizer's trail.
#
#   Rscript dev/rel069-fuzz1.R <lib or "base"> [base lib] [reps]
# Run it with the reference R and with dev/ciharden-out/Rob0.3.32-lapack
# (OPENBLAS_NUM_THREADS=4).
a <- commandArgs(trailingOnly = TRUE)
base <- if (length(a) >= 2) a[2] else "C:/Users/adf44/source/r/rellib-r7"
.libPaths(c(if (!identical(a[1], "base")) a[1], base,
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
reps <- if (length(a) >= 3) as.integer(a[3]) else 1L
suppressMessages({library(testthat); library(frmtmb)})
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
env <- new.env(parent = asNamespace("frmtmb"))
attach(env, name = "fuzzhelpers")
for (h in c("helper-fuzz.R")) {
  sys.source(file.path(root, "tests/testthat", h), envir = env)
}
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "| R:", R.home(), "| threads:",
    Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
plan <- env$fuzz_plan(seed = 20260901L, size = 300L)
i <- which(plan$seed == 20379118L)
stopifnot(length(i) == 1L)
sp <- as.list(plan[i, ])
d <- env$fuzz_data(sp)
cat(env$fuzz_call_text(sp), "\n")
for (r in seq_len(reps)) {
  res <- env$fuzz_fit_one(sp, d)
  if (is.null(res$value)) {
    cat("rep", r, "ERROR:", substr(conditionMessage(res$error %||%
                                     simpleError("?")), 1, 300), "\n")
  } else {
    f <- res$value
    cat("rep", r, "OK code", f$opt$convergence, "objective",
        sprintf("%.10f", f$opt$objective), "warnings:",
        length(res$warnings), "\n")
  }
}

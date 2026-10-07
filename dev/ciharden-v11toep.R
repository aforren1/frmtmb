# test-v11.R "toep matches a hand-rolled reference (glmmTMB when it
# converges)" asserts one more thing with OpenBLAS as LAPACK than with
# the reference: glmmTMB's toep fit gives a finite likelihood there, so
# the comparison against it runs. How close is it to the bound (1e-4)
# where it runs, and to the hand-rolled reference's (1e-3)?
# Usage: Rscript dev/ciharden-v11toep.R <lib or base>
a <- commandArgs(TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1],
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("R:", R.home(), " threads:", Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
root <- "C:/Users/adf44/source/r/frmtmb-wt-ciharden"
e <- new.env()
sys.source(file.path(root, "tests/testthat/helper-reference.R"), envir = e)
dd <- e$sim_ar1_data(seed = 702, n_g = 80, n_t = 4, rho = 0.5)
fit <- suppressWarnings(
  frm(bf(y ~ 1 + toep(tim + 0 | g)) + gaussian(), data = dd))
ref <- suppressWarnings(
  glmmTMB::glmmTMB(y ~ 1 + toep(tim + 0 | g), data = dd, REML = FALSE))
ll_ref <- suppressWarnings(as.numeric(logLik(ref)))
cat(sprintf("frmtmb logLik %.10f glmmTMB %.10f |diff| %.3g (bound 1e-4)\n",
            as.numeric(logLik(fit)), ll_ref,
            abs(as.numeric(logLik(fit)) - ll_ref)))

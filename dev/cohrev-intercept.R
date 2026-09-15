## The one half-sentence in the shipped Rd that carries no measurement:
##
##   "a subject effect shifts both of that subject's conditions
##    together and cancels out of their difference, so it WIDENS THE
##    INTERCEPT'S interval and leaves the contrast's alone."
##
## The second half is measured over 148 replicates (a paired se ratio
## of 0.9965). The first half is not measured anywhere in the lane's
## files, and the TSV records only `coh_condb`. Two fits settle it.
.libPaths(c("C:/Users/adf44/source/r/cohrev-lib",
            "C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
setwd("extensions/frmtmb.coupling/tests/testthat")
test_that <- function(...) invisible(NULL)
source("test-coherence.R", local = TRUE)

se_of <- function(rhs, d) {
  fit <- suppressWarnings(
    frmtmb::frm(cp_bf_coh(rhs), family = cross_wishart(), data = d,
                se = TRUE))
  s <- sqrt(diag(stats::vcov(fit)))
  s[grep("^coh_", names(s))]
}
for (seed in c(2610L, 2611L, 2612L)) {
  d <- cp_cells(seed, n_sub = 16L, n_rep = 6L, sd_idcond = 0.5)
  a <- se_of("cond", d)
  b <- se_of("cond + (1 | id)", d)
  cat(sprintf("seed %d\n  cond            %s\n  cond + (1 | id) %s\n",
              seed,
              paste(sprintf("%s=%.5f", names(a), a), collapse = "  "),
              paste(sprintf("%s=%.5f", names(b), b), collapse = "  ")))
  nm <- intersect(names(a), names(b))
  cat(sprintf("  ratio with/without (1 | id): %s\n",
              paste(sprintf("%s=%.3f", nm, b[nm] / a[nm]),
                    collapse = "  ")))
}

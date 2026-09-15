## Item 2.6 of dev/extension-gaps-plan.md, the replicate sweep.
##
## Phase 0 fit the coupling ladder at ONE seed and found that the
## correct model and the survey's model both cover the condition
## contrast, so the survey's 0.77 was noise. What one seed cannot say is
## whether the correct model covers at the NOMINAL RATE, and what the
## misspecified rungs cost: their intervals are about half as wide, and
## an interval that is too narrow can still cover once.
##
## This script measures coverage AND width, per rung, over replicates.
##
## Two arms:
##   main  the simulator's own truth, sd(id:cond) = 0.2
##   null  the same truth with that one component switched OFF, which
##         makes the `id` rung correct as well. It is the control: if
##         the width gap is the missing term rather than the arithmetic,
##         it has to close here.
##
## Run one chunk:
##   COH_ARM=main COH_FROM=1 COH_TO=20 COH_OUT=dev/coh-recovery.tsv \
##     Rscript dev/coh-recovery.R
##
## Rows are appended, so several processes build one table. Each row
## carries its own seed.

.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
source("dev/coh-sim.R")

## Replicate r of an arm draws seed BASE + r. The two arms use
## different bases so that no data set is shared between them.
COH_SEED_BASE <- c(main = 20260910L, null = 20270910L)

arm <- Sys.getenv("COH_ARM", "main")
from <- as.integer(Sys.getenv("COH_FROM", "1"))
to <- as.integer(Sys.getenv("COH_TO", "1"))
out <- Sys.getenv("COH_OUT", "dev/coh-recovery.tsv")
rungs <- coh_rungs()
sd_ic <- if (identical(arm, "null")) 0 else coupling_truth$sd_idcond

record <- function(...) {
  vals <- list(...)
  fmt <- function(v) {
    if (is.numeric(v)) formatC(v, digits = 8, format = "g") else {
      as.character(v)
    }
  }
  line <- paste(paste0(names(vals), "=", vapply(vals, fmt, character(1))),
                collapse = "\t")
  cat(line, "\n", sep = "", file = out, append = TRUE)
  message("COH ", line)
}

for (r in seq(from, to)) {
  seed <- COH_SEED_BASE[[arm]] + r
  d <- coh_data(seed, sd_idcond = sd_ic)
  for (nm in names(rungs)) {
    res <- try(coh_fit_one(rungs[[nm]], d), silent = TRUE)
    if (inherits(res, "try-error")) {
      record(arm = arm, rep = r, seed = seed, rung = nm,
             ok = FALSE, msg = gsub("[\t\n]", " ",
                                    conditionMessage(attr(res,
                                                          "condition"))))
      next
    }
    record(arm = arm, rep = r, seed = seed, rung = nm, ok = TRUE,
           est = res$est, se = res$se, lo = res$lo, hi = res$hi,
           width = res$hi - res$lo,
           covers = res$lo <= coupling_truth$b_cond &&
             res$hi >= coupling_truth$b_cond,
           sd_id = res$sd_id, sd_idcond = res$sd_idcond,
           loglik = res$loglik, conv = res$conv, maxgrad = res$maxgrad,
           pdhess = res$pdhess, nbadse = res$nbadse, secs = res$secs)
  }
}

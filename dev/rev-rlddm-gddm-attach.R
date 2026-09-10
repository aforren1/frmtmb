# REVIEW, item 1.0b: what ndt_bound_attach() does to gddm(), the fifth
# family in its own package, which it does not refuse.
#
#   Rscript dev/rev-rlddm-gddm-attach.R <lib>
#
# No fit and no seed: this reads the family object the exported seam
# returns. gddm()'s density never reads `ndt_floor` (0 occurrences in
# R/gddm.R), so the fraction the logit produces is what its density
# would use as a time.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
g <- c("a", "a", "a", "b", "b")
bd <- ndt_bound(rt, list(ndt_group = ndt_bound_key(g)))

f0 <- gddm()
cat("gddm() before attach\n")
cat("  ndt link       :", f0[["links"]][["ndt"]][["name"]] %||%
      as.character(f0[["links"]][["ndt"]]), "\n")
cat("  aterm_data     :", is.null(f0[["aterm_data"]]), "(NULL?)\n")
cat("  init ndt       :",
    format(tryCatch(f0[["init_dpars"]][["ndt"]](rt, list()),
                    error = function(e) NA_real_)), "\n")

f1 <- ndt_bound_attach(f0, bd)
cat("\ngddm() after ndt_bound_attach(), which is NOT refused\n")
lk <- f1[["links"]][["ndt"]]
cat("  ndt link       :",
    if (is.character(lk)) lk else lk[["name"]], "\n")
cat("  bound record   :", paste(class(f1[["ndt_bound"]]),
                                collapse = "/"), "\n")
cat("  aterm_data adds:",
    paste(names(f1[["aterm_data"]](rt,
                                   list(ndt_group = ndt_bound_key(g)))),
          collapse = ", "), "\n")
cat("  init ndt       :", format(f1[["init_dpars"]][["ndt"]](rt, list())),
    "\n")
cat("\n  gddm's own sources mention ndt_floor 0 times, so its density\n")
cat("  reads dpars$ndt directly. On a plain logit that is a number in\n")
cat("  (0, 1) SECONDS, not a fraction of the row's own floor.\n")

`%||%` <- function(a, b) if (is.null(a)) b else a

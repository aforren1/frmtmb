# Reviewer: the plausibility table filters brms fits by Rhat but not
# frmtmb fits by convergence. List every model row of the spell pass
# (the pass plausibility.R reads) whose optimizer code is nonzero or
# whose warnings mention convergence or the gradient.
#
#   Rscript dev/vigport-rev-conv.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
for (pass in c("results-spell", "results-keepprior")) {
  res <- list()
  for (f in list.files(file.path(root, "vigport-port-out/r5", pass),
                       "[.]rds$", full.names = TRUE)) res <- c(res, readRDS(f))
  cat("##", pass, "\n")
  n <- 0
  for (r in res) {
    if (!identical(r$kind, "model") || !identical(r$status, "OK")) next
    n <- n + 1
    w <- grep("converg|gradient|Hessian", r$warnings, value = TRUE,
              ignore.case = TRUE)
    cv <- r$fit$conv
    if (length(w) || (!is.null(cv) && !is.na(cv) && cv != 0)) {
      cat(sprintf("%-26s conv=%s  %s\n", r$id, format(cv),
                  substr(paste(w, collapse = " | "), 1, 110)))
    }
  }
  cat("model rows OK:", n, "\n\n")
}

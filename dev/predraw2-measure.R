# Lane predraw2: predict() through the package on each case of a group
# (dev/predraw2-cases.R). Per fit: the draw check's verdict and median
# excess (lane build only), the reweighting's messages and warnings, the
# time of the first predict() (the check and the side law run there) and
# of a second, and a fingerprint of the draws (each draw's mean over the
# rows, and the first row's draws) for a bitwise comparison with the base
# build. 200 draws, seed 5.
# Usage: Rscript dev/predraw2-measure.R <lib or "base"> <group> <out.rds>
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
source("dev/predraw2-cases.R")
group <- args[2]
out <- args[3]
res <- list()
for (id in case_ids(group)) {
  cf <- tryCatch(case_fit(id), error = function(e) NULL)
  fit <- cf$fit
  if (is.null(fit)) {
    cat(id, "fit failed\n")
    next
  }
  conds <- character(0)
  pr <- function() {
    set.seed(5)
    withCallingHandlers(
      tryCatch(predict(fit, ndraws = 200L, summary = FALSE),
               error = function(e) {
                 conds <<- c(conds, paste("ERROR", conditionMessage(e)))
                 NULL
               }),
      warning = function(w) {
        conds <<- c(conds, paste0("W[", class(w)[1], "] ",
                                  substr(conditionMessage(w), 1, 160)))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        conds <<- c(conds, paste0("M[", class(m)[1], "] ",
                                  substr(conditionMessage(m), 1, 160)))
        invokeRestart("muffleMessage")
      })
  }
  t0 <- proc.time()[[3]]
  p <- pr()
  t1 <- proc.time()[[3]]
  p2 <- pr()
  t2 <- proc.time()[[3]]
  chk <- if (is.environment(fit$cache)) fit$cache[["draw_check"]]
  fp <- if (is.null(p)) NULL else if (length(dim(p)) == 2L) {
    list(rm = rowMeans(p), c1 = p[, 1])
  } else {
    list(rm = apply(p, 1, mean), c1 = as.vector(p[, 1, ]))
  }
  res[[id]] <- list(id = id, fire = chk$fire, excess = chk$excess,
                    rel = chk$rel, n = chk$n,
                    conds = unique(conds), t_first = t1 - t0,
                    t_second = t2 - t1, fp = fp,
                    same_twice = identical(p, p2))
  cat(sprintf("%-38s fire %-5s excess %8.3f t1 %.2f t2 %.2f %s\n", id,
              format(chk$fire), chk$excess %||% NA_real_, t1 - t0, t2 - t1,
              paste(substr(unique(conds), 1, 60), collapse = " | ")))
}
saveRDS(res, out)
cat("fits", length(res), "\n")

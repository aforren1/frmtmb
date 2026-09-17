## Reviewer: fit the review models in REAL brms (fresh compiles, pinlib
## ahead of the user library), and save the brmsfits.
##   Rscript dev/brmsnames-rev-brms.R [C1 C2 ...]
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(packageVersion("StanHeaders") == "2.32.10")
suppressMessages(library(brms))
source("dev/brmsnames-rev-data.R")
av <- commandArgs(trailingOnly = TRUE)
ms <- rev_models("brms")
if (length(av)) ms <- ms[av]
dir.create("dev/brmsnames-rev-log", showWarnings = FALSE)
for (nm in names(ms)) {
  out <- sprintf("dev/stan-cache/brmsnames-rev-brms-%s.rds", nm)
  if (file.exists(out)) { cat(nm, "exists\n"); next }
  pr <- if (nm == "C5") c(prior(normal(2, 1), nlpar = "a"),
                          prior(normal(0.7, 1), nlpar = "b")) else NULL
  t0 <- Sys.time()
  fit <- tryCatch(brm(ms[[nm]]$f, data = rev_dd, prior = pr,
                      chains = 2, iter = 60, warmup = 30, refresh = 0,
                      seed = 1, backend = "rstan",
                      init = if (nm == "C5") 0 else "random"),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(nm, "ERROR", conditionMessage(fit), "\n")
    saveRDS(fit, out)
    next
  }
  saveRDS(fit, out)
  cat(nm, "fitted in", format(Sys.time() - t0), "\n")
  cat("  variables:", head(variables(fit), 60), "\n")
}
cat("DONE\n")

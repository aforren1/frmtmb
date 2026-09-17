## Reviewer recheck: real brms 2.23.0 fits of the hostile-name models, so
## names are brms's own output, not a reconstruction. pinlib ahead of the
## user library (StanHeaders 2.32.10 asserted), fresh compile.
##   Rscript dev/brmsnames-rev2-brms.R B1 B2 ...
## Data seed 52 (dev/brmsnames-rev2-data.R), chains 1, iter 60, seed 1.
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(packageVersion("StanHeaders") == "2.32.10",
          packageVersion("brms") == "2.23.0")
suppressMessages(library(brms))
source("dev/brmsnames-rev2-data.R")
d <- rev2_data()
ms <- rev2_models()
keys <- commandArgs(trailingOnly = TRUE)
for (k in keys) {
  M <- ms[[k]]
  f <- rev2_formula(k, M, brms::bf, brms::mvbf)
  fam <- get(M$fam, asNamespace("brms"))()
  out <- sprintf("dev/stan-cache/brmsnames-rev2-brms-%s.rds", k)
  t0 <- proc.time()[[3]]
  res <- tryCatch({
    pri <- if (k == "B5") {
      c(set_prior("normal(2, 1)", nlpar = "a1"),
        set_prior("normal(0, 1)", nlpar = "b2"))
    }
    b <- brm(f, data = d, family = fam, prior = pri, chains = 1, iter = 60,
             warmup = 30, seed = 1, refresh = 0, backend = "rstan",
             init = 0)
    saveRDS(b, out)
    list(ok = TRUE, variables = variables(b))
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  cat("\n==", k, "== elapsed", round(proc.time()[[3]] - t0), "s\n")
  if (res$ok) print(res$variables) else cat("ERROR:", res$msg, "\n")
}

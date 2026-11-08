# The specs of the gated fuzz plan (seed 20260901, 300 rows) that reach
# nested quadrature, run on a given build with this worktree's helper.
#   Rscript dev/nestquad-fuzz-subset.R <lib|base>
a <- commandArgs(TRUE)
lib <- "C:/Users/adf44/source/r/wt-nestquad-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(lib, base, ul))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-nestquad/tests/testthat"
env <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-warnings.R", "helper-fuzz.R")) {
  sys.source(file.path(wt, h), envir = env)
}
plan <- env$fuzz_plan(seed = 20260901L, size = 300L)
hit <- plan$mode == "quadrature" & plan$re == "nested"
sub <- plan[hit, , drop = FALSE]
cat("plan rows:", nrow(plan), " quadrature+nested:", sum(hit), "\n")
for (i in seq_len(nrow(sub))) {
  cat("SPEC", env$fuzz_spec_key(sub[i, ]), "seed", sub$seed[i], "\n")
}
t0 <- proc.time()[[3]]
res <- env$fuzz_run(sub)
sm <- env$fuzz_summary(res)
cat("elapsed", round(proc.time()[[3]] - t0), "s\n")
print(sm$counts)
for (f in sm$triaged) cat(f$class, "|", env$fuzz_format(f), "\n")
cat("END\n")

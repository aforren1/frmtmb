# Fix round 1 copy, which also takes this lane's FUZZ_KNOWN_REFUSAL.
# Which specs of the merged gated plan (cdfs's helper-fuzz.R, read in
# place, read-only) reach the two changed paths, and how do they fare on
# a given build?
#   Rscript dev/fuzzfix-fuzz-subset.R <lib|base> [list]
# "list" prints the subset and stops.
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
cdfs <- "C:/Users/adf44/source/r/frmtmb-wt-cdfs/tests/testthat"
env <- new.env()
for (h in c("helper-warnings.R", "helper-fuzz.R")) {
  sys.source(file.path(cdfs, h), envir = env)
}
# the merged helper: cdfs's, with this lane's fuzz_design()
# "cdfsgen" keeps cdfs's generator, to show what the specs do without
# this lane's fuzz_design(); it also limits the run to the ar1 specs
cdfsgen <- length(a) > 1 && a[2] == "cdfsgen"
if (!cdfsgen) {
  mine <- new.env()
  sys.source("tests/testthat/helper-fuzz.R", envir = mine)
  env$fuzz_design <- mine$fuzz_design
  environment(env$fuzz_design) <- env
  env$FUZZ_KNOWN_REFUSAL <- mine$FUZZ_KNOWN_REFUSAL
}
plan <- env$fuzz_plan(seed = 20260901L, size = 300L)
spike <- c("gaussian", "student", "lognormal", "weibull", "Gamma", "Beta")
hit_w <- plan$re == "ar1" & plan$aterm == "weights" & plan$family %in% spike
hit_q <- plan$mode == "quadrature" & plan$re %in% c("nested")
sub <- plan[hit_w | (hit_q & !cdfsgen), , drop = FALSE]
cat("plan rows:", nrow(plan), " ar1+weights+continuous:", sum(hit_w),
    " quadrature+nested:", sum(hit_q), "\n")
for (i in seq_len(nrow(sub))) {
  cat("SPEC", env$fuzz_spec_key(sub[i, ]), "seed", sub$seed[i], "\n")
}
if (length(a) > 1 && a[2] == "list") quit(save = "no")
t0 <- proc.time()[[3]]
res <- env$fuzz_run(sub)
sm <- env$fuzz_summary(res)
cat("elapsed", round(proc.time()[[3]] - t0), "s\n")
print(sm$counts)
for (f in sm$triaged) cat(f$class, "|", env$fuzz_format(f), "\n")
cat("END\n")

# Reviewer: run fuzz specs of the merged plan (cdfs + fuzzfix trial
# merge) by seed, on the merged build, with the merged helpers.
#   Rscript dev/fuzzfix-rev-merge-fuzz.R <lib> <pkgdir> <seed|subset> ...
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
env <- new.env()
for (h in c("helper-warnings.R", "helper-fuzz.R")) {
  sys.source(file.path(a[2], "tests/testthat", h), envir = env)
}
plan <- env$fuzz_plan(seed = 20260901L, size = 300L)
if (a[3] == "subset") {
  spike <- c("gaussian", "student", "lognormal", "weibull", "Gamma", "Beta")
  sub <- plan[(plan$re == "ar1" & plan$aterm == "weights" &
                 plan$family %in% spike) |
                (plan$mode == "quadrature" & plan$re == "nested"), ,
              drop = FALSE]
  # seeds after "subset" are left out (run them alone, memory capped)
  sub <- sub[!sub$seed %in% as.integer(a[-(1:3)]), , drop = FALSE]
} else {
  sub <- plan[plan$seed %in% as.integer(a[-(1:2)]), , drop = FALSE]
}
cat("specs:", nrow(sub), "\n")
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

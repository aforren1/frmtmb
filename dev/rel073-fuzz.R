# Consolidation item 3: the gated fuzz plan (300 specs, seed 20260901)
# on the merged build, every triaged finding with its class, and the
# three specs the lanes named. From the cdfs reviewer's
# dev/cdfs-rev-fuzz.R, with the library as an argument and this tree's
# test helpers.
#   Rscript dev/rel073-fuzz.R <lib>   (log dev/rel073-log/fuzz.txt)
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_FUZZ = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("lib: frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
testthat::source_test_helpers(
  "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat",
  env = globalenv())
plan <- fuzz_plan(seed = 20260901L, size = 300L)
t0 <- proc.time()[["elapsed"]]
res <- fuzz_run(plan)
sm <- fuzz_summary(res)
cat("specs:", length(plan), "| seconds:", round(proc.time()[["elapsed"]] - t0),
    "\n")
cls <- vapply(sm$triaged, function(f) f$class, "")
print(table(cls))
ids <- vapply(sm$triaged, function(f) f$known_id %||% "", "")
print(table(paste(cls, ids)))
for (f in sm$triaged) {
  cat(f$class, if (!is.null(f$known_id)) paste0("[", f$known_id, "]"), "|",
      fuzz_format(f), "\n")
}
for (s in c(20322452L, 20502220L, 20535438L)) {
  hit <- Filter(function(f) grepl(paste0("seed ", s), fuzz_format(f)),
                sm$triaged)
  i <- which(plan$seed == s)
  cat("SEED", s, "in plan:", length(i) == 1L,
      if (length(i) == 1L) fuzz_spec_key(plan[i, ]), "\n")
  cat("SEED", s, ":", if (length(hit)) {
    paste(vapply(hit, function(f) paste0(f$class, if (!is.null(f$known_id))
      paste0(" [", f$known_id, "]")), ""), collapse = "; ")
  } else "no finding", "\n")
}
cat("END\n")

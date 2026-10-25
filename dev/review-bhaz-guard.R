# Reviewer: guardrail 1 (a) as test-guardrail1.R in a core-only library,
# (b) every probe with ALL extensions loaded in one child.
#   Rscript dev/review-bhaz-guard.R core  > dev/review-bhaz-log/guard-core.txt
#   Rscript dev/review-bhaz-guard.R all   > dev/review-bhaz-log/guard-all.txt
WT <- "C:/Users/adf44/source/r/frmtmb-wt-bhaz"
mode <- commandArgs(trailingOnly = TRUE)[[1L]]
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
if (mode == "core") {
  .libPaths(c(file.path(WT, "dev/review-bhaz-corelib"), user))
} else {
  .libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
              "C:/Users/adf44/source/r/rellib-r12", user))
}
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
cat("lib:", find.package("frmtmb"), " libPaths:", .libPaths(), "\n")
exts <- c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
          "frmtmb.ode", "frmtmb.sample", "frmtmb.spatial", "frmtmb.survival")
cat("extensions installed:",
    exts[nzchar(vapply(exts, function(p) system.file(package = p), ""))], "\n")
f <- file.path(WT, "tests/testthat/test-guardrail1.R")
if (mode == "core") {
  res <- testthat::test_file(f, package = "frmtmb",
                             env = testthat::test_env("frmtmb"),
                             reporter = "silent")
  df <- as.data.frame(res)
  print(df[, c("test", "nb", "passed", "failed", "skipped", "error")])
  quit(save = "no")
}
# all: load the helpers from the test file without running its tests
ex <- parse(f)
keep <- vapply(ex, function(e) is.call(e) && identical(e[[1L]],
                                                         as.name("<-")), NA)
env <- new.env()
for (e in ex[keep]) eval(e, env)
pr <- env$guard_probes()
res0 <- env$guard_run(pr)
resA <- env$guard_run(pr, exts)
cat("child ok core", attr(res0, "ok"), " all", attr(resA, "ok"), "\n")
cls0 <- mapply(env$guard_classify, res0, pr$name)
clsA <- mapply(env$guard_classify, resA, pr$name)
want <- ifelse(startsWith(cls0, "ext:"),
               ifelse(pr$ext == "ok", "core", "refused"), cls0)
bad <- which(clsA != want | is.na(resA))
cat("probes", nrow(pr), " ext homes",
    paste(unique(cls0[startsWith(cls0, "ext:")]), collapse = " "),
    " mismatches with all extensions loaded:", length(bad), "\n")
for (i in bad) cat(" ", pr$kind[i], pr$name[i], pr$code[i], "->", resA[i],
                   "\n")
b <- which(pr$name == "bhaz")
for (i in b) cat("bhaz probe:", pr$code[i], "\n   core:", res0[i],
                 "\n   all :", resA[i], "\n")

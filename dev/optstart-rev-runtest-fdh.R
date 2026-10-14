# Reviewer of lane optstart: run one test file with fd_hessian()'s retake
# branch logged (step 11 of its body, after the branch, before E).
#   Rscript dev/optstart-rev-runtest-fdh.R <lib> <package> <file> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
p <- a[2]
suppressMessages(library(p, character.only = TRUE))
suppressMessages(library(testthat))
cat(a[4], "lib:", dirname(find.package(p)), " frmtmb from",
    dirname(find.package("frmtmb")), " R", R.home(), "\n")
FDH <- new.env()
FDH$log <- list()
ns <- asNamespace("frmtmb")
if ("tau_band" %in% names(formals(ns$fd_hessian))) {
  trace("fd_hessian", at = 11, print = FALSE, where = ns, tracer = quote({
    br <- exists("H2", inherits = FALSE)
    .GlobalEnv$FDH$log[[length(.GlobalEnv$FDH$log) + 1L]] <- data.frame(
      n = n, redo = length(redo), branch = br,
      h2_factors = br && exists("lmin", inherits = FALSE),
      lmin = if (br && exists("lmin", inherits = FALSE)) lmin else NA,
      noise = if (br) norm(abs(H2 - t(H2)) / 2, "2") else NA,
      kept = br && identical(H, H2))
  }))
}
res <- testthat::test_file(a[3], package = p, env = testthat::test_env(p),
                           reporter = "silent")
df <- as.data.frame(res)
for (r in res) {
  for (e in r$results) {
    if (!inherits(e, "expectation_success")) {
      cat("--", class(e)[1], "|", r$test, "\n   ",
          substr(gsub("\n", " ", conditionMessage(e)), 1, 300), "\n")
    }
  }
}
L <- if (length(FDH$log)) do.call(rbind, FDH$log) else NULL
cat(sprintf("FDH calls %d, with a rise above the band %d, branch %d, kept %d\n",
            NROW(L), sum(L$redo > 0), sum(L$branch), sum(L$kept)))
if (!is.null(L) && any(L$branch)) print(L[L$branch, ], digits = 4)
cat(sprintf("RESULT %s pass %d fail %d err %d skip %d warn %d\n",
            basename(a[3]), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))

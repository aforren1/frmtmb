# Lane latperf, fix round: one test file in one process against an
# in-process mutant of R/joint-sparse.R, to see the tests fail on the faults
# the review named. The approach of dev/review-latperf-test.R.
# Usage: Rscript dev/latperf-mutant.R <test file> <mutant>
#   none      the lane build
#   novf      jcs_rowvar() drops the outer part g' Vf g (review m1)
#   norepair  jcs_joint_cov() ignores the SE check's covariance (review m2)
#   perrow    a pair outside the pattern sends its row to one direct solve,
#             the route before the fix of review M1
local({
  .libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
              "C:/Users/adf44/source/r/rellib-r18",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(NOT_CRAN = "true")
  a <- commandArgs(TRUE)
  tf <- a[1]
  mut <- a[2]
  suppressMessages(library(testthat))
  suppressMessages(library(frmtmb))
  cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
      format(packageVersion("frmtmb")), "mutant:", mut, "file:", tf, "\n")
  ns <- asNamespace("frmtmb")
  put <- function(nm, fn) {
    environment(fn) <- ns
    utils::assignInNamespace(nm, fn, ns = "frmtmb")
  }
  src_sub <- function(nm, from, to) {
    txt <- deparse(get(nm, envir = ns), control = "useSource")
    hit <- grepl(from, txt, fixed = TRUE)
    if (sum(hit) != 1L) stop("mutant site not unique in ", nm, ": ", from)
    txt[hit] <- sub(from, to, txt[hit], fixed = TRUE)
    put(nm, eval(parse(text = txt)))
  }
  switch(mut,
    none = NULL,
    novf = src_sub("jcs_rowvar", "v <- rowSums((G %*% S[[\"Vf\"]]) * G)",
                   "v <- 0 * rowSums((G %*% S[[\"Vf\"]]) * G)"),
    norepair = src_sub("jcs_joint_cov", "repair <- !is.null(Vp) &&",
                       "repair <- FALSE && !is.null(Vp) &&"),
    perrow = {
      src_sub("jcs_rowvar_inner",
              "w[out] <- jcs_W_pairs(S, a[out] + 1, b[out] + 1)",
              "w[out] <- 0")
      src_sub("jcs_rowvar_inner", "v[as.integer(rownames(s))] <- s[, 1]",
              paste("v[as.integer(rownames(s))] <- s[, 1];",
                    "bad <- unique(row[e1][out]);",
                    "if (length(bad)) v[bad] <- jcs_rowvar_solve(S,",
                    "Ar[bad, , drop = FALSE], ri)"))
    },
    stop("unknown mutant ", mut))
  path <- file.path("C:/Users/adf44/source/r/frmtmb-wt-latperf/tests/testthat",
                    tf)
  res <- testthat::test_file(path, package = "frmtmb",
                             env = testthat::test_env("frmtmb"),
                             reporter = "silent", stop_on_failure = FALSE)
  df <- as.data.frame(res)
  cat(sprintf("RESULT mutant %s file %s pass %d fail %d err %d skip %d warn %d\n",
              mut, tf, sum(df$passed), sum(df$failed), sum(df$error),
              sum(df$skipped), sum(df$warning)))
  bad <- df[df$failed > 0 | df$error, "test"]
  for (b in bad) cat("  FAILED:", b, "\n")
})

# Review of lane selimit2, item 6: R CMD check run 1 ended abnormally after
# test-gp-noncentered.R recorded failures at lines 481, 527, 528. This runs
# the installed lane build's tests in ONE session, in test_check()'s order,
# from the first file through test-gp-noncentered.R (or every file with
# ALL=1), with R CMD check --as-cran's environment (NOT_CRAN unset,
# _R_CHECK_LIMIT_CORES_, OMP_THREAD_LIMIT=2), on a copy of the test tree.
#   Rscript selimit2-review-session.R <copy of tests/testthat> <tag>
a <- commandArgs(trailingOnly = TRUE)
dir <- a[1]
tag <- a[2]
.libPaths(c("C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.unsetenv("NOT_CRAN")
Sys.setenv("_R_CHECK_LIMIT_CORES_" = "true", OMP_THREAD_LIMIT = "2",
           "_R_CHECK_CRAN_INCOMING_REMOTE_" = "FALSE")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(utils::packageVersion("frmtmb")), "\n")
files <- sort(list.files(dir, "^test-.*\\.R$"))
last <- if (nzchar(Sys.getenv("ALL"))) length(files) else
  match("test-gp-noncentered.R", files)
keep <- files[seq_len(last)]
rx <- NULL
t0 <- proc.time()[["elapsed"]]
res <- testthat::test_dir(dir, package = "frmtmb", load_package = "installed",
                          filter = rx, reporter = "summary",
                          stop_on_failure = FALSE, stop_on_warning = FALSE)
df <- as.data.frame(res)
cat("\nTAG", tag, "files", length(unique(df$file)), "tests", nrow(df),
    "failed", sum(df$failed), "error", sum(df$error), "skipped",
    sum(df$skipped), "elapsed", round(proc.time()[["elapsed"]] - t0), "\n")
bad <- df[df$failed > 0 | df$error, c("file", "test", "failed", "error")]
print(bad)
gp <- df[df$file == "test-gp-noncentered.R", c("test", "nb", "failed",
                                                 "error", "real")]
print(gp)
cat("mem (MB) at end:", round(sum(gc()[, 2])), "\n")

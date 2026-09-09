# Reviewer instrument: the BEFORE arm of lane tmbstan, reconstructed
# without editing the worktree. The suite's testthat directory is
# copied to a scratch tree, skip_sampler() there is replaced by the two
# lines it had at the base commit (the installed check only), and the
# named file is run against the poisoned detector.
#
# usage: Rscript dev/rev-tmbstan-before.R <testfile> <arm>
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
.libPaths(c(LIB, "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
args <- commandArgs(trailingOnly = TRUE)
f <- args[[1L]]
arm <- args[[2L]]

library(testthat)
library(frmtmb)
suppressMessages(library(frmtmb.sample))
root <- "C:/Users/adf44/source/r/frmtmb-wt-tmbstan"
Sys.setenv(FRMTMB_STAN_CACHE = file.path(root, "dev", "stan-cache"))

src <- file.path(root, "extensions/frmtmb.sample/tests/testthat")
dst <- file.path(tempdir(), "rev-before", "testthat")
dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
file.copy(src, dirname(dst), recursive = TRUE)

hp <- file.path(dst, "helper-sampling.R")
txt <- readLines(hp, warn = FALSE)
a <- grep("^skip_sampler <- function[(][)] [{]$", txt)
b <- a + which(txt[(a + 1L):length(txt)] == "}")[[1L]]
stopifnot(length(a) == 1L, txt[[b]] == "}")
cat("## replacing helper lines", a, "to", b, "\n")
old <- c("skip_sampler <- function() {",
         "  testthat::skip_if_not_installed(\"tmbstan\")",
         "  testthat::skip_if_not_installed(\"rstan\")",
         "}")
writeLines(c(txt[seq_len(a - 1L)], old, txt[(b + 1L):length(txt)]), hp)

env <- environment(frmtmb.sample:::tmbstan_build_broken)
if (identical(arm, "broken")) assign("cached", TRUE, envir = env)
cat("## detector:", frmtmb.sample:::tmbstan_build_broken(), "arm:", arm,
    "\n")

setwd(dst)
r <- as.data.frame(test_file(f, package = "frmtmb.sample",
                             reporter = "silent"))
for (i in seq_len(nrow(r))) {
  cat(paste("ROW", paste0("before-", arm), f, r$test[[i]], r$passed[[i]],
            r$failed[[i]], as.integer(r$error[[i]]), r$skipped[[i]],
            sep = "\t"), "\n", sep = "")
}
cat("## TOTAL", paste0("before-", arm), f, nrow(r), sum(r$passed),
    sum(r$failed), sum(r$error), sum(r$skipped), "\n")

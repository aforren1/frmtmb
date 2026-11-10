# Lane support: one block of a generated bin-2 file, run as
# dev/brmsport-bin2-run.R runs a whole file, staged in
# dev/support-bin2-run/<package> so that the shared stage directory is
# not touched.
#
#   Rscript dev/support-bin2-run.R <package> <generated file> <first line>
#     <last line> <record.tsv>
#
# The file's header (every line before its first test_that(): the skip
# and brms's imports) is kept, then lines <first>..<last>, one block.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
src <- readLines(a[2])
first <- as.integer(a[3])
last <- as.integer(a[4])
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p)
tdir <- if (p == "frmtmb") "tests/testthat" else
  file.path("extensions", p, "tests/testthat")
stage <- file.path("dev/support-bin2-run", p)
unlink(stage, recursive = TRUE)
dir.create(stage, recursive = TRUE)
file.copy(list.files(tdir, "^helper.*[.]R$", full.names = TRUE), stage)
file.copy("dev/brmsport-bin2-adapters.R",
          file.path(stage, "helper-brms-bin2-adapters.R"))
f2 <- file.path(stage, "test-support-block.R")
head_end <- grep("^test_that", src)[1] - 1L
writeLines(c(src[seq_len(head_end)], src[first:last]), f2)
rec <- file.path(normalizePath(dirname(a[5]), winslash = "/"),
                 basename(a[5]))
unlink(rec)
Sys.setenv(FRMTMB_BRMSPORT_RECORD = rec)
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(utils::packageVersion("frmtmb")), "\n")
cat("betareg:", format(utils::packageVersion("betareg")), "at",
    dirname(system.file(package = "betareg")), "\n")
r <- as.data.frame(test_file(normalizePath(f2, winslash = "/"), package = p,
                             env = testthat::test_env(p),
                             reporter = "silent", stop_on_failure = FALSE))
cat("RESULT ", p, " pass=", sum(r$passed), " fail=", sum(r$failed),
    " err=", sum(r$error), " skip=", sum(r$skipped), "\n", sep = "")
cat("RECORD\n")
cat(readLines(rec), sep = "\n")

# Release reviewer, 0.73.0: one generated bin-2 file in record mode on a
# chosen library, staged in dev/rel073-rev-bin2-run/<tag>/<package>/ (a
# copy of dev/brmsport-bin2-run.R that writes nothing else).
#   Rscript dev/rel073-rev-bin2.R <lib> <tag> <package> <test file>
# Record: dev/rel073-rev-log/bin2-rec-<tag>-<package>-<file>.tsv
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]; tag <- a[2]; p <- a[3]; f <- a[4]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p,
           FRMTMB_PORT_ROOT = "C:/Users/adf44/source/r/frmtmb-wt-release",
           FRMTMB_STAN_CACHE = file.path(
             "C:/Users/adf44/AppData/Local/Temp/3/claude",
             "c--Users-adf44-source-r-frmtmb",
             "7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad",
             "rev073-stan-cache"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-release")
tdir <- if (p == "frmtmb") "tests/testthat" else
  file.path("extensions", p, "tests/testthat")
stage <- file.path("dev/rel073-rev-bin2-run", tag, p, sub("[.]R$", "", basename(f)))
unlink(stage, recursive = TRUE)
dir.create(stage, recursive = TRUE)
file.copy(list.files(tdir, "^helper.*[.]R$", full.names = TRUE), stage)
file.copy("dev/brmsport-bin2-adapters.R",
          file.path(stage, "helper-brms-bin2-adapters.R"))
file.copy(f, stage)
f2 <- normalizePath(file.path(stage, basename(f)), winslash = "/")
rec <- normalizePath(file.path("dev/rel073-rev-log", paste0(
  "bin2-rec-", tag, "-", p, "-", sub("[.]R$", "", basename(f)), ".tsv")),
  winslash = "/", mustWork = FALSE)
unlink(rec)
Sys.setenv(FRMTMB_BRMSPORT_RECORD = rec)
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(utils::packageVersion("frmtmb")), p,
    format(utils::packageVersion(p)), "\n")
r <- as.data.frame(test_file(f2, package = p, env = testthat::test_env(p),
                             reporter = "silent", stop_on_failure = FALSE))
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " blocks=", nrow(r), "\n", sep = "")

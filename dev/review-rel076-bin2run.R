# Review: one generated bin-2 file against rellib-r14, as
# dev/brmsport-bin2-run.R runs it in tier mode, but staged under the
# review's own scratch so nothing of the release tree is written.
#   Rscript dev/review-rel076-bin2run.R <package> <bin-2 file>
rel <- "C:/Users/adf44/source/r/frmtmb-wt-rel076"
lib <- "C:/Users/adf44/source/r/rellib-r14"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p)
tdir <- file.path(rel, if (p == "frmtmb") "tests/testthat" else
  file.path("extensions", p, "tests/testthat"))
stage <- file.path(rel, "dev/review-rel076-out/bin2stage",
                   sub("[.]R$", "", basename(f)))
unlink(stage, recursive = TRUE)
dir.create(stage, recursive = TRUE)
file.copy(list.files(tdir, "^helper.*[.]R$", full.names = TRUE), stage)
file.copy(file.path(rel, "dev/brmsport-bin2-adapters.R"),
          file.path(stage, "helper-brms-bin2-adapters.R"))
file.copy(f, stage)
f2 <- normalizePath(file.path(stage, basename(f)), winslash = "/")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
port_loaded <- grep("^frmtmb[.]", rownames(utils::installed.packages(lib)),
                    value = TRUE)
for (q in port_loaded) suppressMessages(loadNamespace(q))
cat("extensions loaded:", port_loaded, "\n")
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(utils::packageVersion("frmtmb")), p,
    format(utils::packageVersion(p)), "\n")
r <- as.data.frame(test_file(f2, package = p, env = testthat::test_env(p),
                             reporter = "silent", stop_on_failure = FALSE))
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " warn=", sum(r$warning), " blocks=", nrow(r), "\n", sep = "")
fails <- r[r$failed > 0 | r$error, "test"]
if (length(fails)) cat("FAILING BLOCK:", fails, sep = "\n  ")

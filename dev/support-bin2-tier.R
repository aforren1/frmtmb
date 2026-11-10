# Lane support: one generated bin-2 file in tier mode (its recorded
# verdicts asserted), as dev/brmsport-bin2-run.R runs it, but with a
# `;`-separated library list (the lane's build before rellib-r19, so the
# extensions load) and every failure's message printed. Staged in
# dev/support-bin2-run/<package>.
#
#   Rscript dev/support-bin2-tier.R <package> <generated file>
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
libs <- strsplit(lib, ";")[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
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
file.copy(f, stage)
f2 <- normalizePath(file.path(stage, basename(f)), winslash = "/")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
ext <- unique(unlist(lapply(libs, function(l) {
  grep("^frmtmb[.]", rownames(utils::installed.packages(l)), value = TRUE)
})))
for (q in ext) suppressMessages(loadNamespace(q))
cat("extensions loaded:", ext, "\n")
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
res <- test_file(f2, package = p, env = testthat::test_env(p),
                 reporter = "silent", stop_on_failure = FALSE)
r <- as.data.frame(res)
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), "\n", sep = "")
for (t in res) {
  for (x in t$results) {
    if (inherits(x, c("expectation_failure", "expectation_error"))) {
      cat("FAIL [", t$test, "] ", substr(conditionMessage(x), 1, 300), "\n",
          sep = "")
    }
  }
}

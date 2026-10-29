# One generated bin-2 file, one R process, against a chosen library.
#
#   Rscript dev/brmsport-bin2-run.R <package> <test file> [record.tsv]
#
# As dev/brmsport-run.R, except that the generated file lives under
# dev/brmsport-bin2/ and not in the package's test directory: it is
# copied with that package's helper-*.R files into a scratch directory,
# dev/brmsport-bin2-run/<package>/, and run there, so the harness is
# the package's own and exists in one copy.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p)
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/",
                                               mustWork = FALSE))
}
tdir <- if (p == "frmtmb") "tests/testthat" else
  file.path("extensions", p, "tests/testthat")
stage <- file.path("dev/brmsport-bin2-run", p)
unlink(stage, recursive = TRUE)
dir.create(stage, recursive = TRUE)
helpers <- list.files(tdir, "^helper.*[.]R$", full.names = TRUE)
stopifnot(length(helpers) > 0L)
file.copy(helpers, stage)
# the bin-2P adapters (dev/brmsport-bin2-adapters-spec.R) load as one
# more helper, so the package's own helpers stay untouched
ad <- "dev/brmsport-bin2-adapters.R"
if (file.exists(ad)) {
  file.copy(ad, file.path(stage, "helper-brms-bin2-adapters.R"))
}
file.copy(f, stage)
f2 <- normalizePath(file.path(stage, basename(f)), winslash = "/")
if (length(a) >= 3L) {
  rec <- file.path(normalizePath(dirname(a[3]), winslash = "/"),
                   basename(a[3]))
  unlink(rec)
  Sys.setenv(FRMTMB_BRMSPORT_RECORD = rec)
}
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
suppressMessages(library(p, character.only = TRUE))
# guardrail 4 of dev/core-boundary.md: parity is measured on core plus
# every extension, so with FRMTMB_PORT_EXTENSIONS=all every frmtmb.*
# package of the library under test is loaded before the file runs
port_loaded <- character()
if (identical(Sys.getenv("FRMTMB_PORT_EXTENSIONS"), "all")) {
  port_loaded <- grep("^frmtmb[.]", rownames(utils::installed.packages(lib)),
                      value = TRUE)
  for (q in port_loaded) suppressMessages(loadNamespace(q))
}
cat("extensions loaded:", if (length(port_loaded)) port_loaded else "none",
    "\n")
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(utils::packageVersion("frmtmb")), p,
    format(utils::packageVersion(p)), "\n")
r <- as.data.frame(test_file(f2, package = p, env = testthat::test_env(p),
                             reporter = "silent", stop_on_failure = FALSE))
cat("RESULT ", p, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " blocks=", nrow(r), "\n", sep = "")
fails <- r[r$failed > 0 | r$error, "test"]
if (length(fails)) cat("FAILING BLOCK:", fails, sep = "\n  ")

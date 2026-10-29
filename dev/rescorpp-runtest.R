## Lane rescorpp: run ONE test file of one package, with counts.
##   Rscript dev/rescorpp-runtest.R <lane|base> <pkg> <file>
## <pkg> is frmtmb or an extension's name; <file> is a test file name
## or an absolute path (a lane test file run against the base build).
a <- commandArgs(TRUE)
arm <- a[1]; pkg <- a[2]; tf <- a[3]
base <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-rescorpp-lib", base)
} else base)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-rescorpp"
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev/stan-cache"))
tdir <- if (pkg == "frmtmb") file.path(wt, "tests/testthat") else
  file.path(wt, "extensions", pkg, "tests/testthat")
f <- if (file.exists(tf)) tf else file.path(tdir, tf)
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "|", pkg,
    dirname(find.package(pkg)), format(packageVersion(pkg)), "\n")
t0 <- proc.time()[["elapsed"]]
res <- as.data.frame(testthat::test_file(
  f, package = pkg, env = testthat::test_env(pkg), reporter = "check",
  stop_on_failure = FALSE))
cat("RESULT", basename(f), "pass", sum(res$passed), "fail", sum(res$failed),
    "error", sum(res$error), "skip", sum(res$skipped), "warn",
    sum(res$warning), "seconds", round(proc.time()[["elapsed"]] - t0, 1),
    "\n")

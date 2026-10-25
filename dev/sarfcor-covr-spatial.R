# Fix round 1, F13: frmtmb.spatial's suite in covr's library stack, on
# the lane library. covr::package_coverage() stops on this machine at
# its own install ("lazy loading failed", no message), for the base
# commit's frmtmb.spatial as well. The same R CMD INSTALL line run from
# a shell installs every time; run from inside an R process (covr's, or
# this script's) it failed in most attempts. So the fix round's run did
# covr's two steps by hand, as covr 3.6.5's package_coverage() writes
# them: the install from the shell,
#   R --vanilla CMD INSTALL -l <ip> --example --install-tests
#     --with-keep.source --with-keep.parse.data --no-staged-install
#     --no-multiarch extensions/frmtmb.spatial   (R_LIBS = lane stack)
# then this script with "tests": tools::testInstalledPackage() with the
# install path first and R_LIBS, R_LIBS_USER and R_LIBS_SITE all set to
# the whole stack and R_COVR = "true", the environment the
# absent-package tests must survive. Only covr's instrumentation is
# left out.
# Rscript <this> covr|install|tests
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
mode <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
pkg <- file.path(wt, "extensions/frmtmb.spatial")
ip <- file.path("C:/Users/adf44/AppData/Local/Temp/12/claude",
                "c--Users-adf44-source-r-frmtmb",
                "7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad",
                "covr-package")
if (mode != "tests") unlink(ip, recursive = TRUE)
if (mode == "covr") {
  r <- try(covr::package_coverage(pkg, quiet = FALSE, clean = FALSE,
                                  install_path = ip))
  cat("RESULT covr", if (inherits(r, "try-error")) "ERROR" else "OK", "\n")
  if (!inherits(r, "try-error")) print(r)
}
if (mode == "install") {
  dir.create(ip)
  res <- withr::with_envvar(
    list(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)),
    system2(file.path(R.home("bin"), "R.exe"),
            c("--vanilla", "CMD", "INSTALL", "-l", shQuote(ip), "--example",
              "--install-tests", "--with-keep.source",
              "--with-keep.parse.data", "--no-staged-install",
              "--no-multiarch", shQuote(pkg))))
  cat("INSTALL status", res, "\n")
}
if (mode == "tests") {
  libs <- paste(c(ip, .libPaths()), collapse = .Platform$path.sep)
  withr::with_libpaths(ip, action = "prefix", {
    withr::with_envvar(c(
      R_DEFAULT_PACKAGES = "datasets,utils,grDevices,graphics,stats,methods",
      R_LIBS = libs, R_LIBS_USER = libs, R_LIBS_SITE = libs,
      R_COVR = "true",
      R_TESTS = file.path(R.home("share"), "R", "tests-startup.R")), {
        st <- tools::testInstalledPackage(
          "frmtmb.spatial", outDir = file.path(ip, "frmtmb.spatial"),
          types = "tests", lib.loc = ip)
        cat("TESTS status", st, "\n")
      })
  })
}
f <- list.files(ip, "^testthat[.]Rout", recursive = TRUE, full.names = TRUE)
if (length(f)) {
  x <- readLines(f[1])
  cat(grep("[[] FAIL", x, value = TRUE), sep = "\n")
}

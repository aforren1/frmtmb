# The coordinator's scratchpad compatload.R, copied: test-compat.R and
# test-compat-register.R in one fresh session after preloading the
# extensions named on the command line. TESTISO_DIR names the test
# directory, TESTISO_LIB a library ahead of rellib-r9, and
# TESTISO_ATTACH (non-empty) attaches the preloads instead of only
# loading their namespaces. The last line checks that the two files
# leave frmtmb's registries as they found them.
.libPaths(c(if (nzchar(Sys.getenv("TESTISO_LIB"))) Sys.getenv("TESTISO_LIB"),
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
pre <- commandArgs(TRUE)
suppressMessages(library(frmtmb))
att <- nzchar(Sys.getenv("TESTISO_ATTACH"))
for (p in pre) {
  if (att) library(p, character.only = TRUE) else loadNamespace(p)
}
cat("attach:", att, "\n")
cat("preloaded:", if (length(pre)) pre else "none", "\n")
ns <- asNamespace("frmtmb")
reg <- function() {
  list(ns$frmtmb_compat_contrib$features, ns$frmtmb_compat_contrib$rules,
       ns$frmtmb_aterm_registry$reg, ns$frmtmb_frame_checks$fns,
       ns$frmtmb_prior_defaults$providers,
       ls(ns$frmtmb_notice_state, all.names = TRUE))
}
r0 <- reg()
dir <- Sys.getenv("TESTISO_DIR",
                  "C:/Users/adf44/source/r/frmtmb/tests/testthat")
for (f in c("test-compat.R", "test-compat-register.R")) {
  r <- testthat::test_file(
    file.path(dir, f),
    package = "frmtmb", env = testthat::test_env("frmtmb"),
    reporter = "silent")
  d <- as.data.frame(r)
  cat(f, "pass", sum(d$passed), "fail", sum(d$failed), "err",
      sum(d$error), "\n")
  for (t in r) for (e in t$results) {
    if (inherits(e, "expectation_failure") ||
          inherits(e, "expectation_error")) {
      cat("  [", t$test, "]\n   ", substr(conditionMessage(e), 1, 600),
          "\n")
    }
  }
}
cat("registries unchanged:", identical(reg(), r0), "\n")
r1 <- reg()
cat("changed:", c("features", "rules", "aterms", "frame checks",
                  "prior defaults", "notices")[
                    !mapply(identical, r1, r0)], "\n")
for (i in which(!mapply(identical, r1, r0))) {
  cat("  before", length(r0[[i]]), "after", length(r1[[i]]), "\n")
  if (i == 6L) cat("  notices after:", r1[[i]], "\n")
}
cat("aterms now:", names(ns$frmtmb_aterm_registry$reg), "\n")
cat("frmtmb namespaces:",
    grep("^frmtmb", loadedNamespaces(), value = TRUE), "\n")

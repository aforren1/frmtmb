# Reviewer mutants for lane testiso (2026-10-14). Not part of the lane.
#   Rscript dev/testiso-rev-mutants.R <mutant> <lane|base> [extension ...]
# Each mutant plants what the narrowed test-compat.R could hide, then
# runs test-compat.R from the lane worktree ("lane") or from the main
# checkout at 47cf659d ("base").
#   typo:     a contributed rule whose side is a misspelt expects =
#             name, so it can never match, whatever is loaded
#   method:   a contributed "method" feature that names no function
#   core_hmm: a CORE rule naming hmm, a name frmtmb.sample also waits
#             for, so it sits in the unresolved set with sample loaded
.libPaths(c("C:/Users/adf44/source/r/wt-testiso-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
a <- commandArgs(TRUE)
mut <- a[1]
arm <- a[2]
pre <- a[-(1:2)]
suppressMessages(library(frmtmb))
for (p in pre) loadNamespace(p)
ns <- asNamespace("frmtmb")
dir <- switch(arm,
  lane = "C:/Users/adf44/source/r/frmtmb-wt-testiso/tests/testthat",
  base = "C:/Users/adf44/source/r/frmtmb/tests/testthat")
orig_hand <- ns$compat_hand_rules_tbl
switch(mut,
  none = NULL,
  typo = frmtmb_register_compat(rules = function() {
    b <- compat_rule_builder()
    b$r("gaussian", "hmmm", "works", "mutant: misspelt expects name")
    b$rules()
  }, expects = "hmmm"),
  method = frmtmb_register_compat(features = c(zz_contrib_method = "method")),
  core_hmm = {
    fn <- function() {
      rbind(orig_hand(), data.frame(
        feature_a = "gaussian", feature_b = "hmm", status = "works",
        note = "mutant", override = FALSE, stringsAsFactors = FALSE))
    }
    environment(fn) <- ns
    utils::assignInNamespace("compat_hand_rules_tbl", fn, ns = "frmtmb")
    hc <- ns$compat_hand_cache
    hc$tbl <- NULL
  },
  stop("unknown mutant ", mut))
un <- attr(frm_compat(), "unresolved")
r <- testthat::test_file(file.path(dir, "test-compat.R"), package = "frmtmb",
                         env = testthat::test_env("frmtmb"),
                         reporter = "silent")
d <- as.data.frame(r)
cat("REVMUT", mut, "arm:", arm, "preloaded:",
    if (length(pre)) paste(pre, collapse = ",") else "none",
    "| frm_compat unresolved:", if (length(un)) un else "none",
    "| pass", sum(d$passed), "fail", sum(d$failed), "err", sum(d$error), "\n")
bad <- unique(d$test[d$failed > 0 | d$error])
for (t in bad) cat("  FAILS [", t, "]\n", sep = "")

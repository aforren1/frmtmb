# Mutation checks for the tests lane testiso changed: each mutant breaks
# the behavior one changed test guards, and the test must then fail.
#
#   Rscript dev/testiso-mutants.R <mutant> [extension ...]
#
# The named extensions are loaded first (loadNamespace), so a mutant is
# also seen to fail in the session the change was made for. "none" is
# the control and must fail nothing. TESTISO_DIR is the test directory.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
a <- commandArgs(TRUE)
mut <- a[1]
pre <- a[-1]
suppressMessages(library(frmtmb))
for (p in pre) loadNamespace(p)
ns <- asNamespace("frmtmb")
swap <- function(name, fn) {
  environment(fn) <- ns
  utils::assignInNamespace(name, fn, ns = "frmtmb")
  # the vocabulary is cached against the contributed features, which a
  # mutant of the builder does not change
  cc <- ns$frmtmb_compat_cache
  cc$tbl <- NULL
  cc$features <- NULL
}
dir <- Sys.getenv("TESTISO_DIR",
                  "C:/Users/adf44/source/r/frmtmb-wt-testiso/tests/testthat")
file <- "test-compat.R"
orig_hand <- ns$compat_hand_rules_tbl
orig_build <- ns$compat_features_build
orig_check <- ns$compat_check_rules
switch(mut,
  none = NULL,
  # a contributed rule set whose first row is shadowed by its second:
  # a dead contributed rule that waits for nothing
  contrib_dead = frmtmb_register_compat(rules = function() {
    b <- compat_rule_builder()
    b$r("gaussian", "cens()", "works", "mutant, shadowed")
    b$r("gaussian", "cens()", "refused", "mutant, wins")
    b$rules()
  }),
  # a CORE rule naming a feature nothing supplies: an unresolved side
  # that must not be excused the way a contributed one is
  core_unresolved = swap("compat_hand_rules_tbl", function() {
    rbind(orig_hand(), data.frame(
      feature_a = "gaussian", feature_b = "zz_gone", status = "works",
      note = "mutant", override = FALSE, stringsAsFactors = FALSE))
  }),
  # core loses a special from its vocabulary
  core_special = swap("compat_features_build", function(extra = NULL) {
    out <- orig_build(extra)
    out[out$name != "me()", , drop = FALSE]
  }),
  # core declares a method that does not exist
  core_method = swap("compat_features_build", function(extra = NULL) {
    out <- orig_build(extra)
    rbind(out, data.frame(name = "zz_method", key = "zz_method",
                          kind = "method", stringsAsFactors = FALSE))
  }),
  NULL)
# mutants of test-compat-register.R's guards
reg_mut <- c(
  # the resolver forgets the side it cannot resolve
  unresolved_drop = function() {
    swap("compat_unresolved_sides", function(rules, vocab) character(0))
  },
  # a registration that commits its features before its rules resolve
  commit_early = function() {
    # the mutant reads the original from the global environment, the end
    # of the namespace's own search path
    swap("compat_check_rules", function(df, vocab, expects, arg) {
      cc <- frmtmb_compat_contrib
      cc$features <- c(cc$features, list(c(tfam = "family")))
      orig_check(df, vocab, expects, arg)
    })
  },
  # frmtmb_register_aterm() stops writing the compat feature
  aterm_nofeature = function() {
    swap("compat_new_aterm_feature", function(name) invisible(NULL))
  })
known <- c("none", "contrib_dead", "core_unresolved", "core_special",
           "core_method", "none_reg", names(reg_mut))
if (!mut %in% known) stop("unknown mutant ", mut)
if (mut %in% names(reg_mut)) {
  reg_mut[[mut]]()
  file <- "test-compat-register.R"
}
if (identical(mut, "none_reg")) {
  file <- "test-compat-register.R"
}
r <- testthat::test_file(file.path(dir, file), package = "frmtmb",
                         env = testthat::test_env("frmtmb"),
                         reporter = "silent")
d <- as.data.frame(r)
cat("mutant:", mut, "preloaded:", if (length(pre)) pre else "none",
    "file:", file, "pass", sum(d$passed), "fail", sum(d$failed), "err",
    sum(d$error), "\n")
bad <- unique(d$test[d$failed > 0 | d$error])
for (t in bad) cat("  FAILS [", t, "]\n", sep = "")

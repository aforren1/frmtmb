# Reviewer, claim 8: do the lane's new tests fail BEHAVIOURALLY, not
# only on a missing symbol? Each mutant replaces one lane function in
# the namespace for the session, then runs one test file.
#   Rscript dev/aterms2-rev-09-mutants.R <mutant> <core|sample>
# Log: dev/aterms2-rev-log-09-mutants.txt (driver appends).
args <- commandArgs(TRUE)
mut <- args[[1]]
which_file <- args[[2]]
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  if (which_file == "sample") library(frmtmb.sample)
})
ns <- asNamespace("frmtmb")
put <- function(nm, f) {
  environment(f) <- ns
  assignInNamespace(nm, f, "frmtmb")
}
switch(mut,
  none = NULL,
  # rate() read by nothing after the fit: the mean is mu
  has_rate_false = put("has_rate", function(rspec) FALSE),
  # the shape is not scaled by the exposure
  rate_shape_plain = put("rate_shape", function(shape, aterms) shape),
  # the density ignores the exposure
  rate_mu_plain = put("rate_mu", function(dpars, link, aterms) dpars[["mu"]]),
  # brms's NA rule replaced by "any NA drops the row"
  na_rule_base = put("subset_na_rows", function(spec, mf, resp_cols, exempt) {
    Reduce(`|`, lapply(setdiff(names(mf), exempt), function(cn) {
      v <- mf[[cn]]
      if (is.matrix(v)) rowSums(is.na(v)) > 0 else is.na(v)
    }), rep(FALSE, nrow(mf)))
  }),
  # newdata is not cut to the subset rows
  newdata_unfiltered = put("subset_newdata",
                           function(object, resp, newdata) newdata),
  # the resp check never fires
  resp_check_off = put("subset_resp_check",
                       function(object, resp, what) invisible(NULL)),
  # idx matched by position instead of by value
  idx_by_position = put("mi_idx_rows", function(ent, vn, resp, tgt, mf,
                                                index_vals, sub_rows) {
    if (is.null(ent$idx)) return(NULL)
    rep_len(seq_along(index_vals[[vn]]), nrow(mf))
  }),
  stop("unknown mutant")
)
f <- if (which_file == "core") {
  "C:/Users/adf44/source/r/frmtmb-wt-aterms2/tests/testthat/test-subset-rate.R"
} else {
  "C:/Users/adf44/source/r/frmtmb-wt-aterms2/extensions/frmtmb.sample/tests/testthat/test-subset-rate-draws.R"
}
pkg <- if (which_file == "core") "frmtmb" else "frmtmb.sample"
setwd(dirname(f))
res <- testthat::test_file(f, package = pkg, env = testthat::test_env(pkg),
                           reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("MUTANT %-20s %-6s pass=%d fail=%d error=%d skip=%d\n", mut,
            which_file, sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped)))
bad <- df$test[df$failed > 0 | df$error]
if (length(bad)) cat("   failing blocks:", paste(bad, collapse = " | "), "\n")

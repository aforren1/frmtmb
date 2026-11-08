# Review of lane arstep: my own namespace mutants against the lane's new test
# file, plus the lane's four, rerun.
#   control    no change
#   arsign     gp_ar_left() adds phi * lag instead of subtracting it
#   nomask     gp_ar_left() drops the mask: a group's first rows read row 1
#              (or the row `safe` points at) instead of zero
#   nosettle   gp_second_descent() without obj_settle(): the objective's
#              state is left wherever the second descent ended
#   norestore  gp_start_scan() without the on.exit() restore around reading
#              the alternative start
#   noagree / noalt / noarphi / nogate   the lane's (dev/arstep-mut.R)
# Usage: Rscript review-arstep-mut.R <mutant> [test file basename]
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
a <- commandArgs(TRUE)
mut <- a[1]
tf <- if (length(a) > 1) a[2] else "test-gp-ar-start.R"
rv_libs("lane")
suppressPackageStartupMessages(library(testthat))
ns <- asNamespace("frmtmb")
edit_fn <- function(name, from, to) {
  txt <- deparse(ns[[name]])
  i <- grep(from, txt, fixed = TRUE)
  stopifnot(length(i) >= 1L)
  txt[i] <- gsub(from, to, txt[i], fixed = TRUE)
  rv_put(name, eval(parse(text = txt)))
}
switch(mut,
  control = NULL,
  arsign = edit_fn("gp_ar_left", "out <- out - phi[k]", "out <- out + phi[k]"),
  nomask = {
    orig_gad <<- ns$gp_ar_data
    rv_put("gp_ar_data", function(ar, n) {
      d <- orig_gad(ar, n)
      d$mask <- lapply(d$mask, function(m) m * 0 + 1)
      d$mmask <- lapply(d$mmask, function(m) m * 0 + 1)
      d
    })
  },
  nosettle = edit_fn("gp_second_descent", "obj_settle(obj, opt, best, saved)",
                     "invisible(NULL)"),
  norestore = edit_fn("gp_start_scan", "on.exit(obj_state_restore(obj, saved), add = TRUE)",
                      "invisible(NULL)"),
  noagree = rv_put("gp_marginal_agrees", function(...) TRUE),
  noalt = rv_put("gp_marginal_alts", function(...) list()),
  noarphi = edit_fn("gp_marginal_opt", "c(0, 0.5, 0.9, 1)", "c(0)"),
  nogate = rv_put("gp_marginal_ar", function(frame) {
    if (length(frame[["autocor"]] %||% list())) FALSE else NULL
  }),
  stop("unknown mutant"))
f <- if (file.exists(tf)) tf else file.path(RV_WT, "tests/testthat", tf)
tf <- basename(tf)
res <- as.data.frame(test_file(f, package = "frmtmb", env = test_env("frmtmb"),
                               reporter = "silent"))
cat(sprintf("RMUT %s %s pass %d fail %d error %d warn %d\n", mut, tf,
            sum(res$nb) - sum(res$failed), sum(res$failed), sum(res$error),
            sum(res$warning)))
for (i in which(res$failed > 0 | res$error)) {
  cat(sprintf("RMUTF %s %s\n", mut, res$test[i]))
}
cat("DONE\n")

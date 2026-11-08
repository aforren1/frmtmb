# Lane arstep: namespace mutants against tests/testthat/test-gp-ar-start.R.
# Each mutant replaces one function in frmtmb's namespace, then the file runs.
#   control   no change
#   noagree   gp_marginal_agrees() always TRUE (the agreement guard gone)
#   noalt     gp_marginal_alts() returns nothing (no second start)
#   noarphi   the scan's AR grid holds 0 only (no unit-root cell)
#   nogate    gp_marginal_ar() refuses every AR term (the 0.79.1 route)
# Args: arm mutant
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]; mut <- args[2]
rev_libs(arm)
suppressPackageStartupMessages(library(testthat))
ns <- asNamespace("frmtmb")
put <- function(name, f) {
  environment(f) <- ns
  unlockBinding(name, ns); assign(name, f, envir = ns); lockBinding(name, ns)
}
switch(mut,
  control = NULL,
  noagree = put("gp_marginal_agrees", function(...) TRUE),
  noalt = put("gp_marginal_alts", function(...) list()),
  noarphi = {
    txt <- paste(deparse(ns$gp_marginal_opt), collapse = "\n")
    t2 <- sub("c(0, 0.5, 0.9, 1)", "c(0)", txt, fixed = TRUE)
    stopifnot(t2 != txt)
    put("gp_marginal_opt", eval(parse(text = t2)))
  },
  nogate = put("gp_marginal_ar", function(frame) {
    if (length(frame[["autocor"]] %||% list())) FALSE else NULL
  }),
  # fix round: the second descent without obj_settle() (the reviewer's)
  nosettle = {
    txt <- paste(deparse(ns$gp_second_descent), collapse = "\n")
    t2 <- sub("obj_settle(obj, opt, best, saved)", "NULL", txt, fixed = TRUE)
    stopifnot(t2 != txt)
    put("gp_second_descent", eval(parse(text = t2)))
  },
  # fix round: a constant sigma read as no noise again (the base's reading)
  nosigc = {
    txt <- paste(deparse(ns$gp_marginal_setup), collapse = "\n")
    t2 <- sub("if (!is.null(g$sigc))", "if (FALSE)", txt, fixed = TRUE)
    stopifnot(t2 != txt)
    put("gp_marginal_setup", eval(parse(text = t2)))
  },
  # fix round: the AR marginal on every row again (the first build's B1)
  rows = {
    txt <- paste(deparse(ns$gp_marginal_data), collapse = "\n")
    t2 <- sub("if (all(se2 == se2[1L]) && all(w == w[1L]))",
              "if (is.null(ar) && all(se2 == se2[1L]) && all(w == w[1L]))", txt, fixed = TRUE)
    stopifnot(t2 != txt)
    put("gp_marginal_data", eval(parse(text = t2)))
  },
  stop("unknown mutant"))
res <- as.data.frame(test_file("C:/Users/adf44/source/r/frmtmb-wt-arstep/tests/testthat/test-gp-ar-start.R",
                               package = "frmtmb", env = test_env("frmtmb"), reporter = "silent"))
cat(sprintf("MUT %s %s pass %d fail %d error %d\n", arm, mut, sum(res$nb) - sum(res$failed), sum(res$failed), sum(res$error)))
for (i in which(res$failed > 0 | res$error)) cat(sprintf("MUTF %s %s %s\n", arm, mut, res$test[i]))
cat("DONE\n")

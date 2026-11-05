# Does test-gp-nugget-floor.R kill mutants of the warning? Each mutant rewrites
# one function in the installed namespace, by text, before the file runs.
#
#   Rscript dev/gpgap-mut.R <mutant>      # one mutant, one process
#
# Prints MUT <name> pass= fail= error= and KILLED or SURVIVED.
m <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
ns <- asNamespace("frmtmb")
sub_fn <- function(fname, old, new) {
  f <- get(fname, envir = ns)
  src <- paste(deparse(f, width.cutoff = 200L), collapse = "\n")
  stopifnot(grepl(old, src, fixed = TRUE))
  g <- eval(parse(text = sub(old, new, src, fixed = TRUE)), envir = ns)
  unlockBinding(fname, ns)
  assign(fname, g, envir = ns)
  lockBinding(fname, ns)
}
muts <- list(
  control = function() NULL,
  class_drop = function() sub_fn("fit_assembled",
    "class = \"frmtmb_gp_nugget\"", "class = NULL"),
  no_warning = function() sub_fn("fit_assembled",
    "gp_msg <- gp_nugget_message(fit)", "gp_msg <- NULL"),
  thr_half = function() sub_fn("gp_nugget_message", "sqrt(2/length(info))",
                               "0.5 * sqrt(2/length(info))"),
  thr_double = function() sub_fn("gp_nugget_message", "sqrt(2/length(info))",
                                 "2 * sqrt(2/length(info))"),
  thr_x5 = function() sub_fn("gp_nugget_message", "sqrt(2/length(info))",
                             "5 * sqrt(2/length(info))"),
  thr_x30 = function() sub_fn("gp_nugget_message", "sqrt(2/length(info))",
                              "30 * sqrt(2/length(info))"),
  thr_fixed = function() sub_fn("gp_nugget_message", "sqrt(2/length(info))",
                                "0.01"),
  q_max = function() sub_fn("gp_nugget_message",
    "stats::quantile(nv * info, 0.9, names = FALSE)", "max(nv * info)"),
  q_med = function() sub_fn("gp_nugget_message",
    "stats::quantile(nv * info, 0.9, names = FALSE)", "stats::median(nv * info)"),
  nug_x10 = function() sub_fn("gp_nugget_message",
    "gp_block_nugget(bk) * exp(2 * th[[1L]])",
    "10 * gp_block_nugget(bk) * exp(2 * th[[1L]])"),
  sd1 = function() sub_fn("gp_nugget_message",
    "gp_block_nugget(bk) * exp(2 * th[[1L]])", "gp_block_nugget(bk)"),
  no_trycatch = function() sub_fn("gp_nugget_message",
    "error = function(e) NULL", "error = function(e) stop(e)"),
  no_prior = function() sub_fn("gp_position_info",
    "pmax(H - diag(chol2inv(chol(V))), 0)", "pmax(H, 0)"),
  first_only = function() sub_fn("gp_nugget_message",
    "for (bk in bks) {", "for (bk in bks[1L]) {"),
  one_block_loop = function() sub_fn("gp_nugget_message",
    "for (bk in bks) {", "for (bk in rev(bks)[1L]) {"),
  msg_no_moves = function() sub_fn("gp_nugget_message", "reads low. Where positions repeat",
                                   "moves. Where positions repeat"),
  msg_k_promise = function() sub_fn("gp_nugget_message",
    "but k must be large enough for the ", "and it fits the ")
)
muts[[m]]()
f <- "tests/testthat/test-gp-nugget-floor.R"
Sys.setenv(NOT_CRAN = "true")
res <- suppressWarnings(as.data.frame(
  test_file(f, package = "frmtmb", env = test_env("frmtmb"), reporter = "silent",
            stop_on_failure = FALSE)))
nf <- sum(res$failed) + sum(as.integer(res$error))
cat(sprintf("MUT %-14s pass=%d fail=%d error=%d warn=%d %s\n", m, sum(res$passed),
            sum(res$failed), sum(as.integer(res$error)), sum(res$warning),
            if (nf) "KILLED" else "SURVIVED"))

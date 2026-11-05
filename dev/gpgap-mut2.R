# Does test-reparam.R kill mutants of the sampler's factor order and its start?
# Each rewrites one function in the installed namespace, by text.
#   Rscript dev/gpgap-mut2.R <mutant>
m <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.sample)
})
sub_fn <- function(pkg, fname, old, new) {
  ns <- asNamespace(pkg)
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
  no_order = function() sub_fn("frmtmb.sample", "gp_sample_nugget",
    "gp_farthest_first(bk)", "NULL"),
  no_nugget = function() sub_fn("frmtmb.sample", "gp_sample_nugget",
    "<- gp_sampler_nugget", "<- NULL"),
  first_last = function() sub_fn("frmtmb.sample", "gp_farthest_first",
    "which.max(gap)", "which.min(gap)"),
  start_medoid = function() sub_fn("frmtmb.sample", "gp_farthest_first",
    "which.min(rowSums(d2))", "1L"),
  one_dim = function() sub_fn("frmtmb.sample", "gp_farthest_first",
    "Reduce(`+`, blk[[\"aux_D2\"]])", "blk[[\"aux_D2\"]][[1L]]"),
  no_start_fix = function() sub_fn("frmtmb.sample", "ncp_start_pars",
    "bk[[\"gp_nugget\"]] <- NULL", "NULL"),
  no_random = function() sub_fn("frmtmb.sample", "frm_sample",
    "moved && !init_asked && !from_formula", "FALSE")
)
muts[[m]]()
f <- "extensions/frmtmb.sample/tests/testthat/test-reparam.R"
res <- suppressMessages(suppressWarnings(as.data.frame(
  test_file(f, package = "frmtmb.sample", env = test_env("frmtmb.sample"),
            reporter = "silent", stop_on_failure = FALSE))))
nf <- sum(res$failed) + sum(as.integer(res$error))
cat(sprintf("MUT2 %-14s pass=%d fail=%d error=%d %s\n", m, sum(res$passed),
            sum(res$failed), sum(as.integer(res$error)),
            if (nf) "KILLED" else "SURVIVED"))

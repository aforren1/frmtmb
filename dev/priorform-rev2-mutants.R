# Reviewer recheck round 1: mutants of the round-1 guards, one per
# process, run against the lane test files that should catch them.
#   Rscript dev/priorform-rev2-mutants.R <mutant>
args <- commandArgs(trailingOnly = TRUE); mut <- args[1]
.libPaths(c("C:/Users/adf44/source/r/priorform-lib", "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
root <- "C:/Users/adf44/source/r/frmtmb-wt-priorform"
suppressMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
put <- function(nm, f) { environment(f) <- ns; utils::assignInNamespace(nm, f, ns = "frmtmb") }
edit_fn <- function(nm, from, to) {
  txt <- deparse(get(nm, ns))
  hit <- grepl(from, txt, fixed = TRUE)
  if (sum(hit) != 1L) stop("mutant anchor matched ", sum(hit), " lines in ", nm)
  txt[hit] <- sub(from, to, txt[hit], fixed = TRUE)
  put(nm, eval(parse(text = txt)))
}
switch(mut,
  none = NULL,
  order_identity = put("prior_specificity_order", function(pl) unclass(pl)),
  order_reversed = edit_fn("prior_specificity_order", "order(rank[at])", "order(-rank[at])"),
  order_nogroup = edit_fn("prior_specificity_order", "10L * nz(s, \"group\")", "0L * nz(s, \"group\")"),
  order_nocoef = edit_fn("prior_specificity_order", "100L * nz(s, \"coef\")", "0L * nz(s, \"coef\")"),
  slots_noop = put("check_prior_slots", function(prior) invisible(prior)),
  slots_parens = edit_fn("prior_slot_label", "coef = par_name_bare(s$coef %||% ", "coef = identity(s$coef %||% "),
  slots_density_only = edit_fn("check_prior_slots", "keys <- vapply(unclass(prior), prior_slot_label, \"\")",
    "keys <- vapply(unclass(prior), function(s) if (is.null(s$dist)) paste(runif(1)) else prior_slot_label(s), \"\")"),
  stop("unknown mutant"))
files <- "^(brms-formula-priors|priors-bounds-grcov|setprior|lkj|prior-compat|get-prior-route)$"
setwd(file.path(root, "tests", "testthat"))
res <- test_dir(".", filter = files, package = "frmtmb", reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(res)
cat(sprintf("MUTANT %s pass %d fail %d error %d\n", mut, sum(d$passed), sum(d$failed), sum(d$error)))
bad <- d[d$failed > 0 | d$error, ]
for (i in seq_len(nrow(bad))) cat("   ", bad$file[i], ":", bad$test[i], "\n")

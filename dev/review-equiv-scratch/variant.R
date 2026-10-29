# Measure a VARIANT of one pair's spec, for a control or a diagnosis:
#
#   Rscript dev/equiv-variant.R "<feature_a> x <feature_b>" <edit.R>
#
# <edit.R> holds R lines that modify `spec` (the pair's spec list) before
# it is measured, e.g. `spec$bf <- "bf(y ~ t2(x, z))"; spec$l4 <- TRUE`.
# The edit is read from a file because Rscript on Windows hands its
# arguments to a command line that reads `|` as a pipe. Prints the
# record's verdict and level numbers; writes nothing.
.libPaths(c(strsplit(Sys.getenv("FRMTMB_EQUIV_LIB",
                                "C:/Users/adf44/source/r/rellib-r15"),
                     ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/review-equiv-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
args <- commandArgs(trailingOnly = TRUE)
base <- E$eq_specs[[args[1]]]
if (is.null(base)) stop("no spec for ", args[1])
edit <- paste(readLines(args[2]), collapse = "\n")
mod <- eval(parse(text = paste0("function(spec) {\n", edit, "\nspec\n}")))
environment(mod) <- E
E$eq_specs[["variant"]] <- function() mod(base())
rec <- E$eq_measure("variant", list(status = "untested"))
keep <- c("verdict", "reason", "frm_error", "brms_error", "l1_const",
          "l1_grad", "l1_extra_const", "l1_error", "l2", "l4", "l4_error")
str(rec[intersect(names(rec), keep)], max.level = 3, nchar.max = 300)

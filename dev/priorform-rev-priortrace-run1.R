# Reviewer, lane wt-priorform: blast radius of refusing duplicate prior
# specifications as brms 2.23.0 does. Runs ONE test file (or one Rd
# example set) on the lane build with resolve_priorlist() traced, and
# appends one line per priorlist holding two specifications for the
# same brms slot (class, coef, group, resp, dpar, nlpar).
#   Rscript dev/priorform-rev-priortrace-run1.R <pkgdir> <filter>
args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]; filt <- args[2]
.libPaths(c("C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
root <- "C:/Users/adf44/source/r/frmtmb-wt-priorform"
logf <- file.path(root, "dev", "priorform-rev-priortrace.tsv")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
assign(".rev_label", paste(basename(pkgdir), filt), envir = globalenv())
assign(".rev_logf", logf, envir = globalenv())
tracer <- quote({
  specs <- unclass(pl)
  if (length(specs) > 1L) {
    key <- vapply(specs, function(s) {
      sp <- frmtmb:::spec_spelling(s)
      paste(sp$class, s$coef %||% "", s$group %||% "", s$resp %||% "",
            sp$dpar, s$nlpar %||% "", sep = "|")
    }, "")
    dup <- unique(key[duplicated(key)])
    for (k in dup) {
      ss <- specs[key == k]
      kinds <- vapply(ss, function(s) if (is.null(s$dist)) "bounds" else
        "density", "")
      calls <- vapply(sys.calls(), function(cl) {
        f <- cl[[1L]]
        if (is.name(f)) as.character(f) else if (is.call(f))
          paste(deparse(f), collapse = "") else "?"
      }, "")
      via <- intersect(c("prior_stack", "frm_sample", "validate_prior",
                         "frm", "frm_simulate", "resolve_prior_input",
                         "prior_summary", "update"),
                       sub("^.*::+", "", calls))
      tn <- tryCatch(testthat::get_reporter()$.context %||% "", error = function(e) "")
      cat(get(".rev_label", globalenv()), "\t", k, "\t",
          paste(kinds, collapse = ","), "\t", paste(via, collapse = ","),
          "\t", paste(vapply(ss, function(s) s$prior %||% "", ""),
                      collapse = " ; "), "\n",
          file = get(".rev_logf", globalenv()), append = TRUE, sep = "")
    }
  }
})
suppressMessages(trace("resolve_priorlist", tracer = tracer, where = ns,
                       print = FALSE))
this_pkg <- if (identical(pkgdir, ".")) "frmtmb" else basename(pkgdir)
if (!identical(pkgdir, ".")) {
  suppressMessages(library(this_pkg, character.only = TRUE))
}
setwd(file.path(root, pkgdir, "tests", "testthat"))
res <- testthat::test_dir(".", filter = filt, package = this_pkg,
                          reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(res)
cat(sprintf("PRIORTRACE %s pass %d fail %d error %d skip %d\n", filt,
            sum(d$passed), sum(d$failed), sum(d$error), sum(d$skipped)))

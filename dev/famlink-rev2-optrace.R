## Recheck (rounds 1 and 1b), priority 1: does nlminb_trial_fn() leave
## every optimizer invocation bitwise unchanged, across real test-suite
## fits in core and the extensions?
##
## run_optimizer() is traced on EXIT in the chosen arm, so every
## optimizer call a test file makes (first runs, restarts, recovery
## from best, autoscale pre-fits, importance rounds, frm_allfit()'s
## optimizers) appends one record: optimizer, par, objective,
## convergence, message, iterations, evaluations (function and gradient
## counts, the path signature) and, on the lane, nonfinite_trials. The
## testthat result per file (pass, fail, error, skip, WARNING count) is
## stored too. dev/famlink-rev2-optrace-cmp.R pairs the arms record by
## record.
##
## Usage: Rscript dev/famlink-rev2-optrace.R <lane|base> <pkgdir> <filter>
args <- commandArgs(trailingOnly = TRUE)
ARM <- args[1]; pkgdir <- args[2]; filt <- args[3]
source("dev/famlink-rev-common.R")
suppressMessages(library(testthat))
root <- normalizePath(".")
outdir <- file.path(root, "dev", "famlink-rev2-optrace")
dir.create(outdir, showWarnings = FALSE)
rec <- new.env(); rec$x <- list()
record <- function(res, optimizer) {
  rec$x[[length(rec$x) + 1L]] <- list(
    optimizer = if (is.function(optimizer)) "custom" else optimizer,
    par = res$par, objective = res$objective,
    convergence = res$convergence, message = res$message,
    iterations = res$iterations, evaluations = res$evaluations,
    nonfinite_trials = res$nonfinite_trials)
}
assign("record", record, envir = globalenv())
suppressMessages(trace("run_optimizer", where = asNamespace("frmtmb"),
  exit = quote(tryCatch(record(returnValue(), optimizer),
                        error = function(e) NULL)),
  print = FALSE))
if (!identical(pkgdir, ".")) {
  suppressMessages(library(basename(pkgdir), character.only = TRUE))
}
this_pkg <- if (identical(pkgdir, ".")) "frmtmb" else basename(pkgdir)
setwd(file.path(pkgdir, "tests", "testthat"))
res <- testthat::test_dir(".", filter = filt, package = this_pkg,
                          reporter = "silent", stop_on_failure = FALSE)
d <- as.data.frame(res)
summ <- c(pass = sum(d$passed), fail = sum(d$failed), error = sum(d$error),
          skip = sum(d$skipped), warning = sum(d$warning))
tag <- gsub("[^A-Za-z0-9]+", "", paste0(basename(pkgdir), "-", filt))
saveRDS(list(summary = summ, records = rec$x, tests = d[, c("file", "test",
          "nb", "failed", "skipped", "error", "warning", "passed")]),
        file.path(outdir, paste0(ARM, "-", tag, ".rds")))
cat(sprintf("OPTRACE %s %s %s optimizer calls %d | pass %d fail %d error %d skip %d warning %d\n",
            ARM, pkgdir, filt, length(rec$x), summ[["pass"]], summ[["fail"]],
            summ[["error"]], summ[["skip"]], summ[["warning"]]))

# Lane illse: dev/ciharden-run1.R plus a convergence-warning probe.
# One test file, one R process, for the fragility scan
# (dev/ciharden-scan.sh). The same runner as dev/release/run-tests.R,
# with two differences the scan needs: the library stack is an argument,
# so a lane build and the base build run through one script; and every
# non-passing expectation is printed with its source line, so that two
# configurations can be compared expectation by expectation rather than
# by counts that can cancel.
#
# Usage: Rscript dev/ciharden-run1.R <lib dir or "base"> <package> <file>
local({
# ROUND_BASE lets a later round rerun the scan against its own base
base <- Sys.getenv("ROUND_BASE", "C:/Users/adf44/source/r/rellib-r8")
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1], base, user))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
suppressMessages(library(testthat))
p <- a[2]
f <- a[3]
suppressMessages(library(p, character.only = TRUE))
cat("lib: ", p, " ", as.character(utils::packageVersion(p)), " ",
    dirname(find.package(p)), "; frmtmb ",
    as.character(utils::packageVersion("frmtmb")), " ",
    dirname(find.package("frmtmb")), "\n", sep = "")
# the BLAS actually loaded: a 600 x 600 product is about 0.1 s with the
# reference BLAS and under 0.01 s with OpenBLAS
cat("blas: ", R.home(), " threads=", Sys.getenv("OPENBLAS_NUM_THREADS"),
    " gemm600=", system.time({
      m <- matrix(1, 600, 600); m %*% m
    })[["elapsed"]], "\n", sep = "")
# Lane illse: with ILLSE_CONVSCAN set to a directory, every fit whose
# check_convergence() warns gets one row there: the warning kinds, the
# Newton headroom from the stop (grad_headroom(), over the parameters
# no bound holds) and what a restart of nlminb from the stop gains.
# The objective's state is put back after the probes.
cdir <- Sys.getenv("ILLSE_CONVSCAN")
if (nzchar(cdir)) {
  ns <- asNamespace("frmtmb")
  dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(normalizePath(cdir), paste0(p, "--", basename(f), ".tsv"))
  conv_hook <- function(fit, rv) {
    tryCatch({
      w <- rv$warnings
      if (!length(w)) return(invisible())
      k <- c(code = "did not report convergence",
             grad = "Large maximum absolute gradient",
             nonfin = "gradient at the reported optimum is not finite",
             pd = "Hessian is not positive definite",
             imp = "importance correction")
      kinds <- names(k)[vapply(k, function(s) any(grepl(s, w, fixed = TRUE)),
                               NA)]
      obj <- fit$obj
      pp <- fit$opt$par
      saved <- ns$obj_state_save(obj)
      on.exit(ns$obj_state_restore(obj, saved), add = TRUE)
      g <- suppressWarnings(drop(obj$gr(pp)))
      box <- ns$fit_outer_box(fit)
      act <- ns$grad_bound_active(pp, g, box) %in% TRUE
      gf <- g
      gf[act] <- 0
      head <- suppressWarnings(tryCatch(ns$grad_headroom(fit, gf, act),
                                        error = function(e) NA_real_))
      rs <- if (length(pp) <= 200L) suppressWarnings(tryCatch(
        stats::nlminb(pp, obj$fn, obj$gr, lower = box$lower,
                      upper = box$upper,
                      control = list(iter.max = 300, eval.max = 400)),
        error = function(e) NULL))
      gain <- if (is.null(rs)) NA_real_ else fit$opt$objective - rs$objective
      cat(paste(c(basename(f), paste(kinds, collapse = ","),
                  fit$opt$convergence, length(pp),
                  length(obj$env$random), sum(act),
                  signif(max(abs(g)), 4), signif(head, 4), signif(gain, 4),
                  if (is.null(rs)) NA else rs$convergence,
                  gsub("[\t\n]", " ", substr(fit$opt$message, 1, 60)),
                  substr(gsub("[\t\n ]+", " ",
                              paste(deparse(fit$call), collapse = " ")),
                         1, 120)),
                collapse = "\t"), "\n", file = out, append = TRUE, sep = "")
    }, error = function(e) {
      cat("HOOKERROR\t", basename(f), "\t", conditionMessage(e), "\n",
          file = out, append = TRUE, sep = "")
    })
    invisible()
  }
  suppressMessages(trace("check_convergence", where = ns, print = FALSE,
                         exit = bquote(.(conv_hook)(fit, returnValue()))))
}
res <- NULL
r <- tryCatch({
  res <- test_file(f, package = p, env = testthat::test_env(p),
                   reporter = "silent")
  as.data.frame(res)
}, error = function(e) {
  cat("RESULT ", basename(f), " LOADERROR ", conditionMessage(e), "\n",
      sep = "")
  NULL
})
if (!is.null(res)) {
  for (t in res) {
    for (x in t$results) {
      k <- if (inherits(x, "expectation_failure")) "FAIL" else
        if (inherits(x, "expectation_error")) "ERROR" else
          if (inherits(x, "expectation_skip")) "SKIP" else
            if (inherits(x, "expectation_warning")) "WARN" else NA
      if (!is.na(k)) {
        ln <- if (!is.null(x$srcref)) {
          paste0(basename(utils::getSrcFilename(x$srcref)), ":",
                 x$srcref[1])
        } else "?"
        cat("DETAIL\t", k, "\t", ln, "\t", t$test, "\t",
            substr(gsub("[\r\n\t]+", " | ", conditionMessage(x)), 1, 400),
            "\n", sep = "")
      }
    }
    # a test's pass count can move with the platform when it loops over
    # what a fit returns, so the scan compares it too
    cat("TEST\t", t$test, "\t",
        sum(vapply(t$results, inherits, logical(1),
                   "expectation_success")), "\n", sep = "")
  }
}
if (!is.null(r)) {
  cat("RESULT ", basename(f), " pass=", sum(r$passed), " fail=",
      sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
      " warn=", sum(r$warning), "\n", sep = "")
}
})

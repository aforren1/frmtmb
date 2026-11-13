# Short check of lane quadseam's fix round, item 2: when does the new
# fallback of grad_headroom() / conv_headroom() (a non-finite AD Hessian
# replaced by differences) fire? One test file per process; obj_he_cached()
# is wrapped in this process only and logs each non-finite answer with the
# tape's random-effect count, whether the fit is a quadrature tape (no
# random effects but an integrate spec), and which fit-end callers are on
# the stack.
#   Rscript quadseam-short-hetrace-ext.R <base|lane> <test file> <log> <package>
# (quadseam-short-hetrace.R for an extension package's test file)
a <- commandArgs(TRUE)
arm <- a[1]; tf <- a[2]; logf <- a[3]; pkg <- a[4]
base <- "C:/Users/adf44/source/r/rellib-r20"
lane <- "C:/Users/adf44/source/r/wt-quadseam-lib"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (arm == "lane") c(lane, base, ul) else c(base, ul))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages({library(testthat); library(frmtmb); library(pkg, character.only = TRUE)})
ns <- asNamespace("frmtmb")
orig <- get("obj_he_cached", ns)
watch <- c("grad_headroom", "conv_headroom", "se_check", "cov_from_hessian",
           "newton_polish", "fit_outer_hessian", "grad_verdict",
           "fit_end_checks", "quad_fit", "quad_settle")
wrapped <- function(obj, p) {
  H <- orig(obj, p)
  if (is.numeric(H) && !all(is.finite(H))) {
    nm <- vapply(sys.calls(), function(cl) {
      f <- cl[[1L]]
      if (is.name(f)) as.character(f) else paste(deparse(f), collapse = "")
    }, "")
    who <- intersect(watch, nm)
    quad <- NA
    for (fr in rev(sys.frames())) {
      if (exists("fit", envir = fr, inherits = FALSE)) {
        ft <- get("fit", envir = fr, inherits = FALSE)
        if (is.list(ft) && !is.null(ft$quadrature)) {
          quad <- isTRUE(ft$quadrature)
          break
        }
      }
    }
    if (is.na(quad)) quad <- any(grepl("quad", nm, fixed = TRUE))
    cat(sprintf("NONFINITE_HE file %s quad %s random %d np %d nonfinite %d callers %s\n",
                basename(tf), quad, length(obj$env$random), length(p),
                sum(!is.finite(H)), paste(who, collapse = ",")),
        file = logf, append = TRUE)
  }
  H
}
environment(wrapped) <- environment()
utils::assignInNamespace("obj_he_cached", wrapped, ns = "frmtmb")
cat(sprintf("START %s %s %s\n", arm, basename(tf), dirname(find.package("frmtmb"))),
    file = logf, append = TRUE)
r <- testthat::test_file(tf, package = pkg,
                         env = testthat::test_env(pkg),
                         reporter = "silent")
df <- as.data.frame(r)
cat(sprintf("RESULT %s pass %d fail %d err %d skip %d\n", basename(tf),
            sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped)),
    file = logf, append = TRUE)

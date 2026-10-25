## Records the names of every fit one test file makes, on one arm.
##
##   Rscript dev/brmsnames-1024-suitenames.R base|lane <pkg> <test file>
##
## frm() is traced in frmtmb's namespace, so every fit is caught, the
## ones frmtmb makes internally (update(), refits) included. For each
## fit it records, in call order: the internal coefficient names,
## brms's names of the coefficients (brms_coef_table()), variables(),
## the outer parameter labels (the rows of confint() and
## vcov(full = TRUE)), and the coef column of the prior table. None of
## these computes standard errors, so the hook does not change what a
## test sees. Output: dev/brmsnames-1024-suite/<arm>/<pkg>/<file>.rds
## and a RESULT line with the test counts.
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1L]
pkg <- args[2L]
f <- normalizePath(args[3L], winslash = "/")
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
Sys.setenv(NOT_CRAN = "true")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-brmsnames"
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev/stan-cache"))
}
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(pkg, character.only = TRUE)
})
ns <- asNamespace("frmtmb")
.bn_log <- new.env()
.bn_log$fits <- list()
.bn_rec <- function(rv) {
  if (!inherits(rv, "frmtmb_fit")) return(invisible())
  g <- function(e) {
    tryCatch(suppressWarnings(suppressMessages(e)),
             error = function(err) paste("ERROR:", conditionMessage(err)))
  }
  tab <- g(ns$brms_coef_table(rv))
  rec <- list(
    call = g(paste(deparse(rv$call, width.cutoff = 500L), collapse = " ")),
    internal = if (is.data.frame(tab)) tab$internal else tab,
    brms = if (is.data.frame(tab)) tab$brms else tab,
    variables = g(variables(rv)),
    outer = g(ns$outer_par_labels(rv)),
    prior = g(ns$prior_table(rv$spec, rv$frame, "fit")$coef)
  )
  .bn_log$fits[[length(.bn_log$fits) + 1L]] <- rec
  invisible()
}
assign(".bn_rec", .bn_rec, envir = globalenv())
suppressMessages(trace("frm", where = ns, print = FALSE,
                       exit = quote(.bn_rec(returnValue()))))
# frmtmb.sample's tests reach frm() through the package's own import,
# which the trace above does not reach, so there the draws are recorded
# from frm_sample(): their variables(), beside the names of the fit
if (identical(pkg, "frmtmb.sample")) {
  .bn_rec_draws <- function(rv) {
    if (!inherits(rv, "frmtmb_draws")) return(invisible())
    g <- function(e) {
      tryCatch(suppressWarnings(suppressMessages(e)),
               error = function(err) paste("ERROR:", conditionMessage(err)))
    }
    fit <- g(rv$fit)
    tab <- if (inherits(fit, "frmtmb_fit")) g(ns$brms_coef_table(fit))
    rec <- list(
      call = g(paste(deparse(fit$call, width.cutoff = 500L), collapse = " ")),
      internal = if (is.data.frame(tab)) tab$internal,
      brms = if (is.data.frame(tab)) tab$brms,
      variables = g(posterior::variables(rv)),
      outer = if (inherits(fit, "frmtmb_fit")) g(ns$outer_par_labels(fit)),
      prior = if (inherits(fit, "frmtmb_fit")) {
        g(ns$prior_table(fit$spec, fit$frame, "fit")$coef)
      }
    )
    .bn_log$fits[[length(.bn_log$fits) + 1L]] <- rec
    invisible()
  }
  assign(".bn_rec_draws", .bn_rec_draws, envir = globalenv())
  suppressMessages(trace("frm_sample", where = asNamespace("frmtmb.sample"),
                         print = FALSE,
                         exit = quote(.bn_rec_draws(returnValue()))))
}
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
res <- tryCatch(
  as.data.frame(testthat::test_file(f, package = pkg,
                                    env = testthat::test_env(pkg),
                                    reporter = "silent")),
  error = function(e) e)
cnt <- if (inherits(res, "error")) {
  paste("ABORT", conditionMessage(res))
} else {
  sprintf("pass=%d fail=%d error=%d skip=%d", sum(res$passed),
          sum(res$failed), sum(res$error), sum(res$skipped))
}
od <- file.path(wt, "dev/review-brmsnames-suite", arm, pkg)
dir.create(od, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(file = basename(f), counts = cnt, fits = .bn_log$fits),
        file.path(od, sub("[.]R$", ".rds", basename(f))))
cat("RESULT", basename(f), cnt, "fits=", length(.bn_log$fits), "\n")

# Lane boundfam: one test file through dev/release/run-tests.R, with every
# standard-error analysis and every covariance-interval table recorded.
# As dev/boundfam-run-traced.R, with the limit rule's findings in full.
#
#   Rscript dev/boundfam-run-traced2.R <package> <test file>
#
# FRMTMB_RELLIB picks the build, as for run-tests.R (empty: rellib-r20).
# BOUNDFAM_TRACE names the record, one line per call:
#   COV <k> np=<n> lost=<name:reason,...> se=<sqrt(diag(shown)), %.10g>
#   VCT <k> se_t=<...>
#   LIM <k> <name:kind,...> est= dir= gain= fam= n=   (lane build only)
# A LIM line belongs to the COV line that follows it. The trace is
# attached after frmtmb loads and before the test file runs; it changes
# no value it records.
local({
  LIB <- "C:/Users/adf44/source/r/rellib-r20"
  .libPaths(c(if (nzchar(Sys.getenv("FRMTMB_RELLIB"))) {
    Sys.getenv("FRMTMB_RELLIB")
  }, LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressMessages(library(frmtmb))
  out <- Sys.getenv("BOUNDFAM_TRACE")
  stopifnot(nzchar(out))
  cat("TRACE lib", dirname(find.package("frmtmb")), "\n", file = out)
  ns <- asNamespace("frmtmb")
  rec <- new.env()
  rec$k <- 0L
  assign(".boundfam_rec", rec, envir = globalenv())
  fmt <- function(x) paste(sprintf("%.10g", x), collapse = ",")
  rec$cov <- function(v) {
    rec$k <- rec$k + 1L
    lost <- v$lost
    sh <- v$shown
    se <- if (is.matrix(sh)) suppressWarnings(sqrt(diag(sh))) else NA
    cat("COV ", rec$k, " np=", NROW(sh), " lost=",
        paste(names(lost), lost, sep = ":", collapse = ","),
        " se=", fmt(se), "\n", sep = "", file = out, append = TRUE)
  }
  rec$vct <- function(v) {
    rec$k <- rec$k + 1L
    if (is.data.frame(v)) {
      cat("VCT ", rec$k, " se_t=", fmt(v$se_t), "\n", sep = "", file = out,
          append = TRUE)
    }
  }
  rec$lim <- function(v, fit) {
    if (!any(v)) return(invisible(NULL))
    li <- ns$se_limit_info(fit)
    fam <- tryCatch(fit$spec$responses[[1L]]$family$family,
                    error = function(e) "?")
    cat("LIM ", rec$k, " ",
        paste(names(li$kind), li$kind, sep = ":", collapse = ","),
        " est=", fmt(li$est), " dir=", fmt(li$dir), " gain=", fmt(li$gain),
        " fam=", fam, " n=", fit$frame[["n_obs"]], "\n", sep = "",
        file = out, append = TRUE)
  }
  suppressMessages({
    trace("cov_from_hessian", where = ns, print = FALSE,
          exit = quote(get(".boundfam_rec", globalenv())$cov(returnValue())))
    trace("varcorr_trans_rows", where = ns, print = FALSE,
          exit = quote(get(".boundfam_rec", globalenv())$vct(returnValue())))
    if (exists("se_edge_limit", ns, inherits = FALSE)) {
      trace("se_edge_limit", where = ns, print = FALSE,
            exit = quote(get(".boundfam_rec", globalenv())$lim(returnValue(),
                                                                fit)))
    }
  })
})
source("C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/release/run-tests.R")
rm(".boundfam_rec", envir = globalenv())

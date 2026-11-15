# Lane selimit2: one test file through dev/release/run-tests.R, with every
# standard-error analysis and every covariance-interval table recorded.
# As dev/boundfam-run-traced2.R, on rellib-r21, with the joint and
# one-sided separation fields of the limit rule.
#
#   Rscript dev/selimit2-run-traced.R <package> <test file>
#
# FRMTMB_RELLIB picks the build, as for run-tests.R (empty: rellib-r21).
# SELIMIT2_TRACE names the record, one line per call:
#   COV <k> np=<n> lost=<name:reason,...> se=<sqrt(diag(shown)), %.10g>
#   VCT <k> se_t=<...>
#   LIM <k> <name:kind,...> est= dir= gain= joint= sep= sep1= fam= n=
# A LIM line belongs to the COV line that follows it. The trace is
# attached after frmtmb loads and before the test file runs; it changes
# no value it records.
local({
  LIB <- "C:/Users/adf44/source/r/rellib-r22"
  .libPaths(c(if (nzchar(Sys.getenv("FRMTMB_RELLIB"))) {
    Sys.getenv("FRMTMB_RELLIB")
  }, LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressMessages(library(frmtmb))
  out <- Sys.getenv("SELIMIT2_TRACE")
  stopifnot(nzchar(out))
  cat("TRACE lib", dirname(find.package("frmtmb")), "\n", file = out)
  ns <- asNamespace("frmtmb")
  rec <- new.env()
  rec$k <- 0L
  assign(".selimit2_rec", rec, envir = globalenv())
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
    li <- ns$se_limit_info(fit)
    if (!length(li$kind)) return(invisible(NULL))
    fam <- tryCatch(fit$spec$responses[[1L]]$family$family,
                    error = function(e) "?")
    lg <- function(x) paste(as.integer(x %in% TRUE), collapse = ",")
    cat("LIM ", rec$k, " ",
        paste(names(li$kind), li$kind, sep = ":", collapse = ","),
        " est=", fmt(li$est), " dir=", fmt(li$dir), " gain=", fmt(li$gain),
        " joint=", lg(li$joint), " sep=", lg(li$sep), " sep1=", lg(li$sep1),
        " fam=", fam, " n=", fit$frame[["n_obs"]], "\n", sep = "",
        file = out, append = TRUE)
  }
  suppressMessages({
    trace("cov_from_hessian", where = ns, print = FALSE,
          exit = quote(get(".selimit2_rec", globalenv())$cov(returnValue())))
    trace("varcorr_trans_rows", where = ns, print = FALSE,
          exit = quote(get(".selimit2_rec", globalenv())$vct(returnValue())))
    trace("se_edge_limit", where = ns, print = FALSE,
          exit = quote(get(".selimit2_rec", globalenv())$lim(returnValue(),
                                                              fit)))
  })
})
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/release/run-tests.R")
rm(".selimit2_rec", envir = globalenv())

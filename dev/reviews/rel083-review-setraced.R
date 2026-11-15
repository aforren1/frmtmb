# Release review of 0.83.0: dev/rel083-se-traced.R for any library (a lane
# library ahead of rellib-r22), to attribute the SE records that moved.
# The 0.83.0 consolidation: lane selimit2's standard-error comparison over
# every core test file, keyed by record identity rather than by position.
# As dev/selimit2-run-traced.R, with two changes: the library is an
# argument-free choice of REL083_SE_LIB (rellib-r22 or rellib-r21), and
# each test_that() writes a TEST line, so that a record is named by its
# test and its rank within that test (dev/rel083-se-cmp.R). A file whose
# new tests add fits then shifts no record of another test.
#
#   REL083_SE_LIB=<lib> SELIMIT2_TRACE=<file> Rscript dev/rel083-se-traced.R <package> <test file>
#
# The record, one line per call:
#   TEST <description>
#   COV <k> np=<n> lost=<name:reason,...> se=<sqrt(diag(shown)), %.10g>
#   VCT <k> se_t=<...>
#   LIM <k> <name:kind,...> est= dir= gain= ...
# A LIM line belongs to the COV line that follows it.
local({
  LIB <- Sys.getenv("REL083_SE_LIB")
  LIB2 <- "C:/Users/adf44/source/r/rellib-r22"


  .libPaths(c(LIB, LIB2, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressMessages(library(frmtmb))
  out <- Sys.getenv("SELIMIT2_TRACE")
  stopifnot(nzchar(out))
  cat("TRACE lib", dirname(find.package("frmtmb")), "\n", file = out)
  ns <- asNamespace("frmtmb")
  rec <- new.env()
  rec$k <- 0L
  assign(".selimit2_rec", rec, envir = globalenv())
  fmt <- function(x) paste(sprintf("%.10g", x), collapse = ",")
  rec$test <- function(desc) {
    cat("TEST ", gsub("[\r\n]+", " ", paste(desc, collapse = " ")), "\n",
        sep = "", file = out, append = TRUE)
  }
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
    cat("LIM ", rec$k, " ",
        paste(names(li$kind), li$kind, sep = ":", collapse = ","),
        " est=", fmt(li$est), " fam=", fam, " n=", fit$frame[["n_obs"]],
        "\n", sep = "", file = out, append = TRUE)
  }
  suppressMessages({
    trace("test_that", where = asNamespace("testthat"), print = FALSE,
          tracer = quote(get(".selimit2_rec", globalenv())$test(desc)))
    trace("cov_from_hessian", where = ns, print = FALSE,
          exit = quote(get(".selimit2_rec", globalenv())$cov(returnValue())))
    trace("varcorr_trans_rows", where = ns, print = FALSE,
          exit = quote(get(".selimit2_rec", globalenv())$vct(returnValue())))
    trace("se_edge_limit", where = ns, print = FALSE,
          exit = quote(get(".selimit2_rec", globalenv())$lim(returnValue(),
                                                              fit)))
  })
})
Sys.setenv(FRMTMB_RELLIB = Sys.getenv("REL083_SE_LIB"))
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/release/run-tests.R")
rm(".selimit2_rec", envir = globalenv())

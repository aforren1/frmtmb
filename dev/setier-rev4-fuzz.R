# Reviewer of lane setier, check of punch 2b: are test-fuzz.R's changed
# counts (boundary 119 -> 112, SE warning 131 -> 136) all quadrature
# fits? Runs test-fuzz.R on the lane, logging for every se_report() call
# whether the fit used quadrature, and the conditions that follow. With
# "ungate", the punch-2b quadrature gates are switched off by tracers (the
# fit seen by se_edge_sd(), se_boundary_names() and se_boundary_stop()
# has quadrature FALSE, and se_sd_gain_up() still returns 0 for it, as in
# punch 2).
#   Rscript dev/setier-rev4-fuzz.R gate|ungate <log>
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_FUZZ = "true")
suppressMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
log <- a[2]
cat("", file = log)
suppressMessages({
  trace("se_report", where = ns, print = FALSE, tracer = substitute(
    cat("Q\t", isTRUE(fit$quadrature) || isTRUE(fit$was_q), "\n",
        file = LOG, append = TRUE), list(LOG = log)))
  trace("frm_message", where = ns, print = FALSE, tracer = substitute(
    cat("M\t", substr(paste0(...), 1, 40), "\n", file = LOG,
        append = TRUE), list(LOG = log)))
  trace("frm_warning", where = ns, print = FALSE, tracer = substitute(
    cat("W\t", substr(paste0(...), 1, 40), "\n", file = LOG,
        append = TRUE), list(LOG = log)))
  if (identical(a[1], "ungate")) {
    un <- quote(if (isTRUE(fit$quadrature)) {
      fit$quadrature <- FALSE
      fit$was_q <- TRUE
    })
    for (fn in c("se_edge_sd", "se_boundary_names", "se_boundary_stop",
                 "se_report")) {
      trace(fn, where = ns, print = FALSE, tracer = un)
    }
    trace("se_sd_gain_up", where = ns, print = FALSE,
          tracer = quote(if (isTRUE(fit$was_q)) return(0)))
  }
})
setwd("tests/testthat")
res <- test_file("test-fuzz.R", package = "frmtmb",
                 env = test_env("frmtmb"), reporter = "silent")
print(as.data.frame(res)[, c("test", "nb", "failed", "skipped")])
x <- readLines(log)
q <- NA
tab <- list()
for (l in x) {
  if (startsWith(l, "Q")) q <- grepl("TRUE", l)
  if (grepl("Boundary [(]singular", l)) tab$bnd <- c(tab$bnd, q)
  if (grepl("Standard errors are not available", l)) {
    tab$se <- c(tab$se, q)
  }
}
cat(a[1], ": boundary messages", length(tab$bnd), "(quadrature fits",
    sum(tab$bnd, na.rm = TRUE), ") | SE warnings", length(tab$se),
    "(quadrature fits", sum(tab$se, na.rm = TRUE), ")\n")

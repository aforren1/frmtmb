## Fourth probe: the RTMBp-only helpers, and hard evidence that a tape is
## actually split into parallel terms.
out <- file("dev/rtmbp-probe4.txt", open = "wt")
p <- function(...) cat(..., "\n", file = out, sep = "")
library(RTMBp)
for (nm in c("term_split", "get_term_nodes", "inactivate", "sdreport_xtra")) {
  p("### ", nm, " ###")
  p(paste(deparse(get(nm, envir = asNamespace("RTMBp"))), collapse = "\n"))
  p("")
}

source("dev/bench-rtmbp-data.R")
suppressPackageStartupMessages({library(Matrix); library(lme4)})

## Evidence: TMBad's info string for the tape, with and without autopar.
info <- function(threads, autopar, shape = "B") {
  TMB::openmp(threads, autopar = autopar, DLL = "RTMBp")
  dat <- build_shape(shape)
  obj <- MakeADFun(make_f(shape, dat), dat$parameters,
                   random = dat$random, silent = TRUE)
  i1 <- TMB::TransformADFunObject(obj$env$ADFun, method = "info")
  list(cfg = unlist(TMB::config(DLL = "RTMBp")[c("nthreads", "autopar")]),
       info = i1)
}
for (ap in c(FALSE, TRUE)) {
  r <- try(info(4, ap))
  p("### autopar=", ap, " threads=4 ###")
  p(paste(capture.output(print(r)), collapse = "\n"))
  p("")
}
close(out)
cat("done\n")

# Short check of lane gpfield: the F1 row (sweep:mm_V:bernoulli:40:1) under
# the R copy running, its code and every warning and message, by class
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
cat("R home:", R.home(), "| OPENBLAS_CORETYPE", Sys.getenv("OPENBLAS_CORETYPE"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
r <- gf_row("sweep:mm_V:bernoulli:40:1")
w <- character(0)
f <- withCallingHandlers(frm(r$fml, data = r$data, family = r$family),
  warning = function(e) { w <<- c(w, paste0("W ", class(e)[1], ": ", substr(conditionMessage(e), 1, 160))); invokeRestart("muffleWarning") },
  message = function(e) { w <<- c(w, paste0("M ", class(e)[1], ": ", substr(conditionMessage(e), 1, 160))); invokeRestart("muffleMessage") })
cat("OBF1 conv", f$opt$convergence, "ll", format(as.numeric(logLik(f)), digits = 10), "\n")
writeLines(w)

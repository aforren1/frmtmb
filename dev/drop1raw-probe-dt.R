# Does a fit on a tibble, a data.table or a named list keep its raw
# variables? Arg: base | lane. Log dev/drop1raw-out/probe-dt-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(5)
d <- data.frame(x = rnorm(50), x2 = rexp(50) + 0.1)
d$z <- d$x + log(d$x2) + rnorm(50)
for (kind in c("tibble", "data.table", "list")) {
  if (!requireNamespace(kind, quietly = TRUE) && kind != "list") {
    cat(kind, ": not installed\n"); next
  }
  dd <- switch(kind, tibble = tibble::as_tibble(d),
               data.table = data.table::as.data.table(d), list = as.list(d))
  w <- character(0)
  r <- tryCatch(withCallingHandlers(
    { f <- frm(bf(z ~ x + log(x2)), data = dd)
      paste("raw_vars:", paste(names(f$frame$raw_vars), collapse = ","),
            "| drop1:", tryCatch(paste(rownames(drop1(f)), collapse = ","),
                                 error = function(e) conditionMessage(e))) },
    warning = function(cw) { w <<- c(w, conditionMessage(cw))
                             invokeRestart("muffleWarning") },
    message = function(cm) invokeRestart("muffleMessage")),
    error = function(e) paste("ERROR", conditionMessage(e)))
  cat(kind, ":", r, "\n  warnings:", length(w), paste(substr(w, 1, 90), collapse = " / "), "\n")
}

# Lane setier: test-id-kron.R's merged-dpar fit, printed in full, to see
# why it flips between the boundary message and code 7 across runs.
#   Rscript dev/setier-idk.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ex <- parse("tests/testthat/test-id-kron.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, globalenv())
  }
}
d <- kron_traits()
cs <- character(0)
f <- withCallingHandlers(
  frm(bf(y1 ~ 1 + (1 | q | gr(id, cov = A)),
         sigma ~ 1 + (1 | q | gr(id, cov = A))) + gaussian(),
      data = d$wide, data2 = list(A = d$A),
      control = frmtmb_control(check_olre = "ignore")),
  warning = function(w) {
    cs <<- c(cs, paste("W", substr(conditionMessage(w), 1, 50)))
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cs <<- c(cs, paste("M", substr(conditionMessage(m), 1, 50)))
    invokeRestart("muffleMessage")
  })
o <- f$opt
cat(sprintf("conv %s msg %s obj %.12g iter %s\n", o$convergence, o$message,
            o$objective, o$iterations))
cat("par", sprintf("%.8g", o$par), "\n")
cat("grad", sprintf("%.3g", f$obj$gr(o$par)), "\n")
cat("stop", isTRUE(f$cache$se_boundary_stop), "kind",
    paste(f$cache$se_boundary_kind, collapse = ","), "\n")
cat(cs, sep = "\n")

# Lane wt-conditions, punch round 2: the recheck's R1 and R2 cases, on a
# library arm, printing class and message of each.
#   Rscript dev/conditions-punch2-probe.R <lib> > <log>   (seed 20260917)
.libPaths(c(commandArgs(TRUE)[1L], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260917)
d <- data.frame(y = rnorm(60), x = rnorm(60))
pe <- new.env(); assign("xq", rnorm(60), envir = pe)
de <- new.env(parent = pe)
assign("y", d$y, envir = de); assign("x", d$x, envir = de)
M <- diag(60)
cases <- list(
  env_parent = quote(frm(y ~ x + xq, data = de)),
  sar_data2 = quote(frm(y ~ x + sar(Wm), data = d, data2 = list(Wm = M))),
  fcor_data2 = quote(frm(y ~ x + fcor(Vm), data = d, data2 = list(Vm = M))),
  dot = quote(frm(y ~ ., data = d)))
for (nm in names(cases)) {
  r <- tryCatch({eval(cases[[nm]]); "fits"}, error = function(e)
    paste0(paste(class(e), collapse = "/"), " | ",
           substr(conditionMessage(e), 1, 100)))
  cat(sprintf("%-11s %s\n", nm, r))
}

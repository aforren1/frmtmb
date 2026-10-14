# Lane optstart, item R1: which change of the lane gives fit_imp2
# (brms_missings) code 0. Run on rellib-r8 with "save" first: it stores
# 0.70.0's cold start (par_template()). Then on the lane build: the fit
# from 0.70.0's start (given as `start`, so the least-squares start is
# off and only the mi() units and the restart rule differ), with
# restarts = 0 and 1, and the fit from the lane's start.
#   Rscript dev/optstart-r1-attrib.R <lib or ""> save|fit
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
data("nhanes", package = "mice")
bf2 <- bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) +
  set_rescor(FALSE)
rds <- "dev/optstart-log/r1-start-r8.rds"
if (a[2] == "save") {
  st <- unclass(par_template(bf2, data = nhanes))
  saveRDS(st, rds)
  print(st[c("beta", "betad")])
} else {
  st0 <- readRDS(rds)
  arms <- list(r8_start_restarts0 = list(st0, 0L),
               r8_start = list(st0, 1L),
               lane_start = list(NULL, 1L))
  for (nm in names(arms)) {
    w <- character()
    f <- withCallingHandlers(suppressMessages(frm(bf2, data = nhanes,
      start = arms[[nm]][[1]],
      control = frmtmb_control(restarts = arms[[nm]][[2]], verbose = TRUE))),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    cat(sprintf("%s: code %d %s | logLik %.9f | warnings %d %s\n", nm,
                f$opt$convergence, f$opt$message, as.numeric(logLik(f)),
                length(w), substr(paste(w, collapse = " | "), 1, 80)))
  }
}

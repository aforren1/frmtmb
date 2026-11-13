# Lane boundfam, consolidation item 2: the two student() rows of lane
# gpstop's sweep that move on this lane's build (exp_a n40 s2, pow_a n120
# s2). Fit each with verbose output and print the estimates, the
# log-likelihood, and the gradient and the Hessian's spectrum at the
# stop, so the two builds' points can be compared.
#   Rscript dev/boundfam-gp.R <lib or "base"> <out .rds>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
# the sweep's model builders, read from the merged tree
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/gpstop-fx2-sweep.R")
eval(parse(text = src[45:87]))
keep <- list()
for (r in list(c("exp_a", 40, 2), c("pow_a", 120, 2))) {
  d <- mk(r[1], "student", as.integer(r[2]), as.integer(r[3]))
  w <- character(0)
  m <- character(0)
  f <- withCallingHandlers(
    frm(spell(r[1]), data = d, family = student(),
        control = frmtmb_control(verbose = TRUE)),
    warning = function(e) {
      w <<- c(w, paste0(class(e)[1], ": ", substr(conditionMessage(e), 1, 90)))
      invokeRestart("muffleWarning")
    }, message = function(e) {
      m <<- c(m, substr(conditionMessage(e), 1, 160))
      invokeRestart("muffleMessage")
    })
  nm <- paste0(r[1], "_n", r[2], "_s", r[3])
  cat(sprintf("\n== %s ll %.6f conv %d\n", nm, as.numeric(logLik(f)),
              f$opt$convergence))
  for (x in m) cat("  M:", x, "\n")
  for (x in w) cat("  W:", x, "\n")
  p <- f$opt$par
  cat("  par:", sprintf("%.6g", p), "\n")
  g <- tryCatch(drop(f$obj$gr(p)), error = function(e) NA)
  cat("  grad:", sprintf("%.3g", g), "\n")
  keep[[nm]] <- list(par = p, ll = as.numeric(logLik(f)), data = d)
}
saveRDS(keep, a[2])

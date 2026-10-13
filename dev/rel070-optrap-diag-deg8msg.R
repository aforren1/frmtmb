# Release 0.70.0: the probe's degree-8 raw polynomial (seed 8, x on
# [1, 2]) with every condition's class and full text, because the probe
# cuts each message at 110 characters and the release must show what
# the "illcond" message says.
#   Rscript dev/rel070-optrap-diag-deg8msg.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
for (re in c(FALSE, TRUE)) {
  deg <- 8
  set.seed(deg)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                         if (re) "+ (1 | g)"))
  w <- character()
  f <- withCallingHandlers(frm(fo, data = d),
    warning = function(x) {
      w <<- c(w, paste0("[", paste(class(x), collapse = ","), "] ",
                        conditionMessage(x)))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      w <<- c(w, paste0("[message] ", conditionMessage(x)))
      invokeRestart("muffleMessage")
    })
  cat("\n== degree 8, (1 | g):", re, "| code", f$opt$convergence, "\n")
  cat(w, sep = "\n---\n")
  cat("\n")
}

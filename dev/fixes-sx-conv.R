# Lane fixes, punch round: optimizer verdicts of smooth fits, base vs
# lane, after the diagonal.penalty reparameterization.
#   Rscript dev/fixes-sx-conv.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
forms <- list(y ~ s(x1) + s(x2), y ~ s(x1, bs = "cr", k = 6),
              y ~ s(x1, by = g) + g, y ~ s(x1, by = z), y ~ s(x1, x2),
              y ~ t2(x1, x2), y ~ s(x0) + s(x1) + s(x2) + s(x3))
for (seed in 1:4) {
  set.seed(seed)
  d <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
  d$z <- runif(200)
  d$g <- factor(sample(c("a", "b", "c"), 200, TRUE))
  for (fo in forms) {
    w <- character()
    fit <- withCallingHandlers(frm(bf(fo), data = d), warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
    g <- max(abs(fit$obj$gr(fit$opt$par)))
    cat(sprintf("seed %d %-34s conv %d %-40s logLik %.8f maxgrad %.1e %s\n",
                seed, deparse(fo), fit$opt$convergence,
                substr(fit$opt$message, 1, 40), as.numeric(logLik(fit)), g,
                if (length(w)) paste("WARN:", substr(w[1], 1, 60)) else ""))
  }
}

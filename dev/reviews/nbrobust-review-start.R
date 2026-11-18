# Review of lane nbrobust: a start far along the poisson limit on data
# whose maximum is at a FINITE shape. Does the fit leave the plateau, and
# what does it say? Overdispersed data, seed 5, n = 200, shape 2.
#   Rscript dev/reviews/nbrobust-review-start.R <base|lib>
args <- commandArgs(TRUE)
arm <- args[1]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
catch <- function(expr) {
  w <- character(0)
  mc <- character(0)
  mm <- character(0)
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    mc <<- c(mc, class(x)[1])
    mm <<- c(mm, conditionMessage(x))
    invokeRestart("muffleMessage")
  })
  list(v = v, w = w, mc = mc, mm = mm)
}
set.seed(5)
d0 <- data.frame(x = stats::rnorm(200))
d0$y <- stats::rnbinom(200, size = 2, mu = exp(1 + 0.4 * d0$x))
for (k in c(0, 5, 10)) {
  d <- d0
  d$x <- d$x * (1 + 2e-15 * k)
  for (st in c(NA, 10, 15, 20, 25, 30)) {
    r <- catch(if (is.na(st)) {
      frm(y ~ x, family = negbinomial(), data = d)
    } else {
      frm(y ~ x, family = negbinomial(), data = d, start = list(betad = st))
    })
    f <- r$v
    g <- f$obj$gr(f$opt$par)
    cat(sprintf(paste("k %2d start %4s code %d lsh %8.4f logLik %.10f",
                      "gr(shape) %10.3g\n"),
                k, st, f$opt$convergence, f$opt$par[["betad"]],
                as.numeric(logLik(f)), g[names(f$opt$par) == "betad"][1]))
    if (length(r$mc)) {
      cat("   messages:", paste(unique(r$mc), collapse = ", "), "\n")
      for (mm in unique(r$mm)) cat("   MSG:", substr(gsub("\n", " ", mm), 1, 200), "\n")
    }
    for (w in unique(r$w)) cat("   WARNING:", substr(gsub("\n", " ", w), 1, 200), "\n")
  }
}

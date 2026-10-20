# Re-check of lane illse2, item 6: TMB's ADFun finalizer and the lazy
# load of TMB:::isNullPointer.
#   Rscript dev/illse2-rev-fin.R <rtmb|ordinary|wide> <lane|base> [force]
a <- commandArgs(TRUE)
mode <- a[1]
lib <- if (a[2] == "lane") "C:/Users/adf44/source/r/wt-illse2-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (isTRUE(a[3] == "force")) {
  setHook(packageEvent("TMB", "onLoad"), function(...) {
    invisible(get("isNullPointer", envir = asNamespace("TMB")))
  })
}
cat(mode, a[2], if (isTRUE(a[3] == "force")) "forced" else "", "\n")
if (mode == "rtmb") {
  # no frmtmb: many RTMB tapes left to the collector under memory
  # pressure
  suppressMessages(library(RTMB))
  for (i in 1:300) {
    f <- MakeADFun(function(p) sum(p$x^2), list(x = rnorm(50)),
                   silent = TRUE)
    junk <- numeric(2e5)
  }
  invisible(gc())
  cat("rtmb loop done\n")
} else {
  suppressMessages(library(frmtmb))
  set.seed(1)
  n <- 2000
  for (i in 1:20) {
    d <- data.frame(x = rnorm(n), z = rnorm(n),
                    g = factor(sample(1:30, n, TRUE)))
    d$y <- 1 + 0.5 * d$x + rnorm(30, 0, 0.5)[d$g] + rnorm(n)
    f <- if (mode == "ordinary") {
      frm(y ~ x + z + (1 | g), data = d)
    } else {
      d$c <- rbinom(n, 1, plogis(0.5 * d$x))
      frm(bf(c ~ x + z + (1 | g)) + bernoulli(), data = d)
    }
    s <- summary(f)
  }
  cat("fits done\n")
}

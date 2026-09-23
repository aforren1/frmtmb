# After an interrupted R CMD check (the session was cut at 17:03 while
# `R CMD build` was installing frmtmb to build its vignettes), verify
# the three libraries the lane reads: no hollow package directory, and
# a FIT rather than a version string, as dev/machine-library.md
# requires.  Rscript dev/correct-libcheck.R
libs <- c("C:/Users/adf44/source/r/correct-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
for (lib in libs) {
  d <- list.dirs(lib, recursive = FALSE, full.names = FALSE)
  d <- setdiff(d, grep("^00LOCK", d, value = TRUE))
  hollow <- d[!file.exists(file.path(lib, d, "DESCRIPTION"))]
  cat(sprintf("%-55s dirs %3d hollow %d %s\n", lib, length(d),
              length(hollow),
              if (length(hollow)) paste(hollow, collapse = " ") else ""))
}
.libPaths(libs)
suppressMessages(library(frmtmb))
set.seed(1)
n <- 200
d <- data.frame(x = rnorm(n), g = factor(rep(1:10, each = 20)))
d$y <- 1 + 0.5 * d$x + rnorm(10, 0, 0.7)[d$g] + rnorm(n)
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x + rnorm(10, 0, 0.5)[d$g]))
f1 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
f2 <- frm(bf(cnt ~ x + (1 | g)) + poisson(), data = d)
cat(sprintf("gaussian (1 | g) logLik %.9f\npoisson  (1 | g) logLik %.9f\n",
            as.numeric(logLik(f1)), as.numeric(logLik(f2))))
for (p in c("frmtmb", "frmtmb.sample", "frmtmb.eam", "frmtmb.latent",
            "frmtmb.spline", "frmtmb.coupling", "frmtmb.learn",
            "frmtmb.ode", "brms", "testthat", "bayesplot", "mgcv",
            "rstan", "tmbstan", "RTMB")) {
  ok <- tryCatch({
    suppressMessages(loadNamespace(p))
    paste(format(packageVersion(p)), dirname(find.package(p)))
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("%-16s %s\n", p, ok))
}

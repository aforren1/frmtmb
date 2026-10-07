# Reviewer of lane optima, claim 8a: escape_stationary() does not leave
# the objective at the optimum it keeps. A skew normal (stationary point
# at alpha = 0) with a random intercept: when the escape ran, are the
# reported modes those of the reported optimum?
#   Rscript dev/optima-rev-escape.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
n_run <- 0L
n_off <- 0L
for (s in 1:40) {
  set.seed(s)
  g <- factor(rep(1:15, each = 10))
  u <- rnorm(15, sd = 0.8)
  d <- data.frame(g = g, x = rnorm(150))
  d$y <- 1 + 0.5 * d$x + u[g] + rnorm(150)
  f <- tryCatch(suppressWarnings(frm(bf(y ~ x + (1 | g)),
                                     family = skew_normal(), data = d)),
                error = function(e) NULL)
  if (is.null(f) || is.null(f$opt$stationary_escape)) next
  n_run <- n_run + 1L
  ref <- frmtmb:::solved_par_list(f$obj, f$opt$par)
  dev <- max(abs(f$estimates$b - ref$b))
  if (dev > 1e-6) n_off <- n_off + 1L
  cat(sprintf("seed %2d escape %s | max |b - b at the optimum| %.3g\n", s,
              paste(format(f$opt$stationary_escape, digits = 3),
                    collapse = "/"), dev))
}
cat("escape ran on", n_run, "of 40; modes off by > 1e-6 on", n_off, "\n")

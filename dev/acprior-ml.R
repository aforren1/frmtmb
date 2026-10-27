# frm() with no prior, and with a user prior on each class, before and
# after the lane: the ML and MAP fits must not move. Arg 1: before/after.
# Log: dev/acprior-out/ml-<arm>.txt; compare with diff.
arm <- commandArgs(TRUE)[1]
lib <- c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
         "C:/Users/adf44/source/r/rellib-r13",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(11)
ng <- 30; k <- 5
d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                x = rnorm(ng * k))
e <- as.vector(vapply(seq_len(ng), function(i) {
  as.vector(stats::filter(rnorm(k), 0.5, "recursive"))
}, numeric(k)))
d$y <- 1 + 0.5 * d$x + e
cases <- list(
  ar = list(y ~ x + ar(t, g), set_prior("normal(0.3, 0.2)", class = "ar")),
  arcov = list(y ~ x + ar(t, g, cov = TRUE),
               set_prior("normal(0.3, 0.2)", class = "ar")),
  ar2cov = list(y ~ x + ar(t, g, p = 2, cov = TRUE),
                set_prior("normal(0, 0.5)", class = "ar")),
  macov = list(y ~ x + ma(t, g, cov = TRUE),
               set_prior("normal(0.3, 0.2)", class = "ma")),
  armacov = list(y ~ x + arma(t, g, cov = TRUE),
                 set_prior("normal(0.3, 0.2)", class = "ar")),
  cosy = list(y ~ x + cosy(t, g), set_prior("beta(2, 2)", class = "cosy")),
  unstr = list(y ~ x + unstr(t, g), set_prior("lkj(2)", class = "cortime")))
for (nm in names(cases)) {
  for (wp in c(FALSE, TRUE)) {
    f <- suppressWarnings(frm(cases[[nm]][[1]], data = d,
                              prior = if (wp) cases[[nm]][[2]]))
    cat(sprintf("%-8s prior=%-5s obj %.10f thetaac %s\n", nm, wp,
                f$opt$objective,
                paste(sprintf("%.8f", f$estimates$thetaac), collapse = " ")))
  }
}

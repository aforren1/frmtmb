# Lane optrap: the one-level group with a bound that stops at the
# evaluation limit on the lane (dev/optrap-onelevel.R): does a restart
# in raw units, or from the start with raw units, reach base's optimum?
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
X <- read.delim("dev/optrap-log/onelevel-wt-optrap-lib.tsv")
B <- read.delim("dev/optrap-log/onelevel-rellib-r7.tsv")
bad <- X$seed[X$bound & X$code != 0]
for (s in bad) {
  set.seed(s)
  d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
  d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
  f <- suppressWarnings(frm(bf(y ~ x + (1 | g1)), family = gaussian(),
    data = d3, prior = set_prior("", class = "b", ub = 0.1)))
  o <- ns$optimize_obj(f$obj, f$control, f$cache$bounds, par_units = NULL,
                       start_par = f$opt$par)
  cat(sprintf("seed %d: lane %.6f code %d | raw-unit rerun %.6f code %d | base %.6f\n",
              s, f$opt$objective, f$opt$convergence, o$objective,
              o$convergence, B$obj[B$seed == s & B$bound]))
}

# test-reparam.R's "two routes agree on the correlated model" statistic
# on one build, at iter 1000 (the test's) and 2000, with the bulk ESS of
# each chain's theta.
#   Rscript dev/moparity-rp-iter.R base|lane
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(31)
ng <- 12L
dd <- data.frame(x = stats::rnorm(120),
                 g = factor(rep(seq_len(ng), length.out = 120)))
u <- matrix(stats::rnorm(2 * ng), 2)
dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x + 0.7 * u[1, dd$g] +
                       0.4 * u[2, dd$g] * dd$x, 1)
for (it in c(1000, 2000)) {
  mk <- function(ncp) {
    suppressWarnings(suppressMessages(
      frm_sample(bf(y ~ x + (x | g)), family = gaussian(), data = dd,
                 chains = 1, iter = it, refresh = 0, seed = 33,
                 reparameterize = ncp)))
  }
  cen <- mk(FALSE)$draws
  ncp <- mk(TRUE)$draws
  for (nm in grep("^theta", colnames(cen), value = TRUE)) {
    a <- cen[, nm]
    b <- ncp[, nm]
    cat(sprintf("%s iter %d %s: sd cen %.3f ncp %.3f | |diff| %.3f vs 0.5 max %.3f | ess cen %.0f ncp %.0f\n",
                arm, it, nm, sd(a), sd(b), abs(sd(a) - sd(b)),
                0.5 * max(sd(a), sd(b)), posterior::ess_bulk(a),
                posterior::ess_bulk(b)))
  }
}

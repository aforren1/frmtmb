# Reviewer: test-reparam.R's "the two routes agree on the correlated
# model as well" gate over sampler seeds, at iter 1000 (0.69.0's) and
# 2000 (the lane's), on one build. Per seed: pass/fail of the whole gate
# and the largest sd-statistic ratio |sd(a) - sd(b)| / (0.5 max sd).
#   Rscript dev/moparity-rev-rp-seeds.R base|lane first last
args <- commandArgs(TRUE)
arm <- args[1]
seeds <- as.integer(args[2]):as.integer(args[3])
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
set.seed(31)
ng <- 12L
dd <- data.frame(x = stats::rnorm(120),
                 g = factor(rep(seq_len(ng), length.out = 120)))
u <- matrix(stats::rnorm(2 * ng), 2)
dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x + 0.7 * u[1, dd$g] +
                       0.4 * u[2, dd$g] * dd$x, 1)
gate <- function(cen, ncp) {
  keep <- setdiff(colnames(cen), c("lp__", grep("^b\\[", colnames(cen),
                                                 value = TRUE)))
  fails <- character(0)
  worst <- 0
  worst_nm <- ""
  for (nm in keep) {
    a <- cen[, nm]
    b <- ncp[, nm]
    mcse <- sqrt(stats::var(a) / 50 + stats::var(b) / 50)
    if (!(abs(mean(a) - mean(b)) < 6 * mcse + 1e-8)) {
      fails <- c(fails, paste0(nm, ":mean"))
    }
    r <- abs(stats::sd(a) - stats::sd(b)) / (0.5 * max(stats::sd(a),
                                                      stats::sd(b)))
    if (r > worst) {
      worst <- r
      worst_nm <- nm
    }
    if (!(r < 1)) fails <- c(fails, paste0(nm, ":sd"))
  }
  for (nm in grep("^b\\[", colnames(cen), value = TRUE)) {
    a <- cen[, nm]
    b <- ncp[, nm]
    if (!(abs(mean(a) - mean(b)) < 0.5 * max(stats::sd(a), stats::sd(b)))) {
      fails <- c(fails, paste0(nm, ":bmean"))
    }
  }
  list(fails = fails, worst = worst, worst_nm = worst_nm)
}
for (s in seeds) {
  for (it in c(1000, 2000)) {
    mk <- function(ncp) {
      suppressWarnings(suppressMessages(
        frm_sample(bf(y ~ x + (x | g)), family = gaussian(), data = dd,
                   chains = 1, iter = it, refresh = 0, seed = s,
                   reparameterize = ncp)))
    }
    t0 <- proc.time()[[3]]
    cen <- mk(FALSE)$draws
    ncp <- mk(TRUE)$draws
    g <- gate(cen, ncp)
    cat(sprintf("RES arm %s seed %d iter %d | pass %s | worst sd ratio %.3f (%s) | fails %s | secs %.1f\n",
                arm, s, it, !length(g$fails), g$worst, g$worst_nm,
                if (length(g$fails)) paste(g$fails, collapse = ",") else "-",
                proc.time()[[3]] - t0))
  }
}

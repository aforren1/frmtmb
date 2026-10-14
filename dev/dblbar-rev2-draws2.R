# Reviewer round 2, lane dblbar: are frmtmb.sample draws identical across
# the three spellings of one || term, and across two runs of one
# spelling? set.seed(1) before each call (seed = alone does not fix 2 chains). Data seed 20261013; sampler seed 3, 2 chains of 500.
# Output: dev/dblbar-rev2-log/draws-setseed.txt
.libPaths(c("C:/Users/adf44/source/r/wt-dblbar-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
sink(file.path(wt, "dev/dblbar-rev2-log/draws-setseed.txt"), split = TRUE)
dd <- dblbar_data()
fits <- list(
  dbar = frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd),
  corF = frm(bf(y ~ f + (1 + f | gr(g, cor = FALSE))), family = gaussian(),
             data = dd),
  diag = frm(bf(y ~ f + diag(1 + f | g)), family = gaussian(), data = dd))
smp <- function(fit) { set.seed(1); frm_sample(fit, chains = 2, iter = 500, seed = 3,
                                refresh = 0, cores = 1) }
m <- lapply(fits, function(f) as.matrix(smp(f)))
m$dbar2 <- as.matrix(smp(fits$dbar))
for (k in c("corF", "diag", "dbar2")) {
  a <- m$dbar
  b <- m[[k]]
  cat(sprintf("dbar vs %-5s dims %s/%s colnames %s identical %s max|diff| %.3e\n",
              k, paste(dim(a), collapse = "x"), paste(dim(b), collapse = "x"),
              identical(colnames(a), colnames(b)), identical(a, b),
              if (identical(dim(a), dim(b))) max(abs(unname(a) - unname(b)))
              else NA))
}
cat("identical opt$par:", identical(fits$dbar$opt$par, fits$corF$opt$par),
    identical(fits$dbar$opt$par, fits$diag$opt$par), "\n")
cat("identical obj par:", identical(fits$dbar$obj$par, fits$diag$obj$par),
    "\n")
sink()

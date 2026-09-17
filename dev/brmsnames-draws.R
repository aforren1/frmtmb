## The draws every lane-brmsnames measurement reads, built per ARM so a
## build's own naming is on its own object.
##
##   Rscript dev/brmsnames-draws.R base   # frmtmb 0.58.0 / sample 0.6.0
##   Rscript dev/brmsnames-draws.R lane   # this worktree's build
##
## Construction: the recipe of dev/brmsmatch-measure.R, unchanged, so a
## number here is comparable to that lane's. y ~ x + (1 | g), gaussian,
## n = 120, data seed 9, frm_sample(chains = 4, iter = 1000,
## seed = 20260915). Saved to dev/stan-cache/brmsnames-draws-<arm>.rds.
arm <- commandArgs(trailingOnly = TRUE)[1L]
stopifnot(arm %in% c("base", "lane"))
source("dev/brmsnames-libs.R")
brmsnames_libs(arm)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
cat("frmtmb.sample", format(packageVersion("frmtmb.sample")), "from",
    dirname(system.file(package = "frmtmb.sample")), "\n")
cat("StanHeaders", format(packageVersion("StanHeaders")), "\n")
stopifnot(packageVersion("StanHeaders") == "2.32.10")

set.seed(9)
dd <- data.frame(x = stats::rnorm(120), g = factor(rep(1:6, 20)))
dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x + stats::rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
ds <- q(frm_sample(fit, chains = 4, iter = 1000, refresh = 0,
                   seed = 20260915))
out <- sprintf("dev/stan-cache/brmsnames-draws-%s.rds", arm)
saveRDS(list(ds = ds, fit = fit, data = dd), out)
cat("draws", nrow(ds$draws), "x", ncol(ds$draws), "\n")
cat("names", colnames(ds$draws), "\n")
cat("saved", out, "\n")

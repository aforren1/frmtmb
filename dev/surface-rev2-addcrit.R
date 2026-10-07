# Reviewer re-check of lane surface, B1: add_criterion() overwrite on
# draws, overwrite = FALSE, and loo_compare() after an overwrite.
#   Rscript dev/surface-rev2-addcrit.R > dev/surface-rev-out/rev2-addcrit.txt
source("dev/surface-rev-env.R")
rev_env("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(9)
dd <- data.frame(x = rnorm(60), z = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
sm <- function(f) suppressWarnings(suppressMessages(
  frm_sample(f, family = gaussian(), data = dd, chains = 2, iter = 600,
             refresh = 0, seed = 3)))
ds <- sm(bf(y ~ x)); dz <- sm(bf(y ~ x + z))
dims <- function(l) paste(attr(l, "dims"), collapse = " x ")
sw <- function(e) suppressWarnings(e)
ds2 <- sw(add_criterion(ds, "loo", ndraws = 100))
cat("stored with ndraws = 100:", dims(ds2$criteria$loo), "\n")
ds3 <- sw(add_criterion(ds2, "loo", overwrite = TRUE))
fresh <- sw(loo(ds))
cat("overwrite = TRUE:", dims(ds3$criteria$loo),
    "| estimates identical to fresh loo(ds):",
    identical(ds3$criteria$loo$estimates, fresh$estimates),
    "| identical to the old one:",
    identical(ds3$criteria$loo$estimates, ds2$criteria$loo$estimates), "\n")
ds4 <- sw(add_criterion(ds2, "loo"))
cat("overwrite = FALSE keeps the stored one:",
    identical(ds4$criteria$loo, ds2$criteria$loo), dims(ds4$criteria$loo), "\n")
ds5 <- sw(add_criterion(ds2, c("loo", "waic"), overwrite = TRUE))
cat("overwrite of a set: loo", dims(ds5$criteria$loo), "| waic",
    dims(ds5$criteria$waic), "\n")
dz2 <- sw(add_criterion(dz, "loo"))
cmp_old <- sw(loo_compare(ds2, dz2))
cmp_new <- sw(loo_compare(ds3, dz2))
cmp_ref <- loo::loo_compare(list(ds = fresh, dz = dz2$criteria$loo))
ed <- function(cm) {
  d <- as.data.frame(cm)
  nm <- if ("model" %in% names(d)) d$model else rownames(d)
  setNames(d$elpd_diff, nm)[sort(nm)]
}
cat("loo_compare(ds3, dz2) elpd_diff:", format(ed(cmp_new)),
    "| loo::loo_compare(fresh, dz2):", format(ed(cmp_ref)),
    "| before overwrite (100 draws):", format(ed(cmp_old)), "\n")
cat("names:", names(ed(cmp_new)), "|", names(ed(cmp_ref)), "\n")

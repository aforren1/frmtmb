# Reviewer of lane surface, claim 10: add_criterion() on draws against
# brms 2.23.0's semantics (add_criterion.brmsfit sets
# x$criteria[new_criteria] <- NULL before it computes, so overwrite =
# TRUE recomputes).
#   Rscript dev/surface-rev-addcrit.R > dev/surface-rev-out/addcrit.txt
source("dev/surface-rev-env.R")
rev_env("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(9)
dd <- data.frame(x = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
ds <- suppressWarnings(suppressMessages(
  frm_sample(bf(y ~ x), family = gaussian(), data = dd, chains = 2,
             iter = 600, refresh = 0, seed = 3)))
cat("ndraws(ds):", ndraws(ds), "\n")
dims <- function(l) paste(attr(l, "dims"), collapse = " x ")
ds2 <- suppressWarnings(add_criterion(ds, "loo", ndraws = 100))
cat("add_criterion(ds, 'loo', ndraws = 100): stored dims", dims(ds2$criteria$loo),
    "\n")
ds3 <- suppressWarnings(add_criterion(ds2, "loo", overwrite = TRUE))
cat("add_criterion(ds2, 'loo', overwrite = TRUE): stored dims",
    dims(ds3$criteria$loo), "| identical to the old one (estimates):",
    identical(ds3$criteria$loo$estimates, ds2$criteria$loo$estimates),
    "| fresh loo(ds) dims", dims(suppressWarnings(loo(ds))), "\n")
ds4 <- suppressWarnings(add_criterion(ds2, "waic"))
cat("criteria after adding waic:", names(ds4$criteria), "\n")
cmp <- suppressWarnings(loo_compare(ds2, ds))
print(cmp)
# brms's loo() with no argument but a stored criterion
cat("loo(ds2) dims (stored):", dims(loo(ds2)), "| loo(ds2, ndraws = 200):",
    dims(suppressWarnings(loo(ds2, ndraws = 200))), "\n")
f <- tempfile()
ds5 <- suppressWarnings(add_criterion(ds, "loo", file = f))
cat("file written:", file.exists(paste0(f, ".rds")), "\n")

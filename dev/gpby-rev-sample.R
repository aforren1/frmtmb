# Reviewer: frmtmb.sample's kriging draw. (1) Two rows at one unseen
# position: is the drawn field one value (the lane says the eigen root
# keeps it so)? (2) Rows at observed positions: no draw. (3) The cost of
# posterior_epred() at 300 unseen positions over 300 draws, both arms.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
options(mc.cores = 1)
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d$y <- 0.5 + sin(d$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d)
ds <- suppressWarnings(suppressMessages(
  frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(), chains = 1,
             iter = 600, refresh = 0, seed = 4)))
mm <- as.matrix(ds)
cat("draw columns:", paste(utils::head(colnames(mm), 6), collapse = " "),
    "| sd of column 2", sd(mm[, 2]), "\n")
nd <- data.frame(x = c(d$x[1], 6.7, 6.7, 7.4, 7.4 + 1e-9, d$x[2]))
set.seed(1)
e <- posterior_epred(ds, newdata = nd)
cat(sprintf(paste0("draws %d | rows 2,3 same unseen x: max |diff| %.3e ",
                   "(sd of row 2 %.3e) | rows 4,5 1e-9 apart: max |diff| ",
                   "%.3e\n"), nrow(e), max(abs(e[, 2] - e[, 3])),
            sd(e[, 2]), max(abs(e[, 4] - e[, 5]))))
e_seen <- posterior_epred(ds, newdata = d[c(1, 2), , drop = FALSE])
cat(sprintf("observed rows: epred(newdata) vs epred(newdata alone) max |diff| %.3e\n",
            max(abs(e[, c(1, 6)] - e_seen))))
big <- data.frame(x = seq(6.05, 9, length.out = 300))
tt <- system.time(eb <- posterior_epred(ds, newdata = big, ndraws = 300))
cat(sprintf("TIME %s posterior_epred 300 unseen rows x 300 draws: %.2f s\n",
            arm, tt[["elapsed"]]))
cat("DONE\n")

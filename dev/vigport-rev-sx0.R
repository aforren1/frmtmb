# Reviewer: the sigma ~ s(x0) term of brms_distreg fit_smooth1 on the
# lane's port_seed data: the same column factor as s(x1) and s(x2)?
#
#   Rscript dev/vigport-rev-sx0.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-vigport/dev/brms-port/shim.R")
set.seed(port_seed("brms_distreg.11.1"))
ds <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
bf2 <- bf(y ~ s(x1) + s(x2) + (1 | fac), sigma ~ s(x0) + (1 | fac))
fd <- suppressWarnings(frm(bf2, data = ds))
sdd <- brms::standata(brms::bf(y ~ s(x1) + s(x2) + (1 | fac),
                               sigma ~ s(x0) + (1 | fac)), data = ds)
nm <- names(fd$frame$linpreds)
cat("linpreds:", nm, "\n")
Xs <- fd$frame$linpreds[[grep("sigma", nm)[1]]]$X
cat("frmtmb sigma columns:", colnames(Xs), "\n")
a <- Xs[, grep("s[(]x0[)][.]fx1", colnames(Xs))]
b <- sdd$Xs_sigma[, "sx0_1"]
k <- stats::coef(stats::lm(a ~ 0 + b))
e <- fixef(fd)["sigma_sx0_1", "Estimate"]
cat(sprintf(paste0("s(x0) in sigma: a = %+.4f * b (max resid %.1e); ",
                   "frmtmb %.4f -> brms column %.4f (brms posterior mean ",
                   "-2.245, sd 2.416, lane's brms fit)\n"),
            k, max(abs(a - k * b)), e, e * k))

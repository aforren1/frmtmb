# Reviewer probe 5: mi() with missing predictor values
.libPaths(c("C:/Users/adf44/source/r/wt-hollow-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb.sample"), "\n")
isf <- frmtmb.sample:::draws_newdata_is_fitted
mk <- function(fit) structure(list(stanfit = NULL, draws = matrix(0, 1, 1),
                                   fit = fit), class = "frmtmb_draws")
set.seed(5)
n <- 60
dd <- data.frame(z = rnorm(n))
dd$xm <- dd$z + rnorm(n)
dd$y <- 1 + 0.5 * dd$xm + rnorm(n)
dd$xm[c(3, 10, 20)] <- NA
f <- suppressWarnings(suppressMessages(
  frm(bf(y ~ mi(xm)) + bf(xm | mi() ~ z) + set_rescor(FALSE), data = dd)))
cat("nobs frame:", nrow(f$data), " names:", paste(names(f$data), collapse = ","),
    "\n")
cat("mi: data accepted", isf(mk(f), dd), " fit$data accepted",
    isf(mk(f), f$data), "\n")
d1 <- dd; d1$xm[3] <- 0.5
cat("mi: a missing value filled in, accepted", isf(mk(f), d1), "\n")
d2 <- dd; d2$z[4] <- d2$z[4] + 1
cat("mi: z changed, accepted", isf(mk(f), d2), "\n")

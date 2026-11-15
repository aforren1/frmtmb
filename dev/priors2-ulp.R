# How many ends of a grid miss the round trip through frmtmb's own logit
# and log links, which the fit reads a bound through.
# Rscript dev/priors2-ulp.R > dev/priors2-out/ulp.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
lk <- frmtmb:::get_link("logit")
e1 <- seq(0.01, 0.99, by = 0.01)
r1 <- lk$linkinv(lk$linkfun(e1))
cat("logit: above", sum(r1 > e1), "below", sum(r1 < e1), "of", length(e1),
    "\n")
e2 <- seq(0.05, 5, by = 0.05)
r2 <- exp(log(e2))
cat("log: above", sum(r2 > e2), "below", sum(r2 < e2), "of", length(e2),
    "\n")
cat("0.3 through the logit:", sprintf("%.17g", lk$linkinv(lk$linkfun(0.3))),
    "\n")

# Reviewer round 2: can a location be pinned through an offset-only
# formula? Output: dev/famsimplex-rev2-offset.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
set.seed(31); n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n)); cl <- rbinom(n, 1, 0.4)
d$y <- ifelse(cl == 1, 2 + d$x, -1 + 0.5 * d$z) + rnorm(n, 0, 0.7)
d$o <- -1
r <- function(lab, e) cat(lab, ":", tryCatch({f <- suppressWarnings(e); sprintf("logLik %.6f; coef %s", logLik(f), paste(names(fixef(f)[, 1]), collapse = " "))}, error = function(x) paste("ERROR:", conditionMessage(x))), "\n")
r("mixture mu1 ~ 0 + offset(o)", frm(bf(y ~ x, mu1 ~ 0 + offset(o)), family = mixture(gaussian(), gaussian()), data = d))
dlb <- lba_simulate(300, v = c(2.2, 1.4, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
dlb$o <- 0.8
r("lba v3 ~ 0 + offset(o)", frm(bf(rt | vint(choice) ~ 1, v3 ~ 0 + offset(o)), family = lba(3), data = dlb))

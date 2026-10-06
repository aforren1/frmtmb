# Reviewer of lane ordmix: frmtmb's Links line for a gaussian mixture
# without a theta predictor, base build (pre-existing?)
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(200))
d$y <- d$x + ifelse(runif(200) < 0.4, 3, -1) + rnorm(200)
f <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = d)
cat(grep("Links", capture.output(print(f)), value = TRUE), "\n")

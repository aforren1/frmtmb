# Probe: what summary(), VarCorr(), sigma(), confint(), vcov(full = TRUE)
# return for sigma on an ar(cov = TRUE) fit (rellib-r14).
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-sigmainn/dev/acprior-agree.R")
i0 <- grep("^make_data <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
eval(parse(text = src[i0:i1]))
d <- make_data()$ar
fit <- frm(y ~ x + ar(time = t, gr = g, cov = TRUE), data = d)
s <- summary(fit)
str(s, max.level = 1)
print(s)
print(sigma(fit))
print(VarCorr(fit))
print(confint(fit))
print(names(fit$estimates))
print(fit$estimates)
print(vcov(fit, full = TRUE))
print(fit$obj$par)
print(class(fit$sdr))

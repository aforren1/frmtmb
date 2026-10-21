# Fix round 2: the log shapes the "does not fire" test compares.
.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(testthat)})
env <- test_env("frmtmb")
src <- readLines("tests/testthat/test-weights-unbounded.R")
eval(parse(text = src[1:(grep("^wu_refusal", src) - 1)]), env)
d <- env$wu_data(); G <- Gamma(link = "log")
dw <- d; dw$w <- d$w * 0.7
rw <- frm(bf(y | weights(w) ~ x + ar1(tt + 0 | g)), family = G, data = dw)
ru <- frm(bf(y ~ x + ar1(tt + 0 | g)), family = G, data = d)
se <- function(f) sqrt(diag(vcov(f)))
cat("rw betad", rw$estimates$betad, "ru betad", ru$estimates$betad, "\n")
print(se(rw)); print(se(ru))

# Reviewer: brms 2.23.0's VarCorr() on the lane's negbinomial2 (1 | g)
# fit (read-only rds): does brms report sigma as residual__ as frmtmb does?
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
b <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-out/brmsfit-nb2re.rds")
vc <- VarCorr(b)
print(names(vc))
print(vc$residual__)

## frmtmb's coding of the same mo(k): D, and the fit against mo(k - 1),
## which brms treats as the same model. Data seed 1024.
##   Rscript dev/brmsnames-1024-probe2.R > dev/brmsnames-1024-log/probe2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-brmsnames-lib", "C:/Users/adf44/source/r/rellib-r12", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("dev/brmsnames-1024-cases.R")
d <- brmsnames_data()
f1 <- frm(bf(y ~ mo(k)), data = d, family = gaussian())
d$k0 <- d$k - 1L
f0 <- frm(bf(y ~ mo(k0)), data = d, family = gaussian())
print(logLik(f1)); print(logLik(f0))
print(fixef(f1)); print(fixef(f0))
str(f1$frame$linpreds[[1]]$mo)
print(variables(f1))

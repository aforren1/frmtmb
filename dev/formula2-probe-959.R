# brmsfit-methods:958-959 on the suite's own fixture 3
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
library(testthat)
source("tests/testthat/helper-brms-suite.R")
fit3 <- brms_fixture(3)
up <- tryCatch(update(fit3, bf(~., family = acat())),
               error = function(e) conditionMessage(e))
if (is.character(up)) print(up) else {
  print(family(up)$family)
  print(formula(up))
  print(logLik(up))
  print(is(up, "brmsfit"))
}

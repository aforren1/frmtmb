# The behavioral failures on rellib-r8 that the lane's tests pin, through
# routes that need no new symbol: a brms family object, a family name
# string, brmsfamily(). Argument "before" reads rellib-r8 alone.
# Output: dev/famsimplex-pins-<arm>.log
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (!identical(arm, "before")) "C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
dd <- fs_data_dirichlet(seed = 11)
try1 <- function(lab, expr) {
  r <- tryCatch({ f <- expr; paste("fits; logLik", format(as.numeric(logLik(f)), digits = 10)) },
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(lab, ":", substr(r, 1, 160), "\n")
}
try1("family = brms::dirichlet()", frm(bf(Y ~ x), family = brms::dirichlet(), data = dd))
try1("family = brms::logistic_normal()", frm(bf(Y ~ x), family = brms::logistic_normal(), data = dd))
try1("family = \"dirichlet\"", frm(bf(Y ~ x), family = "dirichlet", data = dd))
try1("brmsfamily(\"dirichlet2\")", frm(bf(Y ~ x), family = brmsfamily("dirichlet2"), data = dd))

# Punch round 1, item 5: what predict() and frm_linpred() do with the
# two unusual newdata the reviewer found (dev/simnewdata-review/rv-edge.R):
# an NA covariate, and a grouping column absent from newdata while an
# object of that name exists in the formula environment. simulate()
# should match them.
#   SIMNEWDATA_LIB=base Rscript dev/simnewdata-edge-predict.R
source("dev/simnewdata-prelude.R")
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(7)
n <- 120
d <- data.frame(x = rnorm(n), g = factor(rep(1:12, 10)))
d$y <- 0.5 + 0.8 * d$x + rnorm(12)[d$g] + rnorm(n, 0, 0.5)
fit <- frm(bf(y ~ x + (1 | g)), data = d)
show <- function(label, expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  cat("--", label, "\n")
  if (is.character(r) && length(r) == 1L) cat("  ", substr(r, 1, 300), "\n")
  else print(r)
  if (length(w)) cat("   warnings:", substr(paste(unique(w), collapse = " | "),
                                            1, 300), "\n")
}
nd_na <- data.frame(x = c(NA, 1), g = factor(c("1", "2")))
set.seed(1)
show("predict, NA covariate", predict(fit, newdata = nd_na, ndraws = 50))
show("frm_linpred, NA covariate", frm_linpred(fit, newdata = nd_na))
if ("newdata" %in% names(formals(frmtmb:::simulate.frmtmb_fit))) {
  show("simulate, NA covariate",
       simulate(fit, nsim = 2, seed = 1, newdata = nd_na))
}
g <- factor(rep("3", 5))
nd_nog <- data.frame(x = c(0, 1))
show("frm_linpred, global g of length 5", frm_linpred(fit, newdata = nd_nog))
set.seed(1)
show("predict, global g of length 5",
     predict(fit, newdata = nd_nog, ndraws = 50))
g <- factor(c("3", "4"))
show("frm_linpred, global g of length 2 (levels 3, 4)",
     frm_linpred(fit, newdata = nd_nog))
show("frm_linpred at g = 3, 4 given in newdata",
     frm_linpred(fit, newdata = data.frame(x = c(0, 1), g = factor(3:4))))
if ("newdata" %in% names(formals(frmtmb:::simulate.frmtmb_fit))) {
  g <- factor(rep("3", 5))
  show("simulate, global g of length 5",
       simulate(fit, nsim = 2, seed = 1, newdata = nd_nog))
  g <- factor(c("3", "4"))
  show("simulate, global g of length 2",
       simulate(fit, nsim = 2, seed = 1, newdata = nd_nog))
}

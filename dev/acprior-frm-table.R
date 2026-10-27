# frmtmb's prior table for each autocorrelation class, and the slots a
# prior takes. Arg 1: "before" (rellib-r13) or "after" (the lane's lib).
# Log: dev/acprior-out/frm-table-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
         "C:/Users/adf44/source/r/rellib-r13",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
set.seed(1)
n <- 40
d <- data.frame(y = rnorm(n), x = rnorm(n), t = rep(1:10, 4),
                g = rep(1:4, each = 10))
cases <- list(
  ar1 = y ~ x + ar(time = t, gr = g, p = 1),
  ar2 = y ~ x + ar(time = t, gr = g, p = 2),
  ma1 = y ~ x + ma(time = t, gr = g),
  arma = y ~ x + arma(time = t, gr = g),
  ar1cov = y ~ x + ar(time = t, gr = g, p = 1, cov = TRUE),
  ar2cov = y ~ x + ar(time = t, gr = g, p = 2, cov = TRUE),
  ma1cov = y ~ x + ma(time = t, gr = g, cov = TRUE),
  armacov = y ~ x + arma(time = t, gr = g, cov = TRUE),
  cosy = y ~ x + cosy(time = t, gr = g),
  unstr = y ~ x + unstr(time = t, gr = g)
)
show <- function(x) {
  df <- as.data.frame(x)
  keep <- intersect(c("prior", "class", "coef", "group", "resp", "dpar",
                      "lb", "ub", "source"), names(df))
  print(df[, keep])
}
for (nm in names(cases)) {
  cat("\n=====", nm, ":", deparse(cases[[nm]]), "\n")
  r <- tryCatch(show(default_prior(cases[[nm]], data = d)),
                error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}
cat("\n===== two slots on one parameter (ar1)\n")
pr <- c(set_prior("normal(0.5, 0.2)", class = "ar"),
        set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1"))
r <- tryCatch(validate_prior(pr, cases$ar1, data = d),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)
r <- tryCatch({
  f <- frm(cases$ar1, data = d, prior = pr)
  length(attr(f$prior, "entries") %||% f$prior_entries %||% list())
}, error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)

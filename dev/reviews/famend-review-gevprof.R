# Review of lane famend, item 3075: the profile interval with parm.
#   Rscript dev/reviews/famend-review-gevprof.R
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
rg <- function(n, mu, s, xi) mu + s * ((-log(runif(n)))^(-xi) - 1) / xi
for (xi0 in c(-0.63, -0.89)) {
  set.seed(4)
  d <- data.frame(y = rg(300, 0, 1, xi0))
  f <- suppressMessages(frm(y ~ 1, family = gen_extreme_value(), data = d))
  for (p in c("(Intercept)", "sigma_(Intercept)", "xi_(Intercept)")) {
    r <- tryCatch(withCallingHandlers(confint(f, parm = p, method = "profile"),
                                      warning = function(w) {
                                        cat("   W:", substr(conditionMessage(w),
                                                            1, 160), "\n")
                                        invokeRestart("muffleWarning")
                                      }),
                  error = function(e) paste("ERR", conditionMessage(e)))
    cat("xi0", xi0, p, ":\n")
    print(r)
  }
}

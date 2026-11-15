# Review of lane famend: the assertions dev/famend-perturb.R does not
# perturb (GEV xi0 -0.89: message, no warning, xi inside (-1 + 1e-3,
# -0.5)), 20 data sets y * (1 + 2e-15 k), under whichever R runs it.
#   Rscript dev/reviews/famend-review-perturb2.R
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "| BLAS",
    extSoftVersion()[["BLAS"]], "| LAPACK", La_library(), "\n")
rg <- function(n, mu, s, xi) mu + s * ((-log(runif(n)))^(-xi) - 1) / xi
set.seed(4)
g89 <- data.frame(y = rg(300, 0, 1, -0.89))
ok <- logical(0)
for (k in 1:20) {
  d <- g89
  d$y <- d$y * (1 + 2e-15 * k)
  w <- list(); m <- list()
  f <- withCallingHandlers(frm(y ~ 1, family = gen_extreme_value(), data = d),
    warning = function(x) { w[[length(w) + 1]] <<- x
      invokeRestart("muffleWarning") },
    message = function(x) { m[[length(m) + 1]] <<- x
      invokeRestart("muffleMessage") })
  xi <- frm_linpred(f, dpar = "xi", type = "response")[1]
  ok[k] <- !length(w) && any(vapply(m, inherits, NA,
                                    "frmtmb_gev_xi_nonregular")) &&
    xi < -0.5 && xi > -1 + 1e-3
}
cat(sprintf("PERTURB gev -0.89 message, xi in band %d of %d\n", sum(ok),
            length(ok)))
cat("END\n")

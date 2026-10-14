# Which GEV and Frechet two-component mixtures converge from frmtmb's
# default start (seed 35, n = 400).
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(35); n <- 400
rg <- function(n, m, s, k) m + s * (rexp(n)^(-k) - 1) / k
d <- data.frame(yg = c(rg(n/2, 0, 1, 0.1), rg(n/2, 6, 1, 0.1)))
d$yf <- c(exp(0.5) / gamma(1 - 1/5) * rexp(n/2)^(-1/5),
          exp(2.5) / gamma(1 - 1/5) * rexp(n/2)^(-1/5))
d$yc <- c(rexp(n/2) + 1, rexp(n/2) * 3 + 8)
for (v in c("yg", "yc")) {
  fm <- withCallingHandlers(frm(as.formula(paste(v, "~ 1")),
            family = mixture(gen_extreme_value(), gen_extreme_value()), data = d),
            warning = function(w) {cat("WARN:", conditionMessage(w), "\n"); invokeRestart("muffleWarning")})
  cat(v, "GEV mixture logLik", as.numeric(logLik(fm)), "conv", fm$opt$convergence, "\n")
  print(fixef(fm, flatten = TRUE))
}
for (v in c("yf", "yc")) {
  fm <- withCallingHandlers(frm(as.formula(paste(v, "~ 1")),
            family = mixture(frechet(), frechet()), data = d),
            warning = function(w) {cat("WARN:", conditionMessage(w), "\n"); invokeRestart("muffleWarning")})
  cat(v, "frechet mixture logLik", as.numeric(logLik(fm)), "conv", fm$opt$convergence, "\n")
  print(fixef(fm, flatten = TRUE))
}

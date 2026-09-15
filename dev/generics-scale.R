# Item 2.5d: measure which SCALE each generic reports, so the contract
# is written from the run rather than from reading the code.
#
# The lognormal fit is the case the plan measured, because it is the
# family where the three candidate scales are all different numbers:
# the linear predictor mu, the median exp(mu), and the mean
# exp(mu + sigma^2/2).
#
#   Rscript dev/generics-scale.R [LIB]
av <- commandArgs(trailingOnly = TRUE)
LIB <- if (length(av)) av[1] else "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))

set.seed(2026)
n <- 400
dd <- data.frame(x = rnorm(n), g = factor(rep(1:20, 20)))
eta <- 8 + 0.4 * dd$x + rnorm(20, 0, 0.3)[dd$g]
dd$y <- exp(rnorm(n, eta, 0.4))
fit <- frm(bf(y ~ x + (1 | g)) + lognormal(), data = dd)

mu <- predict(fit, type = "link")
sg <- sigma(fit)
cat("== lognormal fit, seed 2026, n = 400, dev/generics-scale.R ==\n")
cat(sprintf("sigma(fit)                      %.10f\n", sg))
cat(sprintf("predict(fit)[1]  (default)      %.10f\n", predict(fit)[1]))
cat(sprintf("predict(type='link')[1]         %.10f\n", mu[1]))
cat(sprintf("predict(type='response')[1]     %.10f\n",
            predict(fit, type = "response")[1]))
cat(sprintf("fitted(fit)[1]                  %.10f\n", fitted(fit)[1]))
cat(sprintf("exp(mu)[1]        (median)      %.10f\n", exp(mu[1])))
cat(sprintf("exp(mu+s^2/2)[1]  (mean)        %.10f\n",
            exp(mu[1] + sg^2 / 2)))

id <- function(lbl, a, b) {
  cat(sprintf("%-42s %.3e\n", lbl, max(abs(a - b)) / max(abs(b))))
}
cat("\nrelative residuals at full precision (identities, not fits):\n")
id("predict() vs predict(type='link')", predict(fit), mu)
id("predict(type='response') vs fitted()",
   predict(fit, type = "response"), fitted(fit))
id("fitted() vs exp(mu + sigma^2/2)", fitted(fit), exp(mu + sg^2 / 2))
cat(sprintf("%-42s %.6f\n", "fitted()/exp(mu), the plan's 1.079",
            mean(fitted(fit) / exp(mu))))
cat(sprintf("%-42s %.6f\n", "exp(sigma^2/2)", exp(sg^2 / 2)))

cat("\nresiduals:\n")
r <- residuals(fit)
cat(sprintf("%-42s %.3e\n", "residuals() vs y - fitted()",
            max(abs(r - (dd$y - fitted(fit)))) / max(abs(dd$y))))
cat(sprintf("%-42s %.3e\n", "residuals('pearson') vs resp/sd",
            max(abs(residuals(fit, type = "pearson") -
                      r / sqrt(fitted(fit)^2 *
                                 (exp(sg^2) - 1)))) /
              max(abs(residuals(fit, type = "pearson")))))

cat("\ncoefficient tables, link scale or response scale:\n")
cf <- fixef(fit)
cat(sprintf("fixef()$mu['(Intercept)']       %.10f\n",
            cf$mu[["(Intercept)"]]))
cat(sprintf("fixef()$sigma['(Intercept)']    %.10f\n",
            cf$sigma[["(Intercept)"]]))
cat(sprintf("exp(that)                       %.10f\n",
            exp(cf$sigma[["(Intercept)"]])))
cat(sprintf("sigma(fit)                      %.10f\n", sg))
s <- summary(fit)
cat("\nsummary() coefficient blocks: ",
    paste(names(s$coefficients), collapse = ", "), "\n")
for (nm in names(s$coefficients)) {
  cat(" ", nm, "Estimate[1] =",
      sprintf("%.10f", s$coefficients[[nm]][1, 1]), "\n")
}
vc <- VarCorr(fit)
cat("\nVarCorr()[[1]] (variance of the mu linear predictor) ",
    sprintf("%.10f", vc[[1]][1, 1]), "\n")
cat("sqrt of it                                           ",
    sprintf("%.10f", sqrt(vc[[1]][1, 1])), "\n")

sim <- simulate(fit, nsim = 1, seed = 7)[[1]]
cat("\nsimulate()[1:3] (response scale) ",
    paste(sprintf("%.2f", sim[1:3]), collapse = " "), "\n")
cat("observed y[1:3]                  ",
    paste(sprintf("%.2f", dd$y[1:3]), collapse = " "), "\n")

# a poisson fit, where link and response also differ
dd2 <- data.frame(x = rnorm(300))
dd2$cnt <- rpois(300, exp(0.5 + 0.4 * dd2$x))
f2 <- frm(bf(cnt ~ x) + poisson(), data = dd2)
cat("\n== poisson fit ==\n")
cat(sprintf("predict()[1] (link)             %.10f\n", predict(f2)[1]))
cat(sprintf("fitted()[1]  (response)         %.10f\n", fitted(f2)[1]))
cat(sprintf("exp(link)[1]                    %.10f\n", exp(predict(f2)[1])))
cat(sprintf("sigma(f2) (no dispersion param) %.10f\n", sigma(f2)))
cat("DONE\n")

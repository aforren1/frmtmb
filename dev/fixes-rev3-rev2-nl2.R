# Reviewer of lane fixes, final check: a copy of dev/fixes-rev2-nl2.R that prints the flat warning.
# models in natural units where, at all three of the guard's fixed
# points, one coefficient's column is nonzero and the others underflow
# to exactly zero (a saturated body), so the guard sees a fixed flat
# subspace that the data do not have.
#   Rscript dev/fixes-rev2-nl2.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
run <- function(lab, expr) {
  r <- tryCatch({
    fit <- suppressMessages(withCallingHandlers(expr, warning = function(w) { if (grepl("not identified", conditionMessage(w))) cat("    FLAT WARNING
"); invokeRestart("muffleWarning") }))
    fe <- fixef(fit)
    sprintf("FIT conv=%d est %s se %s", fit$opt$convergence,
            paste(sprintf("%.5g", fe[, "Estimate"]), collapse = " "),
            paste(sprintf("%.3g", fe[, "Est.Error"]), collapse = " "))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 220)))
  cat(sprintf("%-46s %s\n", lab, r))
}
# 1. logistic growth over calendar years, positive scale on the log
set.seed(22)
yr <- seq(1900, 2000, length.out = 120)
d2 <- data.frame(yr = yr, y = 50 / (1 + exp((1950 - yr) / 12)) +
                   rnorm(120, 0, 1.5))
run("Asym / (1 + exp((xmid - yr) / exp(lscal)))",
    frm(bf(y ~ Asym / (1 + exp((xmid - yr) / exp(lscal))), Asym ~ 1,
           xmid ~ 1, lscal ~ 1, nl = TRUE), data = d2,
        start = list(beta = c(50, 1950, log(12)))))
cat("    nls SSlogis:", format(coef(nls(y ~ SSlogis(yr, Asym, xmid, scal),
                                       data = d2)), digits = 6), "\n")
# 2. psychometric function with a lapse rate, stimulus 200 to 400
set.seed(21)
x <- runif(800, 200, 400)
p <- 0.04 * 0.5 + 0.96 * pnorm(x, 300, 25)
d1 <- data.frame(x = x, y = rbinom(800, 1, p))
run("lapse * 0.5 + (1 - lapse) * pnorm(x, m0, exp(ls))",
    frm(bf(y ~ lapse * 0.5 + (1 - lapse) * pnorm(x, m0, exp(ls)),
           lapse ~ 1, m0 ~ 1, ls ~ 1, nl = TRUE),
        family = bernoulli(link = "identity"), data = d1,
        start = list(beta = c(0.05, 300, log(25)))))
# 3. Emax dose-response, dose in mg 100 to 1000, ED50 on the log
set.seed(25)
dose <- rep(c(100, 200, 400, 600, 800, 1000), each = 20)
d3 <- data.frame(dose = dose, y = 2 + 10 * dose^3 / (exp(log(400))^3 +
                                                       dose^3) +
                   rnorm(120, 0, 0.8))
run("e0 + emax * dose^h / (exp(led50)^h + dose^h)",
    frm(bf(y ~ e0 + emax * dose^3 / (exp(led50)^3 + dose^3), e0 ~ 1,
           emax ~ 1, led50 ~ 1, nl = TRUE), data = d3,
        start = list(beta = c(2, 10, log(400)))))

# Reviewer of lane optima, claim 1(d): what summary() and confint()
# report for a simplex after the change of coordinates.
#   Rscript dev/optima-rev-mo-report.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ms <- tryCatch(utils::getFromNamespace("mo_simplex", "frmtmb"),
               error = function(e) function(z) {
                 x <- exp(c(0, z))
                 x / sum(x)
               })
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
dat <- data.frame(income, ls)
# interior simplex (brms_monotonic fit1)
fit1 <- frm(bf(ls ~ mo(income)), data = dat, family = gaussian())
cat("== interior: simplex", format(ms(fit1$estimates$zeta1), digits = 4), "\n")
print(summary(fit1))
print(confint(fit1))
ci <- tryCatch(confint(fit1, parm = c("zeta1_1", "zeta1_2"),
                       method = "profile"), error = function(e) e)
print(ci)
# a face: steps 0, 1, 1
set.seed(2)
d2 <- data.frame(x = sample(0:3, 300, TRUE))
d2$y <- c(0, 0, 1, 2)[d2$x + 1] + rnorm(300)
fit2 <- suppressWarnings(frm(bf(y ~ mo(x)), data = d2, family = gaussian()))
cat("== face: simplex", format(ms(fit2$estimates$zeta1), digits = 4), "\n")
print(confint(fit2))
ci <- tryCatch(confint(fit2, parm = c("zeta1_1", "zeta1_2"),
                       method = "profile"), error = function(e) e)
print(ci)
# a delta-method SE of each weight, from vcov() of the coordinates
J <- function(z, f) {
  h <- 1e-6
  sapply(seq_along(z), function(j) {
    e <- replace(numeric(length(z)), j, h)
    (f(z + e) - f(z - e)) / (2 * h)
  })
}
for (ft in list(fit1, fit2)) {
  z <- ft$estimates$zeta1
  V <- vcov(ft, full = TRUE)
  nm <- grep("^zeta1", rownames(V))
  Jz <- J(z, ms)
  cat("delta-method SE of the simplex:",
      format(sqrt(diag(Jz %*% V[nm, nm] %*% t(Jz))), digits = 3), "\n")
}

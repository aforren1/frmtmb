# Lane fixes, punch round: m7's message and the hypothesis it points
# to, and m4's transformed predictor from the formula environment.
#   Rscript dev/fixes-p1-misc.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(emmeans)})
cat("LIB", find.package("frmtmb"), "\n")
set.seed(62)
n <- 200
d <- data.frame(x = rnorm(n))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
fc <- frm(y ~ x, family = cumulative(), data = d)
cat("m7:", tryCatch(rownames(confint(fc, parm = "Intercept[2]")),
                    error = function(e) conditionMessage(e)), "\n")
h <- tryCatch(hypothesis(fc, "Intercept[2] = 0", method = "profile"),
              error = function(e) conditionMessage(e))
if (is.character(h)) cat("hypothesis:", h, "\n") else
  print(h$hypothesis[, 1:5])
cat("fixef Intercept[2]:", fixef(fc)["Intercept[2]", 1:2], "\n")
# m4: zz is not a column of the data
set.seed(7101)
n <- 150
dd <- data.frame(f = factor(sample(c("a", "b", "c"), n, TRUE)))
zz <- runif(n, 0.5, 4)
dd$y <- rpois(n, exp(0.2 + 0.3 * zz - 0.05 * zz^2 + 0.2 * as.numeric(dd$f)))
fit <- frm(bf(y ~ poly(zz, 2) + f), data = dd, family = poisson())
ref <- glm(y ~ poly(zz, 2) + f, data = dd, family = poisson())
a <- tryCatch(summary(emmeans(fit, "f"))$emmean,
              error = function(e) conditionMessage(e))
b <- summary(emmeans(ref, "f"))$emmean
cat("m4 frmtmb:", a, "\nm4 glm   :", b, "\n")
a2 <- tryCatch(summary(emmeans(fit, ~ zz, at = list(zz = c(1, 3))))$emmean,
               error = function(e) conditionMessage(e))
b2 <- summary(emmeans(ref, ~ zz, at = list(zz = c(1, 3))))$emmean
cat("m4 at zz = 1, 3: frmtmb", a2, " glm", b2, "\n")

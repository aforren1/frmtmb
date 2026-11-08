# What brms, glmmTMB and lme4 call the penalized (log posterior at the
# mode) value, to choose frmtmb's accessor name.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
cat("glmmTMB", format(packageVersion("glmmTMB")), "\n")
print(body(glmmTMB:::logLik.glmmTMB))
pat <- "prior|penal|posterior|lp"
cat("glmmTMB ns:", grep(pat, ls(asNamespace("glmmTMB")), value = TRUE),
    "\n")
cat("brms", format(packageVersion("brms")), "exports:",
    grep(pat, getNamespaceExports("brms"), value = TRUE), "\n")
cat("lme4", format(packageVersion("lme4")), "exports:",
    grep("penal|crit|dev|pwrss", getNamespaceExports("lme4"),
         value = TRUE), "\n")
# glmmTMB with a prior: does logLik() hold the prior?
library(glmmTMB)
set.seed(4101)
n <- 200; ng <- 20
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng), each = n / ng)))
d$y <- 1 + 0.5 * d$x + rnorm(ng)[d$g] + rnorm(n)
m0 <- glmmTMB(y ~ x + (1 | g), data = d)
pr <- data.frame(prior = "normal(0, 100)", class = "fixef", coef = "x")
m1 <- glmmTMB(y ~ x + (1 | g), data = d, priors = pr)
cat(sprintf("glmmTMB logLik none %.10f  normal(0,100) on x %.10f\n",
            as.numeric(logLik(m0)), as.numeric(logLik(m1))))
cat("AIC(m1):", tryCatch(AIC(m1), warning = function(w) conditionMessage(w),
                         error = function(e) conditionMessage(e)), "\n")

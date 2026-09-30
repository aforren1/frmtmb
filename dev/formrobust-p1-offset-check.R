# Punch round 1: the reviewer's five offset cases (formrobust-rev-offset-fix.R) on the lane build, no mutant. Data seed 21.
# mutant (emm_terms keeps the offset attribute). Data seed 21.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(emmeans)})
       
set.seed(21)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 1, 3),
                f = factor(sample(c("a", "b"), n, TRUE)),
                z = runif(n, 0, 1))
d$y <- rpois(n, exp(0.3 + 0.4 * d$x + 0.2 * (d$f == "b")) * d$time)
d$y2 <- rnorm(n, 1 + 0.5 * d$x + d$z, 0.5)
d$y3 <- rnorm(n, 0.5 * d$x, exp(0.2 + log(d$time)))
mx <- mean(d$x); lmt <- log(mean(d$time))
fp <- frm(bf(y ~ x + f + offset(log(time))), data = d, family = poisson())
b <- fixef(fp)[, "Estimate"]
lin <- b[["Intercept"]] + b[["x"]] * mx + c(0, b[["fb"]])
em <- summary(emmeans(fp, ~ f))
cat("pois emmean - (lin + log(mean(time))):", em$emmean - (lin + lmt), "\n")
em1 <- summary(emmeans(fp, ~ f, at = list(time = 1)))
cat("pois emmean at time = 1 - lin:", em1$emmean - lin, "\n")
ep <- summary(emmeans(fp, ~ f, epred = TRUE))
cat("pois epred / exp(lin + lmt):", ep$emmean / exp(lin + lmt), "\n")
fs <- frm(bf(y3 ~ x, sigma ~ 1 + offset(log(time))), data = d)
es <- summary(emmeans(fs, ~ 1, dpar = "sigma"))
cat("sigma emmean - (b_sigma + lmt):",
    es$emmean - (fixef(fs)["sigma_Intercept", 1] + lmt), "\n")
fnl <- frm(bf(y2 ~ a + b * x, a ~ 1 + offset(z), b ~ 1, nl = TRUE),
           data = d)
bn <- fixef(fnl)[, "Estimate"]
en <- summary(emmeans(fnl, ~ 1))
cat("nl emmean - (a + b mx):", en$emmean - (bn[[1]] + bn[[2]] * mx), "\n")

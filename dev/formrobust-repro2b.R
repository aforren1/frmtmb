# Item 2: what the offset contributes to conditional_effects() and to
# emmeans(), against the coefficients. Seed 21, as formrobust-repro2.R.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(emmeans))
set.seed(21)
n <- 80
d <- data.frame(x = rnorm(n), time = runif(n, 1, 3), f = gl(2, 40),
                z = rnorm(n))
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x) * d$time)
f <- frm(bf(yc ~ x + f + offset(log(time))), data = d, family = poisson())
b <- fixef(f)[, "Estimate"]
ce <- conditional_effects(f, effects = "x")[[1]]
i <- 1
eta_no <- b[["Intercept"]] + b[["x"]] * ce$x[i]
cat("CE estimate / exp(eta + log(mean(time))) =",
    format(ce$estimate__[i] / exp(eta_no + log(mean(d$time))), digits = 15),
    "\n")
em <- summary(emmeans(f, ~ f))
lin_no <- b[["Intercept"]] + b[["x"]] * mean(d$x)
cat("emmean[f=1] - (b0 + bx*mean(x)) =",
    format(em$emmean[1] - lin_no, digits = 15), "\n")
cat("log(mean(time)) =", format(log(mean(d$time)), digits = 15), "\n")
em2 <- tryCatch(summary(emmeans(f, ~ f, epred = TRUE)),
                error = function(e) conditionMessage(e))
print(em2)
if (is.data.frame(em2)) {
  cat("epred emmean[f=1] / exp(b0 + bx*mean(x) + log(mean(time))) =",
      format(em2$emmean[1] / exp(lin_no + log(mean(d$time))), digits = 15),
      "\n")
}

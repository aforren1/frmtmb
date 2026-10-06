# Lane fixes, item 1: the grid route (epred = TRUE, re_formula = NULL,
# nonlinear) on a transformed predictor.
#   Rscript dev/fixes-emm4.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(emmeans)})
cat("LIB", find.package("frmtmb"), "\n")
set.seed(20260930)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n, 1, 2), time = runif(n, 1, 5),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(sample(1:8, n, TRUE)))
d$yc <- rpois(n, d$time * exp(0.3 + 0.4 * d$x + 0.2 * d$z))
show <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:",
                                                conditionMessage(e)))
  cat(sprintf("%-40s %s\n", lab,
              if (is.character(r)) r else
                paste(format(r, digits = 8), collapse = " ")))
}
for (fo in list(yc ~ poly(z, 2) + f, yc ~ scale(z) + f + offset(log(time)))) {
  cat("==", deparse(fo), "\n")
  fit <- frm(bf(fo), data = d, family = poisson())
  ref <- glm(fo, data = d, family = poisson())
  show("frm epred", summary(emmeans(fit, "f", epred = TRUE))$emmean)
  show("glm regrid response",
       summary(emmeans(regrid(ref_grid(ref)), "f"))[[2]])
  show("frm epred at z=c(0,3)",
       summary(emmeans(fit, "f", epred = TRUE, at = list(z = c(0, 3))))$emmean)
  show("glm regrid at z=c(0,3)",
       summary(emmeans(regrid(ref_grid(ref, at = list(z = c(0, 3)))),
                       "f"))[[2]])
  show("frm epred SE", summary(emmeans(fit, "f", epred = TRUE))$SE)
  show("glm regrid SE", summary(emmeans(regrid(ref_grid(ref)), "f"))$SE)
}
cat("== random intercept, re_formula = NULL\n")
fo <- yc ~ poly(z, 2) + f + (1 | g)
fit <- frm(bf(fo), data = d, family = poisson())
show("frm re_formula=NULL", summary(emmeans(fit, "f", re_formula = NULL))$emmean)
show("frm re_formula=NA", summary(emmeans(fit, "f"))$emmean)
show("frm grid z", emmeans::ref_grid(fit, re_formula = NULL)@grid$z[1])
show("mean(z)", mean(d$z))
cat("== nonlinear with a transformed covariate\n")
fit <- frm(bf(yc ~ exp(a + b * log(abs(z) + 1)), a ~ f, b ~ 1, nl = TRUE),
           data = d, family = poisson(link = "identity"))
show("frm nl a", summary(emmeans(fit, "f", nlpar = "a"))$emmean)
show("frm nl mu", summary(emmeans(fit, "f", dpar = "mu"))$emmean)
e <- fixef(fit)
cat("by hand mu:", format(exp(e["a_Intercept", 1] + c(0, e["a_fb", 1],
    e["a_fc", 1]) + e["b_Intercept", 1] * log(abs(mean(d$z)) + 1)),
    digits = 8), "\n")

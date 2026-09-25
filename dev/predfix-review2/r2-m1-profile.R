source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# profile confint on an engaged fit vs autoscale = FALSE; the design of
# r2-m1-downstream.R (seed 31, slope column sd 0.03)
set.seed(31)
ng <- 20; per <- 15; n <- ng * per
g <- factor(rep(seq_len(ng), each = per))
xs <- rnorm(n)
y <- 1 + 0.7 * xs + rnorm(ng, 0, 0.8)[g] + rnorm(ng, 0, 0.5)[g] * xs + rnorm(n)
d <- data.frame(y, x = xs * 0.03, g)
a <- frm(y ~ x + (1 + x | g), data = d)
b <- frm(y ~ x + (1 + x | g), data = d, control = frmtmb_control(autoscale = FALSE))
w <- rownames(confint(a))
print(w)
pp <- c("theta_2", "theta_3")

ca <- confint(a, parm = pp, method = "profile")
cb <- confint(b, parm = pp, method = "profile")
print(signif(ca, 7)); print(signif(cb, 7))
cat(sprintf("profile confint max rel diff: %.2e\n", rel(ca, cb)))

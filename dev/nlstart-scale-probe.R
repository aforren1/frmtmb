# Defects 2 and 3: is the cause the start or the scale? Seeds as named.
source("dev/nlstart-common.R")
set.seed(9901)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 20)))
d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4,
            ordered_result = TRUE)
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$xs <- d$x * 1e-10
d$ybig <- 1e6 + 5e4 * d$x + rnorm(n, 0, 1e3)
fm <- bf(ybig ~ a + b * x, a + b ~ 1, nl = TRUE)
cat("sd(ybig)", sd(d$ybig), "\n")
cat("default          :", ll(tt(frm(fm, data = d))), "\n")
cat("a at mean(y)     :", ll(tt(frm(fm, data = d, start = list(
  beta = c("a_(Intercept)" = mean(d$ybig)))))), "\n")
cat("a mean, sigma sd :", ll(tt(frm(fm, data = d, start = list(
  beta = c("a_(Intercept)" = mean(d$ybig)),
  betad = log(sd(d$ybig)))))), "\n")
cat("a, b at LS       :", ll(tt(frm(fm, data = d, start = list(
  beta = c("a_(Intercept)" = 1e6, "b_(Intercept)" = 5e4))))), "\n")
# the same response in units of 1e6: is the body fine at unit scale?
d$ys <- d$ybig / 1e6
cat("y / 1e6          :", ll(tt(frm(bf(ys ~ a + b * x, a + b ~ 1, nl = TRUE),
                                    data = d))),
    " +n*log(1e6) =", format(as.numeric(logLik(frm(bf(ys ~ a + b * x,
      a + b ~ 1, nl = TRUE), data = d))) - n * log(1e6), digits = 10), "\n")
# what does the fit start from?
f0 <- frm(fm, data = d, control = frmtmb_control(optCtrl = list(iter.max = 0)))
print(f0$estimates[c("beta", "betad")])
# defect 2
set.seed(6201); n <- 200
d <- data.frame(x = rnorm(n), g = factor(rep(1:10, each = 20)),
                t = rep(1:20, 10), z = rnorm(n))
d$y <- 0.5 * d$x + rnorm(n)
d$xt <- d$x * 1e-20
cat("tiny: rescaled body a * (xt * 1e20) / 1e20 :",
    ll(tt(frm(bf(y ~ a * (xt * 1e20), a ~ 1, nl = TRUE), data = d))), "\n")

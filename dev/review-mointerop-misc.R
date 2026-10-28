## Reviewer: checks of mechanisms behind findings.
##   Rscript dev/review-mointerop-misc.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
## 1. Gamma(log) HC0: observed-information bread reproduces frmtmb's
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(sample(1:30, n, TRUE)),
                w = runif(n, 0.5, 2), off = log(runif(n, 1, 3)))
d$ypos <- rgamma(n, shape = 2, rate = 2 / exp(0.5 + 0.3 * d$x))
if (identical(arm, "lane")) {
  mg <- glm(ypos ~ x, data = d, family = Gamma(link = "log"))
  X <- model.matrix(mg)
  mu <- fitted(mg)
  r <- d$ypos / mu - 1
  S <- X * r                     # score up to the shape factor
  Hobs <- crossprod(X * sqrt(d$ypos / mu))
  Hexp <- crossprod(X)
  so <- solve(Hobs) %*% crossprod(S) %*% solve(Hobs)
  se <- solve(Hexp) %*% crossprod(S) %*% solve(Hexp)
  fg <- frm(ypos ~ x, data = d, family = Gamma(link = "log"))
  show("Gamma HC0: observed bread, expected bread, frmtmb", rbind(
    observed = sqrt(diag(so)), expected = sqrt(diag(se)),
    frmtmb = sqrt(diag(get_vcov(fg, vcov = "HC0")))[1:2],
    sandwich = sqrt(diag(sandwich::vcovHC(mg, type = "HC0")))))
}
## 2. mo() with 3 categories: HC0 estimates
set.seed(32)
dm <- data.frame(k = sample(1:3, n, TRUE))
dm$y <- 1 + c(0, 0.5, 0.7)[dm$k] + rnorm(n) * dm$k
fm <- frm(y ~ mo(k), data = dm)
show("cell means", tapply(dm$y, dm$k, mean))
show("frm fitted at 1:3", fitted(fm, newdata = data.frame(k = 1:3)))
show("summary simplex", summary(fm))
## 3. a gap fit: posterior predictive predict() and drop1()
set.seed(7)
d6 <- data.frame(x = sample(c(1, 2, 5), 200, TRUE), z = rnorm(200))
d6$y <- 1 + 0.8 * d6$x + 0.3 * d6$z + rnorm(200)
f6 <- suppressWarnings(frm(y ~ mo(x) + z, data = d6))
show("gap fit predict() newdata", predict(f6, newdata = data.frame(
  x = c(1, 5), z = 0), ndraws = 2000))
show("gap fit drop1", drop1(f6))
set.seed(8)
d7 <- data.frame(x = sample(1:4, 200, TRUE), z = rnorm(200))
d7$y <- 1 + 0.8 * d7$x + 0.3 * d7$z + rnorm(200)
f7 <- frm(y ~ mo(x) + z, data = d7)
show("4-category fit drop1", drop1(f7))
show("4-category fit predict()", predict(f7, newdata = data.frame(
  x = c(1, 4), z = 0), ndraws = 2000))
## 4. the RE gap fit fitted() at the empty category
set.seed(23)
n <- 300
d3 <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                 g = factor(sample(1:10, n, TRUE)))
d3$y <- 1 + (0.4 + rnorm(10, 0, 0.2)[d3$g]) * d3$x + 0.5 * d3$z + rnorm(n)
fr <- suppressWarnings(frm(y ~ mo(x) + z + (mo(x) | g), data = d3))
nd4 <- data.frame(x = c(1, 2, 3, 5), z = 0, g = factor(1, levels = 1:10))
show("RE gap fit fitted", fitted(fr, newdata = nd4)[, 1:2])
show("RE gap fit fitted NA", fitted(fr, newdata = nd4, re_formula = NA)[, 1:2])
show("RE gap fit get_vcov", get_vcov(fr))
show("RE gap fit get_coef", get_coef(fr))
## 5. ordered factor with unused levels: D and codes
d4 <- data.frame(xo = factor(c(2, 3, 4, 2, 3, 4), levels = 1:4,
                             ordered = TRUE), y = c(1, 2, 3, 1.5, 2.5, 3.2))
f4 <- show("ordered, empty level 1", frm(y ~ mo(xo), data = d4))
show("D", f4$frame$linpreds[[1]]$mo[[1]][c("D", "codes")])

## Reviewer: gap-aware marginaleffects SEs on cases the lane did not test.
##   Rscript dev/review-mointerop-gap.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
ml_se <- function(m, nd) {
  s <- summary(m)
  n <- length(resid(m))
  predict(m, newdata = nd, se.fit = TRUE)$se.fit *
    sqrt(sum(resid(m)^2) / n) / s$sigma
}
pse <- function(fit, nd, ...) {
  as.data.frame(predictions(fit, newdata = nd, ...))[
    , c("estimate", "std.error")]
}

## Case 1: two mo() terms, each with a gap (seed 21)
set.seed(21)
n <- 300
d <- data.frame(x = sample(c(1, 2, 5), n, TRUE),
                w = sample(c(1, 3, 4), n, TRUE), z = rnorm(n))
d$y <- 1 + 0.4 * d$x + 0.3 * d$w + 0.5 * d$z + rnorm(n)
f <- show("two gap terms fit", frm(y ~ mo(x) + mo(w) + z, data = d))
show("get_coef", get_coef(f))
nd <- expand.grid(x = c(1, 2, 5), w = c(1, 3, 4), z = 0)
m <- lm(y ~ factor(x) + factor(w) + z, data = d)
p <- show("predictions at observed", pse(f, nd))
show("fitted SE at observed", fitted(f, newdata = nd)[, "Est.Error"])
show("lm (additive factors) SE at observed", ml_se(m, nd))
show("logLik frm vs lm(factor + factor)",
     c(frm = as.numeric(logLik(f)), lm = as.numeric(logLik(m))))
nd2 <- data.frame(x = c(1, 3, 4, 5), w = 2, z = 0)
show("predictions with w = 2 (gap of w) and x gaps", pse(f, nd2))
show("avg_comparisons x 1 to 5", avg_comparisons(
  f, variables = list(x = c(1, 5))))
show("lm comparisons x 1 to 5", avg_comparisons(
  m, variables = list(x = c(1, 5))))
show("avg_slopes z", avg_slopes(f, variables = "z"))
show("lm avg_slopes z", avg_slopes(m, variables = "z"))

## set_coef round trip: get_coef -> set_coef -> get_coef
show("set_coef round trip", {
  b <- get_coef(f)
  f2 <- set_coef(f, b)
  p1 <- get_predict(f, newdata = nd)$estimate
  p2 <- get_predict(f2, newdata = nd)$estimate
  list(max_abs_coef = max(abs(get_coef(f2) - b)),
       max_rel_pred = max(abs(p1 - p2) / abs(p1)))
})
show("set_coef perturb cum weight moves predictions only at x >= that cat",
     {
       b <- get_coef(f)
       j <- grep("_cum", names(b))[1]
       b[j] <- b[j] + 1e-4
       f2 <- set_coef(f, b)
       ndx <- data.frame(x = c(1, 2, 3, 4, 5), w = 1, z = 0)
       get_predict(f2, newdata = ndx)$estimate -
         get_predict(f, newdata = ndx)$estimate
     })

## Case 2: mo(x) * z interaction with a gap (seed 22)
set.seed(22)
d2 <- data.frame(x = sample(c(1, 2, 3, 6), n, TRUE), z = rnorm(n))
d2$y <- 1 + 0.4 * d2$x + 0.5 * d2$z + 0.2 * d2$x * d2$z + rnorm(n)
fi <- show("mo(x) * z with a gap", frm(y ~ mo(x) * z, data = d2))
show("get_coef", names(get_coef(fi)))
nd3 <- expand.grid(x = c(1, 2, 3, 4, 6), z = c(-1, 1))
show("predictions", pse(fi, nd3))
show("fitted SE", fitted(fi, newdata = nd3)[, "Est.Error"])
show("avg_slopes z by x", slopes(fi, variables = "z",
                                 newdata = data.frame(x = c(1, 4, 6), z = 0)))

## Case 3: (mo(x) | g) with a gap (seed 23)
set.seed(23)
d3 <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                 g = factor(sample(1:10, n, TRUE)))
d3$y <- 1 + (0.4 + rnorm(10, 0, 0.2)[d3$g]) * d3$x + 0.5 * d3$z + rnorm(n)
fr <- show("(mo(x) | g) with a gap", frm(y ~ mo(x) + z + (mo(x) | g),
                                         data = d3))
show("get_coef", names(get_coef(fr)))
nd4 <- data.frame(x = c(1, 2, 3, 5), z = 0, g = factor(1, levels = 1:10))
show("predictions conditional", pse(fr, nd4))
show("fitted conditional", fitted(fr, newdata = nd4)[, c("Estimate",
                                                         "Est.Error")])
show("predictions re_formula = NA", pse(fr, nd4, re_formula = NA))
show("fitted re_formula = NA", fitted(fr, newdata = nd4, re_formula = NA)[
  , c("Estimate", "Est.Error")])

## Case 4: ordered factor with an unused TOP level (seed 24)
set.seed(24)
d4 <- data.frame(x = sample(1:3, n, TRUE), z = rnorm(n))
d4$y <- 1 + 0.4 * d4$x + 0.5 * d4$z + rnorm(n)
d4$xo <- factor(d4$x, levels = 1:4, ordered = TRUE)
ft <- show("ordered factor, empty top level 4", frm(y ~ mo(xo) + z, data = d4))
show("summary", summary(ft))
nd5 <- data.frame(xo = factor(1:4, levels = 1:4, ordered = TRUE), z = 0)
show("fitted", fitted(ft, newdata = nd5)[, c("Estimate", "Est.Error")])
show("predictions", pse(ft, nd5))

## Case 5: ordered factor with an unused BOTTOM level (seed 25)
d5 <- d4
d5$xo <- factor(d5$x + 1L, levels = 1:4, ordered = TRUE)
fb <- show("ordered factor, empty bottom level 1",
           frm(y ~ mo(xo) + z, data = d5))
show("summary", summary(fb))
show("fitted", fitted(fb, newdata = nd5)[, c("Estimate", "Est.Error")])
show("predictions", pse(fb, nd5))

## Case 6: warning count over a predictions() call on a gap fit
set.seed(7)
d6 <- data.frame(x = sample(c(1, 2, 5), 200, TRUE))
d6$y <- 1 + 0.8 * d6$x + rnorm(200)
f6 <- frm(y ~ mo(x), data = d6)
show("one call with gap rows: warnings", pse(f6, data.frame(x = 1:5)))
show("second call on same fit: warnings", pse(f6, data.frame(x = 1:5)))
f6b <- frm(y ~ mo(x), data = d6)
show("observed rows only on a fresh fit: warnings",
     pse(f6b, data.frame(x = c(1, 2, 5))))
show("hypotheses on the gap fit", hypotheses(f6b, "b2 - b1 = 0"))
show("get_vcov on gap fit", get_vcov(f6b))
show("HC0 on gap fit", pse(f6b, data.frame(x = c(1, 2, 5)), vcov = "HC0"))
show("lm HC0 factor", {
  m6 <- lm(y ~ factor(x), data = d6)
  as.data.frame(predictions(m6, newdata = data.frame(x = c(1, 2, 5)),
                            vcov = "HC0"))[, c("estimate", "std.error")]
})
show("vcov(full = TRUE) on gap fit", pse(f6b, data.frame(x = c(1, 2, 5)),
                                         vcov = vcov(f6b, full = TRUE)))
show("2 * get_vcov on gap fit", pse(f6b, data.frame(x = c(1, 2, 5)),
                                    vcov = 2 * get_vcov(f6b)))

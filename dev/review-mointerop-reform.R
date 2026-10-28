## Reviewer: re_formula and allow_new_levels through marginaleffects.
##   Rscript dev/review-mointerop-reform.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
set.seed(41)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(sample(1:10, n, TRUE)),
                h = factor(sample(1:6, n, TRUE)))
u0 <- rnorm(10, 0, 0.8); u1 <- rnorm(10, 0, 0.5); v0 <- rnorm(6, 0, 0.6)
d$y <- 1 + (0.5 + u1[d$g]) * d$x + u0[d$g] + v0[d$h] + rnorm(n)
d$cnt <- rpois(n, exp(0.2 + (0.3 + u1[d$g] * 0.3) * d$x + u0[d$g] * 0.3))
f <- frm(y ~ x + (1 + x | g) + (1 | h), data = d)
nd <- data.frame(x = c(-1, 0, 1), g = factor(c(1, 2, 3), levels = 1:10),
                 h = factor(c(1, 1, 2), levels = 1:6))
q <- function(expr) suppressWarnings(expr)
show("fitted re_formula = NA", fitted(f, newdata = nd, re_formula = NA)[, 1])
show("predictions re_formula = NA", q(predictions(f, newdata = nd,
                                                  re_formula = NA))$estimate)
show("fitted re_formula = ~ (1 | h)", fitted(f, newdata = nd,
                                            re_formula = ~ (1 | h))[, 1])
show("predictions re_formula = ~ (1 | h)", q(predictions(
  f, newdata = nd, re_formula = ~ (1 | h)))$estimate)
show("fitted re_formula = ~ (1 | g)", fitted(f, newdata = nd,
                                            re_formula = ~ (1 | g))[, 1])
show("predictions re_formula = ~ (1 | g)", q(predictions(
  f, newdata = nd, re_formula = ~ (1 | g)))$estimate)
## slopes of x: conditional uses (0.5 + u1) per group, NA uses the fixed
show("fixef x", fixef(f)["x", 1])
show("ranef g x first 3", ranef(f)$g[1:3, , "x"])
show("slopes(re_formula = NA)", q(slopes(f, variables = "x", newdata = nd,
                                          re_formula = NA))[
                                            , c("estimate", "std.error")])
show("slopes(conditional)", slopes(f, variables = "x", newdata = nd)[
  , c("estimate", "std.error")])
show("avg_comparisons(re_formula = NA)", q(avg_comparisons(
  f, variables = list(x = 1), re_formula = NA)))
show("avg_comparisons(conditional)", avg_comparisons(
  f, variables = list(x = 1)))
## new levels
ndn <- data.frame(x = c(-1, 1), g = factor(c("new1", "new2")),
                  h = factor(c(1, 1), levels = 1:6))
show("fitted new levels, no allow", fitted(f, newdata = ndn)[, 1])
show("predictions new levels, no allow", predictions(f, newdata = ndn)$estimate)
show("fitted new levels allow", fitted(f, newdata = ndn,
                                       allow_new_levels = TRUE)[, 1:2])
show("predictions new levels allow", q(predictions(
  f, newdata = ndn, allow_new_levels = TRUE))[, c("estimate", "std.error")])
show("fitted new levels allow, re NA", fitted(
  f, newdata = ndn, allow_new_levels = TRUE, re_formula = NA)[, 1])
## poisson, response scale: re_formula NA vs fitted
fp <- frm(cnt ~ x + (1 + x | g), data = d, family = poisson())
show("poisson fitted NA", fitted(fp, newdata = nd, re_formula = NA)[, 1:2])
show("poisson predictions NA", q(predictions(fp, newdata = nd,
                                             re_formula = NA))[
                                               , c("estimate", "std.error")])
show("poisson predictions NA type link", q(predictions(
  fp, newdata = nd, re_formula = NA, type = "link"))[
    , c("estimate", "std.error")])
show("poisson fitted NA scale linear", fitted(fp, newdata = nd,
                                              re_formula = NA,
                                              scale = "linear")[, 1:2])
## avg_predictions over the data, re_formula = NA
show("avg_predictions re NA", q(avg_predictions(fp, re_formula = NA)))
show("mean fitted re NA", mean(fitted(fp, re_formula = NA)[, 1]))
## a misspelled argument: is it silently ignored?
show("predictions re.form = NA (lme4 spelling)", q(predictions(
  f, newdata = nd, re.form = NA))$estimate)

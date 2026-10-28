## Reviewer: gap simplex with ordinal extras, a gap in sigma's formula,
## a one-value term in sigma's formula, and a categorical outcome with
## padded new data.
##   Rscript dev/review-mointerop-ord.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
set.seed(61)
n <- 400
d <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                g = factor(sample(1:6, n, TRUE)))
eta <- 0.4 * d$x + 0.5 * d$z
d$yo <- cut(eta + rlogis(n), c(-Inf, 0.5, 1.5, 2.5, Inf),
            labels = FALSE)
fo <- suppressWarnings(frm(yo ~ mo(x) + z, data = d, family = cumulative()))
show("ordinal gap get_coef", names(get_coef(fo)))
nd <- data.frame(x = c(1, 2, 3, 5), z = 0)
p <- show("ordinal gap predictions", as.data.frame(predictions(
  fo, newdata = nd))[, c("rowid", "group", "estimate", "std.error")])
f_se <- show("ordinal fitted SE (category 1)", fitted(fo, newdata = nd)[
  , "Est.Error", 1])
## gap in sigma's formula
d$y <- 1 + 0.5 * d$z + rnorm(n) * exp(0.1 * d$x)
fs <- suppressWarnings(frm(bf(y ~ z, sigma ~ mo(x)), data = d))
show("sigma gap get_coef", names(get_coef(fs)))
show("sigma gap predictions dpar default (mu)", as.data.frame(predictions(
  fs, newdata = nd))[, c("estimate", "std.error")])
## one value in sigma's formula
d$x1 <- 3L
f1 <- suppressWarnings(frm(bf(y ~ z, sigma ~ mo(x1) + x), data = d))
show("sigma one-value get_vcov diag", diag(get_vcov(f1)))
show("sigma one-value predictions", as.data.frame(predictions(
  f1, newdata = data.frame(z = c(0, 1), x = 1, x1 = 3)))[
    , c("estimate", "std.error")])
## categorical outcome, padded new data
d$yc <- factor(sample(c("a", "b", "c"), n, TRUE))
fc <- frm(yc ~ z + g, data = d, family = categorical())
ndc <- data.frame(z = c(-1, 1), g = factor(c(1, 2), levels = 1:6))
pc <- show("categorical predictions on padded new data", as.data.frame(
  predictions(fc, newdata = ndc))[, c("rowid", "group", "estimate")])
show("categorical fitted", fitted(fc, newdata = ndc)[, "Estimate", ])

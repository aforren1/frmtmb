# Short check, item 2 probe: parameter labels of the multivariate design
source("dev/reviews/arspat2-short-common.R")
n <- 40
M <- 0.6^abs(outer(1:n, 1:n, "-"))
set.seed(7)
d <- data.frame(x = rnorm(n))
d$y <- 1 + 0.5 * d$x + as.vector(t(chol(M)) %*% rnorm(n))
d$z <- 0.3 + 0.8 * d$y - 0.4 * d$x + rnorm(n, 0, 0.7)
d$y[c(3, 11, 19, 27, 6, 14, 33)] <- NA
d$z[c(6, 14, 33, 9, 22)] <- NA
fit <- quiet(frm(bf(y | mi() ~ x + fcor(M)) + bf(z | mi() ~ mi(y) + x) +
                   set_rescor(FALSE), data = d, data2 = list(M = M),
                 family = gaussian()))
print(names(fit$opt$par))
print(ns$outer_par_labels(fit))
print(fit$opt$par)
str(lapply(fit$frame$mi_map, function(m) m[c("rows", "unread", "idx")]))
print(fit$estimates$miss)
print(summary(fit))

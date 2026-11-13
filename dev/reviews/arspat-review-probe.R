source("dev/reviews/arspat-review-common.R")
set.seed(3)
n <- 30
d <- data.frame(x = rnorm(n))
M <- 0.5^abs(outer(1:n, 1:n, "-"))
d$y <- 1 + 0.5 * d$x + as.vector(t(chol(M)) %*% rnorm(n))
d$y[c(4, 9)] <- NA
f <- frm(y | mi() ~ x + fcor(M), data = d, data2 = list(M = M),
         family = student())
print(f$opt$par)
print(ns$outer_par_labels(f))
print(ns$outer_par_map(f))
print(fixef(f))

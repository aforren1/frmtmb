# Lane illse: the category route on a cumulative degree-8 raw
# polynomial (seed 1, the test's data). P(y = 1) SE at x = 1.2, 1.5, 1.8
# from fitted() (ord_prob_se(), a coordinate Jacobian times the factor),
# from differences of P along each column of the factor (fit_fd_se()'s
# construction), and from the poly(x, 8) fit of the same likelihood.
#   Rscript dev/illse-cumdiag.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
set.seed(1)
d <- data.frame(x = stats::runif(200, 1, 2))
d$y <- factor(cut(2 * sin(3 * d$x) + stats::rlogis(200),
                  c(-Inf, -1, 0, 1, Inf), labels = FALSE), ordered = TRUE)
for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
f <- frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, data = d,
         family = cumulative())
g <- frm(y ~ poly(x, 8), data = d, family = cumulative())
cat("ll gap", as.numeric(logLik(g) - logLik(f)), "\n")
nd <- data.frame(x = c(1.2, 1.5, 1.8))
for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
a <- fitted(f, newdata = nd)[, "Est.Error", ]
b <- fitted(g, newdata = nd)[, "Est.Error", ]
ns <- asNamespace("frmtmb")
root <- ns$sdr_of(f)$cov_fixed_root
map <- ns$outer_par_map(f)
v0 <- ns$fit_outer_vector(f, map)
pf <- function(fit) fitted(fit, newdata = nd)[, "Estimate", ]
s <- 1e-3
dd <- vapply(seq_len(ncol(root)), function(k) {
  (pf(ns$fit_set_outer(f, v0 + s * root[, k], map)) -
     pf(ns$fit_set_outer(f, v0 - s * root[, k], map))) / (2 * s)
}, matrix(0, 3, 4))
c3 <- sqrt(apply(dd^2, c(1, 2), sum))
cat("fitted() over poly():\n"); print(a / b)
cat("along the factor over poly():\n"); print(c3 / b)

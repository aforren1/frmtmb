# Review of lane sparsese: the "cumulative" design of
# sparsese-review-designs.R, fitted() of new rows on both routes, entry by
# entry (that script's fitted_all read a relative difference of 1).
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(1)
ng <- 30L
n <- 300L
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)),
                h = factor(sample(rep(seq_len(12), length.out = n))),
                f = factor(rep(c("a", "b"), each = 1, length.out = n)))
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
ug <- stats::rnorm(ng, 0, 0.5)
invisible(stats::rnorm(ng, 0, 0.4))
invisible(stats::rnorm(12, 0, 0.3))
eta <- sin(3 * d$x) + ug[d$g] + 0.3 * (d$f == "b")
z <- 0.8 * eta + stats::rlogis(n)
d$y <- factor(cut(z, c(-Inf, 0, 0.8, Inf), labels = c("lo", "mid", "hi")),
              ordered = TRUE)
fo <- y ~ x1 + x2 + x3 + x4 + x5 + (1 | g)
nd <- d[c(1:6, 151:156), ]
out <- list()
for (lim in c(1e6, 20)) {
  options(frmtmb.joint_root_max = lim)
  fit <- q(frm(fo, data = d, family = cumulative(), REML = TRUE))
  invisible(q(frmtmb:::sdr_of(fit)))
  set.seed(3)
  out[[as.character(lim)]] <- q(fitted(fit, newdata = nd))
  set.seed(3)
  out[[paste0("i", lim)]] <- q(fitted(fit))
}
a <- out[["1e+06"]]
b <- out[["20"]]
cat("dims:", dim(a), " dimnames[[2]]:", dimnames(a)[[2]], "\n")
for (s in dimnames(a)[[2]]) {
  x <- a[, s, ]
  y <- b[, s, ]
  cat(sprintf("%-10s max rel %.3e  max abs %.3e\n", s,
              max(abs(x - y) / pmax(abs(x), abs(y), 1e-300)),
              max(abs(x - y))))
}
ai <- out[["i1e+06"]]
bi <- out[["i20"]]
cat("in-sample Est.Error max rel:",
    format(max(abs(ai[, "Est.Error", ] - bi[, "Est.Error", ]) /
                 abs(ai[, "Est.Error", ])), digits = 3), "\n")
cat("dense Est.Error rows 1-3:\n"); print(round(a[1:3, "Est.Error", ], 5))
cat("sparse Est.Error rows 1-3:\n"); print(round(b[1:3, "Est.Error", ], 5))
cat("in-sample Est.Error quantiles, dense:",
    format(stats::quantile(ai[, "Est.Error", ]), digits = 3), "\n")
cat("in-sample Est.Error quantiles, sparse:",
    format(stats::quantile(bi[, "Est.Error", ]), digits = 3), "\n")

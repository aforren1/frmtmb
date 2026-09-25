# Item 2, part 3: does a GENUINE heavy tail cross the one1 > 1/2 line?
# lognormal with sigma 5.5 and 6.5, well identified, row at x = 0.
#   PREDFIX_ARM=base Rscript dev/predfix-absurd3.R > dev/predfix-log/absurd3-base.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
one_share <- function(v) {
  v <- v[is.finite(v)]
  a <- abs(v - stats::median(v))
  max(a) / max(sum(a), .Machine$double.xmin)
}
for (s in c(4.5, 5.5, 6.5)) {
  o <- vapply(1:20, function(seed) {
    set.seed(seed)
    x <- rnorm(200)
    d <- data.frame(x = x, y = exp(rnorm(200, 1 + 0.3 * x, s)))
    f <- suppressWarnings(frm(bf(y ~ x) + lognormal(), data = d))
    set.seed(99)
    dr <- predict(f, newdata = data.frame(x = 0), ndraws = 1000,
                  summary = FALSE)
    c(one_share(dr[, 1]), sigma(f)[1])
  }, c(0, 0))
  cat(sprintf("sigma %.1f: one1 min %.3f median %.3f max %.3f, over 1/2 in %d of 20; fitted sigma %.2f to %.2f\n",
              s, min(o[1, ]), median(o[1, ]), max(o[1, ]), sum(o[1, ] > 0.5),
              min(o[2, ]), max(o[2, ])))
}

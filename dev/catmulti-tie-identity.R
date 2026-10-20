# Lane catmulti, the cumulative tied-threshold fix: in ONE process (so no
# core-type ulp difference, rtmb-pitfalls item 21), the 0.71.0 body of
# ord_cumulative_logpmf() and the lane's, taped and on doubles, at
# random points with distinct thresholds, for every ordinal link with a
# log-odds form. A fit reads only the taped value and gradient, so
# bitwise equality of both there says no fit moves.
# Output: dev/catmulti-out/tie-identity.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
# 0.71.0's body, verbatim but for the namespace prefixes
old <- function(y, eta, tau, lk, disc = 1) {
  q <- lk[["logit_eta"]]
  K1 <- length(tau)
  K <- K1 + 1L
  iK <- as.numeric(y == K)
  i1 <- as.numeric(y == 1)
  out <- i1 * ns$log_inv_logit(q(disc * (tau[1] - eta))) +
    iK * ns$log1m_inv_logit(q(disc * (tau[K1] - eta)))
  ym <- pmin(pmax(y, 2L), K1)
  a <- q(disc * (tau[ym] - eta))
  b <- q(disc * (tau[ym - 1L] - eta))
  out + (1 - i1 - iK) * ns$ord_log_interior(a, b, FALSE)
}
new <- function(y, eta, tau, lk, disc = 1) {
  ns$ord_cumulative_logpmf(y, eta, tau, lk, disc)
}
links <- c("logit", "probit", "probit_approx", "cloglog", "cauchit")
set.seed(20261014)
n <- 40
K1 <- 4L
tot <- 0L
same_v <- 0L
same_g <- 0L
same_d <- 0L
for (lnk in links) {
  lk <- ns$ord_link(lnk, "cumulative")
  if (is.null(lk[["logit_eta"]])) {
    cat(lnk, ": no log-odds form, the path the fix does not touch\n")
    next
  }
  for (rep in 1:50) {
    y <- sample.int(K1 + 1L, n, TRUE)
    x <- stats::rnorm(n)
    f_of <- function(fun) function(p) {
      "c" <- RTMB::ADoverload("c")
      tau <- cumsum(c(p[1], exp(p[2:K1])))
      eta <- p[K1 + 1L] * x
      sum(fun(y, eta, tau, lk, exp(p[K1 + 2L])))
    }
    p0 <- c(stats::rnorm(1), stats::rnorm(K1 - 1L, 0, 0.7),
            stats::rnorm(1), stats::rnorm(1, 0, 0.3))
    to <- RTMB::MakeTape(f_of(old), p0)
    tn <- RTMB::MakeTape(f_of(new), p0)
    tot <- tot + 1L
    same_v <- same_v + identical(to(p0), tn(p0))
    same_g <- same_g + identical(to$jacobian(p0), tn$jacobian(p0))
    same_d <- same_d + identical(f_of(old)(p0), f_of(new)(p0))
  }
  cat(sprintf("%s: %d points\n", lnk, 50L))
}
cat(sprintf("taped value identical %d of %d; taped gradient identical %d of %d; value on doubles identical %d of %d\n",
            same_v, tot, same_g, tot, same_d, tot))

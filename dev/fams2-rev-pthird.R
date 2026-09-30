# Reviewer, claim 2 aside: do other taped CDFs share RTMB::pbeta()'s
# non-finite third derivatives? RTMB::ppois() is poisson()'s lcdf (cens()
# and trunc() rows); RTMB::pgamma(), pbinom(), pnbinom() and pbeta() are
# what a nonlinear formula body gets through nl_rtmb_shadow. Seed 7.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
set.seed(7)
third_ok <- function(f, x0, pts) {
  F <- MakeTape(f, x0)
  T3 <- F$jacfun()$jacfun()$jacfun()
  ok <- vapply(seq_len(nrow(pts)), function(i) {
    v <- F(pts[i, ]); if (!is.finite(v)) return(NA)
    all(is.finite(T3(pts[i, ])))
  }, NA)
  c(points = sum(!is.na(ok)), nonfinite_third = sum(!ok, na.rm = TRUE))
}
N <- 2000
cat("ppois(q, exp(p)) in log lambda, q in 0..30:\n")
for (q in c(0, 3, 10, 30)) {
  print(c(q = q, third_ok(function(p) log(RTMB::ppois(q, exp(p[1]))), 0,
                          matrix(runif(N / 4, -4, 4), ncol = 1))))
}
cat("pgamma(x, shape = exp(p1), scale = exp(p2)):\n")
print(third_ok(function(p) log(RTMB::pgamma(p[3], shape = exp(p[1]), scale = exp(p[2]))),
               c(0, 0, 1), cbind(runif(N, -3, 4), runif(N, -3, 3), runif(N, 0.01, 20))))
cat("pbeta(x, exp(p1), exp(p2)):\n")
print(third_ok(function(p) log(RTMB::pbeta(p[3], exp(p[1]), exp(p[2]))),
               c(0, 0, 0.3), cbind(runif(N, -4, 5), runif(N, -4, 5), runif(N, 0.001, 0.999))))
cat("pnbinom(q = 4, size = exp(p1), mu = exp(p2)):\n")
print(third_ok(function(p) log(RTMB::pnbinom(4, size = exp(p[1]), prob = plogis(p[2]))),
               c(0, 0), cbind(runif(N, -3, 4), runif(N, -3, 4))))
cat("pbinom(q = 4, size = 12, prob = plogis(p1)):\n")
print(third_ok(function(p) log(RTMB::pbinom(4, 12, plogis(p[1]))),
               0, matrix(runif(N / 4, -6, 6), ncol = 1)))

# Review of lane nlstop2, attack 1: two coefficients in one denominator,
# from given starts (the default start fails on both arms).
# REV_ARM=before|after Rscript dev/reviews/nlstop2-review-pole2.R
source("dev/reviews/nlstop2-review-common.R")
mm_data <- function(seed) {
  set.seed(seed * 1000 + 1)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  lat <- 3 * d$xp / (1 + d$xp)
  d$nt <- 10L
  d$y <- rbinom(n, 10, plogis(lat - 2))
  d
}
for (s in c(11, 44)) {
  d <- mm_data(s)
  for (st in list(c(-0.1, -0.2, 1), c(-0.1, -2.5, 1), c(1, 1, 1))) {
    r <- tt(frm(bf(y | trials(nt) ~ a * xp / (b + c * xp), a + b + c ~ 1,
                   nl = TRUE), family = binomial(), data = d,
                start = list(beta = st)))
    f <- r$v
    if (inherits(f, "tterr")) {
      cat("ROW two-coef", s, paste(st, collapse = ","), "ERROR",
          substr(f$msg, 1, 120), "\n")
      next
    }
    cat("ROW two-coef", s, paste(st, collapse = ","), "ll",
        format(llv(f), digits = 12), "code", f$opt$convergence, "beta",
        paste(signif(f$estimates$beta, 6), collapse = " "), "warn",
        paste(unique(r$w), collapse = "+"), "msg",
        paste(unique(r$m), collapse = "+"), "\n")
  }
}
cat("DONE\n")

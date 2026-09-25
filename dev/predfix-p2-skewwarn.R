# Where does the default warn LESS than autoscale = FALSE on the skew
# set? Expected: only where FALSE's warning is about its own stall and
# the default's fit is a verified optimum with a higher logLik.
#   PREDFIX_ARM=lane Rscript dev/predfix-p2-skewwarn.R > dev/predfix-log/p2-skewwarn.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
mk <- function(seed, s, n = 250) {
  set.seed(seed)
  xs <- -abs(rnorm(n)) * 3
  yy <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  data.frame(y = yy, xs = xs * s, xr = xs)
}
fos <- list(alpha0 = bf(y ~ xr, sigma ~ 1, alpha ~ 0 + xs),
            alpha1 = bf(y ~ xr, sigma ~ 1, alpha ~ xs),
            mu = bf(y ~ xs, sigma ~ 1, alpha ~ 1))
run <- function(fo, d, a) {
  w <- character(0)
  f <- withCallingHandlers(frm(fo, family = skew_normal(), data = d,
                               control = frmtmb_control(autoscale = a)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  list(ll = as.numeric(logLik(f)), w = w)
}
n_less <- 0L
for (nm in names(fos)) for (seed in 1:12) for (s in c(5.4e-4, 2e-4, 1e-6)) {
  d <- mk(seed, s)
  a <- run(fos[[nm]], d, NULL)
  b <- run(fos[[nm]], d, FALSE)
  if (length(a$w) < length(b$w)) {
    n_less <- n_less + 1L
    cat(sprintf("%-7s seed %2d s %-7g default %.7f (%d warnings) FALSE %.7f (%d): %s\n",
                nm, seed, s, a$ll, length(a$w), b$ll, length(b$w),
                substr(b$w[1], 1, 60)))
  }
}
cat("fits where the default warns less than FALSE:", n_less, "of 108\n")

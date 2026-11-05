# Item 1 and 2: 256 bits at frm()'s point and at the marginal's best, on the lane's
# noise-free design (seed 20261102, x on 0 to 6, rounded) and on mine. Does the tape's
# ll equal the exact value at its own point? How far short is it?
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
run <- function(lab, y, x) {
  r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x))
  f <- r$fit; g <- rev_gp_read(f)
  s <- cf_setup(y, x)
  m <- cf_max(s, extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 30L)
  ex_f <- tryCatch(cf_ll_mpfr(s, g$beta, g$sigma, g$sd, g$rho), error = function(e) NA)
  ex_m <- tryCatch(cf_ll_mpfr(s, m$beta, m$sigma, m$sd, m$rho), error = function(e) NA)
  fv <- fitted(f)[, 1]
  cat(sprintf(paste0("RMP %s n %d m %d | frm tape ll %.4f exact at its point %.4f (tape-exact %.2e) sigma %.3g sd %.4g rho %.4g conv %d w %d | ",
                     "cf ll %.4f exact %.4f sigma %.3g sd %.4g rho %.4g | exact gap %.3f | max|fit-y| %.2e\n"),
              lab, length(y), s$m, g$ll, ex_f, g$ll - ex_f, g$sigma, g$sd, g$rho, g$conv, length(r$warn), m$ll, ex_m, m$sigma, m$sd, m$rho, ex_f - ex_m, max(abs(fv - y))))
  flush.console()
}
set.seed(20261102)
n <- 60L
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
for (sg in c(1e-4, 1e-5, 3e-6, 1e-6, 1e-7, 1e-8, 0)) run(sprintf("lane-design sin2x s=%g", sg), sin(2 * x) + sg * e0, x)
for (sc in c(1000, 1e-3)) run(sprintf("lane-design noise-free y*%g", sc), sc * sin(2 * x), x)
run("lane-design noise-free x*1000", sin(2 * x), x * 1000)
set.seed(4)
x2 <- sort(runif(25, 0, 10))
for (fn in c("sin", "exp", "x2")) {
  y <- switch(fn, sin = sin(x2), exp = exp(0.8 * x2), x2 = x2^2)
  run(sprintf("n25 %s noise-free", fn), y, x2)
}
cat("DONE\n")

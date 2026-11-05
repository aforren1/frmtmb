# Short check: distinct = 6L in gp_marginal_max(), patched in process, on the other tables: my first 168, my fresh 192, my new-curve 192; then the two test files.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
f <- ns$gp_marginal_max; formals(f)$distinct <- 6L; environment(f) <- ns
unlockBinding("gp_marginal_max", ns); assign("gp_marginal_max", f, envir = ns); lockBinding("gp_marginal_max", ns)
run <- function(lab, fns, sgs, ns_, seeds, seedf) {
  tot <- 0; sh <- 0; sil <- 0; worst <- 0
  for (seed in seeds) for (fn in names(fns)) for (sg in sgs) for (n in ns_) {
    set.seed(seedf(seed, n)); x <- sort(round(runif(n, 0, 10), 4)); y <- fns[[fn]](x) + sg * rnorm(n)
    r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x)); g <- rev_gp_read(r$fit)
    m <- cf_max(cf_setup(y, x), extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 22L)
    gap <- g$ll - m$ll; tot <- tot + 1; worst <- min(worst, gap)
    if (gap < -1e-3) { sh <- sh + 1; if (g$conv == 0 && !length(r$warn)) sil <- sil + 1; cat(sprintf("RD6 %s %s s=%g n %d seed %d gap %.3g conv %d w %d\n", lab, fn, sg, n, seed, gap, g$conv, length(r$warn))) }
    flush.console()
  }
  cat(sprintf("RD6SUM %s rows %d short %d silent %d worst %.3g\n", lab, tot, sh, sil, worst))
}
first <- list(mix = function(x) 2 * sin(3 * x) + x, mix2 = function(x) sin(4 * x) + 0.5 * x, mix3 = function(x) 3 * cos(2 * x) - 2 * x, mix4 = function(x) sin(5 * x) + x^2 / 10, two = function(x) sin(x) + 0.3 * sin(10 * x), sin3 = function(x) sin(3 * x), sin10 = function(x) sin(10 * x))
run("first168", first, c(0.3, 0.05, 1e-2), c(40L, 120L), 1:4, function(s, n) 31000 + s * 7 + n)
fresh <- list(saw = function(x) (x %% 2) - 1 + 0.2 * x, abs3 = function(x) abs(sin(x)) * 3 - 0.5 * x, dbl = function(x) exp(-(x - 3)^2) + 2 * exp(-(x - 7)^2 / 0.1), gauss4 = function(x) sin(7 * x) * exp(-x / 5) + 0.3 * x, poly = function(x) (x - 2) * (x - 5) * (x - 8) / 10 + sin(4 * x), sq = function(x) sign(sin(1.5 * x)) + 0.1 * x, tan1 = function(x) atan(5 * (x - 5)) + cos(9 * x) / 2, big = function(x) 50 + 10 * sin(2.5 * x) + 3 * x)
run("fresh192", fresh, c(1, 0.1, 1e-2, 1e-3), c(30L, 70L, 150L), 1:2, function(s, n) 55000 + s * 11 + n)
new <- list(ramp = function(x) pmin(x, 4) + 0.3 * sin(5 * x), trig2 = function(x) sin(x) * cos(4 * x), gate = function(x) (x > 2 & x < 6) * 2 + 0.1 * x, sqrtw = function(x) 2 * sqrt(x) + 0.5 * cos(6 * x), exp2 = function(x) exp(0.4 * x) + sin(3 * x), kink = function(x) abs(x - 5) + 0.2 * sin(10 * x), beat = function(x) sin(3 * x) + sin(3.7 * x), lin = function(x) -30 + 8 * x + cos(2 * x))
run("new192", new, c(1, 0.1, 1e-2, 1e-3), c(30L, 70L, 150L), 21:22, function(s, n) 73000 + s * 29 + n)
suppressPackageStartupMessages(library(testthat))
for (tf in c("test-gp-start.R", "test-gp-noncentered.R")) {
  r <- as.data.frame(test_file(file.path(REV_WT, "tests/testthat", tf), package = "frmtmb", env = test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE))
  cat(sprintf("RD6T %s pass %d fail %d err %d\n", tf, sum(r$nb - r$failed - r$warning - r$skipped), sum(r$failed), sum(r$error)))
  print(r[r$failed > 0 | r$error, "test"])
}
cat("DONE\n")

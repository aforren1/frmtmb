# Lane gpfield: test-brms-likelihood.R's row 10 (brms_gp_data(): set.seed(23),
# 80 rows, sin(x) + N(0, 0.3) on 0 to 10) and row 10c (gpn_data():
# set.seed(20261102), sin(2x) + 0.3 e on 0 to 6), y ~ gp(x, k = 10): the
# default fit and the profile over the length scale. Args: arm
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs(a[1])
d10 <- local({
  set.seed(23); n <- 80
  d <- data.frame(x = sort(runif(n, 0, 10))); d$y <- sin(d$x) + rnorm(n, 0, 0.3); d
})
d10c <- local({
  set.seed(20261102); n <- 60L
  x <- round(stats::runif(n, 0, 6), 3); z <- round(stats::runif(n, 0, 4), 2)
  e <- stats::rnorm(n)
  data.frame(x = x, y = sin(2 * x) + 0.3 * e)
})
for (nm in c("row10", "row10c")) {
  d <- if (nm == "row10") d10 else d10c
  q <- gf_fit(bf(y ~ gp(x, k = 10)), data = d, family = gaussian())
  f <- q$fit
  cat(nm, "DEFAULT ll", as.numeric(logLik(f)), "theta", f$estimates$theta,
      "| warn", sub(":.*", "", q$warn), "\n")
  th <- which(names(f$opt$par) == "theta")
  for (ls in c(1e-4, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4)) {
    best <- -Inf
    for (s0 in c(0.5, 5)) {
      p0 <- f$opt$par
      p0[th] <- c(log(s0), log(ls))
      fn <- function(q) { p <- p0; p[-th[2]] <- q; f$obj$fn(p) }
      gr <- function(q) { p <- p0; p[-th[2]] <- q; f$obj$gr(p)[-th[2]] }
      o <- tryCatch(nlminb(p0[-th[2]], fn, gr), error = function(e) NULL)
      if (!is.null(o)) best <- max(best, -o$objective)
    }
    cat(sprintf("%s PROF ls %.4g (x %.3g of the spacing 0.25) ll %.4f\n", nm,
                ls, ls / 0.25, best))
  }
}

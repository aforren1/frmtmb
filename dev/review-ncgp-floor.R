# Item 2: where the floor starts, across n, lscale and scale. Args: arm n
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- args[1]; n <- as.integer(args[2])
rev_libs(arm)
for (om in c(0.5, 2, 6)) for (ysc in c(1, 1000, 1e-3)) {
  for (ratio in c(1e-3, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 0)) {
    set.seed(4100 + n)
    x <- sort(runif(n, 0, 10))
    e <- rnorm(n)
    y <- ysc * (sin(om * x) + ratio * e)
    d <- data.frame(y = y, x = x)
    r <- rev_fit(y ~ gp(x), d)
    if (inherits(r$fit, "error")) {
      cat(sprintf("RF n %d om %g ysc %g ratio %g ERROR %s\n", n, om, ysc, ratio,
                  substr(conditionMessage(r$fit), 1, 100)))
      next
    }
    g <- rev_gp_read(r$fit)
    s <- cf_setup(y, x, nug = if (arm == "lane") 1e-12 else 1e-6)
    m <- cf_max(s, extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 28L)
    fv <- fitted(r$fit)
    cat(sprintf(paste0("RF n %d om %g ysc %g ratio %g | frm sigma/ysc %.4g sd/ysc %.4g rho %.4g ",
                       "ll %.4f conv %d nwarn %d | cf sigma/ysc %.4g sd/ysc %.4g rho %.4g ll %.4f | gap %.3e ",
                       "maxresid/ysc %.2e\n"),
                n, om, ysc, ratio, g$sigma / ysc, g$sd / ysc, g$rho, g$ll, g$conv, length(r$warn),
                m$sigma / ysc, m$sd / ysc, m$rho, m$ll, g$ll - m$ll,
                max(abs(y - fv[, 1])) / ysc))
    flush.console()
  }
}
cat("DONE\n")

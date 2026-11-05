# Short check 2, item 4, R4. With the stopgap off on new data, does a new-data
# route read the dense joint covariance? Peak memory of the process per route.
# Usage: Rscript review-latres-r4b.R <route> <n>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
options(warn = -1)
a <- commandArgs(TRUE)
route <- a[[1]]; N <- as.integer(a[[2]])
peak <- function() as.numeric(system(sprintf('powershell -NoProfile -Command "(Get-Process -Id %d).PeakWorkingSet64"', Sys.getpid()), intern = TRUE)[1]) / 1e9
set.seed(1)
d <- data.frame(t = seq_len(N), g = factor(1), x = rnorm(N))
d$y <- rpois(N, exp(.5 + .3 * d$x + as.numeric(arima.sim(list(ar = .5), N, sd = .5))))
f <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)
p0 <- peak()
nd <- d[1:50, ]
t1 <- system.time(r <- switch(route,
  fitted_new = fitted(f, newdata = nd),
  linpred_new = frm_linpred(f, newdata = nd, se.fit = TRUE),
  ce = conditional_effects(f, "x"),
  curve = frm_curve(f, newdata = nd),
  basis_new = frm_lp_basis(f, newdata = nd),
  predictions = marginaleffects::predictions(f, newdata = nd),
  emmeans = emmeans::emmeans(f, ~ x, at = list(x = 0))))[[3]]
cat(sprintf("%-12s n %d: %.1f s, peak %.2f GB before, %.2f GB after; finite SE: %s\n", route, N, t1, p0, peak(),
            tryCatch(if (route == "fitted_new") all(is.finite(r[, 2])) else if (route == "ce") all(is.finite(r$x$lower__)) else NA, error = function(e) NA)))

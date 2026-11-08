# Item 2 probe: a saddle on `disc` beside a ridge on mu (the nlfix review's
# dev/review-nlfix-d10.R block "lost", seed 11001). Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
set.seed(11001); n <- 300
dd <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
invisible(factor(rep(1:15, each = 20))); invisible(rnorm(15, 0, 0.5))
lat <- 1.5 * exp(0.6 * dd$x) + 0.7 * dd$z
dd$yo <- ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
fm4 <- bf(yo ~ 1) + nlf(mu ~ b1 * exp(b2 * x)) + nlf(disc ~ exp(d1 * d2 * z)) +
  lf(b1 ~ 1) + lf(b2 ~ 1) + lf(d1 ~ 1) + lf(d2 ~ 1)
show <- function(lab, r) {
  f <- r$fit
  if (inherits(f, "error")) { cat(lab, "ERROR", conditionMessage(f), "\n"); return() }
  cat(sprintf("%s ll %.4f code %d msg %s\n  beta %s\n  betad %s\n", lab,
              as.numeric(logLik(f)), f$opt$convergence, f$opt$message,
              paste(signif(f$estimates$beta, 5), collapse = " "),
              paste(signif(f$estimates$betad, 5), collapse = " ")))
  for (w in r$warn) cat("  W:", substr(w, 1, 400), "\n")
}
t0 <- proc.time()[[3]]
r0 <- fe_fit(fm4, family = cumulative(), data = dd)
cat(sprintf("default %.1f s\n", proc.time()[[3]] - t0))
show("default", r0)
show("start all moved", fe_fit(fm4, family = cumulative(), data = dd,
  start = list(beta = c("b1_(Intercept)" = 1, "b2_(Intercept)" = 0.5),
               betad = c(0.3, 0.3))))
show("start disc only", fe_fit(fm4, family = cumulative(), data = dd,
  start = list(betad = c(0.3, 0.3))))
show("start mu only", fe_fit(fm4, family = cumulative(), data = dd,
  start = list(beta = c("b1_(Intercept)" = 1, "b2_(Intercept)" = 0.5))))
print(par_template(fm4, family = cumulative(), data = dd))
cat("DONE\n")

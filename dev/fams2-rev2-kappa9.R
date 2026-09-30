# Reviewer, item 3 follow-up: the Beta()-data fits whose kappa stopped
# above 1e-6 without a warning (dev/fams2-rev2-guards.txt 3a), and
# kappa ~ x on seed 31. Is kappa really placed there (logLik above
# Beta(), a usable standard error), or did it run to 0 unseen?
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
sim <- function(n, phi, s) {
  set.seed(s)
  x <- rnorm(n); mu <- plogis(0.2 + 0.5 * x)
  data.frame(x = x, y = rbeta(n, mu * phi, (1 - mu) * phi))
}
q <- function(expr) suppressWarnings(expr)
for (cfg in list(c(400, 5, 9), c(2000, 5, 2), c(2000, 5, 7), c(400, 50, 4),
                 c(400, 50, 5), c(400, 50, 9), c(400, 5, 1))) {
  d <- sim(cfg[1], cfg[2], cfg[3])
  f <- q(frm(y ~ x, family = xbeta(), data = d))
  fb <- frm(y ~ x, family = Beta(), data = d)
  ci <- confint(f)["kappa_(Intercept)", ]
  se <- sqrt(diag(vcov(f, full = TRUE)))[["kappa_(Intercept)"]]
  # the xbeta objective with kappa pinned near 0
  f0 <- q(frm(bf(y ~ x, kappa = 1e-9), family = xbeta(), data = d))
  cat(sprintf(paste0("n %d phi %g seed %d: kappa %.3g (log %.2f, SE %.3g) | logLik xbeta %.6f, ",
                     "kappa held 1e-9 %.6f, Beta %.6f | gain over Beta %.2e\n"),
              cfg[1], cfg[2], cfg[3], exp(ci[["est"]]), ci[["est"]], se,
              as.numeric(logLik(f)), as.numeric(logLik(f0)), as.numeric(logLik(fb)),
              as.numeric(logLik(f)) - as.numeric(logLik(fb))))
}
cat("\nkappa ~ x on seed 31 (n 400, phi 5):\n")
d <- sim(400, 5, 31)
w <- character()
f <- withCallingHandlers(frm(bf(y ~ x, kappa ~ x), family = xbeta(), data = d),
                         warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
cat("warnings:", length(w), if (length(w)) w, "\n")
print(confint(f)[c("kappa_(Intercept)", "kappa_x"), ])
print(sqrt(diag(vcov(f, full = TRUE)))[c("kappa_(Intercept)", "kappa_x")])
cat("logLik", as.numeric(logLik(f)), " Beta", as.numeric(logLik(frm(y ~ x, family = Beta(), data = d))), "\n")

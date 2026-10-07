# Lane optima, item 2: RTMB's pnorm(log.p = TRUE) and the probit
# log-odds built from it, value and derivative far out in the tails.
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- frmtmb:::frmtmb_links$probit$logit_eta
x <- c(1e2, 1e3, 1e4, 3e4, 1e5, 1e6, 1e8, 4.5e9, 1e12, 1e100)
for (s in c(1, -1)) {
  xx <- s * x
  a <- RTMB::MakeTape(function(e) RTMB::pnorm(e, log.p = TRUE), xx)
  b <- RTMB::MakeTape(function(e) q(e), xx)
  cat("x =", format(xx), "\n  log Phi(x):", format(a(xx), digits = 6),
      "\n  d log Phi:", format(diag(a$jacobian(xx)), digits = 6),
      "\n  q(x):", format(b(xx), digits = 6),
      "\n  dq:", format(diag(b$jacobian(xx)), digits = 6), "\n")
}

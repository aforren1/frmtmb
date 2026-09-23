LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
set.seed(8)
n <- 200
xb <- rnorm(n, 0, 1e4)
yb <- 0.001 * xb + RTMBdist::rskewnorm2(n, 0, 1.5, 4)
d3 <- data.frame(y = yb, x = xb)
f <- suppressWarnings(frm(bf(y ~ x, sigma ~ 1, alpha ~ 1),
                          family = skew_normal(), data = d3))
cat("M3 case: alpha", format(f$estimates$betad[["alpha_(Intercept)"]],
                             digits = 8),
    "\n  conv", f$opt$convergence, "'", f$opt$message, "'",
    "\n  max|grad|", format(max(abs(f$obj$gr(f$opt$par))), digits = 3),
    " grad_tol", frmtmb_control()$grad_tol,
    "\n  par_units", if (is.null(f$par_units)) "NULL (autoscale never engaged)"
    else "set",
    "\n  escape", if (is.null(f$opt[["stationary_escape"]])) "did not fire"
    else "fired",
    "\n  logLik", format(as.numeric(logLik(f)), digits = 12),
    " sn", format(sn::selm(y ~ x, family = "SN", data = d3)@logL,
                  digits = 12), "\n")

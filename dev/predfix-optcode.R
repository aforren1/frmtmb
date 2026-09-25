# Item 7, REPORT ONLY: a nonzero optimizer code on a parameter whose MLE
# is at infinity (skew_normal alpha on a half-normal residual). What
# the user sees today, and what a suppression rule keyed on frmtmb's own
# gradient criterion would have hidden. Nothing here changes behavior.
#   PREDFIX_ARM=lane Rscript dev/predfix-optcode.R > dev/predfix-log/optcode.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
fitw <- function(...) {
  w <- character(0)
  f <- withCallingHandlers(frm(...), warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(fit = f, warn = w)
}
gmax <- function(f) max(abs(f$obj$gr(f$opt$par) * (f$par_units %||% 1)))
`%||%` <- function(a, b) if (is.null(a)) b else a

cat("== A. alpha unbounded: y = 1 + 0.5 x + 1.2 |N(0, 1)|, n = 200\n")
cat(sprintf("%-4s %5s %-28s %12s %10s %6s %14s %14s\n", "seed", "conv",
            "message", "alpha", "max|grad|", "warn", "logLik", "sn::selm"))
rows <- list()
wtxt <- NULL
for (seed in 1:20) {
  set.seed(seed)
  n <- 200
  x <- rnorm(n)
  d <- data.frame(x = x, y = 1 + 0.5 * x + 1.2 * abs(rnorm(n)))
  r <- fitw(bf(y ~ x), family = skew_normal(), data = d)
  f <- r$fit
  if (is.null(wtxt) && length(r$warn)) wtxt <- r$warn[1]
  a <- f$estimates$betad[["alpha_(Intercept)"]]
  sn_ll <- tryCatch(sn::selm(y ~ x, family = "SN", data = d)@logL,
                    error = function(e) NA_real_)
  rows[[seed]] <- data.frame(seed = seed, conv = f$opt$convergence,
                             alpha = a, g = gmax(f), warned = length(r$warn),
                             ll = as.numeric(logLik(f)), sn = sn_ll)
  cat(sprintf("%-4d %5d %-28s %12.4g %10.3g %6d %14.6f %14.6f\n", seed,
              f$opt$convergence, substr(f$opt$message, 1, 28), a, gmax(f),
              length(r$warn), as.numeric(logLik(f)), sn_ll))
}
A <- do.call(rbind, rows)
cat("nonzero code:", sum(A$conv != 0), "of", nrow(A),
    "; |alpha| > 1e3:", sum(abs(A$alpha) > 1e3),
    "; nonzero code AND max|grad| < grad_tol:",
    sum(A$conv != 0 & A$g < frmtmb_control()$grad_tol),
    "; logLik - sn (min, max):", paste(signif(range(A$ll - A$sn,
                                                     na.rm = TRUE), 4),
                                       collapse = ", "), "\n")
cat("first warning text:", wtxt, "\n")

cat("\n== B. the reviewer's finite case (dev/skewinit-m3.R)\n")
set.seed(8)
n <- 200
xb <- rnorm(n, 0, 1e4)
yb <- 0.001 * xb + RTMBdist::rskewnorm2(n, 0, 1.5, 4)
d3 <- data.frame(y = yb, x = xb)
r <- fitw(bf(y ~ x, sigma ~ 1, alpha ~ 1), family = skew_normal(), data = d3)
f <- r$fit
cat("alpha", format(f$estimates$betad[["alpha_(Intercept)"]], digits = 8),
    " conv", f$opt$convergence, "'", f$opt$message, "'",
    " max|grad|", format(gmax(f), digits = 3),
    " logLik - sn", format(as.numeric(logLik(f)) -
                             sn::selm(y ~ x, family = "SN", data = d3)@logL,
                           digits = 3),
    " warnings", length(r$warn), "\n")

cat("\n== C. max|grad| below grad_tol is not an optimum certificate\n")
cat("   (0.62.0 default, no autoscale, poisson y ~ 0 + x at 1e-6)\n")
set.seed(23)
xt <- rnorm(250)
dd <- data.frame(y = rpois(250, exp(0.5 + 0.4 * xt)), xs = xt * 1e-6)
f <- frm(bf(y ~ 0 + xs) + poisson(), data = dd,
         control = frmtmb_control(autoscale = FALSE))
gl <- as.numeric(logLik(glm(y ~ 0 + xs, family = poisson(), data = dd)))
cat("conv", f$opt$convergence, " max|grad|", format(gmax(f), digits = 3),
    " grad_tol", frmtmb_control()$grad_tol, " logLik - glm",
    format(as.numeric(logLik(f)) - gl, digits = 8), "\n")

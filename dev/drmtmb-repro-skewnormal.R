# Minimal reproduction, frmtmb only: skew_normal() stalls at alpha = 0
# when the raw response is skewed the other way from the residuals.
# Found by the drmTMB agreement lane; not fixed here (lane changes no R/).
.libPaths(c("C:/Users/adf44/source/r/drmtmb-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
set.seed(1)
n <- 200
xs <- -abs(rnorm(n)) * 3                        # left-skewed covariate
y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + # right-skewed residual
  rnorm(n, 0, 0.3)
d <- data.frame(y = y, xs = xs)

msgs <- character()
f_def <- withCallingHandlers(
  frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(), data = d,
      se = TRUE),
  warning = function(w) { msgs <<- c(msgs, conditionMessage(w))
                          invokeRestart("muffleWarning") },
  message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                          invokeRestart("muffleMessage") })
f_pos <- frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(),
             data = d, start = list(betad = c(0, 2)))
cat("default start: logLik", format(logLik(f_def), digits = 10),
    " alpha", f_def$opt$par[[4]], " convergence", f_def$opt$convergence,
    " message", f_def$opt$message, "\n")
cat("alpha start 2: logLik", format(logLik(f_pos), digits = 10),
    " alpha", f_pos$opt$par[[4]], "\n")
cat("conditions raised by the default fit:", length(msgs), "\n")
if (length(msgs)) print(msgs)
print(summary(f_def)$coefficients)
cat("max |gradient| at the stalled point:",
    max(abs(f_def$obj$gr(f_def$opt$par))), "\n")

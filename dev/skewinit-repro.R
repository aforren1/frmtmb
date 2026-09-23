# Reproduction of the skew_normal() alpha = 0 stall on the BASE build.
# Run with LIB unset so only rellib-r3 supplies frmtmb.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
op <- options(digits = 12)

dead <- as.logical(Sys.getenv("SKEWINIT_DEAD", "FALSE"))
dd <- make_data(1, dead)
cat("\n--- seed 1, dead_draw =", dead, "\n")
cat("skew(raw y)      ", skew(dd$y), "\n")
cat("skew(lm resid)   ", skew(residuals(lm(y ~ xs, data = dd))), "\n")

f0 <- frm_sn(dd)
cat("default logLik   ", as.numeric(logLik(f0)), "\n")
cat("default alpha    ", f0$opt$par[[4]], "\n")
cat("convergence      ", f0$opt$convergence, " '", f0$opt$message, "'\n")
cat("max|grad|        ", max(abs(f0$obj$gr(f0$opt$par))), "\n")
cnt <- 0L
withCallingHandlers(frm_sn(dd),
                    condition = function(c) cnt <<- cnt + 1L)
cat("conditions raised", cnt, "\n")

f1 <- frm_sn(dd, start = list(betad = c(0, 2)))
cat("start(0,2) logLik", as.numeric(logLik(f1)), "\n")
cat("gap              ", as.numeric(logLik(f1)) - as.numeric(logLik(f0)),
    "\n")
s <- sn_ll(dd)
cat("sn::selm logLik  ", s[["ll"]], "\n")
cat("sn::selm alpha   ", s[["alpha"]], "\n")
cat("sn vs frm start  ", s[["ll"]] - as.numeric(logLik(f1)), "\n")

# The standard error at the stall: the detection signal to calibrate.
sd0 <- summary(f0)
print(sd0)
options(op)

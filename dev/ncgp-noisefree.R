# Noise-free y: the case with nothing for sigma to read. What frm() says, on both
# builds. Seed 20261102.
#   Rscript dev/ncgp-noisefree.R <lib or base>
a <- commandArgs(TRUE)
libs <- c("C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (a[1] != "base") libs <- c(a[1], libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261102)
n <- 60L
x <- round(runif(n, 0, 6), 3)
d <- data.frame(x = x, y = sin(2 * x))
w <- character()
t0 <- Sys.time()
f <- withCallingHandlers(frm(bf(y ~ gp(x)), data = d, family = gaussian()),
                         warning = function(c) { w <<- c(w, substr(conditionMessage(c), 1, 110)); invokeRestart("muffleWarning") },
                         message = function(c) { w <<- c(w, substr(conditionMessage(c), 1, 110)); invokeRestart("muffleMessage") })
cat(sprintf("NCGPNF %s ll %.4f conv %d %s sigma %.3g sdgp %.3g rho %.3g %.1fs\n", a[1], -f$opt$objective, f$opt$convergence, f$opt$message,
            sigma(f)[1], exp(f$estimates$theta[1]), exp(f$estimates$theta[2]), as.numeric(Sys.time() - t0, units = "secs")))
print(unique(w))
print(head(frm_linpred(f, newdata = data.frame(x = c(1.2345, 3.3333)))))
print(sin(2 * c(1.2345, 3.3333)))
print(max(abs(fitted(f)[, "Estimate"] - d$y)))

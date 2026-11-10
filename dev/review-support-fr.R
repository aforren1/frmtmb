# Reviewer probe (lane support): Frechet on seed-402 data, ub only, SE
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(402)
n <- 300
e <- rexp(n)
d <- data.frame(y = 2 * e + 0.05, ub = 1e6)
fr0 <- frm(y ~ 1, family = frechet(), data = d)
fr1 <- frm(y | trunc(ub = ub) ~ 1, family = frechet(), data = d)
cat("frechet ub-only dlogLik", as.numeric(logLik(fr1) - logLik(fr0)), "\n")
print(fr0$opt$par, digits = 17); print(as.numeric(logLik(fr0)), digits = 17)
print(fr0$sdr)

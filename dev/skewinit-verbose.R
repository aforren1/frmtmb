# The verbose stage line the escape adds, on a fit that fires it and on
# one that does not.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
dd <- make_data(1, FALSE)
m3 <- skew(dd$y)
bad <- list(betad = c(log(stats::sd(dd$y)), 2 * sign(m3) + 0.5 * m3))
cat("=== does not fire\n")
invisible(frm_sn(dd, verbose = TRUE))
cat("\n=== fires\n")
invisible(frm_sn(dd, start = bad, verbose = TRUE))

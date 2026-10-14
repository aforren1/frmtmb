# Fix round 2: brms's G1 wiener block (tests.distributions.R:212)
# compares brms's dwiener() and rwiener() across its two backends,
# RWiener and rtdists. frmtmb.eam has one implementation, so a shim
# would make d1 and d2 the same call and the assertion hollow. What can
# be checked is that frmtmb.eam's density agrees with both backends at
# the block's own inputs (seed 1234, as brms's block sets it).
#
#   Rscript dev/remeasure-bin2-wiener-g1.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.eam))
cat("frmtmb.eam", format(packageVersion("frmtmb.eam")), "\n")
set.seed(1234)
n <- 10
x <- seq(0.1, 1, 0.1)
alpha <- rexp(n)
tau <- 0.05
beta <- 0.5
delta <- rnorm(n)
resp <- sample(c(0, 1), n, TRUE)
eam <- as.numeric(wiener_lpdf(x - tau, delta, alpha, beta, resp))
d1 <- brms::dwiener(x, alpha, tau, beta, delta, resp, backend = "Rwiener")
d2 <- brms::dwiener(x, alpha, tau, beta, delta, resp, backend = "rtdists")
# on the log scale, where brms's density has not underflowed to 0
ok <- d1 > 0 & d2 > 0
cat(sprintf(paste("rows compared %d of %d; max abs log diff frmtmb.eam",
                  "vs RWiener %.3g, vs rtdists %.3g\n"), sum(ok), n,
            max(abs(eam[ok] - log(d1[ok]))), max(abs(eam[ok] - log(d2[ok])))))

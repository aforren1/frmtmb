# Reviewer 2: the decb grid point where the lane and WienR differ most
# (relative, upper boundary), settled against the Rmpfr reference.
source("dev/phase3b-review2/r2-prelude.R"); r2_lib("lane")
suppressPackageStartupMessages({library(frmtmb.eam); library(WienR); library(Rmpfr)})
source("dev/phase3b-review2/r2-mpfr-lib.R")
ns <- asNamespace("frmtmb.eam")
set.seed(3); n <- 400
g <- data.frame(v = runif(n, -4, 4), a = runif(n, 0.5, 3), w = runif(n, 0.1, 0.9),
                u = exp(runif(n, log(0.02), log(5)))); g$t <- g$u * g$a^2; tau <- 0.2
lu <- ns$ddm_rt_lcdf_b(g$t, g$v, g$a, g$w, 1)
wr <- pWDM(g$t + tau, "upper", g$a, g$v, g$w, tau)$value
i <- order(-abs(exp(lu) - wr) / wr)[1:3]
for (j in i) {
  ref <- lower_images(g$t[j], -g$v[j], g$a[j], 1 - g$w[j])
  cat(sprintf("v %.3f a %.3f w %.3f u %.4f: ref %.15g  lane rel %.2e  WienR rel %.2e\n",
              g$v[j], g$a[j], g$w[j], g$u[j], as.numeric(ref),
              abs(exp(lu[j]) / as.numeric(ref) - 1), abs(wr[j] / as.numeric(ref) - 1)))
}

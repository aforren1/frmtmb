# Reviewer: the negative-binomial far upper tail where log_ibeta_two()
# and stats::pnbinom() disagree; a third reference by summing the mass
# function in log space, and the same for the beta (direct pbeta form)
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
for (cs in list(c(7.5, 1e4, 1000500), c(7.5, 1e6, 100000500), c(7.5, 1e4, 3e5),
                c(7.5, 1e4, 5e5), c(7.5, 1e4, 7e5), c(2, 100, 8000),
                c(1, 10, 5000))) {
  sh <- cs[1]; mu <- cs[2]; q <- cs[3]
  k <- seq(q + 1, q + 2e6, by = 1)
  ref_sum <- lse(dnbinom(k, size = sh, mu = mu, log = TRUE))
  ref_st <- pnbinom(q, size = sh, mu = mu, lower.tail = FALSE, log.p = TRUE)
  got <- ns$nbinom_cdf_log(q, log(mu), log(sh))$lS
  cat(sprintf("shape %g mu %g q %g: lane %.8f | stats %.8f | mass sum %.8f\n",
              sh, mu, q, got, ref_st, ref_sum))
}
cat("== the beta tail behind it: log I_x(a, b) for x near 1 and large a\n")
for (cs in list(c(1e6 + 1, 7.5, 1 - 7.5e-4), c(3e5 + 1, 7.5, 1 - 7.5e-4),
                c(1e5 + 1, 7.5, 1 - 7.5e-4))) {
  a <- cs[1]; b <- cs[2]; x <- cs[3]
  got <- ns$log_ibeta_two(x, 1 - x, a, b)
  cat(sprintf("a %g b %g x %.6g: lF lane %.8f stats %.8f | log_ibeta_half %.8f\n",
              a, b, x, got$lF, pbeta(x, a, b, log.p = TRUE),
              ns$log_ibeta_half(x, a, b)))
}

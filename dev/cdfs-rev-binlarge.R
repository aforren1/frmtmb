# Reviewer: binomial and negbinomial CDF derivatives where
# log_ibeta_half() takes RTMB::pbeta() at whole shapes (large trials,
# large counts), against central differences of stats on doubles
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
for (cs in list(c(1000, 400, 0.4), c(1000, 400, 0.41), c(1e5, 40000, 0.4), c(1e6, 5e5, 0.5),
                c(1e6, 499000, 0.5), c(600, 300, 0.5))) {
  n <- cs[1]; q <- cs[2]; p0 <- cs[3]
  f <- function(e) {
    p <- 1 / (1 + exp(-e))
    r <- ns$binomial_cdf_log(q, n, list(p = p, q = 1 - p))
    c(r$lF, r$lS)
  }
  e0 <- qlogis(p0)
  ad <- RTMB::MakeTape(f, e0)$jacobian(e0)
  h <- 1e-6
  fd <- (c(pbinom(q, n, plogis(e0 + h), log.p = TRUE), pbinom(q, n, plogis(e0 + h), lower.tail = FALSE, log.p = TRUE)) -
         c(pbinom(q, n, plogis(e0 - h), log.p = TRUE), pbinom(q, n, plogis(e0 - h), lower.tail = FALSE, log.p = TRUE))) / (2 * h)
  cat(sprintf("binomial n %g q %g p %g: AD %s FD %s\n", n, q, p0,
              paste(signif(ad, 8), collapse = ","), paste(signif(fd, 8), collapse = ",")))
}

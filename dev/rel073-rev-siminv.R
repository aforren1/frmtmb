# Release reviewer, 0.73.0: sim_trunc_inverse() alone, on negbinomial()
# and brmsfamily("negbinomial2") at one parameter point, 2e5 draws,
# against the exact truncated pmf. And sim_response(), the route
# simulate() takes, at a point where rejection rarely fills a row.
#   Rscript dev/rel073-rev-siminv.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
m <- 200000L
mu <- 30; size <- 4; ub <- 6
p <- dnbinom(0:ub, size, mu = mu); p <- p / sum(p)
cat(sprintf("window mass %.4g\n", pnbinom(ub, size, mu = mu)))
chk <- function(label, y) {
  ob <- tabulate(y + 1L, ub + 1L)
  ex <- p * length(y)
  cat(sprintf(paste0("%-40s out of window %d | mean %.5f exact %.5f (z %.2f)",
                     " | chi-square %.2f on %d df\n"),
              label, sum(y < 0 | y > ub), mean(y), sum((0:ub) * p),
              (mean(y) - sum((0:ub) * p)) /
                sqrt((sum((0:ub)^2 * p) - sum((0:ub) * p)^2) / length(y)),
              sum((ob - ex)^2 / ex), ub))
}
set.seed(11)
fam <- negbinomial()
y <- ns$sim_trunc_inverse(fam, list(mu = rep(mu, m), shape = rep(size, m)),
                          list(), m, seq_len(m), rep(-Inf, m), rep(ub, m))
chk("negbinomial sim_trunc_inverse", y)
fam2 <- brmsfamily("negbinomial2")
y2 <- ns$sim_trunc_inverse(fam2, list(mu = rep(mu, m), sigma = rep(1 / size, m)),
                           list(), m, seq_len(m), rep(-Inf, m), rep(ub, m))
chk("negbinomial2 sim_trunc_inverse", y2)
y3 <- ns$sim_response(fam2, list(mu = rep(mu, m), sigma = rep(1 / size, m)),
                      list(trunc_ub = rep(ub, m)), m)
chk("negbinomial2 sim_response", y3)
# a lower bound too: [8, 12] at mu = 30
lb <- 8; ub2 <- 12
p2 <- dnbinom(lb:ub2, size, mu = mu); p2 <- p2 / sum(p2)
y4 <- ns$sim_trunc_inverse(fam2, list(mu = rep(mu, m), sigma = rep(1 / size, m)),
                           list(), m, seq_len(m), rep(lb, m), rep(ub2, m))
ob <- tabulate(y4 - lb + 1L, ub2 - lb + 1L)
cat(sprintf("negbinomial2 [8, 12] inverse: out %d | chi-square %.2f on %d df\n",
            sum(y4 < lb | y4 > ub2), sum((ob - p2 * m)^2 / (p2 * m)),
            ub2 - lb))
cat("END\n")

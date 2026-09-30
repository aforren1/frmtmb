# The OSA branch of each ordinal density against its data branch, with
# disc, cs() and each threshold structure: with every row kept, the two
# must give the same log density. Seed 20260937.
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260937)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
osa_cls <- methods::getClass("osa", where = asNamespace("RTMB"))
print(osa_cls)
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  for (th in c("flexible", "equidistant", "sum_to_zero")) {
    f <- frm(bf(y ~ x, disc ~ 0 + z), family = get(fam)(threshold = th),
             data = d)
    ff <- family(f)
    dp <- frmtmb:::eval_dpars(f)[[1]]
    ex <- list(tau_raw = f$estimates$tau_raw)
    a <- ff$lpdf(d$y, dp, list(), ex)
    o <- methods::new(osa_cls, x = as.numeric(d$y),
                      keep = matrix(1, n, 1))
    b <- ff$lpdf(o, dp, list(), ex)
    cat(sprintf("%-12s %-12s max|osa - data| / max|data| = %.3g\n", fam, th,
                max(abs(as.numeric(b) - a)) / max(abs(a))))
  }
}

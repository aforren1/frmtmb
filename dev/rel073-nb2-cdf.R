# Consolidation item 1: negbinomial2 in cdfs's CDF tables. On the merged
# build, a censored and a truncated negbinomial2 fit against the same
# negbinomial() fit: the objective at mirrored parameters (log sigma =
# -log shape, the same b) and the fitted log-likelihoods, and what
# frm_compat() says for cens() and trunc() on both families.
#   Rscript dev/rel073-nb2-cdf.R <lib>   (log dev/rel073-log/nb2-cdf.txt)
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "\n")
set.seed(20261021)
n <- 400
d <- data.frame(x = stats::rnorm(n))
d$y <- stats::rnbinom(n, mu = exp(1.2 + 0.4 * d$x), size = 2.5)
q <- stats::quantile(d$y, c(0.15, 0.85))
d$cc <- ifelse(d$y <= q[1], -1, ifelse(d$y >= q[2], 1, 0))
d$y2 <- d$y + 3
d$cc[seq(5, n, by = 25)] <- 2
dt <- d[d$y >= 1, ]
cases <- list(
  cens = list(f = bf(y | cens(cc, y2) ~ x), data = d),
  trunc_lb1 = list(f = bf(y | trunc(lb = 1) ~ x), data = dt),
  trunc_1_40 = list(f = bf(y | trunc(lb = 1, ub = 40) ~ x), data = dt))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  f1 <- frm(cs$f, family = negbinomial(), data = cs$data)
  f2 <- frm(cs$f, family = brmsfamily("negbinomial2"), data = cs$data)
  # no random effects: the outer parameters are beta (mu) and betad,
  # the log shape for negbinomial() and the log sigma for negbinomial2,
  # so the mirrored point negates betad
  p1 <- f1$opt$par
  stopifnot(identical(names(p1), c("beta", "beta", "betad")))
  p2 <- p1
  p2[3] <- -p1[3]
  fn1 <- f1$obj$fn(p1)
  fn2 <- f2$obj$fn(p2)
  cat(sprintf(paste0("%-11s logLik nb %.9f nb2 %.9f diff %.2e | objective at ",
                     "mirrored point nb %.12f nb2 %.12f diff %.2e | ",
                     "sigma %.6f 1/shape %.6f | conv %d %d\n"),
              nm, as.numeric(logLik(f1)), as.numeric(logLik(f2)),
              as.numeric(logLik(f2)) - as.numeric(logLik(f1)), fn1, fn2,
              fn2 - fn1, exp(f2$opt$par[3]), 1 / exp(p1[3]),
              f1$opt$convergence, f2$opt$convergence))
  if (grepl("trunc", nm)) {
    cat(sprintf("%-11s fitted() max relative difference %.2e\n", nm,
                max(abs(fitted(f2)[, 1] / fitted(f1)[, 1] - 1))))
  }
}
tab <- frm_compat()
for (fam in c("negbinomial", "negbinomial2")) {
  for (at in c("cens()", "trunc()")) {
    r <- tab[(tab$feature_a == at & tab$feature_b == fam) |
               (tab$feature_b == at & tab$feature_a == fam), ]
    cat(sprintf("frm_compat %-13s %-8s %s\n", fam, at,
                paste(r$status, collapse = ",")))
  }
}
cat("END\n")

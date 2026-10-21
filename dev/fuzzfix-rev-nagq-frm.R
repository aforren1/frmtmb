# Reviewer: frmtmb quadrature on the data of dev/fuzzfix-rev-nagq.R,
# against GLMMadaptive's nested adaptive quadrature (nAGQ = 15, 25).
.libPaths(c("C:/Users/adf44/source/r/wt-fuzzfix-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(GLMMadaptive)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(31)
ng <- 12; m <- 4
d <- expand.grid(r = seq_len(m), gb = factor(c("i", "ii")),
                 ga = factor(seq_len(ng)))
d$x <- rnorm(nrow(d))
ab <- interaction(d$ga, d$gb, drop = TRUE)
d$y <- rbinom(nrow(d), 1, plogis(0.3 + 0.5 * d$x + rnorm(ng, 0, 0.8)[d$ga] +
                                   rnorm(nlevels(ab), 0, 0.6)[ab]))
fq <- suppressMessages(frm(bf(y ~ x + (1 | ga/gb)), family = bernoulli(),
                           data = d, quadrature = TRUE))
fl <- suppressMessages(frm(bf(y ~ x + (1 | ga/gb)), family = bernoulli(),
                           data = d))
cat(sprintf("frmtmb quad %.8f laplace %.8f\n", as.numeric(logLik(fq)),
            as.numeric(logLik(fl))))
print(fq$estimates$beta); print(VarCorr(fq))
for (k in c(15, 25)) {
  g <- mixed_model(y ~ x, random = ~ 1 | ga/gb, data = d,
                   family = binomial(), nAGQ = k)
  cat(sprintf("GLMMadaptive nAGQ %d logLik %.8f\n", k, as.numeric(logLik(g))))
  print(fixef(g)); print(g$D)
}

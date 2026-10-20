# Lane ndaterm, fix round 1: is the mixture predict(newdata) gap at
# N = 1000 (431 against fitted() 484) the parameter-draw law or the
# simulator? Seed 7 data of dev/ndaterm-r1-probe.R. Log
# dev/ndaterm-out/r1-mixgap.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ndaterm-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n), N = sample(c(2, 50), n, TRUE))
cls <- rbinom(n, 1, 0.5)
d$y <- rbinom(n, d$N, ifelse(cls == 1, plogis(-1 + 0.5 * d$x),
                             plogis(1 + 0.5 * d$x)))
fmx <- frm(bf(y | trials(N) ~ 1, mu1 ~ x, mu2 ~ x),
           family = mixture(binomial(), binomial()), data = d)
nd <- data.frame(x = c(0, 0), N = c(1000, 2))
f <- fitted(fmx, newdata = nd)
cat("fitted", f[, "Estimate"], " Est.Error", f[, "Est.Error"], "\n")
for (pe in c(TRUE, FALSE)) {
  set.seed(3)
  p <- predict(fmx, newdata = nd, summary = FALSE, ndraws = 4000,
               propagate_error = pe)
  cat("propagate_error", pe, "means", colMeans(p), " se",
      apply(p, 2, sd) / sqrt(nrow(p)), "\n")
}
s <- as.matrix(simulate(fmx, nsim = 4000, newdata = nd, seed = 3))
cat("simulate(newdata) means", rowMeans(s), " se",
    apply(s, 1, sd) / sqrt(ncol(s)), "\n")
print(fixef(fmx))

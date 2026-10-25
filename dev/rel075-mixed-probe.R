# Can one model carry a cox() family with a car() field and a sar()
# residual term? sar() takes a gaussian or student response only, so
# the model is multivariate: a cox() response with car() and a gaussian
# response with sar(), over the same rows. Seed 20261024.
#   Rscript dev/rel075-mixed-probe.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.survival)
                  library(frmtmb.spatial)})
cat("lib:", dirname(find.package("frmtmb")),
    format(packageVersion("frmtmb")), "\n")
g <- expand.grid(r = 1:3, c = 1:3)
n <- nrow(g)
W <- 1 * (as.matrix(dist(g, method = "manhattan")) == 1)
dimnames(W) <- list(paste0("L", 1:n), paste0("L", 1:n))
set.seed(20261024)
dd <- data.frame(loc = factor(rep(rownames(W), each = 12),
                              levels = rownames(W)), x = rnorm(12 * n))
u <- rnorm(n, 0, 0.5)
dd$time <- rexp(nrow(dd), exp(-0.5 + 0.6 * dd$x + u[dd$loc]))
dd$out <- rbinom(nrow(dd), 1, 0.25)
# sar() over the rows: a row's neighbours are the rows at the
# neighbouring locations, row-standardized
Wr <- W[dd$loc, dd$loc]
Wr <- Wr / rowSums(Wr)
dd$y <- as.vector(solve(diag(nrow(dd)) - 0.4 * Wr, 1 + 0.5 * dd$x +
                          rnorm(nrow(dd))))
f <- tryCatch(suppressWarnings(suppressMessages(frm(
  mvbf(bf(time | cens(out) ~ x + car(W, gr = loc, type = "escar")) + cox(df = 4),
       bf(y ~ x + sar(Wr)) + gaussian()),
  data = dd, data2 = list(W = W, Wr = Wr)))), error = identity)
if (inherits(f, "error")) {
  cat("MV FIT ERROR:", conditionMessage(f), "\n")
} else {
  cat("MV FIT logLik", format(as.numeric(logLik(f)), digits = 12), "\n")
  print(fixef(f))
}
# the nearest one-model case: a gaussian response with car() and sar()
f2 <- tryCatch(suppressWarnings(suppressMessages(frm(
  bf(y ~ x + car(W, gr = loc, type = "escar") + sar(Wr)),
  family = gaussian(), data = dd, data2 = list(W = W, Wr = Wr)))),
  error = identity)
if (inherits(f2, "error")) {
  cat("CAR+SAR FIT ERROR:", conditionMessage(f2), "\n")
} else {
  cat("CAR+SAR FIT logLik", format(as.numeric(logLik(f2)), digits = 12),
      "\n")
}

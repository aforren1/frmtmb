# Lane wt-resmooth. The statistic the pp_check() test asserts: where the
# data's spread sits among the replicates' spreads, for the three
# group-indexed smooth constructions and for s(x) + (1 | g).
#   RESMOOTH_LIB=base Rscript dev/resmooth-ppcheck.R > dev/resmooth-ppcheck-before.txt
#   Rscript dev/resmooth-ppcheck.R > dev/resmooth-ppcheck-after.txt
base <- identical(Sys.getenv("RESMOOTH_LIB"), "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"), "\n")
set.seed(5)
n <- 300
d <- data.frame(x = runif(n), g = factor(rep(1:10, length.out = n)))
d$y <- sin(2 * pi * d$x) + rnorm(10, 0, 0.5)[d$g] * d$x +
  rnorm(n, 0, 0.3)

probe <- function(lab, f, nsim = 400) {
  fit <- suppressWarnings(frm(f, data = d))
  s <- as.matrix(simulate(fit, nsim = nsim, seed = 3, re_formula = NA))
  sds <- apply(s, 2, stats::sd)
  sy <- stats::sd(fit$frame$y[[1L]])
  # the pp_check reading: is the data's spread in the tail of the
  # replicates' spreads?
  p <- mean(sds >= sy)
  sig <- as.vector(frm_linpred(fit, dpar = "sigma", type = "response",
                               re_formula = NA))
  if (length(sig) == 1L) sig <- rep(sig, nrow(s))
  rsd <- stats::median(apply(s, 1, stats::sd) / sig)
  cat(sprintf(paste0("%-18s sd(y) %.4f | median sd(yrep) %.4f | ",
                     "pooled ratio %.4f | rowsd/sigma %.4f | tol %.4f | ",
                     "P(sd(yrep) >= sd(y)) %.4f | ",
                     "identical(NA, NULL draws) %s\n"),
              lab, sy, stats::median(sds), stats::sd(as.vector(s)) / sy,
              rsd, 6 / sqrt(2 * (nsim - 1)), p,
              identical(simulate(fit, nsim = 5, seed = 1, re_formula = NA),
                        simulate(fit, nsim = 5, seed = 1))))
  pp <- tryCatch(inherits(pp_check(fit, ndraws = 5), "ggplot"),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("    pp_check() returns a ggplot:", format(pp), "\n")
}
probe("s(g, bs = re)", bf(y ~ s(x) + s(g, bs = "re")))
probe("s(x, g, bs = fs)", bf(y ~ s(x, g, bs = "fs", k = 5)))
probe("t2(x, g, cr+re)", bf(y ~ t2(x, g, bs = c("cr", "re"))))
probe("s(x) + (1 | g)", bf(y ~ s(x) + (1 | g)))
probe("s(x)", bf(y ~ s(x)))

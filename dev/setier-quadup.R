# Lane setier, punch round 2: test-quadrature-defects.R's nested beta
# quadrature fit, whose ga sd the upward probe calls short: the
# objective along that sd, quadrature and Laplace.
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(4)
ng <- 20; nt <- 5; n <- ng * nt
d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                gb = factor(rep(c("i", "ii"), length.out = n)),
                x = rnorm(n))
d$ga <- d$g
d$eta <- 0.5 + 0.4 * d$x + rnorm(ng, 0, 0.4)[d$g] +
  rnorm(ng * 2, 0, 0.3)[as.integer(d$ga) * 2 + as.integer(d$gb) - 2]
set.seed(11)
p <- plogis(d$eta)
d$y <- rbeta(n, p * 5, (1 - p) * 5)
for (q in c(TRUE, FALSE)) {
  f <- suppressMessages(suppressWarnings(
    frm(bf(y ~ 1 + x + (1 | ga/gb)) + Beta(), data = d, quadrature = q)))
  pp <- f$opt$par
  cat("quadrature", q, "par", signif(pp, 4), "logLik",
      format(as.numeric(logLik(f)), digits = 10), "\n")
  print(ns$sdr_of(f)$se_lost)
  j <- grep("theta", names(pp))
  f0 <- f$obj$fn(pp)
  for (k in j) for (v in c(pp[k] + c(-2, 2), log(c(0.01, 0.03, 0.1, 0.3, 1)))) {
    qq <- pp
    qq[k] <- v
    cat("  theta", k, "at", signif(v, 3), "dnll",
        signif(f$obj$fn(qq) - f0, 4), "\n")
  }
}

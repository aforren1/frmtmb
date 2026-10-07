# Reviewer of lane setier: the cost of separation_check() at fit end on
# large bernoulli fits, dense and sparse_x designs. Times the check alone
# (minimum of 3 calls) against the fit, and the memory it allocates.
#   Rscript dev/setier-rev-sepcost.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, find.package("frmtmb"), "\n")
tm <- function(expr) {
  gc(reset = TRUE)
  t <- system.time(v <- expr)[["elapsed"]]
  g <- gc()
  list(v = v, t = t, mb = sum(g[, ncol(g)]))
}
cases <- list(
  dense_p20_n1e5 = function() {
    set.seed(1)
    n <- 1e5
    X <- matrix(rnorm(n * 19), n)
    d <- data.frame(X)
    d$y <- rbinom(n, 1, plogis(-2 + X %*% c(3, rep(0.2, 18))))
    list(d = d, f = as.formula(paste("y ~", paste(names(d)[1:19],
                                                  collapse = " + "))),
         ctl = frmtmb_control())
  },
  sparse_f500_n1e5 = function() {
    set.seed(2)
    n <- 1e5
    d <- data.frame(f = factor(sample(500, n, TRUE)), x = rnorm(n))
    d$y <- rbinom(n, 1, plogis(-2 + rnorm(500, 0, 1.5)[d$f] + 2 * d$x))
    list(d = d, f = y ~ f + x, ctl = frmtmb_control(sparse_x = TRUE))
  })
for (nm in names(cases)) {
  cs <- cases[[nm]]()
  r <- tm(suppressWarnings(frm(cs$f, family = bernoulli(), data = cs$d,
                                control = cs$ctl)))
  fit <- r$v
  cat(sprintf("%s: fit %.2f s, max mem %.0f MB, code %d, X class %s\n", nm,
              r$t, r$mb, fit$opt$convergence,
              class(fit$frame$linpreds[[1]]$X)[1]))
  if (exists("separation_check", ns)) {
    ts <- numeric(3)
    for (i in 1:3) {
      s <- tm(ns$separation_check(fit))
      ts[i] <- s$t
    }
    cat(sprintf("   separation_check: min %.2f s (of %s), max mem %.0f MB, named %s\n",
                min(ts), paste(round(ts, 2), collapse = ","), s$mb,
                !is.null(s$v)))
    lp <- fit$frame$linpreds[[1]]
    eta <- as.numeric(lp$X %*% fit$estimates$beta[lp$idx])
    mu <- plogis(eta)
    cat("   obs with a tail <= 1e-3:", sum(pmin(mu, 1 - mu) <= 1e-3), "\n")
  }
}

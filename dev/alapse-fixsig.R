# Lane alapse, fix round 1: with a FIXED quantile and sigma ~ z, is the
# proved fit the global maximum? Multi-start alternating maximization
# (rq.wfit() at weights 1 / sigma_i, then BFGS on sigma's coefficients),
# 15 starts of (sigma intercept, slope) on a grid, against frm()'s fit.
# Design D of dev/alapse-rq.R (200 rows, sd exp(0.4 z)) and a strongly
# heteroscedastic one (sd exp(1.2 z)), quantile 0.5 and 0.2, seeds 1 to 20.
#   Rscript dev/alapse-fixsig.R <base|lane> [seeds]
a <- commandArgs(TRUE)
lib_arg <- a[1]
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:20
source("dev/alapse-common.R")
suppressMessages(library(quantreg))
alt <- function(y, X, Z, p, bs) {
  ll_of <- function(b, bs) ald_ll(y, drop(X %*% b), exp(drop(Z %*% bs)), p)
  cur <- -Inf
  for (it in 1:300) {
    b <- suppressWarnings(quantreg::rq.wfit(X, y, tau = p,
      weights = 1 / exp(drop(Z %*% bs)), method = "br"))$coefficients
    o <- optim(bs, function(t) -ll_of(b, t), method = "BFGS",
               control = list(reltol = 1e-14, maxit = 2000))
    bs <- o$par
    new <- ll_of(b, bs)
    if (abs(new - cur) < 1e-11) { cur <- new; break }
    cur <- new
  }
  cur
}
starts <- expand.grid(i = c(-1, 0, 1), s = c(-2, -1, 0, 1, 2))
res <- NULL
for (des in c("D", "H")) for (p in c(0.5, 0.2)) for (s in seeds) {
  set.seed(4000 + s)
  d <- data.frame(x = rnorm(200), z = rnorm(200))
  k <- if (des == "D") 0.4 else 1.2
  d$y <- 1 + 0.5 * d$x + exp(k * d$z) * rnorm(200)
  X <- cbind(1, d$x)
  Z <- cbind(1, d$z)
  r <- cap(frm(bf(y ~ x, sigma ~ z, quantile = p), family = asym_laplace(),
               data = d))
  f <- r$v
  ll <- as.numeric(logLik(f))
  best <- max(vapply(seq_len(nrow(starts)), function(j) {
    alt(d$y, X, Z, p, c(starts$i[j], starts$s[j]))
  }, 0))
  cat(sprintf("%s p %.1f %2d | code %d cert %-5s | ll %.6f | multistart %.6f | gap %.3e | %s\n",
              des, p, s, f$opt$convergence, isTRUE(f$opt$ald_cert$ok), ll,
              best, best - ll, paste(unique(sub("^frmtmb_", "", r$wc)),
                                     collapse = ",")))
  res <- rbind(res, data.frame(des = des, p = p, gap = best - ll))
}
for (k in split(res, list(res$des, res$p))) {
  cat(sprintf("SUMMARY %s p %.1f: n %d | gap > 1e-6: %d | max gap %.3e\n",
              k$des[1], k$p[1], nrow(k), sum(k$gap > 1e-6), max(k$gap)))
}
cat("DONE\n")

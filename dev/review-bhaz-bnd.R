# Reviewer: what brms 2.23.0's extrapolated basis does past boundary
# knots that sit inside the data, on the lane's probe data (seed 11, g
# drawn before the times, as dev/bhaz-brms-probe.R).
#   Rscript dev/review-bhaz-bnd.R > dev/review-bhaz-log/bnd.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(11); n <- 200; x <- rnorm(n)
g <- factor(sample(c("a", "b", "c"), n, TRUE))
tt <- rexp(n, exp(-0.5 + 0.7 * x)); ct <- rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 g = g)
cat("range(time)", format(range(dd$time), digits = 10), "\n")
for (bk in list(c(0.5, 10), c(0.01, 14.71), c(0, 12), c(0, 14.3))) {
  fo <- eval(bquote(time | cens(cens) + bhaz(Boundary.knots = .(bk)) ~ x))
  sd <- suppressWarnings(standata(fo, data = dd, family = cox()))
  Z <- sd$Zbhaz; C <- sd$Zcbhaz
  lo <- dd$time < bk[1]; hi <- dd$time > bk[2]
  K <- ncol(Z)
  # over all simplexes, the hazard of a row is between its smallest and
  # largest basis entry; it is negative for EVERY simplex when the
  # largest entry is negative, for SOME when the smallest is
  cat(sprintf(paste0("bk %s: below %d above %d | Zbhaz rows: some-s ",
                     "negative %d, every-s negative %d | Zcbhaz rows: ",
                     "some-s negative %d, every-s negative %d | rows with ",
                     "a negative entry %d\n"),
              paste(bk, collapse = ","), sum(lo), sum(hi),
              sum(apply(Z, 1, min) < 0), sum(apply(Z, 1, max) < 0),
              sum(apply(C, 1, min) < 0), sum(apply(C, 1, max) < 0),
              sum(apply(Z < 0, 1, any))))
  # is the cumulative hazard monotone along sorted times for the
  # uniform simplex?
  o <- order(dd$time)
  H <- as.numeric(C %*% rep(1 / K, K))[o]
  cat(sprintf("   uniform simplex: min H %.4g, decreasing steps %d\n",
              min(H), sum(diff(H) < -1e-12)))
}

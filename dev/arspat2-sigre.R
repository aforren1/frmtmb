# A side question of backlog 5289: gaussian() with fcor(), mi() and a
# random effect in sigma. The joint Laplace over (missing values, b) is the
# nested one only when the missing values' Hessian does not move with b,
# and a random effect in sigma moves it. The mi() fit's objective against
# the fit over the observed rows (fcor(M_oo)) at three outer points.
# Data seed 9, n = 60, m = 15, 6 groups. Lane build.
lib <- Sys.getenv("ARSPAT2_LIB", "C:/Users/adf44/source/r/wt-arspat2-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(dirname(system.file(package = "frmtmb"))), "\n")
n <- 60
M <- 0.5^abs(outer(1:n, 1:n, "-"))
set.seed(9)
d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 10)))
d$y <- 1 + 0.5 * d$x + exp(0.6 * rnorm(6)[d$h]) *
  as.vector(t(chol(M)) %*% rnorm(n))
d$y[sort(sample(n, 15))] <- NA
ok <- !is.na(d$y)
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
for (fam in list(gaussian(), student())) {
  f <- quiet(frm(bf(y | mi() ~ x + fcor(M), sigma ~ 1 + (1 | h)),
                 family = fam, data = d, data2 = list(M = M)))
  g <- quiet(frm(bf(y ~ x + fcor(Mo), sigma ~ 1 + (1 | h)), family = fam,
                 data = d[ok, ], data2 = list(Mo = M[ok, ok])))
  lab <- frmtmb:::outer_par_labels(f)
  cat("\n===", fam$family, ": labels", lab, "\n")
  isd <- grep("^theta", names(f$opt$par))
  for (k in 0:2) {
    p <- f$opt$par
    p[isd] <- log(c(0.3, 0.6, 1.2)[k + 1L])
    a <- -f$obj$fn(p)
    b <- -g$obj$fn(p)
    cat(sprintf("  sd(sigma) %.1f: -fn mi %.8f  observed rows %.8f  diff %+.3e\n",
                c(0.3, 0.6, 1.2)[k + 1L], a, b, a - b))
  }
}

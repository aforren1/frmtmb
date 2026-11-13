# escorder seed 6: the stall start ends 5.5e-2 short with no escape on
# both lane arms. Same fit under Laplace, and on the base build.
#   Rscript dev/quadseam-seed6.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(6); n <- 200
xs <- -abs(rnorm(n)) * 3
y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
dd <- data.frame(y = y, xs = xs, g = factor(rep(1:20, length.out = n)))
m3 <- mean((y - mean(y))^3) / sd(y)^3
bad <- list(betad = c(log(sd(y)), 2 * sign(m3) + 0.5 * m3))
fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
for (q in c(FALSE, TRUE)) for (st in list(NULL, bad)) {
  f <- suppressMessages(suppressWarnings(frm(fo, family = skew_normal(),
         data = dd, quadrature = q, start = st)))
  cat(sprintf("quad %-5s start %-7s ll %.9f alpha %.4f esc %s\n", q,
              if (is.null(st)) "default" else "stall", as.numeric(logLik(f)),
              f$estimates$betad[2],
              if (is.null(f$opt$stationary_escape)) "-" else "yes"))
}

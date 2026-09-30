# Punch 2: the worst points of the dev/fams2-p1-sweep.R draw (seed
# 20260930) on the installed build, with where each sits. The reference
# is stats::pbeta(), which the review found within 7e-14 of Rmpfr there.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

f <- frmtmb:::log_ibeta_half
set.seed(20260930); M <- 6000
a <- exp(runif(M, log(1e-3), log(1e7))); b <- exp(runif(M, log(1e-3), log(1e7)))
s <- a + b; m <- (a + 1) / (s + 2); sd <- sqrt(a * b / (s * s * (s + 1)))
u <- runif(M)
x <- ifelse(u < 0.1, m, ifelse(u < 0.6, m + runif(M, -5, 5) * sd, runif(M, 0, 0.5)))
ok <- x > 0 & x < 0.5
d <- data.frame(x, a, b, m, sd)[ok, ]
d$ref <- suppressWarnings(pbeta(d$x, d$a, d$b, log.p = TRUE))
d <- d[is.finite(d$ref) & d$ref > -600, ]
d$ours <- mapply(f, d$x, d$a, d$b)
d$err <- abs(d$ours - d$ref) / pmax(1, abs(d$ref))
d$k <- (d$x - d$m) / d$sd
d$abs <- d$a * d$b / (d$a + d$b)
w <- head(d[order(-d$err), ], 12)
print(signif(w[, c("a", "b", "abs", "k", "ref", "err")], 3))

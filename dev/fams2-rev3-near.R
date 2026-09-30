# Reviewer, punch round 2: x-derivatives near the mean at large shapes.
# Are the third-derivative errors of log_ibeta_half() there RTMB::pbeta()'s
# own (plain log(pbeta) taped), and how do they scale? 256-bit Rmpfr.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(RTMB); library(Rmpfr) })
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-rev3-ibeta.R")
eval(parse(text = src[grep("^prec <- 256", src):(grep("^Tn <- MakeTape", src) - 1)]))
lib <- frmtmb:::log_ibeta_half
mk <- function(f) { F <- MakeTape(function(p) f(p[1], p[2], p[3]), c(0.2, 3, 7)); J <- F$jacfun(); H <- J$jacfun(); list(F = F, J = J, H = H, T = H$jacfun()) }
O <- mk(lib); P <- mk(function(x, a, b) log(RTMB::pbeta(x, a, b)))
xd <- function(M, p) c(M$J(p)[1], M$H(p)[1, 1], M$T(p)[1, 1])
for (ab in list(c(300, 700), c(2e3, 8e3), c(3e4, 7e4), c(1e6, 3e6))) {
  a <- ab[1]; b <- ab[2]; s <- a + b; t1 <- a / s; gap <- b / (s * (s + 1))
  for (nm in c("t1+0.5gap", "t1+0.75gap", "(a+1)/(s+2)", "t1+2gap", "t1+10gap", "t1-10gap")) {
    x <- switch(nm, "t1+0.5gap" = t1 + 0.5 * gap, "t1+0.75gap" = t1 + 0.75 * gap,
                "(a+1)/(s+2)" = (a + 1) / (s + 2), "t1+2gap" = t1 + 2 * gap,
                "t1+10gap" = t1 + 10 * gap, "t1-10gap" = t1 - 10 * gap)
    r <- mp_der(x, a, b)[2:4]
    p <- c(x, a, b)
    eo <- abs(xd(O, p) - r) / abs(r); ep <- abs(xd(P, p) - r) / abs(r)
    cat(sprintf("a %g b %g %-12s rel err d1,d2,d3: package %s | plain log(pbeta) %s\n", a, b, nm,
                paste(sprintf("%.1e", eo), collapse = " "), paste(sprintf("%.1e", ep), collapse = " ")))
  }
}

# Review of lane famend: where the GEV tape log density leaves the exact
# one in the far tails (pre-existing; base and lane alike).
#   Rscript dev/reviews/famend-review-gevtail.R <base|lane>
a <- commandArgs(trailingOnly = TRUE)
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] == "lane") lanelib, base, ulib))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
fam <- gen_extreme_value()
tapev <- function(y, xi) {
  p <- c(0, 0, log1p(xi))
  tp <- RTMB::MakeTape(function(p) {
    fam$lpdf(y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                     .eta_sigma = p[2], .eta_xi = p[3]), list())
  }, p)
  c(tp(p), tp$jacobian(p))
}
ex <- function(x, xi) evd::dgev(x, 0, 1, xi, log = TRUE)
for (xi in c(-0.9, -0.5, -0.2, 0.2, 0.5, 1)) {
  xs <- if (xi < 0) -10^seq(0, 20, by = 0.5) else 10^seq(0, 40, by = 1)
  for (x in xs) {
    v <- tapev(x, xi)
    e <- ex(x, xi)
    flag <- if (is.finite(e) && abs(v[1] - e) > 1e-6 * max(1, abs(e))) " <== off" else ""
    cat(sprintf("xi %5.2f x %10.3g tape %12.6g evd %12.6g grad %s%s\n", xi, x,
                v[1], e, paste(signif(v[-1], 3), collapse = " "), flag))
  }
}
# a fit: does a far row pull the optimum away?
set.seed(3)
rg <- function(n, mu, s, xi) mu + s * ((-log(runif(n)))^(-xi) - 1) / xi
y <- rg(200, 0, 1, -0.3)
r <- tryCatch(frm(y ~ 1, family = gen_extreme_value(), data = data.frame(y = y)),
              error = function(e) conditionMessage(e))
if (!is.character(r)) cat("clean fit: logLik", as.numeric(logLik(r)), "par",
                          signif(r$opt$par, 6), "\n")
cat("\ndone\n")

# Lane optima, item 2 (found on the way): the cloglog log-odds,
# q(eta) = log(1 - exp(-exp(eta))) + exp(eta), against Rmpfr at 256
# bits, value and derivative through RTMB's tape, for the 0.68.1 form
# and the lane's form (clamped at -40).
#   Rscript dev/optima-cloglog.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(Rmpfr)})
old <- function(eta) {
  t <- exp(eta)
  log(-expm1(-t)) + t
}
new <- frmtmb:::frmtmb_links$cloglog$logit_eta
ref <- function(x) {
  X <- mpfr(x, 256)
  t <- exp(X)
  # -expm1(-t), not 1 - exp(-t), which 256 bits round to 0 below -177
  list(q = as.numeric(log(-expm1(-t)) + t),
       dq = as.numeric(t * exp(-t) / (-expm1(-t)) + t))
}
rel <- function(a, b) ifelse(is.finite(a), abs(a / b - 1), Inf)
for (rg in list(c(-40, 6.5), c(-708, -40), c(-745, -708), c(-1000, -745))) {
  x <- seq(rg[1], rg[2], length.out = 1001)
  R <- ref(x)
  cat(sprintf("eta in [%g, %g]:\n", rg[1], rg[2]))
  for (nm in c("old", "new")) {
    f <- get(nm)
    tp <- RTMB::MakeTape(function(e) f(e), x)
    v <- tp(x)
    g <- diag(tp$jacobian(x))
    ev <- rel(v, R$q)
    eg <- rel(g, R$dq)
    cat(sprintf("  %-3s value: non-finite %4d, max rel err %.3g | derivative: non-finite %4d, max rel err %.3g\n",
                nm, sum(!is.finite(v)), max(ev[is.finite(ev)], -Inf),
                sum(!is.finite(g)), max(eg[is.finite(eg)], -Inf)))
  }
}
x <- seq(-40, 6.5, length.out = 1001)
vo <- RTMB::MakeTape(old, x)(x)
vn <- RTMB::MakeTape(new, x)(x)
cat("old vs new on [-40, 6.5]:", sum(vo != vn), "of", length(x),
    "values differ\n")
# the softit's log-odds, log(log1p(exp(eta))), the same way
cat("== softit\n")
old <- function(eta) log(RTMB::logspace_add(0 * eta, eta))
new <- frmtmb:::frmtmb_links$softit$logit_eta
ref <- function(x) {
  X <- mpfr(x, 256)
  s <- log1p(exp(X))
  list(q = as.numeric(log(s)),
       dq = as.numeric(exp(X) / (1 + exp(X)) / s))
}
for (rg in list(c(-40, 30), c(-708, -40), c(-745, -708), c(-1000, -745))) {
  x <- seq(rg[1], rg[2], length.out = 1001)
  R <- ref(x)
  cat(sprintf("eta in [%g, %g]:\n", rg[1], rg[2]))
  for (nm in c("old", "new")) {
    f <- get(nm)
    tp <- RTMB::MakeTape(function(e) f(e), x)
    v <- tp(x)
    g <- diag(tp$jacobian(x))
    ev <- rel(v, R$q)
    eg <- rel(g, R$dq)
    cat(sprintf("  %-3s value: non-finite %4d, max rel err %.3g | derivative: non-finite %4d, max rel err %.3g\n",
                nm, sum(!is.finite(v)), max(ev[is.finite(ev)], -Inf),
                sum(!is.finite(g)), max(eg[is.finite(eg)], -Inf)))
  }
}
x <- seq(-40, 30, length.out = 1001)
cat("old vs new on [-40, 30]:", sum(RTMB::MakeTape(old, x)(x) !=
                                       RTMB::MakeTape(new, x)(x)),
    "of", length(x), "values differ\n")

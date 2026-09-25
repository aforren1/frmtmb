# Reviewer 2, item 1: B1. The tape derivative of the contaminant's log
# density, log survival and distribution function in eta = logit(lambda),
# against Richardson central differences, from lambda near 0 to near 1,
# on rows below ndt (max_ndt above the fastest response), at the window's
# edges, just outside them, and far in the tail.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(RTMB); library(frmtmb); library(frmtmb.eam)})
options(width = 160)
set.seed(48)
d <- ddm_simulate(200, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
cr <- c(0.2, 3)
edge <- c(0.2, 3, 0.2 - 1e-9, 3 + 1e-9, 0.2 + 1e-12, 0.12, 0.15, 0.22, 0.29, 2.99999, 7)
d <- rbind(d, data.frame(rt = edge, upper = rep(0:1, length.out = length(edge))))
o <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5, lambda = 0.1),
         family = wiener(contaminant = TRUE, contaminant_range = cr, max_ndt = 0.5),
         data = d, dry_run = "objective")
fam <- stats::family(o)
cat("window", fam$contaminant_range, "; rows below ndt 0.3:", sum(d$rt < 0.3), "\n")
dp <- function(eta) list(mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5,
                         lambda = 1 / (1 + exp(-eta)), .eta_lambda = eta)
ae <- list(dec = d$upper)
qs <- c(0.5, 1, 2.9, 3, 3.5)   # right-censoring points, one at and past the top
ae_r <- list(dec = rep(0, length(qs)), cens = rep(1, length(qs)))
ae_l <- list(dec = rep(1, length(qs)), cens = rep(-1, length(qs)))
fns <- list(
  lpdf = function(eta) sum(fam$lpdf(d$rt, dp(eta), ae)),
  lccdf = function(eta) sum(fam$lccdf(qs, dp(eta), ae_r)),
  lcdf_left = function(eta) sum(log(fam$lcdf(qs, dp(eta), ae_l))))
rich <- function(f, x, h) {
  D <- function(h) (f(x + h) - f(x - h)) / (2 * h)
  (4 * D(h / 2) - D(h)) / 3
}
for (nm in names(fns)) {
  f <- fns[[nm]]
  tp <- MakeTape(function(e) f(e), 0)
  cat("\n==", nm, "\n")
  cat(sprintf("%6s %12s %16s %16s %16s %10s %10s\n", "eta", "lambda", "value", "tape", "fd h=1e-4",
              "rel", "d2 finite"))
  for (eta in c(-30, -20, -12, -6, -3, -1, 0, 1, 3, 6, 12, 20, 30)) {
    g <- tp$jacobian(eta); f1 <- rich(f, eta, 1e-4); f2 <- rich(f, eta, 1e-6)
    h2 <- tp$jacfun()$jacobian(eta)
    cat(sprintf("%6g %12.4g %16.10g %16.10g %16.10g %10.2e %10s  (fd h=1e-6 %.10g)\n", eta,
                plogis(eta), f(eta), g, f1, abs(g - f1) / max(abs(f1), 1), is.finite(h2), f2))
  }
}

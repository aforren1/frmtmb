# Reviewer 2, item 1: is the staircase gone only in lambda? The value
# lf + m with m = lg - lf rounds lg at the ulp of lf (about 2.4e-7 for
# lf near -1.8e9) on a row below ndt. As ndt, bs or mu move, lf moves by
# many ulps, so the objective should now step in THOSE directions.
# The worker's own construction (test-contaminant.R "smooth in lambda":
# ct_data(400, seed 48), five rows moved below ndt, max_ndt 0.5), with
# the difference taken in ndt, bs and mu instead of lambda.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
set.seed(48)
d <- ddm_simulate(400L, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
d$rt[1:5] <- c(0.12, 0.15, 0.18, 0.2, 0.22)
fit <- suppressWarnings(frm(bf(rt | dec(upper) ~ 1, bias = 0.5, lambda = 0.1),
  family = wiener(contaminant = TRUE, max_ndt = 0.5,
                  contaminant_range = range(d$rt)), data = d, dry_run = "objective"))
fam <- stats::family(fit)
base <- list(mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5, lambda = plogis(-3), .eta_lambda = -3)
obj <- function(nm, x) { dp <- base; dp[[nm]] <- x; sum(fam$lpdf(d$rt, dp, list(dec = d$upper))) }
fd <- function(nm, h) { x0 <- base[[nm]]; (obj(nm, x0 + h) - obj(nm, x0 - h)) / (2 * h) }
cat("rows below ndt:", sum(d$rt < 0.3), "\n")
cat(sprintf("%-4s %16s %16s %16s %12s\n", "par", "fd h=1e-4", "fd h=1e-6", "fd h=1e-8", "|1e-8 - 1e-4|"))
for (nm in c("ndt", "bs", "mu")) {
  a <- fd(nm, 1e-4); b <- fd(nm, 1e-6); c8 <- fd(nm, 1e-8)
  cat(sprintf("%-4s %16.8f %16.8f %16.8f %12.4g\n", nm, a, b, c8, abs(c8 - a)))
}
eta <- function(h) { dp <- base; dp$.eta_lambda <- -3 + h; dp$lambda <- plogis(-3 + h)
  sum(fam$lpdf(d$rt, dp, list(dec = d$upper))) }
cat(sprintf("eta  %16.8f %16.8f %16.8f  (lambda direction, fixed by B1)\n",
            (eta(1e-4) - eta(-1e-4)) / 2e-4, (eta(1e-6) - eta(-1e-6)) / 2e-6,
            (eta(1e-8) - eta(-1e-8)) / 2e-8))
# the per-row value of a below-ndt row as ndt moves by 1e-9 steps
dp <- base; y <- d$rt[1:5]
vals <- sapply(0:6, function(k) { dp$ndt <- 0.3 + k * 1e-9; fam$lpdf(y, dp, list(dec = d$upper[1:5])) })
cat("\nrow values minus their first column, ndt stepped by 1e-9 (true change: 0):\n")
print(signif(vals - vals[, 1], 3))
# a smooth control: rows above ndt only
dp <- base; y2 <- d$rt[d$rt > 0.35][1:5]
vals2 <- sapply(0:6, function(k) { dp$ndt <- 0.3 + k * 1e-9; fam$lpdf(y2, dp, list(dec = d$upper[d$rt > 0.35][1:5])) })
cat("control rows above ndt, same steps:\n"); print(signif(vals2 - vals2[, 1], 3))

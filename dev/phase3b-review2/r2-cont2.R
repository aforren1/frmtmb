# Reviewer 2, item 4 (continued): responses ABOVE the given range, and
# the default window under trunc(ub = ) against a hand-written mixture.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam); library(RWiener)})
set.seed(21)
d <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
hit <- runif(300) < 0.08
d$rt[hit] <- runif(sum(hit), 0.1, 2.5); d$upper[hit] <- rbinom(sum(hit), 1, 0.5)
d <- d[d$rt < 2.5, ]
mk <- function(f) bf(f, bs ~ 1, ndt ~ 1, lambda ~ 1, bias = 0.5)
o <- frm(mk(rt | dec(upper) ~ 1), family = wiener(contaminant = TRUE,
         contaminant_range = c(0, 1.0), max_ndt = 0.5), data = d, dry_run = "objective")
cat("rows above the range top 1.0:", sum(d$rt > 1), "\n")
fam <- stats::family(o)
dp <- list(mu = 0.8, bs = 1.4, ndt = 0.12, bias = 0.5, lambda = 0.05,
           .eta_lambda = qlogis(0.05))
lp <- fam$lpdf(d$rt, dp, list(dec = d$upper))
pl <- frmtmb.eam:::ddm_lpdf_both(d$rt - 0.12, 0.8, 1.4, 0.5, d$upper)
k <- d$rt > 1
cat(sprintf("rows above: log density - (log(0.95) + plain) max abs %.3e; any non-finite %s\n",
            max(abs(lp[k] - (log(0.95) + pl[k]))), any(!is.finite(lp))))
cat("fn at start", o$obj$fn(o$obj$par), " grad finite", all(is.finite(o$obj$gr(o$obj$par))), "\n")

ot <- frm(mk(rt | dec(upper) + trunc(ub = 2.5) ~ 1), family = wiener(contaminant = TRUE),
          data = d, dry_run = "objective")
ft <- stats::family(ot); ub <- ft$ndt_bound$ub; cr <- ft$contaminant_range
cat("default window", cr, " min rt", min(d$rt), "\n")
hand <- function(p) {
  v <- p[1]; a <- exp(p[2]); t0 <- ub * plogis(p[3]); lam <- plogis(p[4])
  resp <- ifelse(d$upper == 1, "upper", "lower")
  f <- mapply(function(y, r) dwiener(y, a, t0, 0.5, v, resp = r), d$rt, resp)
  Fub <- pwiener(2.5, a, t0, 0.5, v, resp = "both")
  -sum(log((1 - lam) * f + lam / (2 * (cr[2] - cr[1]))) - log((1 - lam) * Fub + lam))
}
for (p in list(ot$obj$par, ot$obj$par + c(0.2, 0.1, 0.3, 0.5))) {
  cat(sprintf("trunc(ub) default window: lane %.10f  hand %.10f  rel %.2e\n",
              ot$obj$fn(p), hand(p), abs(ot$obj$fn(p) - hand(p)) / abs(hand(p))))
}

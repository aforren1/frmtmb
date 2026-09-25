# Reviewer 2, item 1, per row: tape d/deta of each row's log density
# against the analytic value where the row's value is dominated by one
# part (below ndt), and against a Richardson difference elsewhere.
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
dp <- function(eta) list(mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5,
                         lambda = 1 / (1 + exp(-eta)), .eta_lambda = eta)
ae <- list(dec = d$upper)
tp <- MakeTape(function(e) fam$lpdf(d$rt, dp(e), ae), 0)
below <- d$rt <= 0.3
inwin <- d$rt >= cr[1] & d$rt <= cr[2]
worst <- c(in_below = 0, out_below = 0, other = 0)
steps <- 0
for (eta in c(-30, -20, -12, -6, -3, -1, 0, 1, 3, 6, 12, 20, 30)) {
  lam <- plogis(eta)
  J <- as.numeric(tp$jacobian(eta))
  H <- as.numeric(tp$jacfun()$jacobian(eta))
  if (any(!is.finite(J)) || any(!is.finite(H))) cat("NON-FINITE at eta", eta, "\n")
  worst["in_below"] <- max(worst["in_below"], abs(J[below & inwin] - (1 - lam)))
  worst["out_below"] <- max(worst["out_below"], abs(J[below & !inwin] - (-lam)))
  f <- function(e) fam$lpdf(d$rt, dp(e), ae)[!below]
  fd <- ((4 * (f(eta + 5e-5) - f(eta - 5e-5)) / 1e-4) - (f(eta + 1e-4) - f(eta - 1e-4)) / 2e-4) / 3
  worst["other"] <- max(worst["other"], abs(J[!below] - fd) / pmax(abs(fd), 1))
  # the staircase: the value of an in-window below-ndt row against
  # log(lambda) + log g, to the bit scale of the answer
  v <- fam$lpdf(d$rt, dp(eta), ae)[below & inwin]
  steps <- max(steps, max(abs(v - (log(lam) - log(2 * diff(cr))))))
}
cat("rows below ndt inside the window:", sum(below & inwin), " outside:", sum(below & !inwin), "\n")
print(worst)
cat("max |value - (log lambda + log g)| on in-window rows below ndt:", steps, "\n")

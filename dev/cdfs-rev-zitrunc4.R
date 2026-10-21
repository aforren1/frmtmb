# Reviewer: zi_cdf_log() sizes its rows by length(q). A trunc() bound
# given as a constant is one value, so the point mass is added to row
# 1 only and every other row's normalizer lacks it. Compare a constant
# bound with the same bound as a data column, and with the hand-rolled
# likelihood, at fixed parameters (the objective at the start).
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
ns <- asNamespace("frmtmb")
set.seed(5)
n <- 60
x <- rnorm(n)
d <- data.frame(x = x, yp = rpois(n, exp(1 + 0.3 * x)), tr = 12)
d$ub30 <- 30
d$ub10 <- 10
d$lb1 <- 1
d$yb <- rbinom(n, 12, 0.3)
d$ub9 <- 9
d$ypos <- pmax(d$yp, 1)
d$yz <- ifelse(runif(n) < 0.2, 0, rbeta(n, 2, 3))
d$ub95 <- 0.95
d$yz <- pmin(d$yz, 0.94)
cmp <- function(lab, f_const, f_col, fam) {
  o1 <- frm(f_const, family = fam, data = d, dry_run = "objective")$obj
  o2 <- frm(f_col, family = fam, data = d, dry_run = "objective")$obj
  p <- o2$par
  cat(sprintf("%-46s -fn: constant bound %.8f | column bound %.8f | diff %.3g\n",
              lab, -o1$fn(p), -o2$fn(p), o2$fn(p) - o1$fn(p)))
}
cmp("zero_inflated_poisson trunc(ub = 30)", yp | trunc(ub = 30) ~ x,
    yp | trunc(ub = ub30) ~ x, zero_inflated_poisson())
cmp("zero_inflated_poisson trunc(ub = 10)", yp | trunc(ub = 10) ~ x,
    yp | trunc(ub = ub10) ~ x, zero_inflated_poisson())
cmp("zero_inflated_negbinomial trunc(ub = 30)", yp | trunc(ub = 30) ~ x,
    yp | trunc(ub = ub30) ~ x, zero_inflated_negbinomial())
cmp("zero_inflated_binomial trunc(ub = 9)", yb | trials(tr) + trunc(ub = 9) ~ x,
    yb | trials(tr) + trunc(ub = ub9) ~ x, zero_inflated_binomial())
cmp("zero_inflated_beta_binomial trunc(ub = 9)",
    yb | trials(tr) + trunc(ub = 9) ~ x,
    yb | trials(tr) + trunc(ub = ub9) ~ x, zero_inflated_beta_binomial())
cmp("zero_inflated_poisson trunc(lb = 1)", ypos | trunc(lb = 1) ~ x,
    ypos | trunc(lb = lb1) ~ x, zero_inflated_poisson())
cmp("zero_inflated_beta trunc(ub = 0.95)", yz | trunc(ub = 0.95) ~ x,
    yz | trunc(ub = ub95) ~ x, zero_inflated_beta())
cmp("hurdle_poisson trunc(ub = 10)", yp | trunc(ub = 10) ~ x,
    yp | trunc(ub = ub10) ~ x, hurdle_poisson())
cmp("poisson trunc(ub = 10)", yp | trunc(ub = 10) ~ x,
    yp | trunc(ub = ub10) ~ x, poisson())
cat("== hand-rolled ZIP trunc(ub = 10) at the column fit's start\n")
o2 <- frm(yp | trunc(ub = ub10) ~ x, family = zero_inflated_poisson(), data = d,
          dry_run = "objective")$obj
p <- o2$par
mu <- exp(p[1] + p[2] * d$x)
zi <- plogis(p[3])
ld <- ifelse(d$yp == 0, log(zi + (1 - zi) * exp(-mu)),
             log1p(-zi) + dpois(d$yp, mu, log = TRUE))
Fub <- zi + (1 - zi) * ppois(10, mu)
cat(sprintf("hand %.8f  column %.8f\n", sum(ld - log(Fub)), -o2$fn(p)))
cat("== the fits\n")
for (f in list(yp | trunc(ub = 10) ~ x, yp | trunc(ub = ub10) ~ x)) {
  r <- tryCatch({
    ft <- suppressWarnings(frm(f, family = zero_inflated_poisson(), data = d))
    sprintf("logLik %.6f  zi %.4f", as.numeric(logLik(ft)),
            plogis(ft$opt$par[3]))
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 70)))
  cat(deparse(f), ":", r, "\n")
}
cat("== simulate() / posterior_predict-style draw under a constant bound\n")
ft <- frm(yp | trunc(ub = ub10) ~ x, family = zero_inflated_poisson(), data = d)
s <- tryCatch(simulate(ft, nsim = 1, seed = 1), error = function(e) conditionMessage(e))
if (!is.character(s)) cat("simulate max:", max(unlist(s)), "\n") else cat(s, "\n")

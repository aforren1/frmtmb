# rel084 release review, alapse follow-up: on a box that binds, is the fit
# at the constrained maximum (fixed-quantile / fixed-sigma fits at the bound
# value), and is the generic gradient warning there the same as on a
# gaussian fit with a binding bound (pre-existing behavior)?
#   REVLIB=r23|r22 Rscript dev/reviews/rel084-review-alapse2.R
source("dev/reviews/rel084-review-common.R")
options(frmtmb.notices = FALSE, width = 160)
ns <- asNamespace("frmtmb")
set.seed(8411)
n <- 120
d <- data.frame(x = rnorm(n))
d$y <- 0.5 + 0.8 * d$x + rexp(n) - 1 + 0.3 * rnorm(n)
al <- asym_laplace()
q0 <- -1.8558
s0 <- -2.1059
ll <- function(r) if (is.null(r$v)) NA else as.numeric(logLik(r$v))
rb <- conds_of(frm(bf(y ~ x, quantile ~ 1), family = al, data = d,
                   prior = set_prior("", class = "Intercept", dpar = "quantile",
                                     lb = q0 - 3, ub = q0 - 0.5)))
pb <- ns$eval_dpars(rb$v)[[1]]$quantile[1]
rf <- conds_of(frm(bf(y ~ x, quantile = pb), family = al, data = d))
cat(sprintf("quantile bound: ll %.8f code %d [%s]; fixed q=%.6f ll %.8f code %d [%s]\n",
            ll(rb), rb$v$opt$convergence, cls_line(rb), pb, ll(rf),
            rf$v$opt$convergence, cls_line(rf)))
cat("  conv_gradient text:", substr(msg_of(rb, "frmtmb_conv_gradient"), 1, 200), "\n")
rs <- conds_of(frm(bf(y ~ x, quantile ~ 1, sigma ~ 1), family = al, data = d,
                   prior = set_prior("", class = "Intercept", dpar = "sigma",
                                     lb = s0 + 0.3, ub = 3)))
sb <- ns$eval_dpars(rs$v)[[1]]$sigma[1]
cat(sprintf("sigma bound: ll %.8f code %d sigma %.6f log %.6f bound %.6f [%s]\n",
            ll(rs), rs$v$opt$convergence, sb, log(sb), s0 + 0.3, cls_line(rs)))
# profile over the quantile at sigma fixed at the bound: grid of fixed-q fits
# cannot fix sigma; instead scan q with sigma's bound and compare
best <- -Inf
for (p in seq(0.05, 0.6, by = 0.01)) {
  r <- conds_of(frm(bf(y ~ x, quantile = p, sigma ~ 1), family = al, data = d,
                    prior = set_prior("", class = "Intercept", dpar = "sigma",
                                      lb = s0 + 0.3, ub = 3)))
  if (!is.null(r$v)) best <- max(best, ll(r))
}
cat(sprintf("  grid of fixed q in [0.05, 0.6] under the sigma box: best ll %.8f\n", best))
# gaussian with a binding bound on sigma: the generic behavior
set.seed(8412)
dg <- data.frame(x = rnorm(n)); dg$y <- 1 + dg$x + rnorm(n)
rg <- conds_of(frm(bf(y ~ x, sigma ~ 1), data = dg,
                   prior = set_prior("", class = "Intercept", dpar = "sigma",
                                     lb = 0.3, ub = 2)))
cat(sprintf("gaussian sigma box binding: code %d [%s]\n", rg$v$opt$convergence, cls_line(rg)))
cat("DONE\n")

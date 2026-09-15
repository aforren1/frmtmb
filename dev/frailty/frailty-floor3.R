# The case frailty-floor.R named and this lane wrongly dismissed: a
# centre with NO EVENTS.
#
# The barrier that keeps a per-centre slope positive is `log(gamma1 + u)`
# in the log DENSITY, and a centre whose rows are all censored
# contributes no density term at all. What it contributes is `-H_i`, and
# for rows at times above 1 the score in `u` is `-sum(x_i H_i) < 0`, so
# the fit PUSHES that centre's slope negative and only the normal prior
# stops it.
#
# `rp_floored()` does not report it, because `R/rp-check.R` tests
# `cens == 0 & detadx <= 0` and that centre has no `cens == 0` row. The
# fit converges clean, warns nothing, passes `rp_floored()` and passes
# `frm_curve()`, and its fitted survival function INCREASES with time.
#
# The construction is the reviewer's, from dev/reviews/20260910-frailty.md
# section 5c, re-derived here so that this lane's findings carry a script
# it owns. `n_nodeath = 0` is the control and is frailty-floor.R's design.
#
# SEEDS 20260910 to 20260915.
source("frailty-common.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})

# frailty-floor.R's design, plus `n_nodeath` centres followed to the
# same administrative time with no deaths in any of them, which is an
# ordinary thing for a small centre in a multi-centre trial.
sim_nodeath <- function(seed, n = 400L, n_centre = 40L, n_nodeath = 5L,
                        sd_u = 0.35, gamma1 = 0.6, scale = 5,
                        beta = 0.6, p_cens = 0.4, floor_true = 0.05) {
  set.seed(seed)
  centre <- rep(seq_len(n_centre), length.out = n)
  repeat {
    u <- stats::rnorm(n_centre, 0, sd_u)
    if (all(gamma1 + u > floor_true)) break
  }
  trt <- stats::rbinom(n, 1L, 0.5)
  g0 <- -gamma1 * log(scale)
  tt <- exp((log(-log(stats::runif(n))) - g0 - beta * trt) /
              (gamma1 + u[centre]))
  tau <- unname(stats::quantile(tt, 1 - p_cens))
  ev <- as.integer(tt <= tau)
  tt <- pmin(tt, tau)
  if (n_nodeath > 0L) {
    j <- centre <= n_nodeath
    tt[j] <- tau
    ev[j] <- 0L
  }
  data.frame(time = tt, event = ev, censored = 1L - ev, trt = trt,
             centre = factor(centre))
}

one <- function(seed, n_nodeath) {
  d <- sim_nodeath(seed, n_nodeath = n_nodeath)
  bk <- range(log(d$time[d$event == 1L]))
  f <- suppressWarnings(frm(
    bf(time | cens(censored) ~ trt, gamma1 ~ (1 | centre)),
    family = royston_parmar(knots = numeric(0), bknots = bk),
    data = d, se = TRUE))
  g1 <- unname(fixef(f)$gamma1[["(Intercept)"]])
  sh <- g1 + as.numeric(frmtmb::ranef(f)[["centre"]])
  r <- rp_floored(f, action = "report")
  refuses <- inherits(try(rp_floored(f), silent = TRUE), "try-error")
  dg <- frmtmb::diagnose(f, quiet = TRUE)
  cat(sprintf(
    paste0("n_nodeath %d seed %d  gamma1 %7.4f  min slope %8.4f",
           "  n<=0 %2d  n_nonmono %d  refuses %-5s  conv %d",
           "  maxgrad %.2e  pdHess %s
"),
    n_nodeath, seed, g1, min(sh), sum(sh <= 0), r$n_nonmonotone,
    refuses, dg$convergence, dg$max_grad, isTRUE(dg$pdHess)))
  invisible(list(fit = f, d = d, g1 = g1, sh = sh))
}

# The sign argument needs `log t > 0` on those rows, so the
# administrative time is reported rather than assumed.
d0 <- sim_nodeath(20260910L, n_nodeath = 5L)
cat("administrative time ", format(max(d0$time), digits = 6),
    ", log of it ", format(log(max(d0$time)), digits = 6),
    ", fraction of the analysed times past 1 ",
    format(mean(d0$time > 1), digits = 4), "

", sep = "")

cat("== the case: five centres with no deaths ==\n")
res <- NULL
for (s in 20260910L + 0:5) {
  r <- one(s, 5L)
  if (s == 20260910L) res <- r
}
cat("\n== the control: frailty-floor.R's own design, no such centre ==\n")
for (s in 20260910L + 0:5) one(s, 0L)

# What a user would see on the worst one. A survival function that
# increases with time is not a survival function.
cat("\n== seed 20260910, centre 1, fitted survival ==\n")
f <- res[["fit"]]
sh1 <- res[["sh"]][1L]
g0 <- unname(fixef(f)$mu[["(Intercept)"]])
tt <- 10^seq(-12, 0.45, length.out = 5)
S <- exp(-exp(g0 + sh1 * log(tt)))
cat("centre 1 fitted slope", format(sh1, digits = 6), "\n")
cat("t   ", paste(formatC(tt, digits = 3, format = "e"),
                  collapse = "  "), "\n")
cat("S(t)", paste(formatC(S, digits = 3, format = "e"),
                  collapse = "  "), "\n")
cat("S increasing in t:", all(diff(S) > 0), "\n")

# and frm_curve() lets it through, which is the documented way to read
# this family
cv <- try(frm_curve(f, newdata = data.frame(
  time = exp(seq(-2, 1, length.out = 5)), trt = 0,
  centre = factor(1, levels = levels(res[["d"]]$centre))),
  dpar = "gamma1", simultaneous = FALSE), silent = TRUE)
cat("frm_curve() refuses:", inherits(cv, "try-error"), "\n")

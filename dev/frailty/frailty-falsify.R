# Every guard in test-frailty.R, run against the case it exists to
# catch. A check that cannot fail is not a check.
#
# SEED 20260910, the same fixtures the test file uses.
source("frailty-common.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})

say <- function(lab, broken_holds, real_holds) {
  cat(sprintf("%-46s broken: %-5s  as shipped: %-5s  %s\n", lab,
              broken_holds, real_holds,
              if (!broken_holds && real_holds) "OK" else "GUARD FAILS OPEN"))
}

d <- frailty_sim(20260910L, n = 800L, n_centre = 20L, sd_b = 0.5)
kn <- rp_knots(log(d$time[d$event == 1L]), 2L)
fit <- frm(bf(time | cens(censored) ~ trt + (1 | centre)),
           family = royston_parmar(knots = kn$ik, bknots = kn$bk),
           data = d, se = TRUE)
rst <- rstpm2::gsm(survival::Surv(time, event) ~ trt, data = d,
                   smooth.formula = ~ rstpm2::nsx(log(time),
                                                  knots = kn$ik,
                                                  Boundary.knots = kn$bk),
                   cluster = d$centre, RandDist = "LogN")
cf <- stats::coef(rst)
sd_f <- sqrt(frmtmb::VarCorr(fit)[[1L]][1L, 1L])
th <- exp(unname(cf[["logtheta"]]))

# 1. theta is a VARIANCE. The broken read takes it for a standard
#    deviation.
say("theta read as an sd instead of a variance",
    abs(sd_f - th) * 50 < abs(sd_f - sqrt(th)),
    abs(sd_f - sqrt(th)) * 50 < abs(sd_f - th))

# 2. the change of basis. The broken arm maps onto a spline with a
#    DIFFERENT interior knot, which spans a different space.
xg <- seq(kn$bk[1L] - 0.5, kn$bk[2L] + 0.5, length.out = 200L)
eta_r <- unname(cf[["(Intercept)"]]) +
  as.vector(rstpm2::nsx(xg, knots = kn$ik, Boundary.knots = kn$bk) %*%
              unname(cf[grep("nsx", names(cf))]))
resid_of <- function(knall) {
  max(abs(stats::lsfit(rp_basis(knall, xg), eta_r,
                       intercept = FALSE)$residuals))
}
bad <- c(kn$bk[1L], kn$ik + 0.6, kn$bk[2L])
say("the two spline spaces differ (knot moved 0.6)",
    resid_of(bad) < 1e-9 * max(abs(eta_r)),
    resid_of(kn$all) < 1e-9 * max(abs(eta_r)))
cat("    residual at the shipped knots", formatC(resid_of(kn$all),
                                                 digits = 3),
    " at the moved knot", formatC(resid_of(bad), digits = 3),
    " eta scale", formatC(max(abs(eta_r)), digits = 3), "\n")

# 3. the Laplace optimum is an optimum of the exact likelihood. The
#    broken arm asks the same question one standard error away from the
#    estimate, where the answer must be no.
beta <- unname(fixef(fit)$mu[["trt"]])
gam <- c(unname(fixef(fit)$mu[["(Intercept)"]]),
         unname(fixef(fit)$gamma1), unname(fixef(fit)$gamma2))
se_b <- summary(fit)[["coefficients"]][["mu"]]["trt", 2L]
at <- function(bb) frailty_exact_ll(d$time, d$event, bb * d$trt,
                                    d$centre, kn$all, gam, sd_f)
peak_at <- function(b0) {
  all(c(at(b0 - 0.25 * se_b), at(b0 + 0.25 * se_b)) < at(b0))
}
say("a point one se off the estimate called an optimum",
    peak_at(beta + se_b), peak_at(beta))

# 4. the random effect on gamma1 has to beat the pooled fit.
dg <- g1_sim(20260910L, n = 800L, n_centre = 20L, sd_u = 0.25)
tr <- attr(dg, "truth")
bk <- range(log(dg$time[dg$event == 1L]))
famg <- royston_parmar(knots = numeric(0), bknots = bk)
on <- suppressWarnings(frm(bf(time | cens(censored) ~ trt,
                             gamma1 ~ (1 | centre)),
                          family = famg, data = dg, se = TRUE))
off <- suppressWarnings(frm(bf(time | cens(censored) ~ trt),
                           family = famg, data = dg, se = TRUE))
g1a <- unname(fixef(on)$gamma1[["(Intercept)"]])
sh_hat <- g1a + as.numeric(frmtmb::ranef(on)[["centre"]])
err_on <- mean(abs(sh_hat - tr$shape))
err_off <- mean(abs(unname(fixef(off)$gamma1[["(Intercept)"]]) -
                      tr$shape))
# broken arm: the component is switched off but its answer is still
# read as if it were on, which is what an inert random effect would do
say("an inert gamma1 block still beats the pooled fit",
    mean(abs(g1a - tr$shape)) < err_off, err_on < err_off)
cat("    per-centre error: on", formatC(err_on, digits = 4),
    " off", formatC(err_off, digits = 4),
    " on-with-deviations-dropped",
    formatC(mean(abs(g1a - tr$shape)), digits = 4), "\n")

# 5. the anchoring. The broken arm compares the correlated block with
#    ITSELF, where the ratio cannot exceed 1.
nev <- sum(dg$event)
at_m <- function(m, form) {
  dd <- dg
  dd$time <- dd$time * m
  f <- suppressWarnings(frm(form, family = royston_parmar(
    knots = numeric(0), bknots = bk + log(m)), data = dd, se = TRUE))
  as.numeric(logLik(f)) + nev * log(m)
}
slope <- bf(time | cens(censored) ~ trt, gamma1 ~ (1 | centre))
both <- bf(time | cens(censored) ~ trt + (1 | c | centre),
           gamma1 ~ (1 | c | centre))
d_slope <- abs(at_m(1, slope) - at_m(12, slope))
d_both <- abs(at_m(1, both) - at_m(12, both))
say("a time-unit-invariant model called anchored",
    d_both > 100 * d_both, d_slope > 100 * d_both)
cat("    slope-only shift", formatC(d_slope, digits = 4),
    " paired-block shift", formatC(d_both, digits = 4),
    " ratio", formatC(d_slope / d_both, digits = 4), "\n")

# The identity, and where it stops.
#
# Three log likelihoods of the SAME model at the SAME parameter point:
#   exact   adaptive quadrature per cluster, written here, rel.tol 1e-12
#   frmtmb  the Laplace approximation RTMB builds
#   rstpm2  adaptive Gauss-Hermite with `nodes` nodes
#
# SEED 20260910. Design: 2000 rows, 40 centres, sd(frailty) 0.5,
# Weibull shape 1.3, beta 0.6, 40 percent censored, RP df = 2.
source("frailty-common.R")

# THE INSTRUMENT THIS WHOLE SCRIPT RESTS ON. `par_of()` builds the
# objective's parameter vector positionally and `setNames()` only
# relabels, so the ordering is an ASSUMPTION. The objective's own names
# are `beta, beta, betad, betad, theta`, duplicated, so they give no
# protection. Check it against a number the package computed by another
# route: the Laplace objective at the fitted estimates must be the fit's
# own logLik().
check_par_order <- function(lap_at_fit, ll) {
  rel <- abs(lap_at_fit - ll) / abs(ll)
  cat(sprintf("par-order check: -obj$fn(estimates) %.10f  logLik %.10f  rel %.2e
",
              lap_at_fit, ll, rel))
  if (rel > 1e-8) {
    stop("the objective's parameter ordering is not what par_of() ",
         "assumes; every offset below would be scored at the wrong ",
         "point", call. = FALSE)
  }
  invisible(TRUE)
}

suppressMessages({
  library(frmtmb); library(frmtmb.spline); library(rstpm2)
  library(survival); library(bbmle)
})

sd_of <- function(fit) sqrt(frmtmb::VarCorr(fit)[[1L]][1L, 1L])

run <- function(seed = 20260910L, df = 2L, sd_b = 0.5, n = 2000L,
                n_centre = 40L) {
  d <- frailty_sim(seed, n = n, n_centre = n_centre, sd_b = sd_b)
  tr <- attr(d, "truth")
  kn <- rp_knots(log(d$time[d$event == 1L]), df)
  fam <- royston_parmar(knots = kn$ik, bknots = kn$bk)
  form <- bf(time | cens(censored) ~ trt + (1 | centre))

  ff <- frm(form, family = fam, data = d, se = TRUE)
  dry <- frm(form, family = fam, data = d, dry_run = "objective")
  rr <- stpm2(Surv(time, event) ~ trt, data = d,
              smooth.formula = ~ nsx(log(time), knots = kn$ik,
                                     Boundary.knots = kn$bk),
              cluster = d$centre, RandDist = "LogN")

  xg <- seq(kn$bk[1L] - 0.5, kn$bk[2L] + 0.5, length.out = 400L)
  mp <- nsx_to_rp(rr, kn$all, xg)

  gam_f <- c(unname(fixef(ff)$mu[["(Intercept)"]]),
             unname(fixef(ff)$gamma1), unname(fixef(ff)$gamma2))
  gam_f <- gam_f[seq_len(df + 1L)]
  beta_f <- unname(fixef(ff)$mu[["trt"]])
  sd_f <- sd_of(ff)

  gam_r <- mp$gam
  beta_r <- unname(coef(rr)[["trt"]])
  sd_r <- sqrt(exp(unname(coef(rr)[["logtheta"]])))

  par_of <- function(gam, beta, sd) {
    p <- c(gam[1L], beta, gam[-1L], log(sd))
    stats::setNames(p, names(dry$obj$par))
  }
  lap <- function(gam, beta, sd) -dry$obj$fn(par_of(gam, beta, sd))
  ex <- function(gam, beta, sd) {
    frailty_exact_ll(d$time, d$event, beta * d$trt, d$centre, kn$all,
                     gam, sd)
  }

  check_par_order(lap(gam_f, beta_f, sd_f), as.numeric(logLik(ff)))
  cat("== design: seed", seed, " df", df, " sd_b", sd_b, " n", n,
      " centres", n_centre, "\n")
  cat("nsx -> RP map residual", format(mp$resid, digits = 3),
      " relative to eta scale", format(mp$resid / mp$scale, digits = 3),
      "\n\n")

  cat("--- estimates ---\n")
  cat(sprintf("%-14s %14s %14s %14s\n", "", "frmtmb", "rstpm2", "truth"))
  cat(sprintf("%-14s %14.8f %14.8f %14.8f\n", "beta(trt)", beta_f,
              beta_r, tr$beta))
  cat(sprintf("%-14s %14.8f %14.8f %14.8f\n", "sd(frailty)", sd_f,
              sd_r, tr$sd_b))
  for (j in seq_along(gam_f)) {
    cat(sprintf("%-14s %14.8f %14.8f %14s\n", paste0("gamma", j - 1L),
                gam_f[j], gam_r[j],
                if (j == 1L) format(tr$gamma0, digits = 8) else
                  if (j == 2L) format(tr$gamma1, digits = 8) else "0"))
  }

  cat("\n--- log likelihood, each package at its OWN optimum ---\n")
  ll_f <- as.numeric(logLik(ff))
  ll_r <- as.numeric(logLik(rr))
  cat(sprintf("frmtmb reported %.6f   rstpm2 reported %.6f   diff %.3e\n",
              ll_f, ll_r, ll_f - ll_r))

  cat("\n--- the same three rules at the SAME parameter point ---\n")
  pts <- list(`frmtmb optimum` = list(gam_f, beta_f, sd_f),
              `rstpm2 optimum` = list(gam_r, beta_r, sd_r))
  for (nm in names(pts)) {
    p <- pts[[nm]]
    e <- ex(p[[1L]], p[[2L]], p[[3L]])
    l <- lap(p[[1L]], p[[2L]], p[[3L]])
    cat(sprintf("%-16s exact %14.8f  laplace %14.8f  (lap-exact) %+.4e\n",
                nm, e, l, l - e))
  }

  cat("\n--- rstpm2's own rule against the exact integral ---\n")
  for (nd in c(3L, 9L, 21L, 51L)) {
    rn <- stpm2(Surv(time, event) ~ trt, data = d,
                smooth.formula = ~ nsx(log(time), knots = kn$ik,
                                       Boundary.knots = kn$bk),
                cluster = d$centre, RandDist = "LogN",
                control = list(nodes = nd))
    mpn <- nsx_to_rp(rn, kn$all, xg)
    bn <- unname(coef(rn)[["trt"]])
    sn <- sqrt(exp(unname(coef(rn)[["logtheta"]])))
    en <- ex(mpn$gam, bn, sn)
    cat(sprintf("nodes %3d  reported %14.8f  exact-at-its-own-point %14.8f",
                nd, as.numeric(logLik(rn)), en))
    cat(sprintf("  (rep-exact) %+.4e  beta %.8f  sd %.8f\n",
                as.numeric(logLik(rn)) - en, bn, sn))
  }

  cat("\n--- who is closer to the exact optimum ---\n")
  ef <- ex(gam_f, beta_f, sd_f)
  er <- ex(gam_r, beta_r, sd_r)
  cat(sprintf("exact ll at frmtmb's point %.8f, at rstpm2's %.8f, ",
              ef, er))
  cat(sprintf("frmtmb - rstpm2 = %+.4e\n", ef - er))
  invisible(NULL)
}

run()

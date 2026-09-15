# Where does the Laplace approximation stop agreeing with the exact
# integral, and does the ESTIMATE follow it?
#
# Cluster size is the axis that matters: the Laplace approximation is
# exact for a Gaussian conditional and a cluster's conditional log
# likelihood approaches one as its event count grows.
#
# Seeds 20260910 + 0..4 per cell. Writes frailty-sweep.tsv.
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

OUT <- "frailty-sweep.tsv"
if (file.exists(OUT)) file.remove(OUT)
rec <- function(...) {
  v <- list(...)
  line <- paste(paste0(names(v), "=", vapply(v, function(z) {
    if (is.numeric(z)) formatC(z, digits = 10, format = "g") else
      as.character(z)
  }, character(1))), collapse = "\t")
  cat(line, "\n", sep = "", file = OUT, append = TRUE)
  message(line)
}

one <- function(cell, seed, n_centre, sd_b, df = 2L, n = 2000L) {
  d <- frailty_sim(seed, n = n, n_centre = n_centre, sd_b = sd_b)
  tr <- attr(d, "truth")
  kn <- rp_knots(log(d$time[d$event == 1L]), df)
  fam <- royston_parmar(knots = kn$ik, bknots = kn$bk)
  form <- bf(time | cens(censored) ~ trt + (1 | centre))
  t0 <- Sys.time()
  ff <- frm(form, family = fam, data = d, se = TRUE)
  t_f <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  rr <- try(stpm2(Surv(time, event) ~ trt, data = d,
                  smooth.formula = ~ nsx(log(time), knots = kn$ik,
                                         Boundary.knots = kn$bk),
                  cluster = d$centre, RandDist = "LogN"), silent = TRUE)
  t_r <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(rr, "try-error")) {
    rec(cell = cell, seed = seed, n_centre = n_centre, sd_b = sd_b,
        rstpm2 = "FAILED")
    return(invisible(NULL))
  }
  xg <- seq(kn$bk[1L] - 0.5, kn$bk[2L] + 0.5, length.out = 400L)
  mp <- nsx_to_rp(rr, kn$all, xg)

  gam_f <- c(unname(fixef(ff)$mu[["(Intercept)"]]),
             unname(fixef(ff)$gamma1), unname(fixef(ff)$gamma2))[
               seq_len(df + 1L)]
  beta_f <- unname(fixef(ff)$mu[["trt"]])
  sd_f <- sqrt(frmtmb::VarCorr(ff)[[1L]][1L, 1L])
  cf <- summary(ff)[["coefficients"]][["mu"]]
  se_beta_f <- unname(cf["trt", 2L])
  cv <- frmtmb::confint_varcorr(ff)
  sd_lo <- cv[["lwr"]][1L]; sd_hi <- cv[["upr"]][1L]
  ci <- suppressWarnings(stats::confint(ff))
  rn <- rownames(ci)
  ib <- grep("mu.*trt|^trt$", rn)[1L]
  b_lo <- ci[ib, 1L]; b_hi <- ci[ib, 2L]

  beta_r <- unname(coef(rr)[["trt"]])
  sd_r <- sqrt(exp(unname(coef(rr)[["logtheta"]])))
  vr <- try(vcov(rr), silent = TRUE)
  se_beta_r <- if (inherits(vr, "try-error")) NA_real_ else
    sqrt(vr["trt", "trt"])
  se_lt_r <- if (inherits(vr, "try-error")) NA_real_ else
    sqrt(vr["logtheta", "logtheta"])

  dry <- frm(form, family = fam, data = d, dry_run = "objective")
  par_of <- function(gam, beta, sd) {
    stats::setNames(c(gam[1L], beta, gam[-1L], log(sd)),
                    names(dry$obj$par))
  }
  ex <- function(gam, beta, sd) {
    frailty_exact_ll(d$time, d$event, beta * d$trt, d$centre, kn$all,
                     gam, sd)
  }
  ex_f <- ex(gam_f, beta_f, sd_f)
  ex_r <- ex(mp$gam, beta_r, sd_r)
  lap_f <- -dry$obj$fn(par_of(gam_f, beta_f, sd_f))
  lap_r <- -dry$obj$fn(par_of(mp$gam, beta_r, sd_r))
  check_par_order(lap_f, as.numeric(logLik(ff)))

  rec(cell = cell, seed = seed, n_centre = n_centre,
      per_centre = n / n_centre,
      ev_per_centre = sum(d$event) / n_centre, sd_b = sd_b, df = df,
      map_resid_rel = mp$resid / mp$scale,
      beta_frm = beta_f, beta_rst = beta_r,
      beta_reldiff = abs(beta_f - beta_r) / abs(beta_r),
      beta_diff_in_se = abs(beta_f - beta_r) / se_beta_f,
      se_beta_frm = se_beta_f, se_beta_rst = se_beta_r,
      sd_frm = sd_f, sd_rst = sd_r,
      sd_reldiff = abs(sd_f - sd_r) / sd_r,
      se_logsd_frm = (log(sd_hi) - log(sd_lo)) / (2 * 1.959964),
      se_logsd_rst = se_lt_r / 2,
      b_cover = as.integer(b_lo <= tr$beta && b_hi >= tr$beta),
      sd_cover = as.integer(sd_lo <= sd_b && sd_hi >= sd_b),
      ll_frm = as.numeric(logLik(ff)), ll_rst = as.numeric(logLik(rr)),
      exact_at_frm = ex_f, exact_at_rst = ex_r,
      lap_err_at_frm = lap_f - ex_f, lap_err_at_rst = lap_r - ex_r,
      ghq_err_at_rst = as.numeric(logLik(rr)) - ex_r,
      exact_gap = ex_f - ex_r, t_frm = t_f, t_rst = t_r)
  invisible(NULL)
}

seeds <- 20260910L + 0:4
for (nc in c(10L, 40L, 200L, 500L)) {
  for (s in seeds) one(paste0("nc", nc), s, nc, 0.5)
}
for (sb in c(0.25, 1.0)) {
  for (s in seeds) one(paste0("sd", sb), s, 40L, sb)
}
cat("done\n")

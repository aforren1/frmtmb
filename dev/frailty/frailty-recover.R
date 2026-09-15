# The recovery table for the plan's realistic RP design:
# 2000 subjects, 40 percent censored, one frailty (40 centres of 50),
# df = 3, sd(frailty) 0.5, beta 0.6. 200 replicates, seeds
# 20260910 + 0..199. Writes frailty-recover.tsv.
#
# The data-generating model is a Weibull PH model with a shared
# log-normal frailty, which is exactly a Royston-Parmar model whose
# spline is linear in log time, so at df = 3 the truth is
# gamma0 = -shape log(scale), gamma1 = shape, gamma2 = gamma3 = 0
# whatever the knots are: a globally linear function has one
# representation in the RP basis and its interior coefficients are 0.
source("frailty-common.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline); library(rstpm2)
  library(survival); library(bbmle)
})

OUT <- "frailty-recover.tsv"
if (file.exists(OUT)) file.remove(OUT)
rec <- function(...) {
  v <- list(...)
  line <- paste(paste0(names(v), "=", vapply(v, function(z) {
    if (is.numeric(z)) formatC(z, digits = 10, format = "g") else
      as.character(z)
  }, character(1))), collapse = "\t")
  cat(line, "\n", sep = "", file = OUT, append = TRUE)
}

DF <- 3L
N <- 2000L
NC <- 40L
SD_B <- 0.5
NREP <- 200L

one <- function(seed) {
  d <- frailty_sim(seed, n = N, n_centre = NC, sd_b = SD_B)
  tr <- attr(d, "truth")
  kn <- rp_knots(log(d$time[d$event == 1L]), DF)
  fam <- royston_parmar(knots = kn$ik, bknots = kn$bk)
  form <- bf(time | cens(censored) ~ trt + (1 | centre))
  ff <- try(suppressWarnings(frm(form, family = fam, data = d,
                                 se = TRUE)), silent = TRUE)
  if (inherits(ff, "try-error")) {
    rec(seed = seed, status = "frm_failed")
    return(invisible(NULL))
  }
  rr <- try(stpm2(Surv(time, event) ~ trt, data = d,
                  smooth.formula = ~ nsx(log(time), knots = kn$ik,
                                         Boundary.knots = kn$bk),
                  cluster = d$centre, RandDist = "LogN"), silent = TRUE)
  ok_r <- !inherits(rr, "try-error")

  gam_f <- c(unname(fixef(ff)$mu[["(Intercept)"]]),
             unname(fixef(ff)$gamma1), unname(fixef(ff)$gamma2),
             unname(fixef(ff)$gamma3))
  beta_f <- unname(fixef(ff)$mu[["trt"]])
  sd_f <- sqrt(frmtmb::VarCorr(ff)[[1L]][1L, 1L])
  ci <- suppressWarnings(stats::confint(ff))
  rn <- rownames(ci)
  cover <- function(pat, truth) {
    i <- grep(pat, rn)[1L]
    if (is.na(i)) return(NA_integer_)
    as.integer(ci[i, 1L] <= truth && ci[i, 2L] >= truth)
  }
  cv <- frmtmb::confint_varcorr(ff)
  sd_cov <- as.integer(cv[["lwr"]][1L] <= SD_B && cv[["upr"]][1L] >= SD_B)
  dg <- frmtmb::diagnose(ff, quiet = TRUE)
  fl <- rp_floored(ff, action = "report")

  ex_f <- frailty_exact_ll(d$time, d$event, beta_f * d$trt, d$centre,
                           kn$all, gam_f, sd_f)
  if (ok_r) {
    xg <- seq(kn$bk[1L] - 0.5, kn$bk[2L] + 0.5, length.out = 400L)
    mp <- nsx_to_rp(rr, kn$all, xg)
    beta_r <- unname(coef(rr)[["trt"]])
    sd_r <- sqrt(exp(unname(coef(rr)[["logtheta"]])))
    ex_r <- frailty_exact_ll(d$time, d$event, beta_r * d$trt, d$centre,
                             kn$all, mp$gam, sd_r)
    ll_r <- as.numeric(logLik(rr))
    vr <- try(vcov(rr), silent = TRUE)
    se_b_r <- if (inherits(vr, "try-error")) NA_real_ else
      sqrt(vr["trt", "trt"])
    map_rel <- mp$resid / mp$scale
  } else {
    beta_r <- sd_r <- ex_r <- ll_r <- se_b_r <- map_rel <- NA_real_
  }

  rec(seed = seed, status = "ok",
      beta_frm = beta_f, beta_rst = beta_r,
      se_beta_frm = summary(ff)[["coefficients"]][["mu"]]["trt", 2L],
      se_beta_rst = se_b_r,
      sd_frm = sd_f, sd_rst = sd_r,
      g0 = gam_f[1L], g1 = gam_f[2L], g2 = gam_f[3L], g3 = gam_f[4L],
      cov_beta = cover("trt", tr$beta),
      cov_g0 = cover("mu.*Intercept|^\\(Intercept\\)$", tr$gamma0),
      cov_g1 = cover("gamma1", tr$gamma1),
      cov_g2 = cover("gamma2", 0),
      cov_g3 = cover("gamma3", 0),
      cov_sd = sd_cov,
      ll_frm = as.numeric(logLik(ff)), ll_rst = ll_r,
      exact_at_frm = ex_f, exact_at_rst = ex_r,
      conv = dg$convergence, maxgrad = dg$max_grad,
      pdHess = as.integer(isTRUE(dg$pdHess)),
      n_nonmono = fl$n_nonmonotone, max_nlogS = fl$max_nlogS,
      map_rel = map_rel)
  invisible(NULL)
}

d0 <- frailty_sim(20260910L, n = N, n_centre = NC, sd_b = SD_B)
kn0 <- rp_knots(log(d0$time[d0$event == 1L]), DF)
f0 <- frm(bf(time | cens(censored) ~ trt + (1 | centre)),
          family = royston_parmar(knots = kn0$ik, bknots = kn0$bk),
          data = d0, se = TRUE)
cat("confint rownames:", paste(rownames(suppressWarnings(confint(f0))),
                               collapse = " | "), "\n")
# the pinned knots ARE the family's own default rule
fd <- frm(bf(time | cens(censored) ~ trt + (1 | centre)),
          family = royston_parmar(df = DF), data = d0, se = TRUE)
knd <- environment(stats::family(fd)[["lpdf"]])$allknots
cat("pinned knots equal default:", isTRUE(all.equal(unname(knd),
                                                    unname(kn0$all))),
    " max abs diff", max(abs(knd - kn0$all)), "\n")

t0 <- Sys.time()
for (i in seq_len(NREP)) {
  one(20260910L + i - 1L)
  if (i %% 20L == 0L) {
    message(i, " done, ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
            " min")
  }
}
cat("done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")),
                     2), "min\n")

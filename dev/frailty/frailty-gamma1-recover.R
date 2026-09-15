# Recovery for a random effect on gamma1: 60 replicates, seeds
# 20260910 + 0..59, n = 2000, 40 centres of 50, sd(gamma1 | centre) 0.2
# against gamma1 = 1.3, beta 0.6, 40 percent censored, df = 1.
#
# df = 1 makes the fitted model EXACTLY the model the data came from, so
# every truth is known: the spline is gamma0 + (gamma1 + u_c) log t and
# the per-centre Weibull shape is gamma1 + u_c.
#
# The second arm has the component switched OFF, because a random
# effect earns its place by beating what the pooled fit already gets,
# not by being estimable. Item 1.0a of the plan learned that the hard
# way on a bounded non-decision time.
#
# Writes frailty-gamma1.tsv.
source("frailty-common.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})

OUT <- "frailty-gamma1.tsv"
if (file.exists(OUT)) file.remove(OUT)
rec <- function(...) {
  v <- list(...)
  line <- paste(paste0(names(v), "=", vapply(v, function(z) {
    if (is.numeric(z)) formatC(z, digits = 10, format = "g") else
      as.character(z)
  }, character(1))), collapse = "\t")
  cat(line, "\n", sep = "", file = OUT, append = TRUE)
}

SD_U <- 0.2
G1 <- 1.3
NREP <- 60L

one <- function(seed) {
  d <- g1_sim(seed, sd_u = SD_U, gamma1 = G1)
  tr <- attr(d, "truth")
  bk <- range(log(d$time[d$event == 1L]))
  fam <- royston_parmar(knots = numeric(0), bknots = bk)
  on_form <- bf(time | cens(censored) ~ trt, gamma1 ~ (1 | centre))
  off_form <- bf(time | cens(censored) ~ trt)
  fa <- try(suppressWarnings(frm(on_form, family = fam, data = d,
                                 se = TRUE)), silent = TRUE)
  fb <- try(suppressWarnings(frm(off_form, family = fam, data = d,
                                 se = TRUE)), silent = TRUE)
  if (inherits(fa, "try-error") || inherits(fb, "try-error")) {
    rec(seed = seed, status = "failed")
    return(invisible(NULL))
  }
  g1a <- unname(fixef(fa)$gamma1[["(Intercept)"]])
  ua <- as.numeric(frmtmb::ranef(fa)[["centre"]])
  sh_hat <- g1a + ua
  sh_true <- tr$shape
  g1b <- unname(fixef(fb)$gamma1[["(Intercept)"]])
  sd_hat <- sqrt(frmtmb::VarCorr(fa)[[1L]][1L, 1L])
  cv <- frmtmb::confint_varcorr(fa)
  ci <- suppressWarnings(stats::confint(fa))
  rn <- rownames(ci)
  cover <- function(pat, truth) {
    i <- grep(pat, rn)[1L]
    if (is.na(i)) return(NA_integer_)
    as.integer(ci[i, 1L] <= truth && ci[i, 2L] >= truth)
  }
  dg <- frmtmb::diagnose(fa, quiet = TRUE)
  fl <- rp_floored(fa, action = "report")
  ex <- frailty_exact_ll_slope(
    d$time, d$event, unname(fixef(fa)$mu[["trt"]]) * d$trt, d$centre,
    c(bk[1L], bk[2L]),
    c(unname(fixef(fa)$mu[["(Intercept)"]]), g1a), sd_hat)
  rec(seed = seed, status = "ok",
      sd_hat = sd_hat, sd_lo = cv[["lwr"]][1L], sd_hi = cv[["upr"]][1L],
      cov_sd = as.integer(cv[["lwr"]][1L] <= SD_U &&
                            cv[["upr"]][1L] >= SD_U),
      g1_on = g1a, g1_off = g1b,
      beta_on = unname(fixef(fa)$mu[["trt"]]),
      beta_off = unname(fixef(fb)$mu[["trt"]]),
      g0_on = unname(fixef(fa)$mu[["(Intercept)"]]),
      cov_beta = cover("trt", tr$beta),
      cov_g1 = cover("gamma1", G1),
      cov_g0 = cover("mu.*Intercept|^\\(Intercept\\)$", tr$gamma0),
      err_on = mean(abs(sh_hat - sh_true)),
      err_off = mean(abs(g1b - sh_true)),
      cor_on = stats::cor(sh_hat, sh_true),
      min_slope = min(sh_hat), margin_sd = min(sh_hat) / sd_hat,
      ll_on = as.numeric(logLik(fa)), ll_off = as.numeric(logLik(fb)),
      exact_on = ex, lap_err = as.numeric(logLik(fa)) - ex,
      n_nonmono = fl$n_nonmonotone,
      conv = dg$convergence, maxgrad = dg$max_grad,
      pdHess = as.integer(isTRUE(dg$pdHess)))
  invisible(NULL)
}

t0 <- Sys.time()
for (i in seq_len(NREP)) {
  one(20260910L + i - 1L)
  if (i %% 10L == 0L) {
    message(i, " done, ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
            " min")
  }
}
cat("done in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    "min\n")

## The item 2.6 simulator, one copy, sourced by every coh-* script.
##
## It draws the SAME data as the Phase 0 tier's
## tests/testthat/test-scale.R at the same seed and the same sizes:
## the random-number order is set.seed, the id deviations, the
## id-by-condition deviations, then one complex Wishart draw per row.
## dev/coh-probe.R checks that byte for byte against the tier.

coupling_truth <- list(n_sub = 40L, n_freq = 60L, n_seg = 8L,
                       b0 = -0.6, b_cond = 0.5, amp = 0.8,
                       sd_id = 0.35, sd_idcond = 0.2,
                       log_s11 = 0.3, log_s22 = 0.1, phase = 0.4)

coupling_fbump <- function(fr) {
  coupling_truth$amp * exp(-0.5 * ((fr - 0.35) / 0.12)^2)
}

coupling_draw <- function(s11, s22, coh, phase, nseg) {
  m <- length(s11)
  a <- sqrt(s11)
  cm <- sqrt(coh * s22)
  b <- sqrt((1 - coh) * s22)
  w11 <- numeric(m); w22 <- numeric(m)
  w12r <- numeric(m); w12i <- numeric(m)
  for (i in seq_len(m)) {
    z1 <- complex(real = stats::rnorm(nseg, 0, sqrt(0.5)),
                  imaginary = stats::rnorm(nseg, 0, sqrt(0.5)))
    z2 <- complex(real = stats::rnorm(nseg, 0, sqrt(0.5)),
                  imaginary = stats::rnorm(nseg, 0, sqrt(0.5)))
    d1 <- a[i] * z1
    d2 <- cm[i] * complex(modulus = 1, argument = -phase[i]) * z1 +
      b[i] * z2
    cr <- sum(d1 * Conj(d2))
    w11[i] <- sum(Mod(d1)^2)
    w22[i] <- sum(Mod(d2)^2)
    w12r[i] <- Re(cr)
    w12i[i] <- Im(cr)
  }
  data.frame(w11 = w11, w22 = w22, w12r = w12r, w12i = w12i, n = nseg)
}

## sd_idcond is an argument because the NULL arm of this study needs
## the same generator with that one component switched off.
coh_data <- function(seed, n_sub = coupling_truth$n_sub,
                     n_freq = coupling_truth$n_freq,
                     sd_idcond = coupling_truth$sd_idcond) {
  tr <- coupling_truth
  set.seed(seed)
  d <- expand.grid(freq = seq_len(n_freq) / n_freq,
                   cond = factor(c("a", "b")),
                   id = factor(seq_len(n_sub)))
  u_id <- stats::rnorm(n_sub, 0, tr$sd_id)
  u_ic <- stats::rnorm(n_sub * 2L, 0, sd_idcond)
  ic <- as.integer(d$id) + n_sub * (as.integer(d$cond) - 1L)
  eta <- tr$b0 + tr$b_cond * (d$cond == "b") + coupling_fbump(d$freq) +
    u_id[as.integer(d$id)] + u_ic[ic]
  w <- coupling_draw(rep(exp(tr$log_s11), nrow(d)),
                     rep(exp(tr$log_s22), nrow(d)),
                     stats::plogis(eta), rep(tr$phase, nrow(d)),
                     tr$n_seg)
  cbind(d, w)
}

coh_bf <- function(rhs) {
  frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
             pow2 ~ 1,
             stats::as.formula(paste0("coh ~ ", rhs)),
             phase ~ 1)
}

## The ladder, by the tier's own names. "full" is the correct model:
## it is the only rung carrying every term the simulator has.
## `idcond` is not a rung of the Phase 0 ladder. It is here to separate
## two explanations of the width gap: "more random effects widen an
## interval" and "only the term the contrast varies within widens it".
coh_rungs <- function(k = 10L) {
  sm <- paste0("cond + s(freq, by = cond, k = ", k, ")")
  list(cond = "cond",
       smooth = sm,
       id = paste(sm, "+ (1 | id)"),
       idcond = paste(sm, "+ (1 | id:cond)"),
       full = paste(sm, "+ (1 | id) + (1 | id:cond)"))
}

## Estimate, standard error and Wald interval for the condition
## contrast, plus what the fit says about itself. The standard error is
## read from the fit rather than reconstructed from the interval, which
## is what item 2.6 asks to capture.
coh_fit_one <- function(rhs, d, level = 0.95) {
  t0 <- Sys.time()
  fit <- suppressWarnings(
    frmtmb::frm(coh_bf(rhs), family = frmtmb.coupling::cross_wishart(),
                data = d, se = TRUE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ci <- suppressWarnings(stats::confint(fit, level = level))
  se_all <- sqrt(diag(stats::vcov(fit)))
  j <- match("coh_condb", rownames(ci))
  k <- match("coh_condb", names(se_all))
  if (is.na(j) || is.na(k)) {
    stop("no coh_condb row: ", paste(rownames(ci), collapse = ","))
  }
  vc <- frmtmb::VarCorr(fit)
  sd_of <- function(pat) {
    i <- which(names(vc) == pat)
    if (length(i)) sqrt(vc[[i[1L]]][1L, 1L]) else NA_real_
  }
  dg <- frmtmb::diagnose(fit, quiet = TRUE)
  list(est = unname(ci[j, "est"]), se = unname(se_all[k]),
       lo = unname(ci[j, "lwr"]), hi = unname(ci[j, "upr"]),
       loglik = as.numeric(stats::logLik(fit)),
       sd_id = sd_of("coh: 1 | id"),
       sd_idcond = sd_of("coh: 1 | id:cond"),
       conv = dg$convergence, maxgrad = dg$max_grad,
       pdhess = isTRUE(dg$pdHess), nbadse = length(dg$bad_se),
       secs = secs)
}

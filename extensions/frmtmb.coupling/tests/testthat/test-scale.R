## Phase 0 of dev/extension-gaps-plan.md: the coupling row.
##
## 40 subjects x 2 conditions x 60 frequencies, a smooth in frequency
## plus random effects. 4,800 rows.
##
## The plan asks for "the four scratchpad benchmarks from the survey"
## promoted to the tier, and quotes 2.4 s, 14.6 s, 52.5 s and 28.9 s.
## Their constructions ARE recorded, in the survey's own session
## transcript rather than in the repository, and the two that matter
## are named there by formula on 4,800 rows:
##
##   52.5 s  coh ~ cond + s(freq, by = cond) + (1 | id)
##   28.9 s  coh ~ cond + s(freq, by = cond) + (1 | id) + (1 | id:cond)
##
## Both are rungs of the ladder below, so the promotion the plan asked
## for is done. The 52.5 s model is the one that produced the survey's
## condition contrast of 0.77 against a truth of 0.5, which is why it
## is a rung of its own here rather than an omission: it is the only
## specification in the set that carries a subject effect WITHOUT the
## subject-by-condition effect the simulator has.
##
## The survey's numbers came from frmtmb 0.53.0 under `load_all()` and
## these come from 0.55.0 installed, so the pairs are not a
## before-and-after on one change and no speedup is claimed from them.
##
## The truth puts all of the structure in the coherence. The two power
## parameters and the phase are constants, so an intercept-only formula
## for them is the correct model and the ladder varies only the
## coherence, which is the quantity the package exists for.
##
## See dev/scale-findings.md for the numbers this produced.

coupling_truth <- list(n_sub = 40L, n_freq = 60L, n_seg = 8L,
                       b0 = -0.6, b_cond = 0.5, amp = 0.8,
                       sd_id = 0.35, sd_idcond = 0.2,
                       log_s11 = 0.3, log_s22 = 0.1, phase = 0.4)

# The smooth in frequency: one bump, so that a spline has something to
# find and a constant is visibly the wrong model.
coupling_fbump <- function(fr) {
  coupling_truth$amp * exp(-0.5 * ((fr - 0.35) / 0.12)^2)
}

# One complex Wishart draw per row, from the spectral matrix the truth
# names. Written here rather than through frm_cross_simulate(), which
# needs a FIT to simulate from and so could not draw from a truth.
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

coupling_scale_data <- function(seed = 20260908L) {
  tr <- coupling_truth
  set.seed(seed)
  n_sub <- if (scale_small()) 6L else tr$n_sub
  n_freq <- if (scale_small()) 15L else tr$n_freq
  d <- expand.grid(freq = seq_len(n_freq) / n_freq,
                   cond = factor(c("a", "b")),
                   id = factor(seq_len(n_sub)))
  u_id <- stats::rnorm(n_sub, 0, tr$sd_id)
  u_ic <- stats::rnorm(n_sub * 2L, 0, tr$sd_idcond)
  ic <- as.integer(d$id) + n_sub * (as.integer(d$cond) - 1L)
  eta <- tr$b0 + tr$b_cond * (d$cond == "b") + coupling_fbump(d$freq) +
    u_id[as.integer(d$id)] + u_ic[ic]
  w <- coupling_draw(rep(exp(tr$log_s11), nrow(d)),
                     rep(exp(tr$log_s22), nrow(d)),
                     stats::plogis(eta), rep(tr$phase, nrow(d)),
                     tr$n_seg)
  cbind(d, w)
}

# The four models of the ladder. Only the coherence carries structure.
coupling_bf <- function(rhs) {
  frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
             pow2 ~ 1,
             stats::as.formula(paste0("coh ~ ", rhs)),
             phase ~ 1)
}

coupling_models <- function() {
  list(c(name = "coh-intercept", rhs = "1"),
       c(name = "coh-cond", rhs = "cond"),
       c(name = "coh-smooth", rhs = "cond + s(freq, by = cond, k = 10)"),
       # the survey's 52.5 s model, and the one its 0.77 came from
       c(name = "coh-id",
         rhs = paste("cond + s(freq, by = cond, k = 10) +",
                     "(1 | id)")),
       # the survey's 28.9 s model, and the one item 2.6 names
       c(name = "coh-full",
         rhs = paste("cond + s(freq, by = cond, k = 10) +",
                     "(1 | id) + (1 | id:cond)")))
}

coupling_one <- function(row, rhs, d) {
  form <- coupling_bf(rhs)
  fam <- cross_wishart()
  scale_mem_reset()
  bd <- scale_build(form, family = fam, data = d)
  g0 <- scale_grad(bd$dry$obj, bd$dry$obj$par)
  ctl <- scale_control(bd$dry$obj, bd$dry$obj$par, g0$calls)
  bd$dry <- NULL

  fit <- NULL
  t_fit <- scale_elapsed(
    fit <- suppressWarnings(frmtmb::frm(form, family = fam, data = d,
                                        se = TRUE)))
  g1 <- scale_grad(fit$obj, fit$opt$par)
  mem <- scale_mem_peak_mb()

  b <- unlist(frmtmb::fixef(fit))
  ci <- suppressWarnings(stats::confint(fit))
  j <- grep("coh_condb", rownames(ci), fixed = TRUE)
  i_c <- if (length(j)) as.numeric(ci[j[1L], 1:2]) else c(NA, NA)
  # Every variance component, BY NAME. A rung whose id components
  # collapsed to zero has silently become the rung below it and its
  # contrast would then say nothing about this model, so the components
  # have to be readable. They are recorded by name rather than by
  # position because an mgcv smooth is a random-effect block too, so
  # position 1 is not the same block from one rung to the next.
  vc <- frmtmb::VarCorr(fit)
  sd_all <- if (length(vc)) {
    paste(names(vc), formatC(vapply(vc, function(m) sqrt(m[1L, 1L]),
                                    numeric(1)),
                             digits = 4, format = "g"),
          sep = "=", collapse = ";")
  } else "none"
  tr <- coupling_truth
  scale_record(
    row, rows = nrow(d), n_par = length(fit$opt$par),
    sd_all = sd_all,
    sd_id_true = tr$sd_id, sd_idcond_true = tr$sd_idcond,
    build_s = bd$build_s, frame_s = bd$frame_s,
    grad_start_s = g0$seconds, grad_start_calls = g0$calls,
    grad_opt_s = g1$seconds, control_ratio = ctl,
    fit_s = t_fit, mem_mb = mem,
    logLik = as.numeric(stats::logLik(fit)),
    coh_cond = unname(b["coh.condb"]) %||% NA_real_,
    coh_cond_true = tr$b_cond,
    coh_cond_lo = i_c[1L], coh_cond_hi = i_c[2L],
    diag = scale_diag(fit))
  list(fit = fit, interval = i_c)
}

for (m in coupling_models()) {
  local({
    mm <- m
    test_that(paste0("the coupling scale row '", mm[["name"]],
                     "' fits and reports its cost"), {
      skip_unless_scale()
      scale_row_on(paste0("coupling-", mm[["name"]]))
      d <- coupling_scale_data()
      r <- coupling_one(paste0("coupling-", mm[["name"]]), mm[["rhs"]], d)
      expect_true(is.finite(as.numeric(stats::logLik(r$fit))))
      if (identical(mm[["name"]], "coh-full") && !scale_small()) {
        # Only the top of the ladder is asserted to recover, and only
        # at the realistic design. The rungs below it OMIT random
        # effects the truth has, so their interval is the wrong width
        # rather than their point estimate being far off. And at the
        # small size this assertion has no power: 6 subjects x 15
        # frequencies gives 0.212 (-0.018, 0.441), which excludes 0.5
        # for want of data rather than because anything is wrong.
        expect_true(r$interval[1L] <= coupling_truth$b_cond &&
                      r$interval[2L] >= coupling_truth$b_cond)
      }
    })
  })
}

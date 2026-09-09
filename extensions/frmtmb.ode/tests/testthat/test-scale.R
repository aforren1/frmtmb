## Phase 0 of dev/extension-gaps-plan.md: the ode row.
##
## 100 subjects x 8 samples, two states, twice-daily dosing for 7 days
## written with `ii`/`addl`, and one `ss` row. 800 rows.
##
## What the row decides, in the plan's own words: whether the segmented
## sensitivity solve is usable at population scale, or whether Phase 5's
## closed-form path (item 5.4) is a prerequisite rather than an option.
##
## Choices this file made where the plan left them open.
##
## "Two compartments" is read as the two-STATE depot and central system
## the package's own vignette fits, which is the one-compartment model
## with first-order absorption of the pharmacokinetic literature. It is
## the system every other test and the vignette use, so the row measures
## the package rather than a new model.
##
## The dosing table is two rows: a steady-state row at t = 0 with
## `ii = 12`, and a row at t = 12 with `ii = 12, addl = 12`, which is 14
## doses over the seven days with exactly one `ss` row, as the plan
## asks. `n_ss` is left at its default of 20, because that default is
## what a user pays.
##
## The starting values are plausible but not the truth (`ka` 0.8 against
## 1.0, `ke` 0.2 against 0.15, `V` 15 against 20), so that the timing is
## not the timing of a fit that starts at its own answer.
##
## See dev/scale-findings.md for the numbers this produced.

ode_truth <- list(ka = 1.0, ke = 0.15, V = 20, sigma = 0.3,
                  sd_lka = 0.3, sd_lke = 0.25, amt = 100, ii = 12,
                  ndose = 14L)

# The depot and central system of the package's vignette, by position.
ode_scale_dyn <- function(t, y, p) {
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}

# Twice daily for seven days: one steady-state row plus twelve more
# doses after the one at t = 12.
ode_scale_events <- function() {
  data.frame(time = c(0, 12), state = "depot", value = ode_truth$amt,
             ii = c(12, 12), addl = c(0L, 12L), ss = c(TRUE, FALSE))
}

# The day-seven profile over one dosing interval, which is what a
# population pharmacokinetic study with eight samples collects. The last
# sample sits exactly on the fourteenth dose, where frm_ode() reads the
# state BEFORE the dose, so it is the trough.
ode_scale_times <- function() {
  144 + c(0.5, 1, 2, 4, 6, 8, 10, 12)
}

# The closed form for this system, used to write the response. The `ss`
# row puts the system at steady state at t = 0 and every later dose is
# the same amount at the same interval, so the profile is PERIODIC in
# the interval and the standard steady-state superposition is exact.
# Writing the response from a solve would make the row a check of the
# solver against itself.
ode_scale_conc <- function(t, ka, ke, V, amt, ii) {
  u <- t %% ii
  amt * ka / (V * (ka - ke)) *
    (exp(-ke * u) / (1 - exp(-ke * ii)) -
       exp(-ka * u) / (1 - exp(-ka * ii)))
}

ode_scale_data <- function(seed = 20260908L,
                           ns = if (scale_small()) 3L else 100L) {
  tr <- ode_truth
  set.seed(seed)
  lka <- log(tr$ka) + stats::rnorm(ns, 0, tr$sd_lka)
  lke <- log(tr$ke) + stats::rnorm(ns, 0, tr$sd_lke)
  tt <- ode_scale_times()
  d <- data.frame(id = factor(rep(seq_len(ns), each = length(tt))),
                  time = rep(tt, times = ns))
  i <- as.integer(d$id)
  mu <- ode_scale_conc(d$time, exp(lka[i]), exp(lke[i]), tr$V, tr$amt,
                       tr$ii)
  d$conc <- mu + stats::rnorm(nrow(d), 0, tr$sigma)
  d
}

# The dosing table is bound in the formula's own environment, which is
# how ?frm_ode says a data.frame reaches a nonlinear body: a name there
# is normally a request for a column, and a data.frame is not something
# a column could hold.
ode_scale_form <- function() {
  ode_doses <- ode_scale_events()
  bf(conc ~ frm_ode(ode_scale_dyn,
                    init = list(0, 0),
                    times = time,
                    parms = list(exp(lka), exp(lke), exp(lV)),
                    group = id,
                    states = c("depot", "central"),
                    output = "central",
                    events = ode_doses),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}

test_that("the ode scale row fits and reports its cost", {
  skip_unless_scale()
  scale_row_on("ode")
  skip_if_not_installed("RTMBode")
  d <- ode_scale_data()
  ns <- nlevels(d$id)
  form <- ode_scale_form()
  st <- list(beta = c(log(0.8), log(0.2), log(15)))
  scale_mem_reset()

  bd <- scale_build(form + gaussian(), data = d, start = st)
  g0 <- scale_grad(bd$dry$obj, bd$dry$obj$par)
  ctl <- scale_control(bd$dry$obj, bd$dry$obj$par, g0$calls)
  bd$dry <- NULL

  fit <- NULL
  t_fit <- scale_elapsed(fit <- frm(form + gaussian(), data = d,
                                    start = st, se = TRUE))
  g1 <- scale_grad(fit$obj, fit$opt$par)
  mem <- scale_mem_peak_mb()

  b <- unlist(fixef(fit))
  ci <- suppressWarnings(stats::confint(fit))
  j <- grep("lke", rownames(ci), fixed = TRUE)
  i_ke <- if (length(j)) as.numeric(ci[j[1L], 1:2]) else c(NA, NA)
  vc <- VarCorr(fit)
  tr <- ode_truth
  scale_record(
    "ode", rows = nrow(d), subjects = ns,
    n_par = length(fit$opt$par),
    n_dose = tr$ndose, n_ss = 20L,
    build_s = bd$build_s, frame_s = bd$frame_s,
    grad_start_s = g0$seconds, grad_start_calls = g0$calls,
    grad_opt_s = g1$seconds, control_ratio = ctl,
    fit_s = t_fit, mem_mb = mem,
    logLik = as.numeric(stats::logLik(fit)),
    ka = exp(unname(b["lka.(Intercept)"])), ka_true = tr$ka,
    ke = exp(unname(b["lke.(Intercept)"])), ke_true = tr$ke,
    ke_lo = exp(i_ke[1L]), ke_hi = exp(i_ke[2L]),
    V = exp(unname(b["lV.(Intercept)"])), V_true = tr$V,
    sd_lka = sqrt(vc[[1L]][1L, 1L]), sd_lka_true = tr$sd_lka,
    sd_lke = sqrt(vc[[2L]][1L, 1L]), sd_lke_true = tr$sd_lke,
    n_failed = frm_ode_failures()$n_groups %||% 0L,
    diag = scale_diag(fit))

  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  # the elimination rate is what the design is powered for, so its own
  # Wald interval covering the simulator's truth is the assertion
  expect_true(i_ke[1L] <= log(tr$ke) && i_ke[2L] >= log(tr$ke))
})

# Shared simulation and ODE helpers for the population-PK ODE probes.
#
# Model: one-compartment, first-order oral absorption, single bolus dose
# D into the depot at t = 0.
#   dA/dt = -ka * A                    A(0) = D
#   dC/dt =  ka * A / V - ke * C       C(0) = 0
# Closed form (ka != ke):
#   C(t) = D*ka / (V*(ka - ke)) * (exp(-ke t) - exp(-ka t))
# The closed form is the simulation truth AND the reference the ODE
# solution is checked against, so an ODE failure is distinguishable from
# a model-specification failure.

PK_TRUTH <- list(
  lka = log(1.0), lke = log(0.2), lV = log(10),
  sd_lka = 0.30, sd_lke = 0.25, sigma = 0.30, dose = 100
)

pk_analytic <- function(t, ka, ke, V, D) {
  D * ka / (V * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))
}

sim_pk <- function(n_id = 12, seed = 2026) {
  set.seed(seed)
  tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
  id <- factor(rep(seq_len(n_id), each = length(tt)))
  time <- rep(tt, n_id)
  b_ka <- rnorm(n_id, 0, PK_TRUTH$sd_lka)
  b_ke <- rnorm(n_id, 0, PK_TRUTH$sd_lke)
  ka <- exp(PK_TRUTH$lka + b_ka)[as.integer(id)]
  ke <- exp(PK_TRUTH$lke + b_ke)[as.integer(id)]
  V <- exp(PK_TRUTH$lV)
  mu <- pk_analytic(time, ka, ke, V, PK_TRUTH$dose)
  data.frame(id = id, time = time, dose = PK_TRUTH$dose,
             conc = mu + rnorm(length(mu), 0, PK_TRUTH$sigma),
             mu_true = mu, b_ka = b_ka[as.integer(id)],
             b_ke = b_ke[as.integer(id)])
}

# Dynamics. Positional indexing only: p[["ka"]] style name lookup on an
# advector is not guaranteed, and c() needs the AD overload because base
# c() strips the advector class (frmtmb gotcha 8).
pk_dyn <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  dA <- -p[1] * y[1]
  dC <- p[1] * y[1] / p[3] - p[2] * y[2]
  list(c(dA, dC))
}

# Row-wise PK solve: one ODE solve per subject over that subject's own
# observation times, scattered back into row order. ka/ke/V arrive as
# ROW-wise vectors (that is all the nl body ever sees); they are constant
# within subject, so the subject's value is read off its first row.
pk_ode <- function(ka, ke, V, time, id, dose) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  ids <- as.integer(id)
  out <- numeric(length(ids))
  for (s in unique(ids)) {
    idx <- which(ids == s)
    ord <- order(time[idx])
    idx <- idx[ord]
    grid <- c(0, time[idx])          # deSolve wants times[1] == t0
    sol <- RTMBode::ode(
      y = c(dose[idx[1]], 0),
      times = grid,
      func = pk_dyn,
      parms = c(ka[idx[1]], ke[idx[1]], V[idx[1]]),
      method = "lsoda", atol = 1e-8, rtol = 1e-8
    )
    out[idx] <- sol[-1, 3]
  }
  out
}

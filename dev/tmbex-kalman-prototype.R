# Kalman-filter node prototype, measured on the `mvrw` TMB example.
#
# mvrw is a multivariate local-level model: a 3-dimensional random walk
# with correlated increments, observed with independent Gaussian noise
# and a flat (improper) prior on the first state. It is linear and
# Gaussian, so the Laplace approximation is EXACT. That gives a
# correctness check no other route offers: the same model marginalized
# two ways - Laplace over the 300 latent states, and a Kalman filter
# that integrates them out in closed form - must return the same
# log-likelihood at every parameter value, not just at the optimum.
#
# The question the prototype answers is whether a Kalman node is worth
# building as an RTMB::ADjoint atomic (a hand-written adjoint) or
# whether plain taped RTMB operations are already good enough.
#
#   Rscript dev/tmbex-kalman-prototype.R

library(RTMB)
say <- function(...) cat(sprintf(...), "\n", sep = "")

## Node count of a taped closure. RTMB's tape printer emits one line per
## node plus a header, which is the only introspection RTMB exposes.
tape_nodes <- function(f, x) {
  F <- MakeTape(f, x)
  length(capture.output(environment(F)$mod$print(0))) - 1L
}

## ------------------------------------------------------------------ data
## mvrw.R's own simulation, inlined (the example plots, which we do not).
set.seed(1)
suppressMessages(library(MASS))
stateDim <- 3
timeSteps <- 100
rho_true <- 0.9
sds_true <- seq(0.5, 2, length = stateDim)
sdObs_true <- rep(1, stateDim)
Sig_true <- rho_true^abs(outer(1:stateDim, 1:stateDim, "-")) *
  (sds_true %o% sds_true)
dstate <- matrix(NA, timeSteps, stateDim)
obs <- dstate
dstate[1, ] <- rnorm(stateDim)
obs[1, ] <- dstate[1, ] + rnorm(stateDim, 0, sdObs_true)
for (i in 2:timeSteps) {
  dstate[i, ] <- dstate[i - 1, ] + mvrnorm(1, rep(0, stateDim), Sig_true)
  obs[i, ] <- dstate[i, ] + rnorm(stateDim, 0, sdObs_true)
}
Y <- t(obs)                            # stateDim x timeSteps, as mvrw.R

trf <- function(x) 2 / (1 + exp(-2 * x)) - 1

## Increment covariance and observation variance from the 7 free
## parameters, in mvrw.R's own order and transforms.
pieces <- function(th) {
  rho <- trf(th[1])
  sds <- exp(th[2:4])
  so <- exp(th[5:7])
  Sigma <- outer(1:stateDim, 1:stateDim,
                 function(i, j) rho^(abs(i - j)) * sds[i] * sds[j])
  list(Sigma = Sigma, R = so^2 * diag(stateDim))
}

## ------------------------------------------------- route A: Laplace (mvrw.R)
joint_nll <- function(p) {
  pc <- pieces(c(p$transf_rho, p$logsds, p$logsdObs))
  -sum(dmvnorm(diff(t(p$u)), 0, pc$Sigma, log = TRUE)) -
    sum(dnorm(Y, p$u, exp(p$logsdObs), log = TRUE))
}
p_lap <- list(transf_rho = 0.1, logsds = rep(0, stateDim),
              logsdObs = rep(0, stateDim), u = Y * 0)
t_tape_lap <- system.time(
  obj_lap <- MakeADFun(joint_nll, p_lap, random = "u", silent = TRUE))[3]
TMB::newtonOption(obj_lap, smartsearch = FALSE)
t_fit_lap <- system.time(
  opt_lap <- nlminb(obj_lap$par, obj_lap$fn, obj_lap$gr))[3]
n_lap <- tape_nodes(function(x) joint_nll(list(
  transf_rho = x[1], logsds = x[2:4], logsdObs = x[5:7],
  u = matrix(x[-(1:7)], stateDim, timeSteps))),
  c(0.1, rep(0, 6), rep(0, stateDim * timeSteps)))

## ---------------------------------------- route B(i): filter in plain RTMB
## Flat prior on the first state: integrating N(y1; u1, R) over u1 gives
## exactly 1, so the filter starts from a1|1 = y1, P1|1 = R and the first
## observation contributes nothing. That is what makes this filter's
## marginal likelihood equal to the Laplace one, constant and all.
kalman_nll <- function(th) {
  pc <- pieces(th)
  a <- Y[, 1]
  P <- pc$R
  nll <- 0
  for (t in 2:timeSteps) {
    P <- P + pc$Sigma                       # predict (transition is I)
    v <- Y[, t] - a                         # innovation
    Fm <- P + pc$R
    nll <- nll - dmvnorm(v, 0, Fm, log = TRUE)
    K <- P %*% solve(Fm)                    # gain
    a <- a + as.vector(K %*% v)
    P <- P - K %*% P
  }
  nll
}
th0 <- c(0.1, rep(0, 6))
t_tape_kal <- system.time(
  obj_kal <- MakeADFun(function(p) kalman_nll(p$th), list(th = th0),
                       silent = TRUE))[3]
t_fit_kal <- system.time(
  opt_kal <- nlminb(obj_kal$par, obj_kal$fn, obj_kal$gr))[3]
n_kal <- tape_nodes(kalman_nll, th0)

## ------------------------------------- route B(ii): the filter as an atomic
## RTMB can turn any taped closure into an atomic function whose adjoint
## is the inner tape's own reverse sweep ($atomic(), as tmb_examples'
## dataeval.R uses it). That is the REGISTER_ATOMIC capability without a
## hand-derived adjoint: the outer tape sees one node.
t_tape_atm <- system.time({
  Fk <- MakeTape(kalman_nll, th0)
  Fk$simplify()
  Fka <- Fk$atomic()
  obj_atm <- MakeADFun(function(p) Fka(p$th), list(th = th0), silent = TRUE)
})[3]
t_fit_atm <- system.time(
  opt_atm <- nlminb(obj_atm$par, obj_atm$fn, obj_atm$gr))[3]
n_atm <- tape_nodes(function(x) Fka(x), th0)

## --------------------------------------------------------- the identity
## Laplace and Kalman are two marginalizations of ONE model, so they must
## agree at every parameter value, not only at the optimum.
grid <- rbind(th0,
              c(0.5, 0.2, -0.3, 0.4, -0.1, 0.2, 0.0),
              c(-0.4, -0.5, 0.6, 0.1, 0.3, -0.2, 0.5),
              opt_lap$par)
ident <- apply(grid, 1, function(th) {
  abs(obj_lap$fn(th) - obj_kal$fn(th))
})

say("=== mvrw: Laplace vs Kalman ===")
say("%-26s %12s %10s %10s %16s", "route", "tape nodes", "tape s", "fit s",
    "logLik")
say("%-26s %12d %10.2f %10.2f %16.8f", "Laplace (mvrw.R, 300 RE)",
    n_lap, t_tape_lap, t_fit_lap, -opt_lap$objective)
say("%-26s %12d %10.2f %10.2f %16.8f", "Kalman, plain RTMB ops",
    n_kal, t_tape_kal, t_fit_kal, -opt_kal$objective)
say("%-26s %12d %10.2f %10.2f %16.8f", "Kalman, MakeTape$atomic()",
    n_atm, t_tape_atm, t_fit_atm, -opt_atm$objective)
say("")
say("identity |nll_laplace - nll_kalman| over 4 parameter vectors: %s",
    paste(sprintf("%.2e", ident), collapse = "  "))
say("|logLik_laplace - logLik_kalman| at each optimum: %.3e",
    abs(opt_lap$objective - opt_kal$objective))
say("max |par_laplace - par_kalman|: %.3e",
    max(abs(opt_lap$par - opt_kal$par)))
say("")
say("sdreport standard errors, Laplace vs Kalman:")
s1 <- sqrt(diag(sdreport(obj_lap)$cov.fixed))
s2 <- sqrt(diag(sdreport(obj_kal)$cov.fixed))
say("  max |se difference| = %.3e", max(abs(s1 - s2)))

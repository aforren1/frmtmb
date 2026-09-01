# Probe E5: the stacked PK differs from the E3/E4 toys in one more way -
# it runs under MakeADFun(random = ), i.e. a Laplace inner problem that
# needs SECOND derivatives through the ODE. Test the stacked solve with
# and without random effects, against the per-subject loop.
library(RTMB)
library(RTMBode)
source("dev/ode/pk-common.R")
d <- sim_pk()
gi <- as.integer(d$id); ng <- nlevels(d$id)
conc <- d$conc; time <- d$time; id <- d$id; dose <- d$dose

pk_dyn_stack <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  m <- length(y) / 2
  h <- 1:m
  A <- y[h]; C <- y[m + h]
  list(c(-p[h] * A, p[h] * A / p[2 * m + h] - p[m + h] * C))
}
pk_ode_stack <- function(ka, ke, V, time, id, dose) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  ids <- as.integer(factor(id)); ns <- max(ids)
  first <- match(seq_len(ns), ids)
  grid <- sort(unique(c(0, time)))
  sol <- RTMBode::ode(y = c(dose[first], rep(0, ns)), times = grid,
                      func = pk_dyn_stack,
                      parms = c(ka[first], ke[first], V[first]),
                      method = "lsoda", atol = 1e-8, rtol = 1e-8)
  out <- numeric(length(time))
  for (s in seq_len(ns)) {
    idx <- which(ids == s)
    out[idx] <- sol[match(time[idx], grid), 1L + ns + s]
  }
  out
}

mk <- function(solver, random) {
  f <- function(p) {
    getAll(p)
    mu <- solver(exp(lka + u_ka[gi]), exp(lke + u_ke[gi]),
                 exp(rep(lV, length(conc))), time, id, dose)
    -sum(dnorm(u_ka, 0, exp(lsd_ka), log = TRUE)) -
      sum(dnorm(u_ke, 0, exp(lsd_ke), log = TRUE)) -
      sum(dnorm(conc, mu, exp(lsigma), log = TRUE))
  }
  MakeADFun(f, list(lka = 0, lke = log(0.25), lV = log(8), lsigma = 0,
                    lsd_ka = log(0.3), lsd_ke = log(0.3),
                    u_ka = numeric(ng), u_ke = numeric(ng)),
            random = random, silent = TRUE)
}
report <- function(label, solver, random) {
  o <- try(mk(solver, random), silent = TRUE)
  if (inherits(o, "try-error")) { cat(label, ": TAPE ERROR\n"); return() }
  v <- suppressWarnings(try(o$fn(o$par), silent = TRUE))
  g <- suppressWarnings(try(o$gr(o$par), silent = TRUE))
  cat(sprintf("%-34s fn = %-12s max|gr| = %s\n", label,
              if (inherits(v, "try-error")) "ERROR" else format(v, digits = 8),
              if (inherits(g, "try-error")) "ERROR" else
                format(max(abs(g)), digits = 6)))
}
report("loop,    random = NULL", pk_ode, NULL)
report("stacked, random = NULL", pk_ode_stack, NULL)
report("loop,    random = u", pk_ode, c("u_ka", "u_ke"))
report("stacked, random = u", pk_ode_stack, c("u_ka", "u_ke"))

cat("\n--- second derivatives through one ode() node ---\n")
# Direct question: is the adjoint ODE node twice differentiable? Laplace
# needs it. One state, one parameter, no frmtmb.
f1 <- function(p) {
  sol <- RTMBode::ode(y = 1, times = c(0, 1, 2), func =
                        function(t, y, q) list(-q * y),
                      parms = exp(p$lr), method = "lsoda")
  sum(sol[, 2]^2)
}
o1 <- MakeADFun(f1, list(lr = log(0.3)), silent = TRUE)
safe <- function(e) suppressWarnings(tryCatch(format(force(e)),
  error = function(err) paste("ERROR:", conditionMessage(err))))
cat("fn:", safe(o1$fn(0)), "\n")
cat("gr:", safe(o1$gr(0)), "\n")
cat("he:", safe(o1$he(0)), "\n")
# and a Laplace approximation over a random effect entering the dynamics
f2 <- function(p) {
  getAll(p)
  sol <- RTMBode::ode(y = 1, times = c(0, 1, 2), func =
                        function(t, y, q) list(-q * y),
                      parms = exp(lr + u), method = "lsoda")
  -dnorm(u, 0, 1, log = TRUE) + sum((sol[, 2] - 0.5)^2)
}
o2 <- try(MakeADFun(f2, list(lr = log(0.3), u = 0), random = "u",
                    silent = TRUE), silent = TRUE)
if (inherits(o2, "try-error")) cat("Laplace tape ERROR\n") else
  cat("Laplace fn:", suppressWarnings(o2$fn(log(0.3))),
      " gr:", suppressWarnings(o2$gr(log(0.3))), "\n")

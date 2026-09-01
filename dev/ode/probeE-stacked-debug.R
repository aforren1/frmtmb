# Probe E: why the stacked (one-solve-for-the-whole-population) variant
# broke. This matters for the design: stacking would collapse N adjoint
# ODE nodes into one.
library(RTMB)
library(RTMBode)
source("dev/ode/pk-common.R")
d <- sim_pk()

pk_dyn_stack <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  ns <- length(y) / 2
  A <- y[1:ns]; C <- y[ns + 1:ns]
  ka <- p[1:ns]; ke <- p[ns + 1:ns]; V <- p[2 * ns + 1:ns]
  list(c(-ka * A, ka * A / V - ke * C))
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
  cat("   sol class:", class(sol)[1], " dim:",
      paste(dim(sol), collapse = "x"), "\n")
  sol[cbind(match(time, grid), 1L + ns + ids)]
}

cat("--- numeric mode ---\n")
ka <- exp(PK_TRUTH$lka + d$b_ka); ke <- exp(PK_TRUTH$lke + d$b_ke)
V <- rep(exp(PK_TRUTH$lV), nrow(d))
num <- try(pk_ode_stack(ka, ke, V, d$time, d$id, d$dose), silent = TRUE)
if (inherits(num, "try-error")) cat("NUMERIC FAILED:", as.character(num)) else
  cat("max |stacked - analytic| =", format(max(abs(num - d$mu_true)),
                                           digits = 3), "\n")

cat("\n--- taped mode: value and gradient at the truth ---\n")
gi <- as.integer(d$id); ng <- nlevels(d$id)
conc <- d$conc; time <- d$time; id <- d$id; dose <- d$dose
mk <- function(solver) {
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
            random = c("u_ka", "u_ke"), silent = TRUE)
}
for (nm in c("loop", "stacked")) {
  solver <- if (nm == "loop") pk_ode else pk_ode_stack
  o <- try(mk(solver), silent = TRUE)
  if (inherits(o, "try-error")) { cat(nm, "TAPE FAILED:", as.character(o));
    next }
  v <- try(o$fn(o$par), silent = TRUE)
  g <- try(o$gr(o$par), silent = TRUE)
  cat(nm, ": fn =", if (inherits(v, "try-error")) "ERROR" else format(v),
      " gr =", if (inherits(g, "try-error")) "ERROR" else
        paste(format(g, digits = 4), collapse = " "), "\n")
}

cat("\n--- is it the matrix indexing or the solve? ---\n")
# same stacked solve, but extract with an explicit per-subject column read
pk_ode_stack2 <- function(ka, ke, V, time, id, dose) {
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
o <- try(mk(pk_ode_stack2), silent = TRUE)
if (inherits(o, "try-error")) cat("stack2 TAPE FAILED:", as.character(o)) else
  cat("stack2: fn =", format(o$fn(o$par)), " gr =",
      paste(format(o$gr(o$par), digits = 4), collapse = " "), "\n")

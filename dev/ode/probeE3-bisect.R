# Probe E3: bisect the stacked-PK AD failure. The minimal 32-state decay
# system tapes fine, so the failure is in some other axis: parms length,
# initial-state magnitude, c() in the dynamics, or the two-block layout.
library(RTMB)
library(RTMBode)
times <- seq(0, 5, length.out = 9)

try_case <- function(label, y0, np, dyn) {
  f <- function(par) {
    sol <- RTMBode::ode(y = y0, times = times, func = dyn,
                        parms = exp(par$lp), method = "lsoda",
                        atol = 1e-8, rtol = 1e-8)
    sum(sol[, -1]^2)
  }
  obj <- try(MakeADFun(f, list(lp = rep(log(0.3), np)), silent = TRUE),
             silent = TRUE)
  v <- if (inherits(obj, "try-error")) NA else
    suppressWarnings(try(obj$fn(obj$par), silent = TRUE))
  cat(sprintf("%-52s fn = %s\n", label,
              if (inherits(v, "try-error") || is.na(v)) "FAIL/NaN" else
                format(v, digits = 8)))
}

ns <- 12
# A: baseline, parms length == states, no c() in dynamics
try_case("A  ns=12, np=12, no c()", rep(1, 2 * ns), 2 * ns,
         function(t, y, p) list(-p * y))
# B: parms LONGER than the state vector
try_case("B  ns=12, np=36 (parms longer than states)", rep(1, 2 * ns),
         3 * ns, function(t, y, p) list(-p[1:(2 * ns)] * y))
# C: c() concatenating two half-vectors in the dynamics
try_case("C  ns=12, np=24, c() of two halves", rep(1, 2 * ns), 2 * ns,
         function(t, y, p) {
           "c" <- RTMB::ADoverload("c")
           h <- 1:ns
           list(c(-p[h] * y[h], -p[ns + h] * y[ns + h]))
         })
# D: large initial state (dose = 100)
try_case("D  ns=12, np=24, y0 = 100", rep(100, 2 * ns), 2 * ns,
         function(t, y, p) list(-p * y))
# E: coupled two-block PK layout, parms length 3*ns
try_case("E  ns=12, np=36, coupled PK layout", c(rep(100, ns), rep(0, ns)),
         3 * ns,
         function(t, y, p) {
           "c" <- RTMB::ADoverload("c")
           h <- 1:ns
           A <- y[h]; C <- y[ns + h]
           ka <- p[h]; ke <- p[ns + h]; V <- p[2 * ns + h]
           list(c(-ka * A, ka * A / V - ke * C))
         })
# F: same as E but ns = 2 (does it depend on size?)
for (nsx in c(1, 2, 3, 4, 6, 8, 12)) {
  local({
    ns <- nsx
    try_case(sprintf("F  coupled PK layout, ns=%d", ns),
             c(rep(100, ns), rep(0, ns)), 3 * ns,
             function(t, y, p) {
               "c" <- RTMB::ADoverload("c")
               h <- 1:ns
               A <- y[h]; C <- y[ns + h]
               ka <- p[h]; ke <- p[ns + h]; V <- p[2 * ns + h]
               list(c(-ka * A, ka * A / V - ke * C))
             })
  })
}
# G: E's stiffness, softer tolerances / different integrator
for (m in c("lsode", "adams", "ode45", "rk4")) {
  ns <- 12
  f <- function(par) {
    sol <- RTMBode::ode(y = c(rep(100, ns), rep(0, ns)), times = times,
                        func = function(t, y, p) {
                          "c" <- RTMB::ADoverload("c")
                          h <- 1:ns
                          A <- y[h]; C <- y[ns + h]
                          ka <- p[h]; ke <- p[ns + h]; V <- p[2 * ns + h]
                          list(c(-ka * A, ka * A / V - ke * C))
                        },
                        parms = exp(par$lp), method = m,
                        atol = 1e-8, rtol = 1e-8)
    sum(sol[, -1]^2)
  }
  obj <- try(MakeADFun(f, list(lp = rep(log(0.3), 3 * ns)), silent = TRUE),
             silent = TRUE)
  v <- if (inherits(obj, "try-error")) NA else
    suppressWarnings(try(obj$fn(obj$par), silent = TRUE))
  cat(sprintf("G  coupled PK, ns=12, method=%-6s          fn = %s\n", m,
              if (inherits(v, "try-error") || is.na(v)) "FAIL/NaN" else
                format(v, digits = 8)))
}

## Events on the RTMBode autodiff path.
##
## Before the patch:
##   - any events$data spelling failed with
##     "too many state variables in 'event'; should be < 0", because the state
##     vector handed to deSolve carried no names;
##   - had that been fixed, 'replace' and 'multiply' would have returned wrong
##     gradients silently (deSolve jumps the state block of the augmented
##     system and leaves the sensitivity block alone).
##
## After the patch: 'add' works and agrees with central differences of the
## numeric solve; 'replace' and 'multiply' are refused by name.
##
## Usage: Rscript 06-events.R <lib> <repro-dir>

args <- commandArgs(TRUE)
if (length(args) >= 1 && nzchar(args[1])) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages({ library(RTMB); library(RTMBode) })

## One-compartment oral PK: depot -> central, first-order in and out.
pk <- function(t, y, p) list(c(-p[1] * y[1], p[1] * y[1] - p[2] * y[2]))
y0 <- c(depot = 100, central = 0)
times <- seq(0, 24, by = 0.5)
dosetimes <- c(6, 12, 18)

events_tab <- function(method)
    list(data = data.frame(var = "depot", time = dosetimes,
                           value = 50, method = method))

## Numeric solve, used for the finite-difference reference.
numsol <- function(par, method) {
    deSolve::ode(y0, times, pk, exp(par), events = events_tab(method),
                 atol = 1e-10, rtol = 1e-10)
}

## Objective through the adjoint path: sum of the central compartment.
adobj <- function(method) {
    f <- function(pl) {
        sol <- RTMBode::ode(y0, times, pk, exp(pl$par),
                            events = events_tab(method),
                            atol = 1e-10, rtol = 1e-10)
        sum(sol[, "central"])
    }
    MakeADFun(f, list(par = log(c(ka = 1.2, ke = 0.35))), silent = TRUE)
}

fd <- function(method) {
    p <- log(c(1.2, 0.35))
    h <- 1e-5
    sapply(seq_along(p), function(i) {
        pp <- pm <- p
        pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
        (sum(numsol(pp, method)[, "central"]) -
         sum(numsol(pm, method)[, "central"])) / (2 * h)
    })
}

cat("== numeric path (unchanged by the patch) ==\n")
s <- numsol(log(c(1.2, 0.35)), "add")
cat("nrow:", nrow(s), " depot at t=6 (pre-dose):", signif(s[s[, 1] == 6, "depot"], 8), "\n")

for (method in c("add", "replace", "multiply")) {
    cat("\n== method =", method, "==\n")
    r <- tryCatch({
        obj <- adobj(method)
        p <- log(c(1.2, 0.35))
        list(value = obj$fn(p), gr = as.numeric(obj$gr(p)))
    }, error = function(e) conditionMessage(e))
    if (is.character(r)) {
        cat("REFUSED:", gsub("\\s+", " ", r), "\n")
    } else {
        ref <- fd(method)
        cat("value  :", signif(r$value, 10), " (numeric solve:",
            signif(sum(numsol(p <- log(c(1.2, 0.35)), method)[, "central"]), 10), ")\n")
        cat("AD grad:", signif(r$gr, 8), "\n")
        cat("FD grad:", signif(ref, 8), "\n")
        cat("max rel err:", signif(max(abs(r$gr - ref) / pmax(abs(ref), 1e-8)), 3), "\n")
    }
}

## The event table with no 'method' column: deSolve's default is 'replace'.
cat("\n== events with no method column (deSolve default = replace) ==\n")
r <- tryCatch({
    f <- function(pl) {
        sol <- RTMBode::ode(y0, times, pk, exp(pl$par),
                            events = list(data = data.frame(var = "depot",
                                                            time = dosetimes,
                                                            value = 50)),
                            atol = 1e-10, rtol = 1e-10)
        sum(sol[, "central"])
    }
    obj <- MakeADFun(f, list(par = log(c(1.2, 0.35))), silent = TRUE)
    obj$fn(log(c(1.2, 0.35)))
}, error = function(e) conditionMessage(e))
cat(if (is.numeric(r)) paste("ACCEPTED (value", r, ")") else paste("REFUSED:", gsub("\\s+", " ", r)), "\n")

## An event function cannot be validated, so it is refused too.
cat("\n== events$func ==\n")
r <- tryCatch({
    f <- function(pl) {
        sol <- RTMBode::ode(y0, times, pk, exp(pl$par),
                            events = list(func = function(t, y, p) y,
                                          time = dosetimes),
                            atol = 1e-10, rtol = 1e-10)
        sum(sol[, "central"])
    }
    obj <- MakeADFun(f, list(par = log(c(1.2, 0.35))), silent = TRUE)
    obj$fn(log(c(1.2, 0.35)))
}, error = function(e) conditionMessage(e))
cat(if (is.numeric(r)) paste("ACCEPTED (value", r, ")") else paste("REFUSED:", gsub("\\s+", " ", r)), "\n")

## No events at all: the named-state change must not disturb ordinary solves.
cat("\n== regression: no events, gradient vs finite differences ==\n")
f <- function(pl) sum(RTMBode::ode(y0, times, pk, exp(pl$par),
                                   atol = 1e-10, rtol = 1e-10)[, "central"])
obj <- MakeADFun(f, list(par = log(c(1.2, 0.35))), silent = TRUE)
p <- log(c(1.2, 0.35))
g <- as.numeric(obj$gr(p))
h <- 1e-5
ref <- sapply(1:2, function(i) {
    pp <- pm <- p; pp[i] <- pp[i] + h; pm[i] <- pm[i] - h
    (sum(deSolve::ode(y0, times, pk, exp(pp), atol = 1e-10, rtol = 1e-10)[, "central"]) -
     sum(deSolve::ode(y0, times, pk, exp(pm), atol = 1e-10, rtol = 1e-10)[, "central"])) / (2 * h)
})
cat("AD:", signif(g, 8), "\nFD:", signif(ref, 8), "\n")
cat("max rel err:", signif(max(abs(g - ref) / abs(ref)), 3), "\n")
cat("colnames of the returned matrix:", colnames(RTMBode::ode(y0, times, pk, exp(p))), "\n")

## RTMBode + tmbstan: minimal reproduction of the warmup-iteration-1 abort.
##
## Usage: Rscript 01-sampler.R <ngroup> <lib> <repro-dir>
##   ngroup : 1 or 2 ODE adjoint nodes in the tape
##   lib    : library path holding the RTMBode build under test ("" = default)
##
## Steps: build the objective, optimize, confirm fn/gr work at the mode, probe
## fn/gr away from the mode (the size of jump Stan's first stepsize search
## makes), then hand the same object to tmbstan starting AT the mode.

args <- commandArgs(TRUE)
ngroup <- as.integer(if (length(args) >= 1) args[1] else 1L)
if (length(args) >= 2 && nzchar(args[2])) .libPaths(c(args[2], .libPaths()))
here <- if (length(args) >= 3) args[3] else "."
source(file.path(here, "lv-model.R"))

cat("RTMBode from:", dirname(system.file(package = "RTMBode")), "\n")
cat("ngroup:", ngroup, "\n\n")

obs <- lv_data(ngroup)
obj <- lv_objective(obs)

cat("== direct evaluation at the start value ==\n")
cat("fn:", obj$fn(obj$par), "\n")
cat("gr:", obj$gr(obj$par), "\n")

opt <- nlminb(obj$par, obj$fn, obj$gr)
cat("\n== optimum ==\n")
cat("objective:", opt$objective, " convergence:", opt$convergence, "\n")
mode <- opt$par
cat("fn(mode):", obj$fn(mode), "\n")
cat("gr(mode):", obj$gr(mode), "\n")

## How far can the parameter move before the solve fails?  Stan's initial
## stepsize is 1 on the unconstrained scale, so warmup iteration 1 routinely
## proposes a point this far away.
cat("\n== direct evaluation away from the mode (no sampler involved) ==\n")
set.seed(42)
for (eps in c(0.1, 0.25, 0.5, 1, 2)) {
    for (rep in 1:3) {
        p <- mode + rnorm(length(mode), sd = eps)
        r <- tryCatch(list(fn = obj$fn(p), gr = max(abs(obj$gr(p)))),
                      error = function(e) conditionMessage(e))
        cat(sprintf("eps=%-5g rep=%d  %s\n", eps, rep,
                    if (is.character(r)) paste("ERROR:", r)
                    else sprintf("fn=%.6g  max|gr|=%.6g", r$fn, r$gr)))
    }
}

cat("\n== tmbstan, init at the mode ==\n")
res <- tryCatch({
    fit <- tmbstan::tmbstan(obj, chains = 1, iter = 30, warmup = 15,
                            init = list(mode), seed = 1)
    n <- tryCatch(nrow(as.matrix(fit)), error = function(e) 0L)
    if (is.null(n) || is.na(n) || n == 0L) {
        cat("SAMPLER PRODUCED NO DRAWS (rstan aborted the chain)\n")
        "no draws"
    } else {
        cat("SAMPLER OK; draws:", n, "\n")
        "ok"
    }
}, error = function(e) {
    cat("SAMPLER ERROR:", conditionMessage(e), "\n")
    "error"
})
cat("RESULT:", res, "\n")

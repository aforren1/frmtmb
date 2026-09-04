## Longer, harder sampler run on the patched build: two ODE adjoint nodes,
## a full 400-iteration chain, and a deliberately over-dispersed init so that
## warmup must survive many failed solves rather than a lucky few.
##
## Usage: Rscript 05-sampler-long.R <lib> <repro-dir>

args <- commandArgs(TRUE)
if (length(args) >= 1 && nzchar(args[1])) .libPaths(c(args[1], .libPaths()))
here <- if (length(args) >= 2) args[2] else "."
source(file.path(here, "lv-model.R"))

obs <- lv_data(2L)
obj <- lv_objective(obs)
mode <- nlminb(obj$par, obj$fn, obj$gr)$par
cat("mode fn:", obj$fn(mode), "\n")

run <- function(label, init, iter) {
    cat("\n==", label, "==\n")
    t0 <- proc.time()[3]
    out <- tryCatch({
        fit <- tmbstan::tmbstan(obj, chains = 1, iter = iter,
                                init = list(init), seed = 11, refresh = 0)
        m <- as.matrix(fit)
        cat("draws:", nrow(m), " params:", ncol(m), "\n")
        cat("posterior mean:", signif(colMeans(m)[1:6], 5), "\n")
        cat("divergences:",
            sum(sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                       function(z) sum(z[, "divergent__"]))), "\n")
        "ok"
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat("outcome:", out, " seconds:", round(proc.time()[3] - t0, 1), "\n")
    out
}

r1 <- run("two ODE nodes, init at the mode, 400 iterations", mode, 400)
## Start well away from the mode: warmup now has to reject many failed solves.
set.seed(99)
r2 <- run("two ODE nodes, over-dispersed init", mode + rnorm(length(mode), sd = 1.5), 400)

cat("\nSUMMARY:", r1, "|", r2, "\n")

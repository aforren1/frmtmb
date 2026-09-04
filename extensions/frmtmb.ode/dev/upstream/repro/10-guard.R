## The workspace guard: what the user sees instead of deSolve's allocator
## message, and proof that the guard does not fire on ordinary models.
##
## Usage: Rscript 10-guard.R <n> <method> <lib>

args <- commandArgs(TRUE)
n <- as.integer(args[1]); method <- args[2]
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages({ library(RTMB); library(RTMBode) })

times <- c(0, 1, 2, 3)
set.seed(1)
obs <- exp(-0.3 * rep(times[-1], each = n)) + rnorm(n * 3, 0, 0.1)
f <- function(p) {
    getAll(p)
    sol <- RTMBode::ode(y = rep(1, n), times = times,
                        func = function(t, y, q) list(-q * y),
                        parms = exp(mu + u), method = method)
    -sum(dnorm(u, 0, 1, log = TRUE)) +
        sum((obs - as.vector(t(sol[-1, -1, drop = FALSE])))^2)
}
o <- MakeADFun(f, list(mu = log(0.3), u = numeric(n)), random = "u", silent = TRUE)
cat(sprintf("n = %d, method = %s\n", n, method))
cat("fn:", format(suppressWarnings(o$fn(log(0.3))), digits = 8), "\n")
w <- NULL
g <- withCallingHandlers(
    tryCatch(o$gr(log(0.3))[1],
             error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
cat("gr:", if (is.character(g)) gsub("\\s+", " ", g) else format(g, digits = 8), "\n")
if (length(w)) for (m in w) cat("WARNING:", gsub("\\s+", " ", m), "\n")

## What actually fails at the state ceiling.
##
## Records the size of every system deSolve is asked to integrate, and the
## error deSolve raises, for a Laplace objective containing one ODE node.
##
## Usage: Rscript 08-ceiling-mechanism.R <n> <method> <lib>

args <- commandArgs(TRUE)
n <- as.integer(args[1])
method <- args[2]
if (length(args) >= 3 && nzchar(args[3])) .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages({ library(RTMB); library(RTMBode) })

seen <- new.env(parent = emptyenv())
seen$neq <- integer(0)
seen$err <- character(0)
trace(deSolve::ode, tracer = quote({
    seen <- get("seen", envir = globalenv())
    seen$neq <- c(seen$neq, length(y))
}), print = FALSE, where = asNamespace("deSolve"))

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
seen$neq <- integer(0)
v <- suppressWarnings(tryCatch(o$fn(log(0.3)), error = function(e) conditionMessage(e)))
fn_neq <- sort(unique(seen$neq))
seen$neq <- integer(0)
g <- suppressWarnings(tryCatch(o$gr(log(0.3))[1], error = function(e) conditionMessage(e)))
gr_neq <- sort(unique(seen$neq))
untrace(deSolve::ode, where = asNamespace("deSolve"))

aug <- function(k) sum(n * (2 * n)^(0:k))
lrw <- function(neq) 22 + neq * max(16, neq + 9)   # deSolve's lsoda work array

cat(sprintf("n = %d, method = %s\n", n, method))
cat("augmented sizes by order:", paste(sprintf("%d:%d", 0:3, sapply(0:3, aug)), collapse = "  "), "\n")
cat("systems integrated during fn():", paste(fn_neq, collapse = ", "), "\n")
cat("systems integrated during gr():", paste(gr_neq, collapse = ", "), "\n")
big <- max(c(fn_neq, gr_neq))
cat(sprintf("largest system: neq = %d -> lsoda lrw = %.4g doubles = %.3g GB\n",
            big, lrw(big), lrw(big) * 8 / 1024^3))
cat("fn:", if (is.character(v)) paste("ERROR:", v) else format(v, digits = 8), "\n")
cat("gr:", if (is.character(g)) paste("ERROR:", g) else format(g, digits = 8), "\n")

## What does a bare deSolve call of that size do?
if (big > 0) {
    y <- rep(1, big)
    r <- tryCatch({
        deSolve::ode(y, c(0, 1), function(t, y, p) list(-y), NULL, method = method)
        "ok"
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat("bare deSolve at neq =", big, ":", r, "\n")
}

## Narrowing: which parameter values make deSolve raise an R *error* (as opposed
## to returning a non-finite solution)?  Stan's first stepsize search multiplies
## a gradient of order 1e4-1e5 by epsilon = 1, so warmup iteration 1 evaluates
## the objective at parameter values many orders of magnitude from the mode.
##
## Usage: Rscript 03-infparm.R <lib> <repro-dir>

args <- commandArgs(TRUE)
if (length(args) >= 1 && nzchar(args[1])) .libPaths(c(args[1], .libPaths()))
here <- if (length(args) >= 2) args[2] else "."
source(file.path(here, "lv-model.R"))

obs <- lv_data(1L)
obj <- lv_objective(obs)
mode <- nlminb(obj$par, obj$fn, obj$gr)$par

classify <- function(p, what) {
    r <- tryCatch(list(v = what(p)), error = function(e) conditionMessage(e))
    if (is.character(r)) paste0("ERROR: ", gsub("\\s+", " ", r))
    else if (any(!is.finite(r$v))) "non-finite" else "finite"
}

cat("== single-coordinate scale sweep (logpar[2] = log growth rate) ==\n")
for (v in c(10, 20, 50, 100, 500, 710, 1000, 1e4, 1e5)) {
    p <- mode; p[2] <- v
    cat(sprintf("logpar[2]=%-8g exp=%-10.3g  fn: %s\n", v, exp(v), classify(p, obj$fn)))
}

cat("\n== all coordinates blown up together ==\n")
for (v in c(100, 710, 1000, 1e5)) {
    p <- mode + v
    cat(sprintf("mode+%-8g  fn: %s\n", v, classify(p, obj$fn)))
}

cat("\n== gradient at the same points ==\n")
for (v in c(100, 1000, 1e5)) {
    p <- mode + v
    cat(sprintf("mode+%-8g  gr: %s\n", v, classify(p, obj$gr)))
}

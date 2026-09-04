## Discriminating experiment: can the *direct* path (no sampler, no tmbstan)
## raise the same R error?  If yes, the tmbstan failure is not a sampler
## interaction at all: it is an unguarded deSolve error inside the adjoint node.
##
## Usage: Rscript 02-direct-error.R <lib> <repro-dir>

args <- commandArgs(TRUE)
if (length(args) >= 1 && nzchar(args[1])) .libPaths(c(args[1], .libPaths()))
here <- if (length(args) >= 2) args[2] else "."
source(file.path(here, "lv-model.R"))

obs <- lv_data(1L)
obj <- lv_objective(obs)
opt <- nlminb(obj$par, obj$fn, obj$gr)
mode <- opt$par
cat("mode fn:", obj$fn(mode), "\n\n")

classify <- function(p, what) {
    r <- tryCatch(list(v = what(p)), error = function(e) conditionMessage(e))
    if (is.character(r)) paste0("ERROR: ", gsub("\\s+", " ", r))
    else if (any(!is.finite(r$v))) "non-finite"
    else "finite"
}

cat("== sweep of |offset| on the direct path ==\n")
set.seed(7)
tab <- list()
for (eps in c(1, 2, 3, 4, 6)) {
    for (rep in 1:6) {
        p <- mode + rnorm(length(mode), sd = eps)
        fnres <- classify(p, obj$fn)
        grres <- if (startsWith(fnres, "ERROR")) "(skipped)" else classify(p, obj$gr)
        cat(sprintf("eps=%-3g rep=%d  fn: %-70s gr: %s\n", eps, rep, fnres, grres))
        tab[[length(tab) + 1]] <- c(eps = eps, err = grepl("ERROR", paste(fnres, grres)))
    }
}
tab <- do.call(rbind, tab)
cat("\nerror rate by eps:\n")
print(tapply(tab[, "err"], tab[, "eps"], mean))

## A hand-built explosive point: huge growth rate, tiny carrying capacity.
cat("\n== a deliberately explosive parameter vector ==\n")
bad <- mode
bad[2] <- log(1e6)   # rGrow
bad[5] <- log(1e-6)  # K
cat("fn:", classify(bad, obj$fn), "\n")
cat("gr:", classify(bad, obj$gr), "\n")

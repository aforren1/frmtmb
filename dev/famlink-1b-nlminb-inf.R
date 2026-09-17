# Round 1b, item 2, step 1: where the NaN enters, and what nlminb does
# with +Inf instead. Lane library, no package change measured here.
#
#   Rscript dev/famlink-1b-nlminb-inf.R > dev/famlink-1b-nlminb-inf-log.txt
#
# Design "band" from dev/famlink-punch-invgauss.R (seed 11, n = 400),
# family = "inverse.gaussian" on its 1/mu^2 default. The taped objective
# is taken with dry_run = "objective" and driven by nlminb directly, so
# the only difference between the arms is the value handed back where
# the objective is NaN.
.libPaths(c("C:/Users/adf44/source/r/famlink-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

set.seed(11)
n <- 400
x <- stats::rnorm(n)
mig <- 1 / sqrt(pmax(0.5 + 0.2 * x, 0.05))
d <- data.frame(y = stats::rgamma(n, shape = 20, rate = 20 / mig), x = x)

obj <- frm(y ~ x, d, family = "inverse.gaussian", dry_run = "objective")
cat("class of dry_run objective:", class(obj), "\n")
if (!is.function(obj$fn)) obj <- obj$obj

# where the NaN enters: the mu intercept on a 1/mu^2 link, pushed below 0
p <- obj$par
cat("start par:", format(p), "\n")
q <- p
q[1] <- -abs(q[1]) - 0.5
cat("fn at eta intercept", format(q[1]), ":", obj$fn(q), "\n")
cat("gr there:", tryCatch(format(obj$gr(q)), error = conditionMessage), "\n")

drive <- function(label, fn) {
  nf <- 0L
  nbad <- 0L
  f2 <- function(par) {
    nf <<- nf + 1L
    v <- fn(par)
    if (!is.finite(v)) nbad <<- nbad + 1L
    v
  }
  warns <- character(0)
  opt <- withCallingHandlers(
    stats::nlminb(obj$par, f2, obj$gr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  cat(sprintf(paste("%-22s objective %.12f conv %d evals %d non-finite %d",
                    "warnings %d %s\n"),
              label, opt$objective, opt$convergence, nf, nbad,
              length(warns), paste(unique(warns), collapse = " | ")))
  opt
}
a <- drive("as taped (NaN)", function(par) obj$fn(par))
b <- drive("NaN mapped to +Inf", function(par) {
  v <- obj$fn(par)
  if (is.nan(v)) Inf else v
})
cat("max |par difference| between the two optima:",
    format(max(abs(a$par - b$par))), "\n")
cat("identical par:", identical(a$par, b$par), " identical objective:",
    identical(a$objective, b$objective), "\n")

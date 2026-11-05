# re_sigma_small seed 1: the fit returned the Laplace estimates with a gradient
# of 5.6e13. What happened, and what is the exact answer?
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
args <- commandArgs(TRUE)
seed <- if (length(args)) as.integer(args[[1]]) else 1L
de <- design_list$re_sigma_small(seed)
fit <- withCallingHandlers(
  frm(de$formula, family = de$family, data = de$data, importance = de$draws,
      verbose = TRUE),
  warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n\n")
    invokeRestart("muffleWarning")
  })
str(fit$importance[c("rounds", "moved", "moves", "capped", "grad", "grads",
                     "start_value", "ess_min", "mcse")])
cat("opt code", fit$opt$convergence, fit$opt$message, "\n")
print(fit$opt$par)
print(fit$obj$gr(fit$opt$par))
lap <- suppressWarnings(frm(de$formula, family = de$family, data = de$data))
print(lap$opt$par)
cat("lap ll", as.numeric(logLik(lap)), " imp ll", as.numeric(logLik(fit)), "\n")
print(fit)
print(tryCatch(summary(fit)$fixed, error = function(e) conditionMessage(e)))
print(tryCatch(vcov(fit), error = function(e) conditionMessage(e)))

# The exact objective by integrate() on the group level, from several starts.
gr <- split(seq_len(nrow(de$data)), de$data$g)
dl <- lapply(gr, function(i) de$data[i, , drop = FALSE])
exact_int <- function(par) {
  tot <- 0
  sd <- exp(par[4])
  for (d in dl) {
    # the integrand in u: two rows, exp(-residual^2 exp(-2(s+u)) / 2 - 2(s+u))
    # is a heavy one-sided shape, so integrate piecewise on a wide window
    f <- function(u) {
      exp(de$ll(par, d, u) + dnorm(u, 0, sd, log = TRUE))
    }
    lo <- -10 * sd - 8
    hi <- 10 * sd + 8
    br <- seq(lo, hi, length.out = 21L)
    v <- 0
    for (i in seq_len(length(br) - 1L)) {
      v <- v + integrate(f, br[i], br[i + 1L], rel.tol = 1e-9, abs.tol = 0,
                         subdivisions = 500L, stop.on.error = FALSE)$value
    }
    tot <- tot - log(v)
  }
  tot
}
cat("exact at fit:", exact_int(unname(fit$opt$par)), "\n")
starts <- list(unname(fit$opt$par),
               c(unname(fit$opt$par[1:3]), 0),
               c(unname(fit$opt$par[1:3]), -0.5),
               c(0.2, 0.5, -0.2, 0))
best <- NULL
for (st in starts) {
  o <- nlminb(st, exact_int, control = list(iter.max = 100, rel.tol = 1e-8))
  cat("start", format(st, digits = 3), "->", format(o$par, digits = 5),
      " nll", format(o$objective, digits = 8), "code", o$convergence, "\n")
  if (is.null(best) || o$objective < best$objective) best <- o
}
cat("exact optimum:", format(best$par, digits = 6), "nll",
    format(best$objective, digits = 9), "\n")
cat("fit minus exact optimum, in the parameter:",
    format(unname(fit$opt$par) - best$par, digits = 4), "\n")
cat("exact nll at fit minus at optimum:",
    format(exact_int(unname(fit$opt$par)) - best$objective, digits = 5), "\n")
H <- numDeriv::hessian(exact_int, best$par)
se <- sqrt(diag(solve(H)))
cat("exact SE:", format(se, digits = 4), "\n")
cat("z:", format((unname(fit$opt$par) - best$par) / se, digits = 3), "\n")

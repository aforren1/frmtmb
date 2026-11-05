# Fix round 1, item 4: a REAL allocation failure inside the refusal's ladder. Run
# under a Windows job object that caps this process's memory, from
# review-impsamp-fr1-alloc.ps1. The rounds are displaced by the lane's own wrapper,
# so the ladder runs. After the refusal the same process must still fit a model.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
displaced <- function(below, code) {
  real_plan <- frmtmb:::imp_plan
  real_opt <- frmtmb:::optimize_obj
  last_n <- 0L
  testthat::with_mocked_bindings(
    code,
    imp_plan = function(lap_obj, frame, lay, par, n_draw, seed) {
      last_n <<- n_draw
      real_plan(lap_obj, frame, lay, par, n_draw, seed)
    },
    optimize_obj = function(obj, control, bounds = list(lower = -Inf, upper = Inf),
                            par_units = NULL, verbose = 0L, start_par = obj$par,
                            tally = NULL) {
      o <- real_opt(obj, control, bounds, par_units, verbose, start_par, tally)
      if (last_n > 0L && last_n < below) o$par[[1L]] <- o$par[[1L]] + 1
      o
    },
    .package = "frmtmb")
}
set.seed(1)
g <- factor(rep(seq_len(61L), each = 4))
x <- rnorm(244)
d <- data.frame(y = rpois(244, exp(0.3 + 0.4 * x + rnorm(61, 0, 0.5)[g])), x, g)
cat("start\n"); flush.console()
msg <- tryCatch(displaced(1e9, frm(bf(y ~ x + (1 | g)), family = poisson(), data = d,
                                   importance = 512L)),
                error = function(e) paste("CAUGHT:", conditionMessage(e)))
cat(if (is.character(msg)) msg else "returned a fit", "\n"); flush.console()
gc()
ok <- tryCatch({
  f <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = d, importance = 64L)
  sprintf("after: importance 64 fits, ll %.4f, rounds %d", as.numeric(logLik(f)), f$importance$rounds)
}, error = function(e) paste("after: ERROR", conditionMessage(e)))
cat(ok, "\n")
l <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = d)
cat("after: Laplace ll", as.numeric(logLik(l)), "\n")
cat("END OF SCRIPT\n")

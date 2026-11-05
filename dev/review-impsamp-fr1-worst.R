# The refusal's worst allowed case under the new bound, no memory cap: 244 rows from
# 512 draws, ladder 1024 to 8192, rounds displaced by the lane's wrapper.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
peak_mb <- function() {
  cmd <- sprintf("(Get-Process -Id %d).PeakWorkingSet64 / 1MB", Sys.getpid())
  as.numeric(system2("powershell", c("-NoProfile", "-Command", shQuote(cmd)), stdout = TRUE))
}
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
t0 <- system.time(f <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = d, importance = 512L))[["elapsed"]]
cat(sprintf("the fit at 512: %.1f s peak %.0f MB\n", t0, peak_mb()))
t1 <- system.time(m <- tryCatch(displaced(1e9, frm(bf(y ~ x + (1 | g)), family = poisson(), data = d,
                                                   importance = 512L)),
                                error = conditionMessage))[["elapsed"]]
cat(sprintf("the refusal: %.1f s, %.1f times; peak %.0f MB\n%s\n", t1, t1 / t0, peak_mb(), m))

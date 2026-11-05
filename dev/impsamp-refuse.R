# M2 and M1 on the final build: the refusal's worst case under the 2e6 cell
# bound, and its text after one doubling and after none. The lane's wrapper
# displaces every round, so no count settles. Usage: Rscript dev/impsamp-refuse.R
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
peak_mb <- function() {
  cmd <- sprintf("(Get-Process -Id %d).PeakWorkingSet64 / 1MB", Sys.getpid())
  as.numeric(system2("powershell", c("-NoProfile", "-Command", shQuote(cmd)),
                     stdout = TRUE))
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
mk <- function(G, n, seed = 1L) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), each = n))
  x <- rnorm(G * n)
  data.frame(y = rpois(G * n, exp(0.3 + 0.4 * x + rnorm(G, 0, 0.5)[g])), x, g)
}
f <- bf(y ~ x + (1 | g))
one <- function(lab, d, N) {
  t0 <- system.time(ok <- frm(f, family = poisson(), data = d, importance = N))
  cat(sprintf("%s: rows %d draws %d: fit %.1f s rounds %d, peak %.0f Mb\n", lab,
              nrow(d), N, t0[["elapsed"]], ok$importance$rounds, peak_mb()))
  t1 <- system.time(msg <- tryCatch(displaced(1e9, frm(
    f, family = poisson(), data = d, importance = N)),
    error = conditionMessage))[["elapsed"]]
  cat(sprintf("  refusal %.1f s, %.1f times the fit, peak %.0f Mb\n  %s\n", t1,
              t1 / t0[["elapsed"]], peak_mb(), msg))
}
one("the worst the bound allows, 244 rows", mk(61L, 4L), 512L)
one("the review's 2400 rows, no ladder", mk(100L, 24L), 512L)
one("2400 rows at 64 draws", mk(100L, 24L), 64L)
one("20000 rows, draws over the bound", mk(200L, 100L), 500L)

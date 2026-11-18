# Lane nlband: Vm * x / (K + x) on y = 2 x / (x - 2) + noise with the
# band data's hole, seeds 2 and 3: the hole start at K -2 ran and lost.
# From where does a run in the hole reach the maximum there?
# NLBAND_ARM=after Rscript dev/nlband-dbg2.R
source("dev/nlband-common.R")
fm <- bf(y ~ Vm * x / (K + x), Vm + K ~ 1, nl = TRUE)
for (s in 1:3) {
  d <- band_data(s)
  d$y <- 2 * d$x / (d$x - 2) + rnorm(nrow(d), 0, 0.2)
  r <- tt(frm(fm, data = d))
  cat("SEED", s, "default ll", llv(r), "beta", r$estimates$beta,
      "msg", msum(r), "warn", wsum(r), "\n")
  off <- tt(frm(fm, data = d, control = frmtmb_control(nl_starts = FALSE)))
  cat("  nl_starts FALSE ll", llv(off), "beta", off$estimates$beta, "\n")
  for (vm in c(r$estimates$beta[[1L]], 100, 10, 2, 1, 0.1, 0, -1)) {
    g <- tt(frm(fm, data = d, start = list(beta = c(vm, -2))))
    cat(sprintf("  start Vm %g K -2: ll %.4f code %s beta %s warn %s\n", vm,
                llv(g), if (iserr(g)) "ERR" else g$opt$convergence,
                paste(signif(g$estimates$beta, 5), collapse = " "), wsum(g)))
  }
  # the tape's runs, as the escape makes them, from the fit
  f <- r
  p <- f$opt$par
  for (vm in c(p[[1L]], 1, 0)) {
    q <- p
    q[1:2] <- c(vm, -2)
    o <- tryCatch(ns$optimize_obj(f$obj, f$control, ns$fit_outer_box(f),
                                  f$par_units, start_par = q),
                  error = function(e) conditionMessage(e))
    if (is.character(o)) {
      cat("  run Vm", vm, "ERROR", o, "\n")
    } else {
      cat(sprintf("  run Vm %g: objective %.4f code %s par %s\n", vm,
                  o$objective, o$convergence,
                  paste(signif(o$par, 5), collapse = " ")))
    }
  }
}
